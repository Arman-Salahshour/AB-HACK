#!/usr/bin/env python3
"""Complex MD at one pH: Fv + HER2 domain IV, histidine states drawn from bound-state pKa.

Each histidine is protonated (HIP) with probability 1/(1+10^(pH - pKa_bound)), else
neutral (HIE), using the propka pKa of that His *in the complex*. The draw is seeded,
so replicas at the same pH sample different protonation microstates. The two pH points
therefore differ only where the predicted pKa says they should. Fixed-state MD is a
sampling approximation to constant-pH, not a replacement for it.

Writes <out>.dcd, <out>_top.pdb, <out>.json (His states, timing). Analyse with
scripts/analyze_ph_md.py. CUDA required; no CPU fallback.

    python scripts/run_ph_md.py --pdb cplx.pdb --out runs/phmd/ctrl_pH6.8_r1 --ph 6.8 --ns 1 --seed 1
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import tempfile
import time
from pathlib import Path

from mdtraj.reporters import DCDReporter
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit, CustomExternalForce
from openmm.app import PME, ForceField, HBonds, Modeller, PDBFile, Simulation, StateDataReporter
from pdbfixer import PDBFixer

from abhack.phscore import _apo, _propka

ANTIGEN = ("C", "D")


def his_states(pdb: Path, ph: float, rng: random.Random) -> tuple[dict, dict]:
    """Return ({(chain,resnum): 'HIP'|'HIE'}, {(chain,resnum): bound pKa}) for every His."""
    with tempfile.TemporaryDirectory(prefix="phmd_") as d:
        wd = Path(d)
        (wd / "c.pdb").write_text(Path(pdb).read_text())
        pka = _propka(wd / "c.pdb", wd)
    states, pk = {}, {}
    for (chain, num, name), (pka_val, _) in pka.items():
        if name != "HIS":
            continue
        p_prot = 1.0 / (1.0 + 10.0 ** (ph - pka_val))
        states[(chain, num)] = "HIP" if rng.random() < p_prot else "HIE"
        pk[(chain, num)] = pka_val
    return states, pk


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", required=True, help="complex, chains A,B (Fv) + C,D (HER2), heavy atoms")
    p.add_argument("--out", required=True)
    p.add_argument("--ph", type=float, required=True)
    p.add_argument("--ns", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--temp", type=float, default=310)
    p.add_argument("--equil-ps", type=float, default=100)
    p.add_argument("--platform", default="CUDA")
    a = p.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    rng = random.Random(f"{a.seed}-{a.ph}")
    states, pkas = his_states(Path(a.pdb), a.ph, rng)

    fixer = PDBFixer(filename=a.pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.topology.createDisulfideBonds(modeller.positions)

    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    variants = []
    for res in modeller.topology.residues():
        variants.append(states.get((res.chain.id, int(res.id)))) if res.name == "HIS" else variants.append(None)
    modeller.addHydrogens(ff, pH=7.0, variants=variants)
    n_protein = modeller.topology.getNumAtoms()
    chains = {c.id: [r for r in c.residues()] for c in modeller.topology.chains()}
    with open(f"{a.out}_top.pdb", "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    modeller.addSolvent(ff, model="tip3p", padding=1.0 * unit.nanometer,
                        ionicStrength=0.15 * unit.molar, neutralize=True)
    system = ff.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1 * unit.nanometer,
                             constraints=HBonds, hydrogenMass=4 * unit.amu)

    # Weak restraint on HER2 CA atoms far from the Fv so the cropped target cannot unravel,
    # while the epitope and interface stay free.
    import numpy as np
    pos = np.array(modeller.positions.value_in_unit(unit.nanometer))
    fv_idx = [at.index for at in modeller.topology.atoms() if at.residue.chain.id not in ANTIGEN and at.element.symbol != "H"]
    fv_xyz = pos[fv_idx]
    restr = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restr.addGlobalParameter("k", 100.0)
    for q in ("x0", "y0", "z0"):
        restr.addPerParticleParameter(q)
    n_restr = 0
    for at in modeller.topology.atoms():
        if at.residue.chain.id in ANTIGEN and at.name == "CA":
            if np.min(np.linalg.norm(fv_xyz - pos[at.index], axis=1)) > 1.2:  # > 12 A from the Fv
                restr.addParticle(at.index, pos[at.index])
                n_restr += 1
    system.addForce(restr)
    system.addForce(MonteCarloBarostat(1 * unit.bar, a.temp * unit.kelvin, 25))

    dt = 4.0
    integ = LangevinMiddleIntegrator(a.temp * unit.kelvin, 1 / unit.picosecond, dt * unit.femtosecond)
    integ.setRandomNumberSeed(a.seed)
    try:
        plat = Platform.getPlatformByName(a.platform)
    except Exception as e:  # noqa: BLE001
        sys.exit(f"OpenMM {a.platform!r} unavailable: {e}; refusing CPU fallback")
    sim = Simulation(modeller.topology, system, integ, plat)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=5000)
    sim.context.setVelocitiesToTemperature(a.temp * unit.kelvin, a.seed)
    sim.step(int(a.equil_ps * 1000 / dt))
    stride = max(1, int(10 * 1000 / dt))  # 10 ps frames
    sim.reporters.extend([
        DCDReporter(f"{a.out}.dcd", stride, atomSubset=list(range(n_protein))),
        StateDataReporter(f"{a.out}.log", stride, step=True, time=True, potentialEnergy=True,
                          temperature=True, speed=True),
    ])
    t0 = time.time()
    sim.step(int(a.ns * 1e6 / dt))
    meta = {
        "pdb": os.path.abspath(a.pdb), "ph": a.ph, "seed": a.seed, "ns": a.ns, "temp_K": a.temp,
        "his_states": {f"{c}{n}": s for (c, n), s in states.items()},
        "his_pka_bound": {f"{c}{n}": v for (c, n), v in pkas.items()},
        "n_her2_ca_restrained": n_restr, "forcefield": "amber14-all + tip3p",
        "wall_seconds": round(time.time() - t0, 1), "n_total_atoms": modeller.topology.getNumAtoms(),
    }
    json.dump(meta, open(f"{a.out}.json", "w"), indent=1)
    print(json.dumps({k: meta[k] for k in ("ph", "seed", "wall_seconds", "n_total_atoms")}))


if __name__ == "__main__":
    main()
