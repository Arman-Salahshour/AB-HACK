#!/usr/bin/env python3
"""CUDA-required apo Fv/VHH thermal screen; relative stability only, never Tm."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

from mdtraj.reporters import DCDReporter
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
from openmm.app import PME, ForceField, HBonds, Modeller, PDBFile, Simulation, StateDataReporter
from pdbfixer import PDBFixer


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--design-id", required=True)
    p.add_argument("--binder-chains", default="A,B")
    p.add_argument("--format", choices=("fv", "vhh"), default="fv")
    p.add_argument("--temp", type=float, default=400)
    p.add_argument("--ns", type=float, default=5)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--equil-ps", type=float, default=200)
    p.add_argument("--padding", type=float, default=1.0)
    p.add_argument("--frame-ps", type=float, default=10)
    p.add_argument("--dt-fs", type=float, default=0)
    p.add_argument("--platform", default="CUDA")
    p.add_argument("--allow-no-disulfide", action="store_true")
    a = p.parse_args()
    chains = [x.strip() for x in a.binder_chains.split(",") if x.strip()]
    expected = 2 if a.format == "fv" else 1
    if len(chains) != expected or a.temp < 310:
        p.error("use 2 Fv chains (or 1 VHH chain) and temperature >=310 K")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    dt = a.dt_fs or (4 if a.temp <= 320 else 3)
    fixer = PDBFixer(filename=a.pdb)
    all_chains = list(fixer.topology.chains())
    available = {c.id for c in all_chains}
    if set(chains) - available:
        sys.exit(f"missing binder chains {set(chains) - available}; available={sorted(available)}")
    fixer.removeChains(
        [i for i, c in enumerate(all_chains) if c.id not in chains]
    )  # strips HER2 C,D
    if {c.id for c in fixer.topology.chains()} != set(chains):
        sys.exit("binder-chain extraction failed")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.topology.createDisulfideBonds(modeller.positions)
    disulfides = sum(1 for x, y in modeller.topology.bonds() if x.name == y.name == "SG")
    if disulfides < expected and not a.allow_no_disulfide:
        sys.exit(f"need >= {expected} disulfide(s), found {disulfides}")
    protein_atoms = list(range(modeller.topology.getNumAtoms()))
    with open(f"{a.out}_top.pdb", "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    ff = ForceField("charmm36.xml", "charmm36/water.xml")
    modeller.addSolvent(
        ff,
        model="tip3p",
        padding=a.padding * unit.nanometer,
        ionicStrength=0.15 * unit.molar,
        neutralize=True,
    )
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * unit.nanometer,
        constraints=HBonds,
        hydrogenMass=4 * unit.amu,
    )
    barostat = system.addForce(MonteCarloBarostat(1 * unit.bar, 310 * unit.kelvin, 25))
    integrator = LangevinMiddleIntegrator(
        310 * unit.kelvin, 1 / unit.picosecond, dt * unit.femtosecond
    )
    integrator.setRandomNumberSeed(a.seed)
    try:
        platform = Platform.getPlatformByName(a.platform)
    except Exception as e:  # noqa: BLE001 - OpenMM backend errors are implementation-specific.
        sys.exit(f"OpenMM {a.platform!r} unavailable: {e}; refusing CPU fallback")
    sim = Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=5000)
    sim.context.setVelocitiesToTemperature(310 * unit.kelvin, a.seed)
    sim.step(int(a.equil_ps * 1000 / dt))
    system.removeForce(barostat)
    sim.context.reinitialize(preserveState=True)
    for n in range(1, 11):
        integrator.setTemperature((310 + (a.temp - 310) * n / 10) * unit.kelvin)
        sim.step(int(100000 / dt / 10))
    integrator.setTemperature(a.temp * unit.kelvin)
    stride = max(1, int(a.frame_ps * 1000 / dt))
    sim.reporters.extend(
        [
            DCDReporter(f"{a.out}.dcd", stride, atomSubset=protein_atoms),
            StateDataReporter(
                f"{a.out}.log",
                stride,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=True,
                speed=True,
            ),
        ]
    )
    started = time.time()
    sim.step(int(a.ns * 1e6 / dt))
    wall = time.time() - started
    meta = {
        "design_id": a.design_id,
        "source_pdb": os.path.abspath(a.pdb),
        "binder_chains": chains,
        "format": a.format,
        "temp_K": a.temp,
        "ns": a.ns,
        "seed": a.seed,
        "dt_fs": dt,
        "n_disulfides": disulfides,
        "n_protein_atoms": len(protein_atoms),
        "n_total_atoms": system.getNumParticles(),
        "platform": platform.getName(),
        "wall_seconds": round(wall, 1),
        "ns_per_day": round(a.ns * 86400 / wall, 1) if wall else None,
    }
    with open(f"{a.out}_meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
