#!/usr/bin/env python3
"""Build the proof-of-concept MD inputs: top-ranked design, its framework twin, trastuzumab reference.

Reads the top row of designs/best50.csv, finds its Boltz-2 refold (falling back to the designed
cif), and writes PDBs with chains A=VL, B=VH, C/D=HER2 into runs/thermal-poc/inputs and
runs/phmd/inputs. The twin applies VH S85N+R87T (abhack.tickets) to the parent, checking the
wild-type residues first. The reference is the trastuzumab Fv + HER2 domain IV from 1N8Z.

    python scripts/poc_inputs.py [--rank 1]
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import gemmi
import polars as pl
from openmm.app import PDBFile
from pdbfixer import PDBFixer

from abhack.constants import HER2_TEMPLATE_CIF
from abhack.tickets import FRAMEWORK_MUTATIONS

THREE = {"S": "SER", "R": "ARG", "N": "ASN", "T": "THR"}


def _to_pdb(cif: str, out: Path) -> Path:
    st = gemmi.read_structure(cif)
    st.remove_ligands_and_waters()
    out.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(out))
    return out


def _refold_cif(row: dict) -> str:
    did = row["design_id"].removesuffix("-ft")
    run = Path(row["cif"]).parts[1]  # runs/<run>/design/...
    hits = sorted(glob.glob(f"runs/{run}-gated/boltz2/{did}/seed_*/{did}_seed_*_sample_1_model.cif")
                  + glob.glob(f"runs/{run}/boltz2/{did}/seed_*/{did}_seed_*_sample_1_model.cif"))
    return hits[0] if hits else row["cif"]


def _twin(parent_pdb: Path, out: Path) -> Path:
    fx = PDBFixer(filename=str(parent_pdb))
    vh = next(c for c in fx.topology.chains() if c.id == "B")
    res = {int(r.id): r.name for r in vh.residues()}
    muts = []
    for pos, wt, mut in FRAMEWORK_MUTATIONS:
        if res.get(pos) != {"S": "SER", "R": "ARG"}[wt]:
            raise ValueError(f"VH {pos} is {res.get(pos)}, expected {wt}: structure numbering does not match")
        muts.append(f"{THREE[wt]}-{pos}-{THREE[mut]}")
    fx.applyMutations(muts, "B")
    fx.findMissingResidues()
    fx.missingResidues = {}
    fx.findMissingAtoms()
    fx.addMissingAtoms()
    with open(out, "w") as f:
        PDBFile.writeFile(fx.topology, fx.positions, f, keepIds=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank", type=int, default=1)
    ap.add_argument("--best", default="designs/best50.csv")
    a = ap.parse_args()
    row = pl.read_csv(a.best, infer_schema_length=10_000).row(a.rank - 1, named=True)
    did = row["design_id"]
    cif = _refold_cif(row)
    print(f"rank {a.rank}: {did}  CDRH3 {row['cdrh3']}  structure {cif}")

    th = Path("runs/thermal-poc/inputs")
    ph = Path("runs/phmd/inputs")
    parent = _to_pdb(cif, th / f"{did}.pdb")
    twin = _twin(parent, th / f"{did}-ft.pdb")
    ref = _to_pdb(str(HER2_TEMPLATE_CIF), th / "trastuzumab_ref.pdb")
    for src, name in ((parent, did), (twin, f"{did}-ft")):
        (ph / f"{name}.pdb").parent.mkdir(parents=True, exist_ok=True)
        (ph / f"{name}.pdb").write_text(src.read_text())
    print("wrote", sorted(p.name for p in th.glob("*.pdb")))
    Path("runs/thermal-poc/design_id.txt").write_text(did + "\n")


if __name__ == "__main__":
    main()
