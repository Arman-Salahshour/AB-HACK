#!/usr/bin/env python
"""Structural stability and interface quality from the structures we already have.

Computed per design, from the Boltz-2 refold and the RFD4 generation metrics:

* **CDRH3 pLDDT** — per-residue confidence over the redesigned loop, from the
  refolded structure's B-factor column. A loop the folding model is unsure of is
  a loop that moves.
* **loop vs framework pLDDT** — the gap says whether uncertainty is local to the
  design or general.
* **buried surface area** — how much surface the interface removes from solvent
  (Shrake-Rupley), the usual proxy for interface size.
* **interface contacts** — heavy-atom pairs under 4.5 A, plus RFD4's own counts
  of hydrogen bonds, salt bridges and hydrophobic contacts.
* **secondary structure content** — from RFD4's metrics.

This is not a thermostability prediction. It measures how well-determined and how
extensive the designed interface is, which is what our data supports.

    python scripts/structural_stability.py --run runs/opt-epitope-s21 -o designs/s21-stability.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import gemmi
import polars as pl

CDRH3 = range(99, 109)          # chain B, trastuzumab numbering
ANTIBODY = ("A", "B")
ANTIGEN = ("C", "D")
CONTACT = 4.5


def plddt(pdb_or_cif: Path) -> dict:
    st = gemmi.read_structure(str(pdb_or_cif))
    st.setup_entities()
    if "B" not in [c.name for c in st[0]]:
        return {}
    ch = st[0]["B"]
    loop = [a.b_iso for r in ch if r.seqid.num in CDRH3 for a in r]
    fw = [a.b_iso for r in ch if r.seqid.num not in CDRH3 for a in r]
    if not loop or not fw:
        return {}
    return {
        "cdrh3_plddt": round(sum(loop) / len(loop), 2),
        "framework_plddt": round(sum(fw) / len(fw), 2),
        "plddt_gap": round(sum(fw) / len(fw) - sum(loop) / len(loop), 2),
        "cdrh3_plddt_min": round(min(loop), 1),
    }


def interface(cif: Path) -> dict:
    """Buried surface area and heavy-atom contacts across the interface."""
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.remove_hydrogens()
    model = st[0]
    names = [c.name for c in model]

    from Bio.PDB import ShrakeRupley as _SR
    from Bio.PDB.MMCIFParser import MMCIFParser

    def sasa(keep: tuple[str, ...]) -> float:
        """Solvent-accessible area of the named chains, via Biopython."""
        bio = MMCIFParser(QUIET=True).get_structure("s", str(cif))
        mdl = bio[0]
        for ch in [c.id for c in mdl if c.id not in keep]:
            mdl.detach_child(ch)
        sr = _SR()
        sr.compute(mdl, level="C")
        return float(sum(c.sasa for c in mdl))

    ab = tuple(c for c in ANTIBODY if c in names)
    ag = tuple(c for c in ANTIGEN if c in names)
    if not ab or not ag:
        return {}
    bsa = (sasa(ab) + sasa(ag) - sasa(ab + ag)) / 2

    ab_atoms = [a for cn in ab for r in model[cn] for a in r]
    n_contacts = 0
    loop_contacts = 0
    for cn in ag:
        for res in model[cn]:
            for atom in res:
                for b in ab_atoms:
                    if atom.pos.dist(b.pos) <= CONTACT:
                        n_contacts += 1
    for cn in ag:
        for res in model[cn]:
            for atom in res:
                for r2 in model["B"]:
                    if r2.seqid.num not in CDRH3:
                        continue
                    for b in r2:
                        if atom.pos.dist(b.pos) <= CONTACT:
                            loop_contacts += 1
    return {
        "buried_sasa_A2": round(bsa, 1),
        "interface_atom_contacts": n_contacts,
        "cdrh3_atom_contacts": loop_contacts,
        "cdrh3_contact_share": round(loop_contacts / n_contacts, 3) if n_contacts else 0.0,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("--limit", type=int, default=12)
    args = ap.parse_args()

    run = Path(args.run)
    metrics = next(run.glob("**/inference_metrics/per_example_metrics.csv"), None)
    m = pl.read_csv(metrics, infer_schema_length=10_000) if metrics else None

    rows = []
    for folded in sorted(run.glob("boltz2/*/**/*_model.cif"))[: args.limit]:
        did = folded.parts[folded.parts.index("boltz2") + 1]
        designed = next(run.glob(f"**/inference_samples/{did}/{did}_generated.cif"), None)
        row = {"design_id": did, **plddt(folded)}
        if designed:
            row |= interface(designed)
        if m is not None:
            mm = m.filter(pl.col("example_id") == did)
            if mm.height:
                for src, dst in (
                    ("interchain_interactions.hbond_atom_count_coi", "hbonds"),
                    ("interchain_interactions.saltbridge_residue_count_coi", "salt_bridges"),
                    ("interchain_interactions.hydrophobic_atom_count_coi", "hydrophobic"),
                    ("sse_distribution.sheet_fraction_coi", "sheet_fraction"),
                    ("sse_distribution.loop_fraction_coi", "loop_fraction"),
                ):
                    if src in mm.columns:
                        row[dst] = mm[src][0]
        rows.append(row)

    df = pl.DataFrame(rows)
    pl.Config.set_tbl_cols(15)
    print(df.to_pandas().to_string(index=False))
    if args.output:
        df.write_csv(args.output)
        print(f"\n-> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
