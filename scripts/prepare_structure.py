#!/usr/bin/env python
"""Turn any structure file into an mmCIF the RFD4-Proteina module will accept.

The module takes mmCIF only, and rejects a file whose polymer records are
incomplete — with nothing more than `error[input-data]: Invalid RFD4-Proteina
design input`, so the cause has to be found by comparison. A PDB file carries no
`label_*` records at all, so a naive gemmi conversion writes rows like

    ATOM 1 N N . CYS Cxp C  .  ... 489 C 1     <- label_asym_id, entity, seq_id

where a usable file has

    ATOM 1 N N . CYS C   3  489 ... 489 C 1

This module writes the second form: hydrogens, waters, ligands and alternative
conformations removed, subchains named after their chain, `label_seq_id` taken
from the author numbering, and entities numbered and mapped to those subchains.

    python scripts/prepare_structure.py in.pdb out.cif --chain C --range 489-607
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import gemmi


def prepare(
    src: Path | str,
    dest: Path | str,
    *,
    chains: list[str] | None = None,
    residue_range: tuple[int, int] | None = None,
    rename_chain: str | None = None,
) -> gemmi.Structure:
    """Clean a structure and write it as a module-ready mmCIF."""
    st = gemmi.read_structure(str(src))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_ligands_and_waters()
    st.setup_entities()

    if chains or residue_range:
        lo, hi = residue_range or (None, None)
        parts = []
        for name in chains or [c.name for c in st[0]]:
            span = f"{lo}-{hi}" if residue_range else "*"
            parts.append(f"/1/{name}/{span}")
        st = gemmi.Selection(" or ".join(parts)).copy_structure_selection(st)

    st.setup_entities()
    if rename_chain and len(st[0]) == 1:
        st[0][0].name = rename_chain

    # Restore the polymer bookkeeping a PDB cannot carry.
    for chain in st[0]:
        for res in chain:
            res.subchain = chain.name
            res.label_seq = res.seqid.num
    for i, entity in enumerate(st.entities, start=1):
        entity.name = str(i)
        entity.subchains = [s[0] for s in entity.subchains]

    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    st.make_mmcif_document().write_file(str(dest))
    return st


def main() -> int:  # pragma: no cover - thin CLI
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("src")
    ap.add_argument("dest")
    ap.add_argument("--chain", action="append", dest="chains")
    ap.add_argument("--range", dest="residue_range", help="e.g. 489-607")
    ap.add_argument("--rename-chain")
    args = ap.parse_args()

    rng = None
    if args.residue_range:
        lo, _, hi = args.residue_range.partition("-")
        rng = (int(lo), int(hi))

    st = prepare(
        args.src, args.dest, chains=args.chains, residue_range=rng,
        rename_chain=args.rename_chain,
    )
    for chain in st[0]:
        nums = [r.seqid.num for r in chain]
        atoms = sum(len(r) for r in chain)
        print(f"chain {chain.name}: {len(nums)} res {nums[0]}-{nums[-1]}, {atoms} atoms")
    print(f"-> {args.dest}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
