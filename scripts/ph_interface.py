#!/usr/bin/env python3
"""Map the trastuzumab paratope against HER2 acidic residues for His-switch design.

A histidine switch needs a protonated His (+) to pair with an Asp/Glu (-) on the
antigen, so a candidate position is a paratope residue that contacts HER2 or whose
side chain sits within reach of a HER2 carboxylate.
Reference: her2_fv_1n8z.cif (chains A=VL, B=VH, C/D=HER2).

    python scripts/ph_interface.py [out.csv]
"""

from __future__ import annotations

import sys

import gemmi
import polars as pl

from abhack.epitope import ANTIGEN_CHAINS, REFERENCE_CIF

CONTACT = 4.5  # A, heavy-atom contact
REACH = 7.0  # A, CB to carboxylate oxygen: a His could plausibly reach
ACID_ATOMS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}


def main(out: str = "designs/ph-interface.csv") -> None:
    st = gemmi.read_structure(str(REFERENCE_CIF))
    st.setup_entities()
    m = st[0]
    ns = gemmi.NeighborSearch(m, st.cell, 6).populate()
    antigen = [c for c in m if c.name in ANTIGEN_CHAINS]
    acids = [
        (c.name, r.seqid.num, r.name, a.pos)
        for c in antigen
        for r in c
        if r.name in ACID_ATOMS
        for a in r
        if a.name in ACID_ATOMS[r.name]
    ]
    her2_his = [(c.name, r.seqid.num) for c in antigen for r in c if r.name == "HIS"]
    rows = []
    for ch in m:
        if ch.name not in ("A", "B"):
            continue
        for r in ch:
            contacts: set[tuple[str, int, str]] = set()
            mind = 99.0
            for a in r:
                for mk in ns.find_atoms(a.pos, "\0", radius=CONTACT):
                    cra = mk.to_cra(m)
                    if cra.chain.name in ANTIGEN_CHAINS:
                        contacts.add((cra.chain.name, cra.residue.seqid.num, cra.residue.name))
                        mind = min(mind, a.pos.dist(cra.atom.pos))
            ref = next((a for a in r if a.name == "CB"), next((a for a in r if a.name == "CA"), None))
            if ref is None:
                continue
            best = min(((ref.pos.dist(p), c, n, nm) for c, n, nm, p in acids), default=(99.0, "", 0, ""))
            d_acid, ac_ch, ac_n, ac_nm = best
            if contacts or d_acid <= REACH:
                rows.append(
                    {
                        "chain": ch.name,
                        "resnum": r.seqid.num,
                        "resname": r.name,
                        "min_dist_to_her2": round(mind, 2) if mind < 99 else None,
                        "n_her2_contacts": len(contacts),
                        "contacts": ";".join(f"{c}{n}{nm}" for c, n, nm in sorted(contacts)),
                        "nearest_acid": f"{ac_ch}{ac_n}{ac_nm}",
                        "cb_to_acid_A": round(d_acid, 2),
                    }
                )
    df = pl.DataFrame(rows).sort(["chain", "resnum"])
    df.write_csv(out)
    print(f"{len(df)} candidate paratope positions; {len(acids)} HER2 carboxylate atoms; HER2 His: {her2_his}")
    with pl.Config(tbl_rows=80, tbl_cols=-1, fmt_str_lengths=60):
        print(df.filter(pl.col("n_her2_contacts") > 0))


if __name__ == "__main__":
    main(*sys.argv[1:])
