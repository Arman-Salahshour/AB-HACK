#!/usr/bin/env python3
"""Cross-architecture check of the top designs against HER2 with AF3-family models (JapanFold).

For each design and each model (Mol-1, hosted OpenFold3; unguided, MSA on, no template):

* ``iptm`` and the HER2-VH / HER2-VL pair ipTM from the model's own results;
* interface PAE: mean PAE over the Fv-to-HER2 token blocks, both directions (Angstrom);
* pose error: HER2 is superposed on the trastuzumab crystal (1N8Z), then the RMSD of the VH framework and
  VL C-alpha atoms (CDRs excluded) to the crystal Fv. Small = the same binding mode as trastuzumab;
* epitope overlap: Jaccard index between the HER2 residues the predicted Fv touches (< 4.5 A) and the ones
  trastuzumab touches in the crystal.

A similar predicted orientation is not proof of epitope competition; this is a consistency check.

    python scripts/af3_check.py runs/japanfold -o designs/af3-check.csv
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import gemmi
import numpy as np
import polars as pl

from abhack.epitope import REFERENCE_CIF

CONTACT = 4.5
HER2_OFFSET = 488  # predicted chain A residue p is HER2 residue 488 + p (domain IV, 489-607)
CDR_VH = set(range(26, 36)) | set(range(50, 66)) | set(range(95, 109))  # generous CDR windows, 1-based
CDR_VL = set(range(24, 35)) | set(range(50, 57)) | set(range(89, 98))


def _atoms(st, chains):
    return [(c.name, r, a) for c in st[0] if c.name in chains for r in c for a in r]


def crystal():
    st = gemmi.read_structure(str(REFERENCE_CIF))
    st.setup_entities()
    her2 = {}  # resseq -> CA pos
    for c in st[0]:
        if c.name in ("C", "D"):
            for r in c:
                ca = r.find_atom("CA", "*")
                if ca:
                    her2[r.seqid.num] = ca.pos
    vh = {r.seqid.num: r.find_atom("CA", "*").pos for r in st[0]["B"] if r.find_atom("CA", "*")}
    vl = {r.seqid.num: r.find_atom("CA", "*").pos for r in st[0]["A"] if r.find_atom("CA", "*")}
    ns = gemmi.NeighborSearch(st[0], st.cell, 6).populate()
    contacts = set()
    for ch in ("A", "B"):
        for r in st[0][ch]:
            for a in r:
                for m in ns.find_atoms(a.pos, "\0", radius=CONTACT):
                    cra = m.to_cra(st[0])
                    if cra.chain.name in ("C", "D"):
                        contacts.add(cra.residue.seqid.num)
    return her2, vh, vl, contacts


def pose_and_epitope(pdb: str, ref) -> tuple[float | None, float | None, int]:
    her2_ref, vh_ref, vl_ref, epi_ref = ref
    st = gemmi.read_structure(pdb)
    m = st[0]
    pa, pb = {}, {}
    for ch, store in (("H", {}), ("L", {})):
        store.update({r.seqid.num: r.find_atom("CA", "*").pos for r in m[ch] if r.find_atom("CA", "*")})
        (pa if ch == "H" else pb).update(store)
    her2_pred = {488 + r.seqid.num: r.find_atom("CA", "*").pos for r in m["A"] if r.find_atom("CA", "*")}
    common = sorted(set(her2_pred) & set(her2_ref))
    if len(common) < 20:
        return None, None, 0
    sup = gemmi.superpose_positions([her2_ref[i] for i in common], [her2_pred[i] for i in common])
    tr = sup.transform
    sq = []
    for store, ref_store, cdr in ((pa, vh_ref, CDR_VH), (pb, vl_ref, CDR_VL)):
        for i, p in store.items():
            if i in ref_store and i not in cdr:
                q = tr.apply(p)
                r_ = ref_store[i]
                sq.append((q.x - r_.x) ** 2 + (q.y - r_.y) ** 2 + (q.z - r_.z) ** 2)
    rmsd = float(np.sqrt(np.mean(sq))) if sq else None
    ns = gemmi.NeighborSearch(m, st.cell, 6).populate()
    touched = set()
    for ch in ("H", "L"):
        for r in m[ch]:
            for a in r:
                for mk in ns.find_atoms(a.pos, "\0", radius=CONTACT):
                    cra = mk.to_cra(m)
                    if cra.chain.name == "A":
                        touched.add(488 + cra.residue.seqid.num)
    jac = len(touched & epi_ref) / len(touched | epi_ref) if (touched | epi_ref) else 0.0
    return rmsd, jac, len(touched)


def interface_pae(npz: str, n_her2=119, n_vh=120) -> float:
    pae = np.load(npz)["pae"]
    her2 = slice(0, n_her2)
    fv = slice(n_her2, pae.shape[0])
    return float((pae[her2, fv].mean() + pae[fv, her2].mean()) / 2)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="runs/japanfold")
    ap.add_argument("-o", "--out", default="designs/af3-check.csv")
    ap.add_argument("--models", nargs="+", default=["mol1", "openfold3"])
    a = ap.parse_args()
    ref = crystal()
    rows = []
    for model in a.models:
        res = json.load(open(Path(a.root) / model / "results.json"))
        for r in res:
            did = r["id"]
            pdb = Path(a.root) / model / f"{did}.pdb"
            npz = Path(a.root) / model / f"{did}_pae.npz"
            pair = r.get("pair_chains_iptm") or {}
            rmsd, jac, ntouch = pose_and_epitope(str(pdb), ref) if pdb.exists() else (None, None, 0)
            rows.append({
                "design_id": did, "model": model, "iptm": r.get("iptm"), "plddt": r.get("complex_plddt"),
                "her2_vh_iptm": (pair.get("0") or {}).get("1"), "her2_vl_iptm": (pair.get("0") or {}).get("2"),
                "interface_pae": interface_pae(str(npz)) if npz.exists() else None,
                "pose_rmsd_A": rmsd, "epitope_jaccard": jac, "her2_residues_touched": ntouch,
            })
    df = pl.DataFrame(rows, infer_schema_length=None)
    df.write_csv(a.out)
    order = {d: i for i, d in enumerate(pl.read_csv("designs/best50.csv", infer_schema_length=10_000)["design_id"].to_list())}
    wide = (df.with_columns(rank=pl.col("design_id").replace_strict(order, default=999))
              .sort("rank", "model"))
    with pl.Config(tbl_rows=40, tbl_cols=-1, fmt_str_lengths=20, float_precision=2):
        print(wide.select("rank", "design_id", "model", "iptm", "her2_vh_iptm", "her2_vl_iptm",
                          "interface_pae", "pose_rmsd_A", "epitope_jaccard"))


if __name__ == "__main__":
    main()
