#!/usr/bin/env python3
"""Interface stability from the two-pH complex MD.

Per trajectory: fraction of the initial Fv-HER2 heavy-atom contacts (< 4.5 A) still present
in the last half (Q_if), interface Fv-HER2 minimum distance, and His-carboxylate salt-bridge
occupancy (His ND1/NE2 to Asp/Glu O within 4 A) for every Fv histidine. Then pH 6.8 minus
pH 7.2 per design: positive dQ means the interface holds better at the lower pH.

    python scripts/analyze_ph_md.py "runs/phmd/*_pH*_r*" -o runs/phmd/ph_md_summary.csv
"""

from __future__ import annotations

import argparse
import glob
import json
import re
import sys
import warnings

import mdtraj as md
import numpy as np
import polars as pl

warnings.filterwarnings("ignore")
ANTIGEN = ("C", "D")


def analyse(stem: str) -> dict:
    meta = json.load(open(f"{stem}.json"))
    traj = md.load(f"{stem}.dcd", top=f"{stem}_top.pdb")
    top = traj.topology
    fv = [a.index for a in top.atoms if a.residue.chain.chain_id not in ANTIGEN and a.element.symbol != "H"]
    ag = [a.index for a in top.atoms if a.residue.chain.chain_id in ANTIGEN and a.element.symbol != "H"]
    pairs = np.array([(i, j) for i in fv for j in ag])
    d0 = md.compute_distances(traj[0], pairs)[0]
    native = pairs[d0 < 0.45]
    half = max(1, traj.n_frames // 2)
    d = md.compute_distances(traj[half:], native)
    q = float((d < 0.55).mean())  # 1.0 A tolerance on 4.5 A contacts
    dmin = md.compute_distances(traj[half:], pairs[:: max(1, len(pairs) // 40000)]).min(axis=1).mean()
    sb = {}
    for res in top.residues:
        if res.name != "HIS" or res.chain.chain_id in ANTIGEN:
            continue
        n_atoms = [a.index for a in res.atoms if a.name in ("ND1", "NE2")]
        acid = [a.index for a in top.atoms if a.residue.chain.chain_id in ANTIGEN
                and a.residue.name in ("ASP", "GLU") and a.name in ("OD1", "OD2", "OE1", "OE2")]
        if not n_atoms or not acid:
            continue
        pp = np.array([(i, j) for i in n_atoms for j in acid])
        occ = float((md.compute_distances(traj[half:], pp).min(axis=1) < 0.40).mean())
        sb[f"{res.chain.chain_id}{res.resSeq}"] = round(occ, 3)
    m = re.search(r"([^/]+)_pH([\d.]+)_r(\d+)$", stem)
    return {"design": m.group(1), "ph": float(m.group(2)), "rep": int(m.group(3)), "ns": meta["ns"],
            "n_native_contacts": int(len(native)), "q_interface": round(q, 4),
            "mean_min_dist_nm": round(float(dmin), 3), "his_acid_occupancy": json.dumps(sb),
            "n_his_protonated": sum(v == "HIP" for v in meta["his_states"].values())}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pattern")
    ap.add_argument("-o", "--out", default="runs/phmd/ph_md_summary.csv")
    a = ap.parse_args()
    stems = sorted({p[:-5] for p in glob.glob(a.pattern + ".json")})
    if not stems:
        sys.exit("no trajectories matched")
    df = pl.DataFrame([analyse(s) for s in stems])
    df.write_csv(a.out)
    with pl.Config(tbl_rows=60, tbl_cols=-1):
        print(df)
    piv = (df.group_by("design", "ph").agg(pl.col("q_interface").mean().alias("q"), pl.len().alias("n"))
             .pivot(on="ph", index="design", values="q"))
    cols = {str(c) for c in piv.columns}
    if {"6.8", "7.2"} <= cols:
        piv = piv.with_columns((pl.col("6.8") - pl.col("7.2")).alias("dQ_6.8_minus_7.2"))
    print(piv)


if __name__ == "__main__":
    main()
