#!/usr/bin/env python3
"""Score protein-only thermal trajectories using native contacts and DSSP retention."""

from __future__ import annotations

import argparse
import glob
import json
import os

import mdtraj as md
import numpy as np
import pandas as pd


def q_native(traj):
    heavy = traj.topology.select_atom_indices("heavy")
    pairs = np.array(
        [
            (a, b)
            for i, a in enumerate(heavy)
            for b in heavy[i + 1 :]
            if abs(traj.topology.atom(a).residue.index - traj.topology.atom(b).residue.index) >= 4
        ]
    )
    if not len(pairs):
        return np.ones(traj.n_frames)
    d0 = md.compute_distances(traj[0], pairs)[0]
    keep = d0 < 0.45
    if not keep.any():
        return np.ones(traj.n_frames)
    d = md.compute_distances(traj, pairs[keep])
    return (1 / (1 + np.exp(50 * (d - 1.8 * d0[keep])))).mean(axis=1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--runs", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--tail-frac", type=float, default=0.2)
    a = p.parse_args()
    rows = []
    for dcd in sorted(glob.glob(f"{a.runs}.dcd")):
        prefix = dcd[:-4]
        top = f"{prefix}_top.pdb"
        meta_path = f"{prefix}_meta.json"
        if not os.path.exists(top) or not os.path.exists(meta_path):
            print(f"skip {prefix}: missing topology/metadata")
            continue
        traj = md.load(dcd, top=top)
        if traj.n_frames < 5:
            print(f"skip {prefix}: <5 frames")
            continue
        with open(meta_path) as f:
            meta = json.load(f)
        traj.superpose(traj, 0)
        q = q_native(traj)
        dssp = md.compute_dssp(traj, simplified=True)
        structured = np.isin(dssp[0], ["H", "E"])
        ss = (
            np.ones(traj.n_frames)
            if not structured.any()
            else (dssp[:, structured] == dssp[0, structured]).mean(axis=1)
        )
        tail = max(1, int(traj.n_frames * a.tail_frac))
        bb = traj.topology.select("backbone")
        rg = md.compute_rg(traj)
        rmsd = md.rmsd(traj, traj, 0, atom_indices=bb)
        rows.append(
            {
                **meta,
                "q_final": float(q[-tail:].mean()),
                "q_min": float(q.min()),
                "rg_ratio": float(rg[-tail:].mean() / rg[0]),
                "ss_keep": float(ss[-tail:].mean()),
                "rmsd_final_nm": float(rmsd[-tail:].mean()),
                "unfolded": bool(q[-tail:].mean() < 0.6),
            }
        )
    per = pd.DataFrame(rows)
    if per.empty:
        raise SystemExit("no complete trajectories to score")
    per.to_csv(a.out.replace(".csv", "_per_run.csv"), index=False)
    score = per.groupby(["design_id", "temp_K"], as_index=False).agg(
        n_replicas=("seed", "count"),
        q_mean=("q_final", "mean"),
        q_sd=("q_final", "std"),
        q_worst=("q_final", "min"),
        ss_keep=("ss_keep", "mean"),
        rg_ratio=("rg_ratio", "mean"),
        rmsd_nm=("rmsd_final_nm", "mean"),
        n_unfolded=("unfolded", "sum"),
    )
    score["q_sd"] = score.q_sd.fillna(0)
    score["thermal_score"] = 0.6 * score.q_worst + 0.3 * score.ss_keep + 0.1 * score.q_mean
    score.sort_values("thermal_score", ascending=False).to_csv(a.out, index=False)
    print(score.sort_values("thermal_score", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
