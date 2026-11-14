#!/usr/bin/env python3
"""Build the hand-in CSV (odyssey.csv) from designs/final-30.csv - the same file the results page is built from.

Columns: sequence_id,kind,heavy_sequence,light_sequence,design_description,track. Row order is the rank order
(d01 = rank 1). design_description says what the design is and why it sits where it does.

Then cross-checks every row against docs/final-ranking.html (rank, CDRH3, track) so the file and the page cannot drift.

    python scripts/make_odyssey_csv.py
"""

from __future__ import annotations

import csv
import re
import shutil
from pathlib import Path

import polars as pl

SRC = Path("designs/final-30.csv")
OUT = Path("designs/odyssey.csv")
HOME_COPY = Path("/home/hackathon/odyssey.csv")
PAGE = Path("docs/final-ranking.html")


def describe(r: dict, sid: dict[int, str]) -> str:
    rank = r["rank"]
    if r["track"] == "denovo":
        m = re.search(r"ipTM ([\d.]+), interface PAE (\d+) A", r["why"] or "")
        ev = f"AF2-IG ipTM {m.group(1)}, interface PAE {m.group(2)} A, below the credibility bar" if m else "no independent structure check passed"
        return f"de novo VHH against HER2, rank {rank}; {ev}; included last as a lottery ticket"
    base = (f"CDRH3 {r['cdrh3']} on trastuzumab framework; binder probability {r['p_high']:.2f}, predicted KD {r['kd_pred_nM']:.0f} nM, "
            f"Boltz-2 interface PAE {r['ipae_boltz2']:.2f} A, Mol-1 pose {r['mol1_pose_A']:.1f} A from trastuzumab, "
            f"pH 6.5/7.4 selectivity {r['ph_sel']:+.2f}, "
            f"{int(r['liab_cdrh3'])} CDRH3 liabilit{'y' if int(r['liab_cdrh3']) == 1 else 'ies'}")
    if r["design_id"].endswith("-ft"):
        parent = next(k for k, v in sid.items() if v == r["design_id"][:-3])
        return (f"{base}; framework twin of {sid_label(parent)}: adds VH S85N+R87T "
                f"(published: higher Tm with HER2 binding retained)")
    return f"{base}; rank {rank} by binding evidence"


def sid_label(rank: int) -> str:
    return f"d{rank:02d}"


def main() -> None:
    df = pl.read_csv(SRC, infer_schema_length=10_000).sort("rank")
    # parents by design_id, for twin descriptions (a twin immediately follows its parent)
    sid = {r["rank"]: r["design_id"] for r in df.to_dicts()}
    rows = []
    for r in df.to_dicts():
        rows.append({
            "sequence_id": sid_label(r["rank"]),
            "kind": r["kind"],
            "heavy_sequence": r["vh"],
            "light_sequence": r["vl"] or "",
            "design_description": describe(r, sid),
            "track": r["track"],
        })
    cols = ["sequence_id", "kind", "heavy_sequence", "light_sequence", "design_description", "track"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    shutil.copy(OUT, HOME_COPY)
    print(f"wrote {OUT} and {HOME_COPY}: {len(rows)} rows from {SRC}")

    # cross-check against the page
    html = PAGE.read_text()
    trs = re.findall(r'<tr class="[^"]*"><td class="rk">(\d+)</td><td class="id"><span class="m">([^<]+)</span><br><span class="tag [^"]*">([^<]+)</span></td><td class="m">([^<]*)</td>', html)
    assert len(trs) == len(rows), f"page has {len(trs)} rows, csv {len(rows)}"
    bad = []
    for (pr, pid, ptag, pcdr3), r, src in zip(trs, rows, df.to_dicts(), strict=True):
        ptrack = "denovo" if ptag == "de novo" else "optimization"
        if not (int(pr) == src["rank"] and pid == src["design_id"] and ptrack == r["track"] and pcdr3 == (src["cdrh3"] or "")):
            bad.append((pr, pid))
        # the CDRH3 in the page must be the one inside the submitted heavy chain
        if r["kind"] == "paired" and pcdr3 not in r["heavy_sequence"]:
            bad.append((pr, "cdrh3 not in heavy_sequence"))
    print("page <-> csv cross-check:", "all 30 rows match (rank, design, track, CDRH3 in sequence)" if not bad else f"MISMATCH {bad}")


if __name__ == "__main__":
    main()
