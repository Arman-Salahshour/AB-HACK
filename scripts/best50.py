#!/usr/bin/env python3
"""Funnel every generated optimisation design down to the best N under the full filter stack.

    gates (novelty, interface PAE) -> pH gate (propka selectivity >= 0) -> composite -> top N

De-duplicates on the (VH, VL) pair first: design_id repeats across runs, and the same sequence
can come back from different seeds. Writes designs/best50.csv and prints the count at each stage.
Does not touch designs/submission.csv.

    python scripts/best50.py -n 50
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import polars as pl

from abhack import ranking

TABLES = ["designs/opt-s11.csv", "designs/opt-s12.csv", "designs/opt-s13.csv"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=50)
    ap.add_argument("--ph-scores", default="designs/ph-scores.csv")
    ap.add_argument("--ph-md", default=None, help="runs/phmd/ph_md_summary.csv when real MD results exist")
    ap.add_argument("--rescore", action="store_true", help="recompute pH scores even if the file exists")
    ap.add_argument("-o", "--out", default="designs/best50.csv")
    a = ap.parse_args()

    if a.rescore or not Path(a.ph_scores).exists():
        subprocess.run([sys.executable, "scripts/ph_score_pool.py", *TABLES, "-o", a.ph_scores], check=True)

    raw = pl.concat([pl.read_csv(t, infer_schema_length=10_000).with_columns(source=pl.lit(t)) for t in TABLES],
                    how="diagonal_relaxed")
    uniq = raw.unique(subset=["vh_sequence", "vl_sequence"], keep="first", maintain_order=True)
    print(f"generated rows          {raw.height}")
    print(f"unique sequence pairs   {uniq.height}")

    # stage counts, applying each filter in turn on a copy
    stages = [("unique", uniq.height)]
    tmp = uniq
    nov = tmp.filter(pl.col("passes").cast(pl.String).str.to_lowercase() != "false")
    stages.append(("novelty (undecided kept)", nov.height))
    ip = pl.col("min_ipae_boltz2").cast(pl.Float64, strict=False)
    pae = nov.filter(ip.is_null() | (ip <= ranking.IPAE_GATE))
    stages.append((f"interface PAE <= {ranking.IPAE_GATE} (unrefolded kept)", pae.height))
    ph = ranking.add_ph_scores(pae, a.ph_scores, a.ph_md)
    scored = ph.filter(pl.col("ph_log10_sel").is_not_null()).height
    gated = ranking.gate(ph, max_ipae=None, ph_min=ranking.PH_MIN_SELECTIVITY)
    stages.append((f"pH gate >= {ranking.PH_MIN_SELECTIVITY} (scored {scored}, unscored kept)", gated.height))

    # rank the survivors through the standard path so the composite is identical to production
    tmp_in = Path("designs/_best50_pool.csv")
    uniq.write_csv(tmp_in)
    try:
        ranked = ranking.rank_table([tmp_in], ph_scores=a.ph_scores, ph_md=a.ph_md)
    finally:
        tmp_in.unlink(missing_ok=True)
    best = ranked.head(a.n)
    stages.append((f"best {a.n} by composite", best.height))
    best.write_csv(a.out)

    print()
    for name, n in stages:
        print(f"  {n:4d}  {name}")
    show = [c for c in ("design_id", "cdrh3", "composite", "min_ipae_boltz2", "p_high", "kd_pred_nM",
                        "ph_log10_sel", "liab_cdrh3", "source") if c in best.columns]
    with pl.Config(tbl_rows=12, tbl_cols=-1, fmt_str_lengths=40):
        print(best.select(show).head(10))
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
