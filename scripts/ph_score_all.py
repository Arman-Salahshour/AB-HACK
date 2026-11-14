#!/usr/bin/env python3
"""Score every pH-switch candidate with propka proton linkage; writes designs/ph-scored.csv."""
import sys
from multiprocessing import Pool
import polars as pl
from abhack.phscore import score_variant

def work(r):
    try:
        s = score_variant(r["vh_sequence"], r["vl_sequence"])
        return {**r, **{f"ph_{k}": v for k, v in s.items()}, "ph_error": ""}
    except Exception as e:  # noqa: BLE001 - keep going, record the failure
        return {**r, "ph_error": f"{type(e).__name__}: {e}"[:200]}

if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "designs/ph-candidates.csv"
    out = sys.argv[2] if len(sys.argv) > 2 else "designs/ph-scored.csv"
    rows = pl.read_csv(src).to_dicts()
    with Pool(12) as p:
        res = p.map(work, rows, chunksize=4)
    df = pl.DataFrame(res, infer_schema_length=None)
    df.write_csv(out)
    ok = df.filter(pl.col("ph_error") == "")
    print(f"{ok.height}/{df.height} scored -> {out}; errors: {df.height - ok.height}")
    print(ok.select("ph_log10_sel").describe())
    top = ok.sort("ph_log10_sel", descending=True).head(10)
    with pl.Config(tbl_rows=20, tbl_cols=-1, fmt_str_lengths=70):
        print(top.select("design_id", "cdrh3", "ph_his_sites", "p_high", "ph_log10_sel", "ph_his_summary"))
