#!/usr/bin/env python3
"""Score pH selectivity (6.5 vs 7.4) for every design in the given tables.

Each design is scored on its own designed complex (its ``cif`` column, else looked up in
designs/*.csv or runs/boltzgen). Output is one row per design_id for
``abhack.ranking rank --ph-scores``. A design with no structure gets an empty score, which
the ranking treats as undecided (kept, no credit), never as a failure.

    python scripts/ph_score_pool.py designs/opt-s11.csv designs/opt-s12.csv designs/denovo-ranked.csv \
        -o designs/ph-scores.csv

Twins (``-ft``) reuse the parent's structure; their two framework mutations are far from the
interface, so the score is the parent's by construction.
"""

from __future__ import annotations

import argparse
import glob
import tempfile
from multiprocessing import Pool
from pathlib import Path

import gemmi
import polars as pl

from abhack.phscore import score_variant


def _lookup(did: str) -> str | None:
    for f in glob.glob("designs/*.csv"):
        try:
            d = pl.read_csv(f, infer_schema_length=10_000)
        except Exception:  # noqa: BLE001
            continue
        if {"design_id", "cif"} <= set(d.columns):
            hit = d.filter((pl.col("design_id") == did) & pl.col("cif").is_not_null())
            if hit.height:
                return hit["cif"][0]
    hit = glob.glob(f"runs/boltzgen/design/rank*_{did}.cif")
    return hit[0] if hit else None


def _work(item: tuple[str, str, str, str | None]) -> dict:
    did, vh, vl, cif = item
    cif = cif or _lookup(did.removesuffix("-ft"))
    out = {"design_id": did, "vh_sequence": vh, "vl_sequence": vl, "ph_log10_sel": None, "ph_his_summary": "", "ph_top_groups": "", "ph_error": ""}
    if not cif or not Path(cif).exists():
        out["ph_error"] = "no structure"
        return out
    try:
        st = gemmi.read_structure(cif)
        st.remove_ligands_and_waters()
        antigen = ("A",) if "shard" in did else ("C", "D")
        with tempfile.TemporaryDirectory() as d:
            pdb = Path(d) / "x.pdb"
            st.write_pdb(str(pdb))
            s = score_variant("", "", pdb=pdb, antigen=antigen)
        out.update(ph_log10_sel=s["log10_sel"], ph_his_summary=s["his_summary"], ph_top_groups=s["top_groups"])
    except Exception as e:  # noqa: BLE001 - record and continue
        out["ph_error"] = f"{type(e).__name__}: {e}"[:160]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("tables", nargs="+")
    ap.add_argument("-o", "--out", default="designs/ph-scores.csv")
    ap.add_argument("-j", type=int, default=12)
    a = ap.parse_args()
    # design_id repeats across runs (her2_cdrh3_28 exists in several), so the key is the sequence pair.
    items: dict[tuple[str, str], tuple] = {}
    for t in a.tables:
        df = pl.read_csv(t, infer_schema_length=10_000)
        cif = df["cif"] if "cif" in df.columns else pl.Series([None] * df.height)
        vl = df["vl_sequence"] if "vl_sequence" in df.columns else pl.Series([None] * df.height)
        for did, vh_, vl_, c in zip(df["design_id"], df["vh_sequence"], vl, cif, strict=True):
            items.setdefault((vh_ or "", vl_ or ""), (did, vh_ or "", vl_ or "", c))
    with Pool(a.j) as p:
        rows = p.map(_work, list(items.values()), chunksize=4)
    out = pl.DataFrame(rows, infer_schema_length=None)
    out.write_csv(a.out)
    ok = out.filter(pl.col("ph_log10_sel").is_not_null())
    print(f"{ok.height}/{out.height} scored -> {a.out}; unscored: {out.height - ok.height}")


if __name__ == "__main__":
    main()
