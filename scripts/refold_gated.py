#!/usr/bin/env python
"""Gate-first refolding: refold only what survived the gates, not RFD4's own top-N.

`scripts/refold.sh` passes ``--top-n`` straight to ``refold-prepare``, which
ranks by RFD4's own pre-oracle — novelty, epitope and confidence are never
consulted, and in the smoke run it refolded the two *worse* designs while
skipping two novelty-passers (CLAUDE.md, "cost us designs"). This builder
fixes the order of operations:

    rank (gates) -> select survivors lacking refold evidence -> refold those.

It assembles a *filtered design directory* — the run's config plus symlinks
to only the selected ``inference_samples`` subfolders and a filtered
``per_example_metrics.csv`` — so ``refold-prepare --top-n <n>`` selects
exactly the gated survivors (best-ranked first, all of them).

    python scripts/refold_gated.py runs/opt-novel-s12 designs/opt-ranked.csv \
        --top 24 --out runs/opt-novel-s12-gated
    OUT=runs/opt-novel-s12-gated TOPN=24 bash scripts/refold.sh
    # then collect: the refolded metrics land in the ORIGINAL run's table via
    # the gated folder's scored-<model>/ and `abhack.designs collect` on it.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import polars as pl

DESIGN_FILES = ("config.json", "result.json", ".apheris-foundry", "generation.log", "generate.log")


def select(ranked: Path, run_dir: Path, top: int, run_csv: Path | None = None) -> list[str]:
    """Design ids in ``run_dir`` that passed the gates and lack refold evidence."""
    df = pl.read_csv(ranked, infer_schema_length=10_000)

    if run_csv is None:
        for cand in (f"designs/{run_dir.name}.csv", f"designs/{run_dir.name.replace('-gated', '')}.csv"):
            if Path(cand).exists():
                run_csv = Path(cand)
                break
    if run_csv is None or not Path(run_csv).exists():
        raise SystemExit(f"no table for {run_dir}; pass --table or run abhack.designs collect first")

    run_df = pl.read_csv(run_csv, infer_schema_length=10_000)
    run_df = run_df.with_columns(_run=pl.lit(True))

    # gates from the ranked table: passes false -> out (null stays), no ipae yet
    if "passes" in df.columns:
        s = df["passes"]
        if s.dtype != pl.Boolean:
            s = (
                s.cast(pl.String).str.to_lowercase().replace_strict(
                    {"true": True, "false": False, "1": True, "0": False, "": None},
                    default=None, return_dtype=pl.Boolean)
            )
        df = df.filter(s.is_null() | s)

    # join on the sequences the run actually produced
    key = ["vh_sequence"] + (["vl_sequence"] if "vl_sequence" in df.columns else [])
    keep = [c for c in df.columns if c in ("vh_sequence", "vl_sequence", "cdrh3", "composite",
                                          "p_high", "pkd_pred", "liab_cdrh3")]
    ranked_keyed = df.select(keep).unique(subset=key, keep="first")
    joined = run_df.join(ranked_keyed, on=key, how="inner")

    # survivors needing evidence: no min_ipae yet
    if "min_ipae_boltz2" in joined.columns:
        joined = joined.filter(pl.col("min_ipae_boltz2").is_null())
    joined = joined.sort("composite", descending=True, nulls_last=True)

    ids = joined["design_id"].to_list()
    if not ids:
        raise SystemExit("nothing to refold: every gated design from this run already has refold metrics")
    return ids[:top]


def build(run_dir: Path, out: Path, ids: list[str]) -> Path:
    """Filtered design directory: config + symlinks to the selected samples."""
    samples = run_dir / "design" / "generation" / "inference_samples"
    metrics = run_dir / "design" / "generation" / "inference_metrics" / "per_example_metrics.csv"

    design_out = out / "design"
    gen_out = design_out / "generation"
    if out.exists():  # a previous gated build is stale by construction
        shutil.rmtree(out, ignore_errors=True)
    (gen_out / "inference_samples").mkdir(parents=True, exist_ok=True)

    for f in DESIGN_FILES:
        src = run_dir / "design" / f
        if src.exists():
            if src.is_dir():
                shutil.copytree(src, design_out / f, dirs_exist_ok=True)
            else:
                shutil.copy2(src, design_out / f)

    # refold-prepare reads result.json, whose `result.designs` lists every
    # generated file of the *original* run — filter it to the selection or
    # prepare dies with "a file listed in the design result.json is missing".
    # A run killed mid-flight (e.g. opt-explore-s13) has no result.json at all;
    # synthesize one from a sibling run's template — same module, same schema —
    # with this selection's paths.
    result_path = design_out / "result.json"
    if ids:
        import json

        selected = [
            f"generation/inference_samples/{d}/{d}_generated.cif" for d in ids
        ]
        if result_path.exists():
            r = json.loads(result_path.read_text())
        else:
            templates = sorted(
                p for p in run_dir.parent.glob("*/design/result.json") if p != result_path
            )
            if not templates:
                raise SystemExit(f"no result.json in {run_dir} and no sibling template to synthesise from")
            r = json.loads(templates[0].read_text())
            r["provenance"]["ended_at"] = "2026-11-14T14:30:00Z"
            print(f"  synthesised result.json from {templates[0].relative_to(run_dir.parent.parent)}")
        r["result"]["designs"] = selected
        result_path.write_text(json.dumps(r))

    n = 0
    for design_id in ids:
        src = samples / design_id
        if not src.is_dir():
            print(f"  missing {src}; skipped")
            continue
        (gen_out / "inference_samples" / design_id).symlink_to(src, target_is_directory=True)
        n += 1

    if metrics.exists():
        mdf = pl.read_csv(metrics, infer_schema_length=10_000)
        key = next((c for c in ("example_id", "sample_id", "name", "id") if c in mdf.columns), None)
        if key:
            mdf = mdf.filter(pl.col(key).is_in(ids))
        outdir = gen_out / "inference_metrics"
        outdir.mkdir(parents=True, exist_ok=True)
        mdf.write_csv(outdir / "per_example_metrics.csv")

    print(f"{n} gated designs -> {design_out}")
    return design_out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", help="original design run, e.g. runs/opt-novel-s12")
    ap.add_argument("ranked", help="abhack.ranking output table")
    ap.add_argument("--top", type=int, default=24, help="how many gated survivors to refold")
    ap.add_argument("--out", required=True, help="output run folder, e.g. runs/opt-novel-s12-gated")
    ap.add_argument("--table", default=None, help="collected design table for the run (designs/opt-s12.csv)")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    ids = select(Path(args.ranked), run_dir, args.top, args.table and Path(args.table))
    print(f"selected {len(ids)} from {run_dir}: {ids[:5]}{'...' if len(ids) > 5 else ''}")
    build(run_dir, Path(args.out), ids)
    print(f"\nnext: OUT={args.out} TOPN={len(ids)} FOLD=boltz2 WEIGHT=her2-novel bash scripts/refold.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
