"""Build the ranked submission table `abguard` validates.

Columns: sequence_id, kind (vhh/scfv/paired), heavy_sequence, light_sequence,
design_description. Up to 20 rows, best first.

Usage::

    python -m abhack.submission build designs/*.csv -o designs/submission.csv
    bash scripts/validate.sh designs/submission.csv
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from abhack.constants import MAX_TOTAL_LENGTH, SUBMISSION_COLUMNS

MAX_DESIGNS = 20


def to_submission(
    df: pl.DataFrame,
    *,
    kind: str = "paired",
    prefix: str = "D",
    description: str | None = None,
) -> pl.DataFrame:
    """Map a design table (vh_sequence/vl_sequence/...) to the submission schema."""
    if "vh_sequence" not in df.columns:
        raise ValueError("design table needs a vh_sequence column")
    light = df["vl_sequence"] if "vl_sequence" in df.columns else pl.Series([None] * df.height)

    def describe(row: dict) -> str:
        if description:
            return description
        bits = [f"{row.get('weight', 'rfd4')}"]
        if row.get("track"):
            bits.append(str(row["track"]))
        if row.get("cdrh3"):
            bits.append(f"CDRH3={row['cdrh3']}")
        for col in ("min_ipae_boltz2", "iptm_boltz2", "binder_plddt_boltz2"):
            if row.get(col) is not None:
                bits.append(f"{col}={row[col]:.3g}")
        return "; ".join(bits)

    rows = df.to_dicts()
    return pl.DataFrame(
        {
            "sequence_id": [f"{prefix}{i:02d}" for i in range(1, df.height + 1)],
            "kind": [kind] * df.height,
            "heavy_sequence": df["vh_sequence"],
            "light_sequence": light,
            "design_description": [describe(r) for r in rows],
        }
    ).select(SUBMISSION_COLUMNS)


def check_lengths(df: pl.DataFrame) -> pl.DataFrame:
    """Flag rows breaking the 250 aa total-length filter."""
    return df.with_columns(
        total_length=pl.col("heavy_sequence").str.len_chars()
        + pl.col("light_sequence").fill_null("").str.len_chars()
    ).with_columns(length_ok=pl.col("total_length") <= MAX_TOTAL_LENGTH)


def build(
    tables: list[Path | str],
    *,
    rank_by: str = "min_ipae_boltz2",
    ascending: bool = True,
    only_passing: bool = True,
    limit: int = MAX_DESIGNS,
) -> pl.DataFrame:
    """Concatenate design tables, drop duplicates and failures, rank and cut to 20."""
    frames = []
    for path in tables:
        df = pl.read_csv(path, infer_schema_length=10_000)
        frames.append(df.with_columns(source=pl.lit(str(path))))
    df = pl.concat(frames, how="diagonal_relaxed")

    if only_passing and "passes" in df.columns:
        df = df.filter(pl.col("passes"))
    if "cdrh3" in df.columns:
        df = df.unique(subset=["cdrh3"], keep="first")
    if rank_by in df.columns:
        df = df.sort(rank_by, descending=not ascending, nulls_last=True)
    return df.head(limit)


def main() -> None:  # pragma: no cover - thin CLI
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("command", choices=["build"])
    ap.add_argument("tables", nargs="+")
    ap.add_argument("-o", "--output", default="designs/submission.csv")
    ap.add_argument("--kind", default="paired", choices=["vhh", "scfv", "paired"])
    ap.add_argument("--rank-by", default="min_ipae_boltz2")
    ap.add_argument("--descending", action="store_true", help="higher is better")
    ap.add_argument("--limit", type=int, default=MAX_DESIGNS)
    args = ap.parse_args()

    ranked = build(
        args.tables,
        rank_by=args.rank_by,
        ascending=not args.descending,
        limit=args.limit,
    )
    sub = check_lengths(to_submission(ranked, kind=args.kind))
    bad = sub.filter(~pl.col("length_ok"))
    if bad.height:
        print(f"WARNING: {bad.height} designs exceed {MAX_TOTAL_LENGTH} aa")
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    sub.select(SUBMISSION_COLUMNS).write_csv(out)
    print(f"{sub.height} designs -> {out}\nNow run: bash scripts/validate.sh {out}")


if __name__ == "__main__":  # pragma: no cover
    main()
