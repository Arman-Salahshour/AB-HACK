"""Build the ranked submission table `abguard` validates.

Columns: sequence_id, kind (vhh/scfv/paired), heavy_sequence, light_sequence,
design_description. Up to 20 rows, best first.

Usage::

    python -m abhack.submission build designs/*.csv -o designs/submission.csv
    bash scripts/validate.sh designs/submission.csv

Three rules this module enforces, learned the hard way (see CLAUDE.md):

* a **null** ``passes`` means "the local check could not decide", not "failed" — a de
  novo VHH has no trastuzumab CDRH3 window, so only `abguard` can grade its CDR3.
  Dropping nulls silently deletes the whole de novo track;
* designs are deduplicated on their **sequences**, not on ``cdrh3``, which is null for
  every VHH and would collapse them all to one row;
* ``kind`` is derived **per row** from the design's own ``format``, because a mixed de
  novo + optimisation submission otherwise mislabels one track and fails abguard's
  per-partition schema.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from abhack.constants import MAX_TOTAL_LENGTH, SUBMISSION_COLUMNS
from abhack.novelty import unpaired_cysteines

MAX_DESIGNS = 20

#: `format` as `abhack.designs.parse_design` emits it -> abguard's `kind` partition.
FORMAT_TO_KIND = {"paired": "paired", "vhh": "vhh", "scfv": "scfv"}


def _as_bool(df: pl.DataFrame, column: str) -> pl.Series:
    """``column`` as a nullable Boolean, whatever the CSV round-trip made of it.

    An all-empty column reads back as Null or String, so a plain ``filter`` on it
    either raises or silently drops everything.
    """
    s = df[column]
    if s.dtype == pl.Boolean:
        return s
    if s.dtype == pl.Null:
        return pl.Series(column, [None] * df.height, dtype=pl.Boolean)
    return (
        s.cast(pl.String)
        .str.to_lowercase()
        .replace_strict(
            {"true": True, "false": False, "1": True, "0": False, "": None},
            default=None,
            return_dtype=pl.Boolean,
        )
    )


def resolve_kind(df: pl.DataFrame, *, default: str = "paired") -> pl.Series:
    """One abguard ``kind`` per row, from the design's own ``format`` column."""
    if "format" not in df.columns:
        return pl.Series("kind", [default] * df.height, dtype=pl.String)
    fmt = df["format"].cast(pl.String).str.to_lowercase()
    return pl.Series(
        "kind",
        [FORMAT_TO_KIND.get(f, default) for f in fmt.to_list()],
        dtype=pl.String,
    )


def to_submission(
    df: pl.DataFrame,
    *,
    kind: str | None = None,
    prefix: str = "D",
    description: str | None = None,
) -> pl.DataFrame:
    """Map a design table (vh_sequence/vl_sequence/...) to the submission schema.

    ``kind`` forces one value for every row; left as None (the default) it is derived
    per row from ``format``.
    """
    if "vh_sequence" not in df.columns:
        raise ValueError("design table needs a vh_sequence column")
    light = df["vl_sequence"] if "vl_sequence" in df.columns else pl.Series([None] * df.height)
    kinds = (
        pl.Series("kind", [kind] * df.height, dtype=pl.String)
        if kind is not None
        else resolve_kind(df)
    )

    def describe(row: dict) -> str:
        if description:
            return description
        bits = [str(row.get("weight") or "rfd4")]
        if row.get("framework_ticket"):
            bits.append("fw S85N+R87T (Ng mAbs 2025: +Tm, binding retained)")
        if row.get("track"):
            bits.append(str(row["track"]))
        if row.get("cdrh3"):
            bits.append(f"CDRH3={row['cdrh3']}")
        elif row.get("cdr3_imgt"):
            bits.append(f"CDR3={row['cdr3_imgt']}")
        for col in ("min_ipae_boltz2", "binder_plddt_boltz2", "epitope_fraction",
                    "p_high", "af2ig_iptm", "af2ig_ipae"):
            if row.get(col) is not None:
                bits.append(f"{col}={row[col]:.3g}")
        return "; ".join(bits)

    rows = df.to_dicts()
    return pl.DataFrame(
        {
            "sequence_id": [f"{prefix}{i:02d}" for i in range(1, df.height + 1)],
            "kind": kinds,
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


def enforce_hard_filters(df: pl.DataFrame, *, verbose: bool = True) -> pl.DataFrame:
    """Drop rows abguard would reject on format-independent grounds.

    The 250-residue cap and the unpaired-cysteine rule were both already computed
    upstream and then discarded; a design that breaks either is never sent to the lab,
    so there is no reason to spend a submission slot on it.
    """
    before = df.height
    heavy = df["vh_sequence"] if "vh_sequence" in df.columns else None
    if heavy is None:
        return df

    light = (
        df["vl_sequence"].cast(pl.String).fill_null("")
        if "vl_sequence" in df.columns
        else pl.Series([""] * df.height)
    )
    total = heavy.cast(pl.String).fill_null("").str.len_chars() + light.str.len_chars()
    # abguard catches an unpaired cysteine anywhere in the submitted polypeptide.
    cys_ok = [
        not (unpaired_cysteines(h or "") or unpaired_cysteines(l or ""))
        for h, l in zip(heavy.to_list(), light.to_list())
    ]
    keep = (total <= MAX_TOTAL_LENGTH) & pl.Series(cys_ok)
    out = df.filter(keep)
    if verbose and out.height < before:
        print(f"hard filters: dropped {before - out.height} of {before} (length/cysteine)")
    return out


def build(
    tables: list[Path | str],
    *,
    rank_by: str = "min_ipae_boltz2",
    ascending: bool = True,
    only_passing: bool = True,
    limit: int = MAX_DESIGNS,
    hard_filters: bool = True,
) -> pl.DataFrame:
    """Concatenate design tables, drop duplicates and failures, rank and cut to 20.

    Raises if ``rank_by`` is absent from the concatenated table: a silent no-op sort
    returns concatenation order, which looks exactly like a ranking and is not one.
    """
    frames = []
    for path in tables:
        df = pl.read_csv(path, infer_schema_length=10_000)
        frames.append(df.with_columns(source=pl.lit(str(path))))
    df = pl.concat(frames, how="diagonal_relaxed")

    if only_passing and "passes" in df.columns:
        passes = _as_bool(df, "passes")
        # Keep nulls: "undecided locally, abguard is the authority".
        df = df.filter(passes.is_null() | passes)

    subset = [c for c in ("kind", "format", "vh_sequence", "vl_sequence") if c in df.columns]
    if subset:
        df = df.unique(subset=subset, keep="first")
    elif "cdrh3" in df.columns:
        df = df.unique(subset=["cdrh3"], keep="first")

    if hard_filters:
        df = enforce_hard_filters(df)

    if rank_by not in df.columns:
        numeric = [c for c, t in df.schema.items() if t.is_numeric()]
        raise ValueError(
            f"rank_by={rank_by!r} is not a column of the concatenated table, so the "
            f"result would be in concatenation order rather than ranked. "
            f"Numeric columns available: {sorted(numeric)}"
        )
    df = df.sort(rank_by, descending=not ascending, nulls_last=True)
    return df.head(limit)


def main() -> None:  # pragma: no cover - thin CLI
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("command", choices=["build"])
    ap.add_argument("tables", nargs="+")
    ap.add_argument("-o", "--output", default="designs/submission.csv")
    ap.add_argument(
        "--kind",
        default=None,
        choices=["vhh", "scfv", "paired"],
        help="force one kind for every row (default: per-row, from the design's format)",
    )
    ap.add_argument("--rank-by", default="min_ipae_boltz2")
    ap.add_argument("--descending", action="store_true", help="higher is better")
    ap.add_argument("--limit", type=int, default=MAX_DESIGNS)
    ap.add_argument("--keep-failing", action="store_true", help="do not drop passes=False")
    args = ap.parse_args()

    ranked = build(
        args.tables,
        rank_by=args.rank_by,
        ascending=not args.descending,
        only_passing=not args.keep_failing,
        limit=args.limit,
    )
    sub = check_lengths(to_submission(ranked, kind=args.kind))
    bad = sub.filter(~pl.col("length_ok"))
    if bad.height:
        raise SystemExit(
            f"{bad.height} design(s) exceed {MAX_TOTAL_LENGTH} aa after ranking: "
            f"{bad['sequence_id'].to_list()}"
        )
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    sub.select(SUBMISSION_COLUMNS).write_csv(out)
    counts = sub.group_by("kind").len().sort("kind")
    print(f"{sub.height} designs -> {out}")
    print(f"by kind: {dict(zip(counts['kind'].to_list(), counts['len'].to_list()))}")
    print(f"Now run: bash scripts/validate.sh {out}")


if __name__ == "__main__":  # pragma: no cover
    main()
