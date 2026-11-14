"""The one ranking path: gates, then a composite, then a quota merge.

This replaces the two parallel rankers the repo grew (`abhack.evaluate`, whose
`composite` never reached the submission, and `abhack.submission`, which
re-ranked on its own). `abhack.submission` now only converts a table ranked
here into the abguard schema; it never invents an order.

Ordering within a track, and why (see reports/calibration.md and CLAUDE.md):

* **gates first** — the local novelty verdict (``passes`` false -> out; null
  means *undecided, abguard decides* and stays), and
  :data:`IPAE_GATE` on the measured interface PAE. A design that fails a gate
  is never ranked, let alone refolded;
* **`min_ipae`** (lower better) is the primary evidence — the organisers'
  own criterion, and the only folding metric the calibration left standing;
  ``iptm_boltz2`` is deliberately absent: at the conventional cutoff it
  admitted 81.7 % of known non-binders (1.10x enrichment);
* ``binder_plddt`` (weight 0.5) — fold quality of the designed chain;
* ``p_high`` / ``pkd_pred`` — the models trained on this system (FACS 0.91
  accuracy, SPR R^2 0.56), preferred over any generic proxy;
* ``liab_cdrh3`` (weight 0.5, inverted) — developability ordering only.

AF2-IG columns from JapanFold, when present, join the composite for the
de novo track: ``af2ig_iptm`` (up) and ``af2ig_ipae`` (down).

Usage::

    python -m abhack.ranking rank designs/opt-s11.csv designs/opt-s12.csv \
        -o designs/opt-ranked.csv
    python -m abhack.ranking merge designs/opt-ranked.csv designs/denovo-ranked.csv \
        --opt-slots 12 --denovo-slots 8 -o designs/final-ranked.csv
    python -m abhack.submission build designs/final-ranked.csv \
        --rank-by submission_order --descending
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

from abhack.liabilities import add_liabilities
from abhack.model import MODELS_DIR, AffinityRegressor, BinderClassifier

#: interface-PAE pass mark (A) — organisers' vendored antibody criterion.
IPAE_GATE = 2.0

#: (column, higher_is_better, weight) for the composite. Order = tie-break order.
COMPONENTS: tuple[tuple[str, bool, float], ...] = (
    ("min_ipae_boltz2", False, 1.0),
    ("binder_plddt_boltz2", True, 0.5),
    ("p_high", True, 1.0),
    ("pkd_pred", True, 1.0),
    ("liab_cdrh3", False, 0.5),
    ("af2ig_ipae", False, 1.0),
    ("af2ig_iptm", True, 1.0),
)


def _as_bool(df: pl.DataFrame, column: str) -> pl.Series:
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


def add_model_scores(df: pl.DataFrame) -> pl.DataFrame:
    """Append ``p_high`` and ``pkd_pred`` when the trained models exist.

    Moved from the retired `abhack.evaluate`; identical behaviour. Only rows
    with a 10-residue ``cdrh3`` window (the FACS/SPR featurisation domain)
    are scored; VHH rows are left null and rank on their folding evidence.
    """
    if "cdrh3" not in df.columns:
        return df
    seqs = [s for s in df["cdrh3"].to_list()]
    mask = [bool(s) for s in seqs]
    usable = [s for s in seqs if s]
    if not usable:
        return df

    binder_path = MODELS_DIR / "binder.joblib"
    if binder_path.exists():
        p = BinderClassifier.load(binder_path).p_high(usable)
        df = df.with_columns(p_high=pl.Series(_scatter(p, mask)))
    affinity_path = MODELS_DIR / "affinity.joblib"
    if affinity_path.exists():
        model = AffinityRegressor.load(affinity_path)
        df = df.with_columns(
            pkd_pred=pl.Series(_scatter(model.predict(usable), mask)),
            kd_pred_nM=pl.Series(_scatter(model.predict_kd_nM(usable), mask)),
        )
    return df


def _scatter(values, mask: list[bool]) -> list[float | None]:
    it = iter(values)
    return [float(next(it)) if m else None for m in mask]


def _zscore(s: pl.Series, higher_is_better: bool) -> pl.Series:
    """Z-scored, direction-corrected, with unmeasured rows neutral at 0.

    A design without a refold yet (``min_ipae`` null) must still rank on its
    model scores; a null would propagate through the sum and bury it at the
    bottom regardless of everything else.
    """
    mean, std = s.mean(), s.std()
    if std in (None, 0) or mean is None:
        return pl.Series([0.0] * len(s))
    z = (s - mean) / std
    z = z.fill_null(0.0)
    return (z if higher_is_better else -z).fill_null(0.0)


def gate(df: pl.DataFrame, *, max_ipae: float | None = IPAE_GATE) -> pl.DataFrame:
    """Drop what the gates reject. Nulls are *undecided*, never dropped.

    * novelty: ``passes`` false -> out (the local check is a strict subset of
      abguard's, so anything it fails is dead); null stays;
    * interface PAE: a *measured* ``min_ipae_boltz2`` above ``max_ipae`` is out.
      Unmeasured designs stay: they may simply not have been refolded yet.
    """
    before = df.height
    if "passes" in df.columns:
        passes = _as_bool(df, "passes")
        df = df.filter(passes.is_null() | passes)
    if max_ipae is not None and "min_ipae_boltz2" in df.columns:
        ipae = pl.col("min_ipae_boltz2").cast(pl.Float64, strict=False)
        df = df.filter(ipae.is_null() | (ipae <= max_ipae))
    dropped = before - df.height
    if dropped:
        print(f"gates: dropped {dropped} of {before} (novelty/ipae)")
    return df


def composite_score(df: pl.DataFrame) -> pl.DataFrame:
    """Weighted mean z-score over the components present, then sort by it."""
    parts, used = [], []
    for col, higher, weight in COMPONENTS:
        if col in df.columns and df[col].dtype.is_numeric():
            parts.append(_zscore(df[col], higher) * weight)
            used.append(col)
    if not parts:
        raise ValueError(
            f"no rankable columns found; have {df.columns}. Rank on something "
            "real or admit the table is unranked — a silent concatenation order "
            "looks exactly like a ranking and is not one."
        )
    total = parts[0]
    for p in parts[1:]:
        total = total + p
    print(f"composite over: {', '.join(used)}")
    return df.with_columns(composite=(total / len(parts))).sort(
        "composite", descending=True, nulls_last=True
    )


def rank_table(tables: list[Path | str], *, max_ipae: float | None = IPAE_GATE) -> pl.DataFrame:
    """Collect design tables -> gate -> score -> one ranked table."""
    frames = [
        pl.read_csv(t, infer_schema_length=10_000).with_columns(source=pl.lit(str(t)))
        for t in tables
    ]
    df = pl.concat(frames, how="diagonal_relaxed")
    df = gate(df, max_ipae=max_ipae)
    df = add_liabilities(add_model_scores(df))
    return composite_score(df)


def merge(
    opt: Path | str,
    denovo: Path | str,
    *,
    opt_slots: int = 12,
    denovo_slots: int = 8,
) -> pl.DataFrame:
    """Concatenate the two tracks into one ranked list with ``submission_order``.

    Optimisation designs come first: the track has the stronger evidence
    (templated refolds, models trained on 525k measured variants). Slot
    counts are ceilings — a track never pads with designs it cannot rank,
    and a slot one track cannot fill spills to the other.
    """
    ranked = {}
    for path, track in ((opt, "optimisation"), (denovo, "de_novo")):
        df = pl.read_csv(path, infer_schema_length=10_000)
        if "composite" not in df.columns:
            raise ValueError(f"{path}: no composite column — rank it with `rank` first")
        ranked[track] = df.sort("composite", descending=True)

    takes = {
        "optimisation": ranked["optimisation"].head(opt_slots),
        "de_novo": ranked["de_novo"].head(denovo_slots),
    }
    for track, other, slots in (
        ("optimisation", "de_novo", opt_slots),
        ("de_novo", "optimisation", denovo_slots),
    ):
        spare = slots - takes[track].height
        pool = ranked[other]
        if spare > 0 and pool.height > len(takes[other]):
            takes[other] = pl.concat(
                [takes[other], pool.slice(len(takes[other]), spare)]
            )

    out = pl.concat(
        [
            takes["optimisation"].with_columns(track=pl.lit("optimisation")),
            takes["de_novo"].with_columns(track=pl.lit("de_novo")),
        ],
        how="diagonal_relaxed",
    )
    return out.with_row_index("submission_order", offset=1).with_columns(
        submission_order=pl.col("submission_order").cast(pl.Int64)
    )


def main() -> None:  # pragma: no cover - thin CLI
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="command", required=True)

    r = sub.add_parser("rank", help="gate + score + rank design tables")
    r.add_argument("tables", nargs="+")
    r.add_argument("-o", "--output", required=True)
    r.add_argument("--max-ipae", type=float, default=IPAE_GATE)
    r.add_argument("--no-ipae-gate", action="store_true")

    m = sub.add_parser("merge", help="quota-merge two ranked tables")
    m.add_argument("opt")
    m.add_argument("denovo")
    m.add_argument("-o", "--output", required=True)
    m.add_argument("--opt-slots", type=int, default=12)
    m.add_argument("--denovo-slots", type=int, default=8)

    args = ap.parse_args()
    if args.command == "rank":
        df = rank_table(
            args.tables,
            max_ipae=None if args.no_ipae_gate else args.max_ipae,
        )
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(args.output)
        show = [c for c in ("design_id", "track", "cdrh3", "composite", "passes",
                            "min_ipae_boltz2", "p_high", "kd_pred_nM", "liab_cdrh3")
                if c in df.columns]
        print(df.select(show).head(20))
        print(f"{df.height} ranked -> {args.output}")
    else:
        df = merge(args.opt, args.denovo, opt_slots=args.opt_slots, denovo_slots=args.denovo_slots)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(args.output)
        print(df.select([c for c in ("submission_order", "design_id", "track", "composite") if c in df.columns]))
        print(f"{df.height} merged -> {args.output}")


if __name__ == "__main__":  # pragma: no cover
    main()
