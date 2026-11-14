"""Score and rank collected designs.

Combines what we can measure in silico into one ordering:

* ``min_ipae`` (lower better), ``iptm`` and ``binder_plddt`` from refolding —
  how confidently a folding model reproduces the designed interface;
* ``p_high`` from the FACS classifier and ``pkd_pred`` from the SPR regressor,
  if the models in ``models/`` have been trained;
* the novelty verdict, which is a gate, not a score.

Ranks are compared only within one folding model (see data/docs/18-vm.md).

    python -m abhack.evaluate designs/novel-s1.csv -o designs/novel-s1.ranked.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

from abhack.model import MODELS_DIR, AffinityRegressor, BinderClassifier

#: (column suffix, higher_is_better) for the refold metrics we rank on.
REFOLD_METRICS = {"min_ipae": False, "iptm": True, "binder_plddt": True}


def add_model_scores(df: pl.DataFrame) -> pl.DataFrame:
    """Append ``p_high`` and ``pkd_pred`` when the trained models exist."""
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
    mean, std = s.mean(), s.std()
    if std in (None, 0) or mean is None:
        return pl.Series([0.0] * len(s))
    z = (s - mean) / std
    return z if higher_is_better else -z


def composite_score(df: pl.DataFrame, model: str = "boltz2", weights: dict | None = None):
    """Mean z-score over the available refold metrics and model predictions."""
    weights = weights or {"min_ipae": 1.0, "iptm": 1.0, "binder_plddt": 0.5, "p_high": 1.0,
                          "pkd_pred": 1.0}
    parts, used = [], []
    for metric, higher in REFOLD_METRICS.items():
        col = f"{metric}_{model}"
        if col in df.columns and df[col].dtype.is_numeric():
            parts.append(_zscore(df[col], higher) * weights.get(metric, 1.0))
            used.append(col)
    for col, higher in (("p_high", True), ("pkd_pred", True)):
        if col in df.columns and df[col].dtype.is_numeric():
            parts.append(_zscore(df[col], higher) * weights.get(col, 1.0))
            used.append(col)
    if not parts:
        raise ValueError(f"no rankable columns found; have {df.columns}")
    total = parts[0]
    for p in parts[1:]:
        total = total + p
    print(f"composite over: {', '.join(used)}")
    return df.with_columns(composite=total / len(parts)).sort("composite", descending=True)


def main() -> None:  # pragma: no cover - thin CLI
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("table")
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("--model", default="boltz2", help="folding model suffix to rank on")
    ap.add_argument("--only-passing", action="store_true")
    args = ap.parse_args()

    df = pl.read_csv(args.table, infer_schema_length=10_000)
    if args.only_passing and "passes" in df.columns:
        df = df.filter(pl.col("passes"))
    df = composite_score(add_model_scores(df), model=args.model)

    show = [c for c in ("design_id", "cdrh3", "track", "passes", "composite",
                        f"min_ipae_{args.model}", f"iptm_{args.model}", "p_high", "kd_pred_nM")
            if c in df.columns]
    print(df.select(show).head(20))
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(args.output)
        print(f"{df.height} rows -> {args.output}")


if __name__ == "__main__":  # pragma: no cover
    main()
