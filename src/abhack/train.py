"""Train the CDRH3 models.

    python -m abhack.train binder      # FACS class classifier  -> models/binder.joblib
    python -m abhack.train affinity    # SPR pKD regressor      -> models/affinity.joblib
    python -m abhack.train all

The SPR set is tiny (190 KDs), so the affinity model is reported with repeated
K-fold cross-validation rather than a single hold-out.
"""

from __future__ import annotations

import argparse

import numpy as np
from sklearn.model_selection import RepeatedKFold, cross_val_score, train_test_split

from abhack.data import load_facs, load_spr_quantified
from abhack.model import MODELS_DIR, AffinityRegressor, BinderClassifier


def train_binder(kind: str = "gbt", sample: int | None = None, seed: int = 0) -> BinderClassifier:
    df = load_facs()
    if sample:
        df = df.sample(n=min(sample, df.height), seed=seed)
    seqs = df["cdrh3"].to_list()
    labels = df["binding_class"].to_list()
    s_tr, s_te, y_tr, y_te = train_test_split(
        seqs, labels, test_size=0.1, random_state=seed, stratify=labels
    )
    model = BinderClassifier(kind=kind).fit(s_tr, y_tr)
    acc = (model.predict(s_te) == np.array(y_te)).mean()
    print(f"binder[{kind}] n={len(s_tr):,}  held-out accuracy {acc:.3f}  (3 classes, chance 0.33)")
    path = model.save(MODELS_DIR / "binder.joblib")
    print(f"  -> {path}")
    return model


def train_affinity(kind: str = "gbt", seed: int = 0) -> AffinityRegressor:
    df = load_spr_quantified()
    seqs, pkd = df["cdrh3"].to_list(), df["pkd"].to_numpy()
    model = AffinityRegressor(kind=kind)
    X = model._features(seqs)
    cv = RepeatedKFold(n_splits=5, n_repeats=3, random_state=seed)
    est = model.__class__(kind=kind).fit(seqs, pkd).estimator
    scores = cross_val_score(est, X, pkd, cv=cv, scoring="r2")
    print(
        f"affinity[{kind}] n={len(seqs)}  CV R2 {scores.mean():.3f} +/- {scores.std():.3f}"
        f"  (pKD range {pkd.min():.2f}-{pkd.max():.2f})"
    )
    model.fit(seqs, pkd)
    path = model.save(MODELS_DIR / "affinity.joblib")
    print(f"  -> {path}")
    return model


def main() -> None:  # pragma: no cover - thin CLI
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("what", choices=["binder", "affinity", "all"])
    ap.add_argument("--kind", default="gbt")
    ap.add_argument("--sample", type=int, default=None, help="subsample FACS rows")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.what in ("binder", "all"):
        train_binder(args.kind, args.sample, args.seed)
    if args.what in ("affinity", "all"):
        train_affinity("ridge" if args.kind == "logreg" else args.kind, args.seed)


if __name__ == "__main__":  # pragma: no cover
    main()
