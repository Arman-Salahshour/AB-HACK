"""Sequence models over the trastuzumab CDRH3 window.

Two jobs, two labels:

* :class:`BinderClassifier` — FACS classes (high/medium/low), 524k rows. Says
  *whether* a variant binds.
* :class:`AffinityRegressor` — SPR pKD, 190 rows. Says *how tightly*, and is the
  one that matters for the affinity axis of the judging.

Both take raw CDRH3 strings and featurise internally, so a trained model can be
applied straight to a design table.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, Ridge

from abhack.features import featurise

MODELS_DIR = Path(__file__).resolve().parents[2] / "models"


@dataclass
class _Base:
    kind: str = "gbt"
    estimator: object | None = field(default=None, repr=False)
    classes_: list | None = None

    def _features(self, sequences: list[str]) -> np.ndarray:
        return featurise(list(sequences))

    def save(self, path: Path | str) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(self, path)
        return path

    @staticmethod
    def load(path: Path | str):
        return joblib.load(Path(path))


@dataclass
class BinderClassifier(_Base):
    """FACS class predictor. ``kind`` is ``gbt`` or ``logreg``."""

    def fit(self, sequences: list[str], labels: list[str]) -> BinderClassifier:
        X = self._features(sequences)
        if self.kind == "logreg":
            est = LogisticRegression(max_iter=2000, n_jobs=-1, C=1.0)
        else:
            est = HistGradientBoostingClassifier(
                max_iter=300, learning_rate=0.1, early_stopping=True, random_state=0
            )
        est.fit(X, labels)
        self.estimator, self.classes_ = est, list(est.classes_)
        return self

    def predict(self, sequences: list[str]) -> np.ndarray:
        return self.estimator.predict(self._features(sequences))

    def predict_proba(self, sequences: list[str]) -> np.ndarray:
        return self.estimator.predict_proba(self._features(sequences))

    def p_high(self, sequences: list[str]) -> np.ndarray:
        """Probability of the ``high`` FACS class — the useful scalar for ranking."""
        proba = self.predict_proba(sequences)
        return proba[:, self.classes_.index("high")]


@dataclass
class AffinityRegressor(_Base):
    """pKD = -log10(KD [M]) regressor. ``kind`` is ``gbt`` or ``ridge``."""

    def fit(self, sequences: list[str], pkd: np.ndarray) -> AffinityRegressor:
        X = self._features(sequences)
        if self.kind == "ridge":
            est = Ridge(alpha=1.0)
        else:
            est = HistGradientBoostingRegressor(
                max_iter=400, learning_rate=0.05, early_stopping=False, random_state=0
            )
        est.fit(X, np.asarray(pkd, dtype=float))
        self.estimator = est
        return self

    def predict(self, sequences: list[str]) -> np.ndarray:
        return self.estimator.predict(self._features(sequences))

    def predict_kd_nM(self, sequences: list[str]) -> np.ndarray:
        return 10 ** (-self.predict(sequences)) * 1e9
