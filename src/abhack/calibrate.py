"""Calibrate a confidence filter against known binders and known non-binders.

Why this module exists
----------------------
`abhack.evaluate` ranks designs on **within-batch z-scores**. A z-score says how a
design compares with the others in its own batch, so a batch of uniformly bad designs
still produces a confident-looking top 20. Nothing in the pipeline had an *absolute*
boundary, and nothing had ever been checked against data where the answer is known.

The Smorodina benchmark supplies exactly that check. `iptm_confidence_scores.zip` holds
one 91x91 matrix per predictor, indexed and columned by PDB id:

* the **diagonal** is a real, cognate nanobody-antigen pair;
* every **off-diagonal** cell is a shuffled pairing — the same nanobody against somebody
  else's antigen, which is a non-binder by construction.

So each matrix is ~91 positives and ~8,100 negatives, scored by the same model we would
use to rank designs. That is enough to ask the only question that matters: *at a cutoff
that keeps most real binders, how many non-binders come with them?*

The headline result (see :func:`report`) is that **Boltz-2 ipTM barely separates the two
populations at all**, which matters because `scripts/refold.sh` folds with Boltz-2 and
`submission.build` used to rank on its output.

Two honesty constraints, both load-bearing
------------------------------------------
1. The negatives are **shuffled pairings, scored computationally** — not measured
   non-binders. The organisers' own README is explicit that the SPR "no reportable fit"
   rows must *not* be relabelled as measured true negatives, and the 454 SPR negatives
   referred to in the task description do not exist in the released data.
2. Shuffled pairs are an *easier* negative than a real design: a design has been
   optimised against this target, so it is engineered to look plausible. The false
   positive rate measured here is therefore a **lower bound** on the rate we would see
   on our own candidates. A cutoff that cannot separate shuffled pairs certainly cannot
   separate designs; one that can is necessary, not sufficient.

Usage::

    python -m abhack.calibrate report
    python -m abhack.calibrate report --fpr 0.01 -o reports/calibration.md
"""

from __future__ import annotations

import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

#: Where the released matrices live, relative to the repo root.
IPTM_ZIP = Path("data/raw/zenodo-18390239-ab-ag-complexes/iptm_confidence_scores.zip")

#: Members of that zip, by the name we use for the predictor.
IPTM_MEMBERS = {
    "af3": "iptm_confidence_scores/af3_matrix_clean.csv",
    "boltz2": "iptm_confidence_scores/boltz2_matrix_clean.csv",
    "chai1": "iptm_confidence_scores/chai1_matrix_clean.csv",
}

#: Which predictor each local folding model's ipTM should be read against. `refold.sh`
#: defaults to Boltz-2, whose calibration is the worst of the three.
MODEL_TO_PREDICTOR = {
    "boltz2": "boltz2",
    "boltz1": "boltz2",
    "chai1": "chai1",
}


@dataclass(frozen=True)
class Calibration:
    """What a predictor's ipTM is worth, measured on known pairs."""

    predictor: str
    real: np.ndarray
    shuffled: np.ndarray

    @property
    def n_real(self) -> int:
        return int(self.real.size)

    @property
    def n_shuffled(self) -> int:
        return int(self.shuffled.size)

    @property
    def auroc(self) -> float:
        """P(a real pair outscores a shuffled one), ties counted as a half.

        0.5 is coin-flip. Computed by rank-sum rather than by thresholding so it does
        not depend on a grid.
        """
        both = np.concatenate([self.real, self.shuffled])
        order = both.argsort(kind="mergesort")
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, both.size + 1, dtype=float)
        # average ranks within tie groups, so ties contribute 0.5
        _, inverse, counts = np.unique(both, return_inverse=True, return_counts=True)
        sums = np.zeros(counts.size)
        np.add.at(sums, inverse, ranks)
        ranks = (sums / counts)[inverse]
        r_real = ranks[: self.real.size].sum()
        n_r, n_s = self.real.size, self.shuffled.size
        return float((r_real - n_r * (n_r + 1) / 2) / (n_r * n_s))

    def tpr_at(self, threshold: float) -> float:
        """Fraction of real pairs kept by ``score >= threshold`` (sensitivity)."""
        return float((self.real >= threshold).mean())

    def fpr_at(self, threshold: float) -> float:
        """Fraction of shuffled pairs that also survive (the false-positive rate)."""
        return float((self.shuffled >= threshold).mean())

    def threshold_at_fpr(self, target_fpr: float) -> float:
        """Lowest cutoff whose false-positive rate is <= ``target_fpr``.

        Uses the empirical quantile of the negatives, so the cutoff is defined by the
        non-binders rather than by the positives — the direction that matters when the
        cost of a false positive is a wasted submission slot.
        """
        if not 0.0 < target_fpr <= 1.0:
            raise ValueError(f"target_fpr must be in (0, 1], got {target_fpr}")
        return float(np.quantile(self.shuffled, 1.0 - target_fpr))

    def enrichment_at(self, threshold: float) -> float:
        """TPR / FPR: how much more often a real pair survives than a shuffled one.

        1.0 means the cutoff carries no information. This is the number to quote,
        because it is invariant to how many candidates are screened.
        """
        fpr = self.fpr_at(threshold)
        if fpr == 0.0:
            return float("inf")
        return self.tpr_at(threshold) / fpr


def load_matrix(predictor: str, *, zip_path: Path | str = IPTM_ZIP) -> tuple[np.ndarray, list[str]]:
    """The raw ipTM matrix and its PDB ids, read straight out of the zip."""
    member = IPTM_MEMBERS[predictor]
    with zipfile.ZipFile(zip_path) as z:
        text = z.read(member).decode()
    lines = [ln for ln in text.splitlines() if ln.strip()]
    header = lines[0].split(",")[1:]
    ids, rows = [], []
    for line in lines[1:]:
        parts = line.split(",")
        ids.append(parts[0])
        rows.append([np.nan if p == "" or p.upper() == "NA" else float(p) for p in parts[1:]])
    matrix = np.asarray(rows, dtype=float)
    if [i.strip() for i in ids] != [h.strip() for h in header]:
        raise ValueError(f"{predictor}: row and column ids differ, cannot take a diagonal")
    return matrix, ids


def calibrate(predictor: str, *, zip_path: Path | str = IPTM_ZIP) -> Calibration:
    """Split a predictor's matrix into its real diagonal and shuffled off-diagonal."""
    matrix, _ = load_matrix(predictor, zip_path=zip_path)
    real = np.diag(matrix)
    off = ~np.eye(matrix.shape[0], dtype=bool)
    shuffled = matrix[off]
    return Calibration(
        predictor=predictor,
        real=real[np.isfinite(real)],
        shuffled=shuffled[np.isfinite(shuffled)],
    )


def calibrate_all(*, zip_path: Path | str = IPTM_ZIP) -> dict[str, Calibration]:
    return {name: calibrate(name, zip_path=zip_path) for name in IPTM_MEMBERS}


def report(target_fpr: float = 0.05, *, zip_path: Path | str = IPTM_ZIP) -> str:
    """A markdown summary of what each predictor's ipTM is actually worth."""
    cals = calibrate_all(zip_path=zip_path)
    out: list[str] = []
    out.append("# ipTM calibration against known cognate / shuffled pairs\n")
    out.append(
        "Source: Smorodina et al., `iptm_confidence_scores.zip`. Diagonal = real cognate "
        "pair, off-diagonal = shuffled pairing (non-binder by construction).\n"
    )
    out.append(
        "Negatives are **computational** shuffled pairings, not measured non-binders, and "
        "are an easier negative than an optimised design — so every FPR below is a "
        "*lower bound* on what we would see on our own candidates.\n"
    )

    out.append("## Separation\n")
    out.append("| predictor | n real | n shuffled | median real | median shuffled | AUROC |")
    out.append("|---|---|---|---|---|---|")
    for name, c in cals.items():
        out.append(
            f"| {name} | {c.n_real} | {c.n_shuffled} | {np.median(c.real):.3f} | "
            f"{np.median(c.shuffled):.3f} | **{c.auroc:.3f}** |"
        )

    out.append(f"\n## Cutoff at a {target_fpr:.0%} false-positive rate\n")
    out.append("| predictor | cutoff | real kept (TPR) | shuffled kept (FPR) | enrichment |")
    out.append("|---|---|---|---|---|")
    for name, c in cals.items():
        t = c.threshold_at_fpr(target_fpr)
        out.append(
            f"| {name} | ipTM >= {t:.3f} | {c.tpr_at(t):.1%} | {c.fpr_at(t):.1%} | "
            f"{c.enrichment_at(t):.2f}x |"
        )

    out.append("\n## What a 0.8 cutoff buys, the number people actually use\n")
    out.append("| predictor | real kept | shuffled kept | enrichment |")
    out.append("|---|---|---|---|")
    for name, c in cals.items():
        out.append(
            f"| {name} | {c.tpr_at(0.8):.1%} | {c.fpr_at(0.8):.1%} | {c.enrichment_at(0.8):.2f}x |"
        )

    b = cals["boltz2"]
    a = cals["af3"]
    out.append("\n## Conclusion\n")
    out.append(
        f"- **Boltz-2 ipTM must not be used to rank designs.** AUROC {b.auroc:.3f}; "
        f"{b.fpr_at(0.8):.0%} of known non-binders score >= 0.8, against "
        f"{b.tpr_at(0.8):.0%} of real pairs — an enrichment of {b.enrichment_at(0.8):.2f}x. "
        f"`scripts/refold.sh` folds with Boltz-2, so its `iptm` column carries almost no "
        f"evidence of binding."
    )
    out.append(
        f"- AF3 ipTM does separate (AUROC {a.auroc:.3f}), but only "
        f"{a.tpr_at(0.8):.0%} of real pairs reach 0.8, so a 0.8 cutoff on AF3 would also "
        f"discard most true binders. AF3 is not available locally."
    )
    out.append(
        "- Therefore rank on **interface PAE** against the organisers' vendored antibody "
        "criterion (`min_ipae <= 2.0` A, `data/docs/06-evaluate.md`) rather than on ipTM, "
        "and treat an orthogonal predictor (AF2-IG) as the confirmatory check. "
        "Note this benchmark contains ipTM only, so the PAE threshold is vendored, "
        "not calibrated here."
    )
    return "\n".join(out) + "\n"


def main() -> None:  # pragma: no cover - thin CLI
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("command", choices=["report"])
    ap.add_argument("--fpr", type=float, default=0.05, help="target false-positive rate")
    ap.add_argument("--zip", dest="zip_path", default=str(IPTM_ZIP))
    ap.add_argument("-o", "--output", default=None)
    args = ap.parse_args()

    text = report(args.fpr, zip_path=args.zip_path)
    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
        print(f"-> {out}")
    print(text)


if __name__ == "__main__":  # pragma: no cover
    main()
