"""Check the calibration maths on cases with a known answer.

The AUROC here is hand-rolled (no sklearn in this env), and a wrong AUROC would make a
useless predictor look usable — the exact error this module exists to prevent. So it is
pinned against cases whose value can be computed by hand.
"""

from __future__ import annotations

import numpy as np
import pytest

from abhack.calibrate import Calibration, calibrate, calibrate_all


def _cal(real, shuffled) -> Calibration:
    return Calibration("test", np.asarray(real, dtype=float), np.asarray(shuffled, dtype=float))


def test_auroc_perfect_separation():
    assert _cal([0.9, 0.8, 0.7], [0.1, 0.2, 0.3]).auroc == pytest.approx(1.0)


def test_auroc_perfectly_inverted():
    assert _cal([0.1, 0.2, 0.3], [0.9, 0.8, 0.7]).auroc == pytest.approx(0.0)


def test_auroc_all_ties_is_a_coin_flip():
    """Every score identical must be 0.5, not 0.0 or 1.0 — the tie correction."""
    assert _cal([0.5] * 4, [0.5] * 6).auroc == pytest.approx(0.5)


def test_auroc_known_fraction():
    """real=[2,4], shuffled=[1,3,5]: pairs won = (2>1) + (4>1) + (4>3) = 3 of 6."""
    assert _cal([2, 4], [1, 3, 5]).auroc == pytest.approx(3 / 6)


def test_auroc_half_credit_for_one_tie():
    """real=[2], shuffled=[1,2]: one win, one tie -> (1 + 0.5) / 2."""
    assert _cal([2], [1, 2]).auroc == pytest.approx(0.75)


def test_tpr_and_fpr_are_inclusive_at_the_threshold():
    c = _cal([0.5, 0.9], [0.1, 0.5])
    assert c.tpr_at(0.5) == pytest.approx(1.0)
    assert c.fpr_at(0.5) == pytest.approx(0.5)


def test_threshold_at_fpr_controls_the_negatives():
    rng = np.random.default_rng(0)
    c = _cal(rng.uniform(0.4, 1.0, 200), rng.uniform(0.0, 0.6, 2000))
    for target in (0.5, 0.1, 0.05, 0.01):
        t = c.threshold_at_fpr(target)
        assert c.fpr_at(t) <= target + 1e-9, f"FPR exceeded its budget at {target}"


def test_threshold_at_fpr_is_monotonic():
    rng = np.random.default_rng(1)
    c = _cal(rng.uniform(0.4, 1.0, 100), rng.uniform(0.0, 0.8, 1000))
    assert c.threshold_at_fpr(0.01) >= c.threshold_at_fpr(0.1) >= c.threshold_at_fpr(0.5)


def test_threshold_at_fpr_rejects_nonsense():
    c = _cal([0.9], [0.1])
    with pytest.raises(ValueError):
        c.threshold_at_fpr(0.0)
    with pytest.raises(ValueError):
        c.threshold_at_fpr(1.5)


def test_enrichment_of_an_uninformative_cutoff_is_one():
    """If real and shuffled survive a cutoff equally often, it carries no information."""
    c = _cal([0.9, 0.1], [0.9, 0.1])
    assert c.enrichment_at(0.5) == pytest.approx(1.0)


# --- the real data, as a regression on the conclusions we are acting on ---


def test_real_matrices_load_square_and_labelled():
    for name in ("af3", "boltz2", "chai1"):
        c = calibrate(name)
        assert c.n_real >= 89, f"{name}: expected ~91 cognate pairs, got {c.n_real}"
        assert c.n_shuffled > 7000, f"{name}: expected ~8100 shuffled pairs"
        assert np.isfinite(c.real).all() and np.isfinite(c.shuffled).all()


def test_boltz2_iptm_is_near_uninformative_at_the_conventional_cutoff():
    """The finding the ranking change rests on. If this breaks, re-read the data."""
    c = calibrate("boltz2")
    assert c.enrichment_at(0.8) < 1.3, "Boltz-2 ipTM@0.8 should carry almost no signal"
    assert c.fpr_at(0.8) > 0.5, "most shuffled pairs should clear 0.8 under Boltz-2"


def test_af3_separates_better_than_boltz2_which_beats_chai1():
    cals = calibrate_all()
    assert cals["af3"].auroc > cals["boltz2"].auroc > cals["chai1"].auroc
    assert cals["af3"].auroc > 0.7
