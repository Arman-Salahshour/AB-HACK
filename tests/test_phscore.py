import math

import polars as pl
import pytest

from abhack.constants import TRASTUZUMAB_VH, TRASTUZUMAB_VL
from abhack.phscore import PH_HIGH, PH_LOW, mutations, selectivity
from abhack.ranking import add_ph_scores


def test_unchanged_pka_gives_zero_selectivity():
    groups = {("A", 91, "HIS"): (5.0, 6.5), ("B", 35, "HIS"): (6.8, 6.5)}
    assert selectivity(groups, dict(groups))["log10_sel"] == pytest.approx(0.0)


def test_raised_bound_pka_prefers_low_ph():
    free = {("A", 30, "HIS"): (6.0, 6.5)}
    bound = {("A", 30, "HIS"): (7.5, 6.5)}
    r = selectivity(bound, free, PH_LOW, PH_HIGH)
    expect = (math.log10(1 + 10 ** (7.5 - 6.5)) - math.log10(1 + 10 ** (6.0 - 6.5))) - (
        math.log10(1 + 10 ** (7.5 - 7.4)) - math.log10(1 + 10 ** (6.0 - 7.4))
    )
    assert r["log10_sel"] == pytest.approx(expect)
    assert r["log10_sel"] > 0


def test_lowered_bound_pka_prefers_high_ph():
    free = {("A", 30, "HIS"): (6.5, 6.5)}
    bound = {("A", 30, "HIS"): (4.0, 6.5)}
    assert selectivity(bound, free)["log10_sel"] < 0


def test_antigen_only_groups_are_ignored():
    bound = {("C", 497, "HIS"): (7.0, 6.5)}
    assert selectivity(bound, {})["log10_sel"] == 0.0


def test_mutations_reports_differences_and_rejects_indels():
    vh = TRASTUZUMAB_VH[:32] + "H" + TRASTUZUMAB_VH[33:]
    assert mutations(vh, TRASTUZUMAB_VL) == [("B", 33, "Y", "H")]
    with pytest.raises(ValueError):
        mutations(TRASTUZUMAB_VH + "A", TRASTUZUMAB_VL)


def test_add_ph_scores_joins_by_design_id(tmp_path):
    scores = tmp_path / "ph.csv"
    pl.DataFrame({"design_id": ["a", "b"], "ph_log10_sel": [0.3, -0.1]}).write_csv(scores)
    out = add_ph_scores(pl.DataFrame({"design_id": ["a", "b", "c"]}), scores)
    assert out["ph_log10_sel"].to_list() == [0.3, -0.1, None]
    dup = tmp_path / "dup.csv"
    pl.DataFrame({"design_id": ["a", "a"], "ph_log10_sel": [0.1, 0.2]}).write_csv(dup)
    with pytest.raises(ValueError):
        add_ph_scores(pl.DataFrame({"design_id": ["a"]}), dup)
