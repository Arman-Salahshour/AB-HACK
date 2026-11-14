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


def test_add_ph_scores_joins_on_sequence_pair_when_design_ids_repeat(tmp_path):
    # design_id repeats across runs; the same id with different sequences must get different scores
    scores = tmp_path / "ph.csv"
    pl.DataFrame({"design_id": ["x", "x"], "vh_sequence": ["AAA", "CCC"], "vl_sequence": ["L", "L"],
                  "ph_log10_sel": [0.4, -0.2]}).write_csv(scores)
    df = pl.DataFrame({"design_id": ["x", "x", "x"], "vh_sequence": ["AAA", "CCC", "GGG"],
                       "vl_sequence": ["L", "L", None]})
    out = add_ph_scores(df, scores)
    assert out["ph_log10_sel"].to_list() == [0.4, -0.2, None]


def test_ph_gate_drops_measured_negatives_and_keeps_unscored():
    from abhack.ranking import gate

    df = pl.DataFrame({"design_id": ["a", "b", "c"], "ph_log10_sel": [0.5, -0.1, None]})
    assert gate(df, max_ipae=None, ph_min=0.0)["design_id"].to_list() == ["a", "c"]
    assert gate(df, max_ipae=None, ph_min=None).height == 3


def test_ph_credit_needs_propka_and_md_agreement(tmp_path):
    scores = tmp_path / "ph.csv"
    pl.DataFrame({"design_id": ["agree", "propka_only", "md_only", "ctrl"],
                  "ph_log10_sel": [0.5, 0.5, 0.0, 0.0]}).write_csv(scores)
    md = tmp_path / "md.csv"
    rows = []
    for name, q65, q74 in (("agree", 0.97, 0.90), ("propka_only", 0.90, 0.90), ("md_only", 0.97, 0.90),
                           ("ctrl", 0.95, 0.95)):
        for rep, d in ((1, 0.0), (2, 0.01)):
            rows += [{"design": name, "ph": 6.5, "rep": rep, "q_interface": q65 + d},
                     {"design": name, "ph": 7.4, "rep": rep, "q_interface": q74 + d}]
    pl.DataFrame(rows).write_csv(md)
    df = pl.DataFrame({"design_id": ["agree", "propka_only", "md_only", "ctrl"]})
    out = add_ph_scores(df, scores, md).sort("design_id")
    credit = dict(zip(out["design_id"], out["ph_credit"], strict=True))
    assert credit["agree"] == pytest.approx(0.5)
    assert credit["propka_only"] == 0.0 and credit["md_only"] == 0.0 and credit["ctrl"] == 0.0
