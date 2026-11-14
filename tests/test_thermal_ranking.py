import polars as pl
import pytest

from abhack.ranking import add_thermal_scores


def test_thermal_join_and_control_gate(tmp_path):
    hot, control = tmp_path / "hot.csv", tmp_path / "control.csv"
    pl.DataFrame(
        {"design_id": ["D1", "D2"], "thermal_score": [0.9, 0.7], "q_worst": [0.9, 0.7]}
    ).write_csv(hot)
    pl.DataFrame({"design_id": ["D1", "D2"], "q_worst": [0.95, 0.92]}).write_csv(control)
    out = add_thermal_scores(pl.DataFrame({"design_id": ["D1", "D2"]}), hot, control, require=True)
    assert out["thermal_score"].to_list() == [0.9, 0.7]


def test_thermal_require_rejects_unscreened_rows(tmp_path):
    hot = tmp_path / "hot.csv"
    pl.DataFrame({"design_id": ["D1"], "thermal_score": [0.9]}).write_csv(hot)
    with pytest.raises(ValueError, match="lack thermal"):
        add_thermal_scores(pl.DataFrame({"design_id": ["D1", "D2"]}), hot, require=True)
