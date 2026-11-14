"""Tests for the single ranking path: liabilities, gates, composite, merge.

The invariants these protect, in one place:

* a table with nothing rankable must **raise**, never fall back to
  concatenation order (that looks exactly like a ranking and is not one);
* the novelty gate drops only ``passes=False`` — null means *undecided,
  abguard decides* and must survive, or the de novo track disappears;
* the ipae gate drops only designs it **measured** as bad; unmeasured
  designs stay (they may just not be refolded yet);
* merge gives a track only what it can rank and spills spare slots.
"""

import polars as pl
import pytest

from abhack.liabilities import count_motifs, liability_penalty, penalty
from abhack.ranking import IPAE_GATE, composite_score, gate, merge

# --- liabilities ---------------------------------------------------------------

def test_count_motifs_sees_the_calibrated_tags():
    # CDRH3 probe calibrated against `abguard liability-report` 2026-11-14.
    c = count_motifs("WNGWDSMDP")
    assert c["deamd_high"] == 1  # NG
    assert c["isom"] == 1        # DS
    assert c["frag_high"] == 1   # DP
    assert c["met_ox"] == 1       # M
    assert c["trp_ox"] == 2       # W


def test_count_motifs_clean_loop():
    assert all(n == 0 for n in count_motifs("AYSYFYTYV").values())


def test_ngly_sequon_but_not_proline_blocked():
    assert count_motifs("NNTS")["ngly"] == 1     # N-T-S sequon
    assert count_motifs("NPTS")["ngly"] == 0    # proline blocks
    assert count_motifs("NGTS")["ngly"] == 1


def test_penalty_weights_and_detail():
    p, detail = liability_penalty("NGDP")
    # NG -> deamd_high 3.0, DP -> frag_high 3.0; no DG/DS/DT, so no isom.
    assert p == pytest.approx(6.0)
    assert "NG" in detail and "DP" in detail
    assert liability_penalty(None) == (0.0, "")


def test_penalty_zero_for_clean():
    assert penalty(count_motifs("AYSYFYTYV")) == 0.0


# --- gates ---------------------------------------------------------------------

def _table(passes, ipae):
    return pl.DataFrame(
        {
            "design_id": [f"d{i}" for i in range(len(passes))],
            "passes": passes,
            "min_ipae_boltz2": ipae,
        }
    )


def test_gate_keeps_null_and_true_drops_false():
    out = gate(_table([True, None, False], [0.3, None, 0.2]))
    assert out["design_id"].to_list() == ["d0", "d1"]


def test_gate_ipae_drops_only_measured_failures():
    out = gate(_table([True] * 3, [0.3, IPAE_GATE, 2.5]))
    assert out["design_id"].to_list() == ["d0", "d1"]


def test_gate_ipae_nulls_survive():
    out = gate(_table([True, True], [None, 0.4]))
    assert out.height == 2


# --- composite -----------------------------------------------------------------

def test_composite_raises_when_nothing_rankable():
    with pytest.raises(ValueError, match="no rankable"):
        composite_score(pl.DataFrame({"design_id": ["d0"]}))


def test_composite_prefers_low_ipae_and_high_plddt():
    df = pl.DataFrame(
        {
            "design_id": ["good", "bad"],
            "min_ipae_boltz2": [0.3, 1.9],
            "binder_plddt_boltz2": [97.0, 90.0],
        }
    )
    out = composite_score(df)
    assert out["design_id"][0] == "good"


def test_composite_penalises_liabilities():
    df = pl.DataFrame(
        {
            "design_id": ["clean", "dirty"],
            "min_ipae_boltz2": [0.4, 0.4],
            "liab_cdrh3": [0.0, 6.0],
        }
    )
    out = composite_score(df)
    assert out["design_id"][0] == "clean"


# --- merge ---------------------------------------------------------------------

def _ranked(n, track_prefix):
    return pl.DataFrame(
        {
            "design_id": [f"{track_prefix}{i}" for i in range(n)],
            "composite": [float(n - i) for i in range(n)],
        }
    )


def test_merge_orders_opt_first_then_denovo(tmp_path):
    opt = tmp_path / "opt.csv"
    dn = tmp_path / "dn.csv"
    _ranked(4, "o").write_csv(opt)
    _ranked(3, "n").write_csv(dn)
    out = merge(opt, dn, opt_slots=2, denovo_slots=2)
    assert out["track"].to_list() == ["optimisation"] * 2 + ["de_novo"] * 2
    assert out["submission_order"].to_list() == [1, 2, 3, 4]


def test_merge_spills_unused_slots(tmp_path):
    opt = tmp_path / "opt.csv"
    dn = tmp_path / "dn.csv"
    _ranked(5, "o").write_csv(opt)
    _ranked(1, "n").write_csv(dn)
    out = merge(opt, dn, opt_slots=4, denovo_slots=4)
    # de novo filled 1 of its 4 slots; opt takes the spare up to its pool:
    # 4 + 1 spilled = 5 optimisation, plus the 1 de novo, 6 of the 8 slots.
    assert out.height == 6
    assert out["track"].to_list().count("optimisation") == 5
    assert out["design_id"].to_list().count("o4") == 1  # the spilled design


def test_merge_never_pads_beyond_evidence(tmp_path):
    opt = tmp_path / "opt.csv"
    dn = tmp_path / "dn.csv"
    _ranked(1, "o").write_csv(opt)
    _ranked(1, "n").write_csv(dn)
    out = merge(opt, dn, opt_slots=8, denovo_slots=8)
    assert out.height == 2
