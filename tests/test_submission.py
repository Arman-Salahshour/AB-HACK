"""The submission path's failure modes, each of which silently lost designs.

Every test here corresponds to a bug that produced a *plausible-looking* submission
rather than an error, which is why they are worth pinning.
"""

from __future__ import annotations

import polars as pl
import pytest

from abhack.submission import build, enforce_hard_filters, resolve_kind, to_submission

# A VHH framework long enough to look like a domain; content is irrelevant to these
# tests, only length, cysteine parity and uniqueness matter.
VHH_A = "EVQLVESGGGLVQPGGSLRLSCAASGSNITRYADSVKGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAARLAGLYSTLESDYDYWGQGTQVTVS"
VHH_B = "EVQLVESGGGLVQPGGSLRLSCAASGLTFTRYADSVKGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAARTYSNWWESDYDYWGQGTQVTVS"
TZ_VH = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
TZ_VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"


def _write(tmp_path, rows: list[dict], name: str = "t.csv") -> str:
    path = tmp_path / name
    pl.DataFrame(rows).write_csv(path)
    return str(path)


def test_null_passes_survives_to_submission(tmp_path):
    """A de novo VHH cannot be graded locally, so `passes` is null — not False.

    Filtering `pl.col("passes")` drops nulls, which deleted the entire de novo track.
    """
    table = _write(
        tmp_path,
        [
            {"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None, "passes": None, "s": 1.0},
            {"format": "vhh", "vh_sequence": VHH_B, "vl_sequence": None, "passes": None, "s": 2.0},
        ],
    )
    out = build([table], rank_by="s", ascending=False)
    assert out.height == 2, "null `passes` must be kept for abguard to decide"


def test_explicit_failure_is_still_dropped(tmp_path):
    table = _write(
        tmp_path,
        [
            {"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None, "passes": False, "s": 1.0},
            {"format": "vhh", "vh_sequence": VHH_B, "vl_sequence": None, "passes": True, "s": 2.0},
        ],
    )
    out = build([table], rank_by="s", ascending=False)
    assert out.height == 1
    assert out["vh_sequence"][0] == VHH_B


def test_distinct_vhhs_do_not_collapse(tmp_path):
    """Deduplicating on `cdrh3` treated every null as the same value."""
    table = _write(
        tmp_path,
        [
            {"format": "vhh", "cdrh3": None, "vh_sequence": VHH_A, "vl_sequence": None, "s": 1.0},
            {"format": "vhh", "cdrh3": None, "vh_sequence": VHH_B, "vl_sequence": None, "s": 2.0},
        ],
    )
    out = build([table], rank_by="s", ascending=False)
    assert out.height == 2, "two different VHHs must not dedupe to one"


def test_identical_sequences_do_dedupe(tmp_path):
    table = _write(
        tmp_path,
        [
            {"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None, "s": 1.0},
            {"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None, "s": 2.0},
        ],
    )
    out = build([table], rank_by="s", ascending=False)
    assert out.height == 1


def test_kind_is_per_row(tmp_path):
    """A mixed submission stamped one global `kind` and failed abguard's schema."""
    df = pl.DataFrame(
        [
            {"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None},
            {"format": "paired", "vh_sequence": TZ_VH, "vl_sequence": TZ_VL},
        ]
    )
    assert resolve_kind(df).to_list() == ["vhh", "paired"]
    sub = to_submission(df)
    assert sub["kind"].to_list() == ["vhh", "paired"]


def test_kind_can_still_be_forced(tmp_path):
    df = pl.DataFrame([{"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None}])
    assert to_submission(df, kind="scfv")["kind"].to_list() == ["scfv"]


def test_over_length_design_is_rejected():
    """The 250 aa cap was warn-only, so over-length rows shipped."""
    long_vh = "A" * 200
    long_vl = "G" * 60  # 260 total
    df = pl.DataFrame(
        [
            {"format": "paired", "vh_sequence": long_vh, "vl_sequence": long_vl},
            {"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None},
        ]
    )
    out = enforce_hard_filters(df, verbose=False)
    assert out.height == 1
    assert out["vh_sequence"][0] == VHH_A


def test_unpaired_cysteine_is_rejected():
    odd_cys = VHH_A.replace("C", "", 1)  # VHH_A has 2 C; removing one leaves 1
    assert odd_cys.count("C") % 2 == 1
    df = pl.DataFrame(
        [
            {"format": "vhh", "vh_sequence": odd_cys, "vl_sequence": None},
            {"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None},
        ]
    )
    out = enforce_hard_filters(df, verbose=False)
    assert out.height == 1
    assert out["vh_sequence"][0] == VHH_A


def test_unpaired_cysteine_in_light_chain_is_rejected():
    """The old check looked at VH only and never at VL."""
    odd_vl = TZ_VL.replace("C", "", 1)
    df = pl.DataFrame([{"format": "paired", "vh_sequence": TZ_VH, "vl_sequence": odd_vl}])
    assert enforce_hard_filters(df, verbose=False).height == 0


def test_missing_rank_by_raises(tmp_path):
    """A silent no-op sort returns concatenation order and looks like a ranking."""
    table = _write(
        tmp_path,
        [{"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None, "s": 1.0}],
    )
    with pytest.raises(ValueError, match="rank_by"):
        build([table], rank_by="min_ipae_boltz2")


def test_string_passes_column_is_parsed(tmp_path):
    """An all-empty `passes` round-trips through CSV as Null/String, not Boolean."""
    table = _write(
        tmp_path,
        [
            {"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None, "passes": "true", "s": 1.0},
            {"format": "vhh", "vh_sequence": VHH_B, "vl_sequence": None, "passes": "false", "s": 2.0},
        ],
    )
    out = build([table], rank_by="s", ascending=False)
    assert out.height == 1
    assert out["vh_sequence"][0] == VHH_A


def test_ranking_order_is_respected(tmp_path):
    table = _write(
        tmp_path,
        [
            {"format": "vhh", "vh_sequence": VHH_A, "vl_sequence": None, "s": 1.0},
            {"format": "vhh", "vh_sequence": VHH_B, "vl_sequence": None, "s": 9.0},
        ],
    )
    asc = build([table], rank_by="s", ascending=True)
    assert asc["s"].to_list() == [1.0, 9.0]
    desc = build([table], rank_by="s", ascending=False)
    assert desc["s"].to_list() == [9.0, 1.0]
