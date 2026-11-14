"""Tests for the thermostability ticket designs.

The protected properties, all of which come from the mAbs 2025 evidence:

* the mutations only ever apply to the **trastuzumab** framework (wild-type
  check) — any other framework must raise, not mutate;
* the CDRH3 window is untouched (novelty gate, paratope, model scores);
* cysteines are unchanged (the unpaired-cysteine hard filter).
"""

import polars as pl
import pytest

from abhack.constants import CDRH3_WINDOW_START, TRASTUZUMAB_VH, TRASTUZUMAB_VL
from abhack.tickets import (
    FRAMEWORK_MUTATIONS,
    TICKET_TAG,
    apply_framework_mutations,
    build_tickets,
)


def test_mutations_hit_the_paper_positions():
    out = apply_framework_mutations(TRASTUZUMAB_VH)
    assert TRASTUZUMAB_VH[84] == "S" and out[84] == "N"
    assert TRASTUZUMAB_VH[86] == "R" and out[86] == "T"
    assert out[80:94] == "LQMNNLTAEDTAVY"


def test_mutations_leave_cdrh3_and_cysteines_alone():
    out = apply_framework_mutations(TRASTUZUMAB_VH)
    assert out[CDRH3_WINDOW_START:] == TRASTUZUMAB_VH[CDRH3_WINDOW_START:]
    assert out.count("C") == TRASTUZUMAB_VH.count("C")


def test_sequon_created_is_the_paper_one():
    out = apply_framework_mutations(TRASTUZUMAB_VH)
    assert "NLT" in out  # the new N-glycosylation site, a known trade-off


def test_foreign_framework_raises():
    other = list(TRASTUZUMAB_VH)
    other[84] = "A"  # not S at 85
    with pytest.raises(ValueError, match="position 85"):
        apply_framework_mutations("".join(other))


def test_build_tickets_twins_only_the_vh():
    ranked = pl.DataFrame(
        {
            "design_id": ["d1", "d2"],
            "cdrh3": ["AYSYFYTYVA", "WSSAVGYAMD"],
            "vh_sequence": [TRASTUZUMAB_VH] * 2,
            "vl_sequence": [TRASTUZUMAB_VL] * 2,
            "format": ["paired", "paired"],
            "composite": [1.0, 0.5],
        }
    )
    out = build_tickets(ranked, n=4)
    assert out.height == 2
    assert out["design_id"].to_list() == ["d1-ft", "d2-ft"]
    assert (out["framework_ticket"]).all()
    for r in out.to_dicts():
        assert r["cdrh3"] in ("AYSYFYTYVA", "WSSAVGYAMD")  # parent's loop kept
        assert r["vh_sequence"] != TRASTUZUMAB_VH  # framework mutated


def test_build_tickets_skips_vhh_rows():
    ranked = pl.DataFrame(
        {
            "design_id": ["v1", "d1"],
            "cdrh3": [None, "AYSYFYTYVA"],
            "vh_sequence": [TRASTUZUMAB_VH[:30] + "A" * 90, TRASTUZUMAB_VH],
            "vl_sequence": [None, TRASTUZUMAB_VL],
            "format": ["vhh", "paired"],
            "composite": [2.0, 1.0],
        }
    )
    out = build_tickets(ranked, n=4)
    assert out["design_id"].to_list() == ["d1-ft"]


def test_ticket_tag_documents_the_source():
    assert "S85N+R87T" in TICKET_TAG
    assert len(FRAMEWORK_MUTATIONS) == 2
