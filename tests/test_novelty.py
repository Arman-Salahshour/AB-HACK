from abhack.constants import TRASTUZUMAB_CDRH3_WINDOW as TZ
from abhack.constants import TRASTUZUMAB_VH, TRASTUZUMAB_VL
from abhack.designs import extract_cdrh3
from abhack.novelty import (
    classify,
    hamming,
    identity,
    levenshtein,
    min_mutations_to_provided,
    mutations_from_trastuzumab,
    unpaired_cysteines,
)
from abhack.numbering import TRASTUZUMAB_CDR3_IMGT, cdr3, imgt_to_window, window_to_imgt


def test_distances():
    assert hamming("ABC", "ABD") == 1
    assert hamming("ABC", "ABCD") == 4  # length mismatch is a full miss
    assert levenshtein("ABC", "ABCD") == 1
    assert identity(TZ, TZ) == 1.0
    assert identity("WGGDGFYAMD", "WGGDGFYAMA") == 0.9


def test_mutations_from_trastuzumab():
    assert mutations_from_trastuzumab(TZ) == 0
    assert mutations_from_trastuzumab("AGGDGFYAMD") == 1
    assert mutations_from_trastuzumab("AAADGFYAMD") == 3


def test_min_mutations_to_provided():
    provided = {"WGGDGFYAMD", "AAGDGFYAMD"}
    assert min_mutations_to_provided("WGGDGFYAMD", provided) == 0
    assert min_mutations_to_provided("AAGDGFYAMA", provided) == 1


def test_classify_rejects_trastuzumab_itself():
    v = classify(TZ, frozenset({TZ}))
    assert v.track == "optimisation"
    assert not v.passes


def test_classify_needs_three_mutations_from_every_provided():
    # 4 changes from trastuzumab, but 2 from a provided variant -> abguard rejects it.
    provided = frozenset({TZ, "QGNDGFYAMD"})
    v = classify("QGNDGFVAMV", provided)
    assert v.min_mutations_to_provided == 2
    assert not v.passes

    # 4 changes from everything provided -> passes.
    v = classify("QGNDGFVAMV", frozenset({TZ}))
    assert v.min_mutations_to_provided == 4
    assert v.passes


def test_cdrh3_identity_rule_is_toggleable():
    # 3 of 10 changed -> identity exactly 70 %, which the task description rejects
    # but the shipped abguard accepts.
    three = "AAADGFYAMD"
    assert not classify(three, frozenset({TZ})).passes
    assert classify(three, frozenset({TZ}), enforce_cdrh3_identity=False).passes


def test_de_novo_track_needs_only_two_mutations():
    v = classify("QGNDGFVAMV", frozenset({TZ}), domain_identity_to_known=0.5)
    assert v.track == "de_novo"
    assert v.required_mutations == 2


def test_unpaired_cysteines():
    assert unpaired_cysteines("ACDEC") is False
    assert unpaired_cysteines("ACDE") is True


def test_extract_cdrh3_from_framework():
    assert extract_cdrh3(TRASTUZUMAB_VH) == TZ
    assert extract_cdrh3("NO ANCHORS HERE") is None


def test_immunum_agrees_with_our_anchors():
    """abguard numbers with immunum; the IMGT CDR3 adds the SR/Y flanks."""
    assert cdr3(TRASTUZUMAB_VH) == TRASTUZUMAB_CDR3_IMGT == "SRWGGDGFYAMDY"
    assert window_to_imgt(TZ) == TRASTUZUMAB_CDR3_IMGT
    assert imgt_to_window(TRASTUZUMAB_CDR3_IMGT) == TZ
    assert cdr3(TRASTUZUMAB_VL) == "QQHYTTPPT"
