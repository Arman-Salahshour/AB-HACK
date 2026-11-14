"""IMGT numbering with `immunum` — the same library `abguard` uses.

`abguard` compares the **IMGT CDR3**, which for trastuzumab is ``SRWGGDGFYAMDY``
(13 residues): the 10-residue window the hackathon tables call ``cdrh3``, plus the
``SR`` and ``Y`` flanks. Use :func:`cdr3` whenever a number has to agree with
`abguard`, and the plain window when working with the provided tables.
"""

from __future__ import annotations

from functools import lru_cache

import immunum

from abhack.constants import TRASTUZUMAB_CDRH3_WINDOW

SCHEME = "imgt"
CHAINS = ["H", "K", "L"]

#: Flanks IMGT adds around the 10-residue window on the trastuzumab framework.
CDR3_PREFIX = "SR"
CDR3_SUFFIX = "Y"
TRASTUZUMAB_CDR3_IMGT = CDR3_PREFIX + TRASTUZUMAB_CDRH3_WINDOW + CDR3_SUFFIX


@lru_cache(maxsize=1)
def _annotator() -> immunum.Annotator:
    return immunum.Annotator(chains=CHAINS, scheme=SCHEME)


def segment(sequence: str):
    """Framework/CDR segmentation of one variable domain."""
    return _annotator().segment(sequence)


def cdr3(sequence: str) -> str | None:
    """The IMGT CDR3 of a variable domain, or None if it cannot be numbered."""
    seg = segment(sequence)
    if seg.error or not seg.cdr3:
        return None
    return seg.cdr3


def numbering_confidence(sequence: str) -> float:
    """Confidence that the sequence is an antibody domain (hard filter: >= 0.1)."""
    result = _annotator().number(sequence)
    return float(getattr(result, "confidence", 0.0) or 0.0)


def window_to_imgt(window: str) -> str:
    """Add the trastuzumab framework flanks to a 10-residue CDRH3 window."""
    return f"{CDR3_PREFIX}{window}{CDR3_SUFFIX}"


def imgt_to_window(cdr3_imgt: str) -> str:
    """Strip the trastuzumab framework flanks, if present."""
    s = cdr3_imgt
    s = s.removeprefix(CDR3_PREFIX)
    s = s.removesuffix(CDR3_SUFFIX)
    return s
