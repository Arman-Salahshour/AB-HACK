"""Local liability penalty for ranking — an approximation, never a gate.

`abguard liability-report` prints one aggregate table per submission file, with
no per-design breakdown, so it cannot be joined into ranking directly. This
module reimplements the motif side of the same LAP tags sequence-locally.

What is counted, and why only the redesigned loop:

* every optimisation design shares trastuzumab's framework, CDRH1/2 and the
  whole VL, so the liabilities that *differ* between designs are exactly those
  introduced by the new CDRH3 window — the penalty below counts motifs there;
* abguard itself counts CDR-region motifs (calibrated on a two-design probe,
  2026-11-14: a CDRH3 carrying NG/DS/DP/M moved exactly the DeAmdH/Isom/FragH/
  MetOx counts, framework motifs did not), so this mirrors the authority's
  own focus;
* the authoritative check stays `abguard liability-report` at validation time
  (`scripts/validate.sh`). A design is never dropped for liabilities here —
  the penalty only orders designs within a track.

Tags and weights (severity per abguard's report):

===============  ============================  ========  ====
tag              motif                         severity  weight
===============  ============================  ========  ====
``deamd_high``   ``NG``                         high      3.0
``isom``         ``DG``/``DS``/``DT``           high      3.0
``frag_high``    ``DP``                         high      3.0
``ngly``         ``N[^P][ST]``                 high      3.0
``deamd_medium`` ``NS``/``NT``                 medium    1.0
``trp_ox``       ``W``                         medium    1.0
``met_ox``       ``M``                         medium    1.0
``intbind``      ``RGD``/``LDV``               low       0.25
===============  ============================  ========  =====

Usage::

    from abhack.liabilities import add_liabilities
    df = add_liabilities(df)   # -> liab_cdrh3, liab_detail columns
"""

from __future__ import annotations

import re

import polars as pl

#: tag -> (motifs, severity weight). Order matters only for the detail string.
LIABILITY_MOTIFS: dict[str, tuple[tuple[str, ...], float]] = {
    "deamd_high": (("NG",), 3.0),
    "isom": (("DG", "DS", "DT"), 3.0),
    "frag_high": (("DP",), 3.0),
    "deamd_medium": (("NS", "NT"), 1.0),
    "trp_ox": (("W",), 1.0),
    "met_ox": (("M",), 1.0),
    "intbind": (("RGD", "LDV"), 0.25),
}

#: N-linked glycosylation sequon, weighted like the other high-severity tags.
_NGLY = re.compile("N[^P][ST]")


def count_motifs(sequence: str) -> dict[str, int]:
    """Occurrences of every liability motif in ``sequence``."""
    if not sequence:
        return {tag: 0 for tag in LIABILITY_MOTIFS} | {"ngly": 0}
    out: dict[str, int] = {}
    for tag, (motifs, _w) in LIABILITY_MOTIFS.items():
        out[tag] = sum(sequence.count(m) for m in motifs)
    out["ngly"] = len(_NGLY.findall(sequence))
    return out


def penalty(counts: dict[str, int]) -> float:
    """Severity-weighted total. ``ngly`` joins the high-severity tags."""
    total = sum(
        w * counts.get(tag, 0)
        for tag, (_motifs, w) in LIABILITY_MOTIFS.items()
    )
    return total + 3.0 * counts.get("ngly", 0)


def liability_penalty(cdrh3: str | None) -> tuple[float, str]:
    """``(penalty, human-readable detail)`` for one redesigned loop."""
    if not cdrh3:
        return 0.0, ""
    counts = count_motifs(cdrh3)
    detail = ";".join(
        f"{tag}={n}" for tag, n in counts.items() if n
    ) or "clean"
    return round(penalty(counts), 2), f"{detail} [{cdrh3}]"


def add_liabilities(df: pl.DataFrame) -> pl.DataFrame:
    """Append ``liab_cdrh3`` (penalty) and ``liab_detail`` per design.

    VHH rows have no ``cdrh3`` window; their IMGT CDR3 (``cdr3_imgt``) is
    scored instead, so de novo tables can use the same column.
    """
    seqs = (
        df["cdrh3"]
        if "cdrh3" in df.columns
        else df["cdr3_imgt"]
        if "cdr3_imgt" in df.columns
        else None
    )
    if seqs is None:
        return df.with_columns(
            liab_cdrh3=pl.lit(0.0, dtype=pl.Float64),
            liab_detail=pl.lit("", dtype=pl.String),
        )
    pairs = [liability_penalty(s) for s in seqs.to_list()]
    return df.with_columns(
        liab_cdrh3=pl.Series("liab_cdrh3", [p for p, _ in pairs], dtype=pl.Float64),
        liab_detail=pl.Series("liab_detail", [d for _, d in pairs], dtype=pl.String),
    )
