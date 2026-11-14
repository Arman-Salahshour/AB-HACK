"""Histidine-switch candidates: bind HER2 at tumour pH (6.8), weaker at 7.2.

Mechanism. A histidine (pKa ~6.0 free) gains a positive charge as pH falls. If that
His sits against a HER2 carboxylate (Glu558, Asp560, Asp570, Asp596, Glu598 in
1N8Z) the salt bridge forms when protonated and is lost when neutral, so binding
weakens as pH rises. Free His is only ~14 % protonated at 6.8 vs ~6 % at 7.2
(about 2.5x), so a useful switch needs the environment to lift the His pKa near 7.
That lift is what the pKa stage must estimate; this module only proposes candidates.

Constraints it respects:

* CDRH3 must be novel: >= 3 Levenshtein mutations from every provided loop and
  < 70 % window identity to trastuzumab (>= 4 of the 10 positions changed). A
  His scan outside CDRH3 alone would fail the novelty gate, so every candidate
  carries a changed CDRH3 and the His can live in CDRH3, CDR1/2 or VL.
* Positions come from ``scripts/ph_interface.py`` (paratope residues contacting
  or within 7 A of a HER2 carboxylate). Sequence numbering is sequential, 1-based.

The binder classifier (``p_high``) is trained on FACS at the screen's own pH, so
it says nothing about pH dependence; it only stops the search leaving the
trastuzumab-binding region of sequence space.

    python -m abhack.phswitch -n 400 -o designs/ph-candidates.csv
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import polars as pl

from abhack.constants import (
    CDRH3_WINDOW_END,
    CDRH3_WINDOW_START,
    TRASTUZUMAB_CDRH3_WINDOW,
    TRASTUZUMAB_VH,
    TRASTUZUMAB_VL,
)
from abhack.novelty import classify

#: (chain, 1-based position, wild-type) paratope positions facing a HER2 acid.
#: CDRH3 positions are handled separately through the 10-residue window.
FRAMEWORK_HIS_SITES: tuple[tuple[str, int, str], ...] = (
    ("VH", 33, "Y"),
    ("VH", 50, "R"),
    ("VH", 52, "Y"),
    ("VH", 57, "Y"),
    ("VH", 58, "T"),
    ("VH", 59, "R"),
    ("VL", 30, "N"),
    ("VL", 31, "T"),
    ("VL", 94, "T"),
)
#: CDRH3 window indices (0-based) that contact HER2 / sit near Asp570: W99 G103 F104 Y105.
CDRH3_HIS_INDEX: tuple[int, ...] = (0, 4, 5, 6)
#: Indices that stay free to mutate for novelty without touching the contact core.
CDRH3_FREE_INDEX: tuple[int, ...] = (1, 2, 3, 7, 8, 9)
#: Conservative substitutions for the free positions.
FREE_CHOICES = "SGATNDYVLIRKEQ"


def _check_wt() -> None:
    for chain, pos, wt in FRAMEWORK_HIS_SITES:
        seq = TRASTUZUMAB_VH if chain == "VH" else TRASTUZUMAB_VL
        if seq[pos - 1] != wt:
            raise ValueError(f"{chain}{pos} is {seq[pos - 1]!r}, expected {wt!r}")


def _splice(vh: str, window: str) -> str:
    return vh[:CDRH3_WINDOW_START] + window + vh[CDRH3_WINDOW_END:]


def _apply(seq: str, edits: dict[int, str]) -> str:
    out = list(seq)
    for pos, aa in edits.items():
        out[pos - 1] = aa
    return "".join(out)


def enumerate_candidates(
    n: int = 400, seed: int = 1, provided: frozenset[str] | None = None
) -> pl.DataFrame:
    """Sample His-bearing designs that clear the local novelty gate."""
    _check_wt()
    rng = random.Random(seed)
    provided = provided if provided is not None else frozenset()
    seen: set[tuple[str, str]] = set()
    rows: list[dict] = []
    tries = 0
    while len(rows) < n and tries < n * 200:
        tries += 1
        win = list(TRASTUZUMAB_CDRH3_WINDOW)
        # 1-3 contact-core histidines in CDRH3.
        n_core = rng.choice((1, 1, 2, 3))
        core = rng.sample(CDRH3_HIS_INDEX, n_core)
        for i in core:
            win[i] = "H"
        # enough extra changes to reach >= 4 of 10 positions, plus up to 2 more.
        need = max(0, 4 - n_core)
        for i in rng.sample(CDRH3_FREE_INDEX, need + rng.choice((0, 0, 1, 2))):
            win[i] = rng.choice(FREE_CHOICES)
        window = "".join(win)
        # 0-2 histidines elsewhere in the paratope.
        sites = rng.sample(FRAMEWORK_HIS_SITES, rng.choice((0, 0, 1, 2)))
        vh_e = {p: "H" for c, p, _ in sites if c == "VH"}
        vl_e = {p: "H" for c, p, _ in sites if c == "VL"}
        vh = _apply(_splice(TRASTUZUMAB_VH, window), vh_e)
        vl = _apply(TRASTUZUMAB_VL, vl_e)
        key = (vh, vl)
        if key in seen:
            continue
        verdict = classify(window, provided)
        if not verdict.passes:
            continue
        seen.add(key)
        h_cdrh3 = window.count("H")
        rows.append(
            {
                "design_id": f"phsw_{len(rows):04d}",
                "cdrh3": window,
                "vh_sequence": vh,
                "vl_sequence": vl,
                "format": "paired",
                "ph_his_cdrh3": h_cdrh3,
                "ph_his_sites": ",".join(f"{c}{p}H" for c, p, _ in sorted(sites, key=lambda s: (s[0], s[1]))),
                "ph_his_total": h_cdrh3 + len(sites),
                "mutations_from_trastuzumab_cdrh3": verdict.mutations_from_trastuzumab,
                "min_mutations_to_provided": verdict.min_mutations_to_provided,
            }
        )
    return pl.DataFrame(rows)


def add_binder_score(df: pl.DataFrame) -> pl.DataFrame:
    """Append ``p_high`` from the FACS classifier when it is trained."""
    from abhack.constants import ROOT
    from abhack.model import BinderClassifier

    path = ROOT / "models" / "binder.joblib"
    if not path.exists() or df.is_empty():
        return df
    p = BinderClassifier.load(path).p_high(df["cdrh3"].to_list())
    return df.with_columns(p_high=pl.Series(p))


def main() -> None:  # pragma: no cover - thin CLI
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-n", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("-o", "--out", type=Path, default=Path("designs/ph-candidates.csv"))
    a = ap.parse_args()
    from abhack.data import provided_cdrh3

    df = add_binder_score(enumerate_candidates(a.n, a.seed, provided_cdrh3()))
    df.write_csv(a.out)
    print(f"{df.height} candidates -> {a.out}")
    if "p_high" in df.columns:
        print(df.select("p_high").describe())
    print(df.group_by("ph_his_total").len().sort("ph_his_total"))


if __name__ == "__main__":
    main()
