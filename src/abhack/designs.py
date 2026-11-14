"""Read RFD4-Proteina design runs into a flat table.

A design run writes::

    <run>/generation/inference_samples/<name>_<i>/<name>_<i>_generated.cif
    <run>/generation/inference_metrics/per_example_metrics.csv

and, after refolding::

    <run>/../scored-<model>/refold_metrics.csv

Usage::

    python -m abhack.designs collect runs/novel-s1 -o designs/novel-s1.csv
"""

from __future__ import annotations

import re
from pathlib import Path

import gemmi
import polars as pl

from abhack.constants import (
    MIN_CDRH3_MUTATIONS_DE_NOVO,
    NOVELTY_IDENTITY_THRESHOLD,
    TRASTUZUMAB_VH,
    TRASTUZUMAB_VL,
)
from abhack.novelty import classify, identity, unpaired_cysteines
from abhack.numbering import segment

#: Framework residues flanking the redesigned CDRH3 window in trastuzumab's VH.
CDRH3_PREFIX = "YYCSR"
CDRH3_SUFFIX = "YWGQG"
_CDRH3_RE = re.compile(f"{CDRH3_PREFIX}(.*?){CDRH3_SUFFIX}")


def chain_sequences(cif_path: Path | str) -> dict[str, str]:
    """One-letter sequence per polymer chain of a (m)mCIF file."""
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    out: dict[str, str] = {}
    for chain in st[0]:
        poly = chain.get_polymer()
        seq = gemmi.one_letter_code(poly.extract_sequence())
        if seq.strip("X"):
            out[chain.name] = seq.upper()
    return out


def extract_cdrh3(vh: str) -> str | None:
    """The redesigned window between the fixed framework anchors."""
    m = _CDRH3_RE.search(vh)
    return m.group(1) if m else None


def _looks_like(seq: str, reference: str, *, anchor: int = 25) -> bool:
    """Cheap test: does ``seq`` start like ``reference``'s framework 1?"""
    return bool(seq) and seq[:anchor] == reference[:anchor]


def parse_design(cif_path: Path) -> dict:
    """One row per generated structure: chains, CDRs and basic liabilities.

    Handles both shapes this repo designs:

    * **paired** — a CDRH3 redesign on trastuzumab's framework, where ``cdrh3``
      is the 10-residue window the hackathon tables use and the novelty gate
      compares against;
    * **vhh** — a de novo single domain, which has no such window. Its loops are
      read with `immunum`, so ``cdr3_imgt`` is filled and ``cdrh3`` is null.

    The antigen chains of a complex are ignored: only chains `immunum` numbers
    as an antibody domain are considered, longest first.
    """
    chains = chain_sequences(cif_path)

    # Antibody domains only — this drops HER2 chains C/D from a de novo complex.
    domains: list[tuple[str, str, object]] = []
    for name, seq in chains.items():
        if len(seq) < 60:
            continue
        seg = segment(seq)
        if seg and not seg.error and seg.cdr3:
            domains.append((name, seq, seg))
    domains.sort(key=lambda d: len(d[1]), reverse=True)

    vh = next((s for _, s, _ in domains if _looks_like(s, TRASTUZUMAB_VH)), None)
    vl = next((s for _, s, _ in domains if _looks_like(s, TRASTUZUMAB_VL)), None)
    if vh is None and domains:  # a de novo VHH, or a framework that moved
        vh = domains[0][1]
    seg = next((g for _, s, g in domains if s == vh), None)

    cdrh3 = extract_cdrh3(vh) if vh else None  # trastuzumab window, or None
    fmt = "paired" if vl else ("vhh" if vh else "unknown")
    return {
        "design_id": cif_path.parent.name,
        "cif": str(cif_path),
        "format": fmt,
        "cdrh3": cdrh3,
        "cdr1": seg.cdr1 if seg else None,
        "cdr2": seg.cdr2 if seg else None,
        "cdr3_imgt": seg.cdr3 if seg else None,
        "vh_sequence": vh,
        "vl_sequence": vl,
        "vh_length": len(vh) if vh else None,
        "total_length": (len(vh) if vh else 0) + (len(vl) if vl else 0),
        "identity_to_trastuzumab_vh": round(identity(vh, TRASTUZUMAB_VH), 4) if vh else None,
        "unpaired_cys_vh": unpaired_cysteines(vh) if vh else None,
        "n_chains": len(chains),
    }


def collect(run_dir: Path | str) -> pl.DataFrame:
    """Every design in a run, joined with its generation and refold metrics."""
    run_dir = Path(run_dir)
    samples = sorted(run_dir.glob("**/inference_samples/*/*_generated.cif"))
    if not samples:
        raise FileNotFoundError(f"no *_generated.cif under {run_dir}")

    df = pl.DataFrame([parse_design(p) for p in samples])

    gen_metrics = next(run_dir.glob("**/inference_metrics/per_example_metrics.csv"), None)
    if gen_metrics is not None:
        mdf = pl.read_csv(gen_metrics, infer_schema_length=10_000)
        key = next((c for c in ("example_id", "sample_id", "name", "id") if c in mdf.columns), None)
        if key is not None:
            df = df.join(mdf.rename({key: "design_id"}), on="design_id", how="left")

    for scored in sorted(run_dir.glob("**/refold_metrics.csv")):
        model = scored.parent.name.replace("scored-", "")
        sdf = pl.read_csv(scored, infer_schema_length=10_000)
        # refold_metrics.csv keys designs by `query_id`; `design` is a path.
        key = next(
            (c for c in ("query_id", "design_id", "example_id", "name", "id") if c in sdf.columns),
            None,
        )
        if key is None:
            continue
        keep = [
            c
            for c in ("min_ipae", "iptm", "binder_plddt", "binder_ptm", "backbone_rmsd", "pass")
            if c in sdf.columns
        ]
        sdf = sdf.select([key, *keep]).rename(
            {key: "design_id", **{c: f"{c}_{model}" for c in keep}}
        )
        df = df.join(sdf, on="design_id", how="left")

    return df


def add_novelty(df: pl.DataFrame, provided: frozenset[str] | None = None) -> pl.DataFrame:
    """Append the local novelty verdict for every row.

    A trastuzumab-window CDRH3 is graded against the provided data. A de novo
    VHH has no such window: its track is decided by domain identity instead, and
    the 2-mutation bar is against every known CDR3, which only `abguard` can
    check — so it is marked ``de_novo`` and left for `scripts/validate.sh`.
    """
    if provided is None:
        from abhack.data import provided_cdrh3

        provided = provided_cdrh3()

    ident = (
        df["identity_to_trastuzumab_vh"].to_list()
        if "identity_to_trastuzumab_vh" in df.columns
        else [None] * df.height
    )
    verdicts = []
    for cdrh3, dom_id in zip(df["cdrh3"].to_list(), ident):
        if cdrh3:
            verdicts.append(classify(cdrh3, provided).as_dict())
            continue
        de_novo = dom_id is not None and dom_id < NOVELTY_IDENTITY_THRESHOLD
        verdicts.append(
            {
                "cdrh3": None,
                "track": "de_novo" if de_novo else "unknown",
                "identity_to_trastuzumab_cdrh3": None,
                "mutations_from_trastuzumab": None,
                "min_mutations_to_provided": None,
                "required_mutations": MIN_CDRH3_MUTATIONS_DE_NOVO if de_novo else None,
                "passes": None,  # abguard decides; see scripts/validate.sh
                "reason": (
                    f"de novo domain (identity {dom_id:.2f} < 0.70); "
                    "CDR3 novelty must be checked with abguard"
                    if de_novo
                    else "could not classify: no CDRH3 window and no domain identity"
                ),
            }
        )
    cols = pl.DataFrame(verdicts).drop("cdrh3", strict=False)
    return pl.concat([df, cols], how="horizontal")


def main() -> None:  # pragma: no cover - thin CLI
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("command", choices=["collect"])
    ap.add_argument("run_dir")
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("--no-novelty", action="store_true")
    args = ap.parse_args()

    df = collect(args.run_dir)
    if not args.no_novelty:
        df = add_novelty(df)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(args.output)
        print(f"{df.height} designs -> {args.output}")
    else:
        print(df)


if __name__ == "__main__":  # pragma: no cover
    main()
