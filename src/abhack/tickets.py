"""Thermostability "ticket" designs: paper-validated framework mutations.

Every optimisation design inherits trastuzumab's framework, so its
thermostability is the baseline's — it can tie that axis but not win it,
and S-tier (all three axes improved) is out of reach. The fix costs two
framework residues the novelty gate never looks at:

    VH S85N + R87T  (sequential numbering from the VH N-terminus)

from Ng et al., mAbs 2025, "Tuning antibody stability and function by
rational designs of framework mutations" (bioRxiv 2025.03.13.642419), where
this exact double mutant on trastuzumab:

* raised the melting temperature (DSF) vs WT — the thermostability axis;
* sits in FW3, solvent-exposed, on the face opposite the CDRs — 12 residues
  from the redesigned CDRH3 window, which is untouched here;
* retained HER2 binding (FACS on SKBR3 comparable to WT) and ADCC.

Their other candidates are deliberately not used: R50S+R59N raised Tm but
*lost* binding, and A40R+K43D had low production yield.

The double mutant creates an ``NLT`` N-glycosylation sequon. That is a
liability flag, not a hard filter, and it is why this is a *ticket* design
submitted alongside — never instead of — its parent: Adaptyv's cell-free
expression does not glycosylate, and the Tm benefit was measured on the
protein actually produced. The parent design keeps its slot either way.

Usage::

    python -m abhack.tickets designs/opt-ranked.csv -n 4 -o designs/opt-tickets.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

#: (1-based VH position, wild-type residue, mutant) — the paper's construct.
FRAMEWORK_MUTATIONS: tuple[tuple[int, str, str], ...] = (
    (85, "S", "N"),
    (87, "R", "T"),
)

#: What to record in the design_description of a ticket.
TICKET_TAG = "fw S85N+R87T (Ng mAbs 2025: +Tm, HER2 binding retained)"


def apply_framework_mutations(vh: str) -> str:
    """Apply the validated stabilising pair, verifying the wild-type first.

    Raises ValueError on a mismatch so a framework drift upstream can never
    be silently "stabilised" — the evidence is for these two residues on
    this framework, nothing else.
    """
    for pos, wt, mut in FRAMEWORK_MUTATIONS:
        i = pos - 1
        if vh[i] != wt:
            raise ValueError(
                f"VH position {pos} is {vh[i]!r}, expected {wt!r}: not the "
                "trastuzumab framework the mAbs 2025 mutations were validated on"
            )
    out = list(vh)
    for pos, wt, mut in FRAMEWORK_MUTATIONS:
        out[pos - 1] = mut
    return "".join(out)


def build_tickets(ranked: pl.DataFrame, *, n: int = 4) -> pl.DataFrame:
    """Framework-mutated twins of the top ``n`` ranked optimisation designs.

    The CDRH3 — and therefore every sequence model score and the novelty
    verdict — is the parent's; only ``vh_sequence`` changes. ``design_id``
    gains a ``-ft`` suffix so the twin never collides with its parent in a
    join or a submission table.
    """
    if "cdrh3" not in ranked.columns or "vh_sequence" not in ranked.columns:
        raise ValueError("ranked table needs vh_sequence and cdrh3 columns")

    # Only paired (trastuzumab-framework) designs can carry the ticket.
    fmt = (
        ranked["format"].cast(pl.String).str.to_lowercase()
        if "format" in ranked.columns
        else pl.Series(["paired"] * ranked.height)
    )
    parents = ranked.with_columns(_fmt=fmt).filter(pl.col("_fmt") == "paired")
    if parents.height < n:
        n = parents.height
        print(f"tickets: only {n} paired designs available")
    take = parents.head(n)

    rows = []
    for r in take.to_dicts():
        vh = r["vh_sequence"]
        try:
            mutated = apply_framework_mutations(vh)
        except ValueError as e:
            print(f"tickets: skipping {r.get('design_id')}: {e}")
            continue
        r = dict(r)
        r["design_id"] = f"{r['design_id']}-ft"
        r["vh_sequence"] = mutated
        r["framework_ticket"] = True
        rows.append(r)
    if not rows:
        raise ValueError("no tickets could be built")
    return pl.DataFrame(rows)


def main() -> None:  # pragma: no cover - thin CLI
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("ranked", help="ranked design table (abhack.ranking output)")
    ap.add_argument("-n", type=int, default=4, help="how many top designs to twin")
    ap.add_argument("-o", "--output", default="designs/opt-tickets.csv")
    args = ap.parse_args()

    df = pl.read_csv(args.ranked, infer_schema_length=10_000)
    tickets = build_tickets(df, n=args.n)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    tickets.write_csv(args.output)
    show = [c for c in ("design_id", "cdrh3", "framework_ticket") if c in tickets.columns]
    print(tickets.select(show))
    print(f"{tickets.height} tickets -> {args.output}")


if __name__ == "__main__":  # pragma: no cover
    main()
