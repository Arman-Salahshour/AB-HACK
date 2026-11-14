"""A thermostability / developability proxy, and why it is only a proxy.

Thermostability is one of the three judged axes and nothing in this repo scored it.
There is **no provided thermostability data**, so nothing here is fitted or validated —
this module computes transparent sequence and liability features and combines them into
an explicit *ordering*, never a predicted Tm. Say so in the pitch; claiming a Tm from
these features would be the kind of overclaim the judgment criterion punishes.

What is actually known about the quantity we are proxying
---------------------------------------------------------
For the optimisation track the headroom is small, and that is measured, not guessed.
McConnell et al. (mAbs 2014, PMID 25517312) grafted CDRs onto an engineered stable
framework across ten therapeutics: adalimumab gained 9.9 °C, denosumab 7.0 °C,
cetuximab 6.9 °C — and **trastuzumab 0.8 °C**. Trastuzumab's framework is already
near-optimal, so a CDRH3-only redesign on that framework should sit at baseline.

Consequence for ranking: on a fixed framework the only lever is the CDRH3 loop itself,
so that is what this module scores. Ng et al. (mAbs 2025, PMID 40653772) do design
stabilising *framework* mutations for trastuzumab specifically, which is the route to a
real gain — but note they also found a light-chain framework mutation that kept antigen
binding and lost effector function, so framework changes are not free of consequence.

Two components
--------------
``liability_burden`` — from ``abguard team-validate``, which already computes LAP
sequence liabilities per design and reports the region each one sits in. A liability in
a CDR is weighted above one in the framework: trastuzumab itself carries several
framework liabilities, so those are not disqualifying, while a *new* one introduced into
the redesigned loop is a change we chose.

``cdrh3_features`` — length, net charge, GRAVY hydrophobicity, aromatic and Gly/Pro
content of the designed loop. Longer, more flexible, more charged loops are the ones
that tend to cost stability; this is directional prior knowledge, not a fitted model.

Usage::

    python -m abhack.stability report designs/opt-s11-submission.csv
    python -m abhack.stability annotate designs/opt-s11.csv -o designs/opt-s11.stab.csv
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

import polars as pl

#: Kyte-Doolittle hydrophobicity.
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
AROMATIC = set("FWY")
#: Net charge at neutral pH, counting His as neutral.
POSITIVE, NEGATIVE = set("KR"), set("DE")

#: How much a liability counts, by abguard's own severity label.
SEVERITY_WEIGHT = {"high": 3.0, "medium": 1.5, "low": 0.5}

#: The region we actually redesigned. Measured on the seed-11 shortlist, abguard reports
#: cdr1 and cdr2 liabilities at exactly 2 per design and fr4 at exactly 1 — they are
#: trastuzumab's own and identical across every design, so they carry no information and
#: including them only dilutes the term. cdr3 varied (50 occurrences across 61 designs),
#: because cdr3 is the loop we changed. Weight it; report the rest for context only.
DESIGNED_REGIONS = {"cdr3"}
DESIGNED_MULTIPLIER = 4.0


def sequence_features(seq: str | None) -> dict[str, float | None]:
    """Transparent per-loop features. Empty/None sequence gives all-None."""
    if not seq:
        return {k: None for k in
                ("length", "net_charge", "gravy", "frac_aromatic", "frac_gly_pro")}
    n = len(seq)
    return {
        "length": float(n),
        "net_charge": float(sum(c in POSITIVE for c in seq) - sum(c in NEGATIVE for c in seq)),
        "gravy": round(sum(KD.get(c, 0.0) for c in seq) / n, 4),
        "frac_aromatic": round(sum(c in AROMATIC for c in seq) / n, 4),
        "frac_gly_pro": round(sum(c in "GP" for c in seq) / n, 4),
    }


def _abguard_bundle() -> str | None:
    """abguard is a uv tool, so its reference bundle is not importable from here."""
    base = Path(os.path.expanduser("~/.local/share/uv/tools/abguard/lib"))
    for p in sorted(base.glob("python*/site-packages/abguard/reference_data")):
        if p.is_dir():
            return str(p)
    return None


def liabilities_from_abguard(
    submission_csv: Path | str, *, timeout: int = 900
) -> pl.DataFrame:
    """Per-design liability and novelty columns, straight from `abguard team-validate`.

    abguard already computes all of this and the pipeline threw it away. Returns one row
    per ``sequence_id``; raises if abguard is unavailable rather than inventing numbers.

    Do **not** pass ``--drop`` here: it still writes a row per surviving design but
    leaves every liability column null, which silently turns the whole developability
    signal into zeros. Without it abguard reports a non-compliant table on stderr, exits
    0 anyway and writes the full annotated set, which is what we want.
    """
    ref = _abguard_bundle()
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "tv.json"
        cmd = ["abguard", "team-validate"]
        if ref:
            cmd += ["--reference-db", ref]
        # --all-liabilities is what populates n_liabilities / liability_tags / the
        # per-liability region detail. Without it every liability column comes back null.
        cmd += ["--all-liabilities", "-o", str(out), str(submission_csv)]
        subprocess.run(cmd, capture_output=True, timeout=timeout, check=False)
        # abguard writes the path it was given; tolerate a suffix swap.
        written = out if out.exists() else next(Path(tmp).glob("*"), None)
        if written is None:
            raise RuntimeError(f"abguard produced no output for {submission_csv}")
        rows = json.loads(Path(written).read_text())

    recs = []
    for r in rows:
        liabs = r.get("liabilities") or []
        weighted = 0.0
        designed_hits = 0
        inherited_hits = 0
        for liab in liabs:
            w = SEVERITY_WEIGHT.get(str(liab.get("severity", "")).lower(), 1.0)
            if str(liab.get("region", "")).lower() in DESIGNED_REGIONS:
                w *= DESIGNED_MULTIPLIER
                designed_hits += 1
            else:
                inherited_hits += 1
            weighted += w
        recs.append(
            {
                "sequence_id": r.get("sequence_id"),
                "is_valid": r.get("is_valid"),
                "n_liabilities": r.get("n_liabilities"),
                "n_high_severity": r.get("n_high_severity"),
                # the loop we designed: the only part that varies, so the only part that ranks
                "n_designed_liabilities": designed_hits,
                # trastuzumab's own, kept for context — do not rank on these
                "n_inherited_liabilities": inherited_hits,
                "liability_burden": round(weighted, 3),
                "liability_tags": ";".join(r.get("liability_tags") or []),
                # abguard's own novelty verdict — more authoritative than ours.
                "abguard_novelty_class": r.get("novelty_class"),
                "abguard_cdrh3_min_distance": r.get("cdrh3_min_distance"),
                "abguard_cdrh3_nearest_source": r.get("cdrh3_nearest_source"),
                "abguard_is_novel_enough": r.get("is_novel_enough"),
                "abguard_identity_to_trastuzumab": r.get("identity_to_trastuzumab"),
            }
        )
    return pl.DataFrame(recs)


def add_cdrh3_features(df: pl.DataFrame, column: str = "cdrh3") -> pl.DataFrame:
    """Append `cdrh3_*` sequence features for whichever loop column exists."""
    col = column if column in df.columns else ("cdr3_imgt" if "cdr3_imgt" in df.columns else None)
    if col is None:
        raise ValueError("need a `cdrh3` or `cdr3_imgt` column to score the loop")
    feats = [sequence_features(s) for s in df[col].to_list()]
    block = pl.DataFrame(feats).rename(lambda c: f"cdrh3_{c}")
    return pl.concat([df, block], how="horizontal")


def stability_proxy(df: pl.DataFrame) -> pl.DataFrame:
    """A single `stability_proxy` ordering column. Lower is better. NOT a Tm.

    Deliberately a small, readable sum of standardised directional terms rather than a
    fitted model, because there is nothing to fit it on. Every term's sign is a stated
    prior, so a reader can disagree with a specific term instead of with a black box.
    """
    terms, used = [], []

    def z(col: str, *, higher_is_worse: bool) -> None:
        if col not in df.columns:
            return
        s = df[col].cast(pl.Float64, strict=False)
        if s.null_count() == s.len():
            return
        mean, std = s.mean(), s.std()
        if std in (None, 0.0):
            return
        zs = (s - mean) / std
        terms.append(zs if higher_is_worse else -zs)
        used.append(f"{'+' if higher_is_worse else '-'}{col}")

    # Chemical liabilities we introduced: strictly bad.
    z("liability_burden", higher_is_worse=True)
    z("n_designed_liabilities", higher_is_worse=True)
    # Loop geometry: longer and more flexible loops cost stability.
    z("cdrh3_length", higher_is_worse=True)
    z("cdrh3_frac_gly_pro", higher_is_worse=True)
    # Charge extremes and exposed hydrophobicity both hurt; aromatics pack well.
    z("cdrh3_net_charge", higher_is_worse=True)
    z("cdrh3_gravy", higher_is_worse=True)
    z("cdrh3_frac_aromatic", higher_is_worse=False)
    # Fold confidence from the refold, where present.
    for c in df.columns:
        if c.startswith("binder_plddt"):
            z(c, higher_is_worse=False)
            break

    if not terms:
        return df.with_columns(
            stability_proxy=pl.lit(None, dtype=pl.Float64),
            stability_terms=pl.lit("no usable columns"),
        )
    total = terms[0]
    for t in terms[1:]:
        total = total + t
    return df.with_columns(
        stability_proxy=(total / len(terms)).round(4),
        stability_terms=pl.lit(" ".join(used)),
    )


def annotate(
    designs_csv: Path | str, submission_csv: Path | str | None = None
) -> pl.DataFrame:
    """Design table + CDRH3 features + (optionally) abguard liabilities + the proxy."""
    df = pl.read_csv(designs_csv, infer_schema_length=10_000)
    df = add_cdrh3_features(df)
    if submission_csv is not None:
        liab = liabilities_from_abguard(submission_csv)
        if "sequence_id" in df.columns:
            df = df.join(liab, on="sequence_id", how="left")
        else:
            # No shared key: attach the burden distribution as a report, not a join.
            print(
                f"note: {designs_csv} has no sequence_id, so liabilities were read but "
                f"not joined ({liab.height} rows available)"
            )
    return stability_proxy(df)


def main() -> None:  # pragma: no cover - thin CLI
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("command", choices=["report", "annotate"])
    ap.add_argument("designs")
    ap.add_argument("--submission", default=None, help="table to run abguard against")
    ap.add_argument("-o", "--output", default=None)
    args = ap.parse_args()

    df = annotate(args.designs, args.submission)
    cols = [c for c in df.columns if c.startswith(("cdrh3_", "liability", "stability", "n_"))]
    keep = [c for c in ("design_id", "sequence_id", "cdrh3") if c in df.columns] + cols
    view = df.select(keep)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(args.output)
        print(f"{df.height} rows -> {args.output}")
    print(view.head(12))
    if "stability_proxy" in df.columns and df["stability_proxy"].null_count() < df.height:
        print(f"\nterms used: {df['stability_terms'][0]}")
        print("stability_proxy: lower is better. This is an ORDERING, not a Tm.")
        print("On a fixed trastuzumab framework expect ~baseline thermostability")
        print("(McConnell 2014: trastuzumab gained only 0.8 C from stable-framework grafting).")


if __name__ == "__main__":  # pragma: no cover
    main()
