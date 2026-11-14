#!/usr/bin/env python
"""Assemble the de novo evidence table: BoltzGen + RFD4 VHHs, AF2-IG-scored.

Every de novo candidate this repo can field, one table, ranked by the only
orthogonal evidence that exists for the track: AF2-IG ipTM / interface PAE
from the designed coordinates (plus BoltzGen's own refold confidence for its
designs). Both generators failed the credibility bar outright — best AF2-IG
ipTM 0.46 / PAE 16 A against a real-binder bar of ipTM >= 0.6 / PAE <= 3 A —
so these are explicitly ranked *lottery tickets* for a thin track, and the
table says so in a column rather than pretending otherwise.

Kappa-like BoltzGen outputs (a lone light-chain domain is not a VHH and
would fail abguard's format check) are dropped here rather than at
validation time.

    python scripts/build_denovo_table.py -o designs/denovo-ranked.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

from abhack.liabilities import add_liabilities
from abhack.numbering import segment

#: heavy-chain frameworks start like this; kappa/lambda do not (cheap filter,
#: cross-checked with immunum below).
_KAPPA_LIKE = ("DIQMTQ", "GDIVMTQ", "EIVLTQ", "AIRMTQ", "SYELTQ")


def looks_heavy(seq: str) -> bool:
    if any(seq.startswith(p) for p in _KAPPA_LIKE):
        return False
    seg = segment(seq)
    return bool(seg and not seg.error and seg.cdr3)


def build(boltzgen: Path, rfd4_tables: list[Path], af2ig: Path) -> pl.DataFrame:
    bg = pl.read_csv(boltzgen, infer_schema_length=10_000)
    keep_ok = [looks_heavy(s) for s in bg["vh_sequence"].to_list()]
    dropped = bg.height - sum(keep_ok)
    if dropped:
        print(f"boltzgen: dropped {dropped} non-heavy (kappa-like) design(s)")
    bg = bg.filter(pl.Series(keep_ok)).with_columns(
        design_id=pl.col("boltzgen_id"),
        vl_sequence=pl.lit(None, dtype=pl.String),
        format=pl.lit("vhh"),
        source=pl.lit("boltzgen"),
        cdrh3=pl.lit(None, dtype=pl.String),
        passes=pl.lit(None, dtype=pl.Boolean),
        generator_iptm=pl.col("boltzgen_iptm"),
        generator_ipae=pl.col("boltzgen_ipae"),
    )

    frames = [bg]
    scored = pl.read_csv(af2ig, infer_schema_length=10_000)
    for table in rfd4_tables:
        run_name = Path(table).stem  # denovo-vhh-s1 -> matches af2ig's `run`
        df = pl.read_csv(table, infer_schema_length=10_000)
        df = (df.filter(pl.col("format") == "vhh")
                .with_columns(run=pl.lit(run_name))
                .join(scored, on=["run", "design_id"], how="inner")
                .with_columns(
                    vl_sequence=pl.lit(None, dtype=pl.String),
                    cdrh3=pl.lit(None, dtype=pl.String),
                    source=pl.col("run"),
                    generator_iptm=pl.lit(None, dtype=pl.Float64),
                    generator_ipae=pl.lit(None, dtype=pl.Float64),
                ))
        frames.append(df)

    common = ["design_id", "vh_sequence", "vl_sequence", "format", "cdrh3", "cdr3_imgt",
              "source", "af2ig_iptm", "af2ig_ipae", "af2ig_plddt", "generator_iptm",
              "generator_ipae", "passes"]
    out = pl.concat([f.select([c for c in common if c in f.columns]) for f in frames],
                    how="diagonal_relaxed")
    return add_liabilities(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--boltzgen", default="designs/boltzgen-scored.csv")
    ap.add_argument("--rfd4", nargs="*", default=["designs/denovo-vhh-s1.csv",
                                                  "designs/denovo-vhh-s2.csv"])
    ap.add_argument("--af2ig", default="designs/af2ig-denovo.csv")
    ap.add_argument("-o", "--output", default="designs/denovo-ranked.csv")
    args = ap.parse_args()

    df = build(Path(args.boltzgen), [Path(p) for p in args.rfd4], Path(args.af2ig))
    # evidence order: AF2-IG first, BoltzGen's own confidence as the tiebreak
    df = df.sort(["af2ig_iptm", "generator_iptm"], descending=[True, True], nulls_last=True)
    df.write_csv(args.output)
    print(df.select(["design_id", "source", "af2ig_iptm", "af2ig_ipae", "generator_iptm", "liab_cdrh3"]))
    print(f"{df.height} de novo candidates -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
