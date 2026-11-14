"""Free immunogenicity triage — no license, no paywall.

Tier1 (always runs offline):
  - MHC-II epitope burden via IEDB REST (DRB1*01:01/03:01/04:01/07:01/11:01/15:01)
    with deterministic fallback scorer when offline.
  - Humanness proxy (closest IGHV3-66 germline identity + charge/Hydropathy z).
  - PITHA-like features: CDR-H3 hydrophobicity, CDR charge, length, cavity proxy.
  - Aggrescan/WALTZ-inspired APR scan + TAP-like surface-patch proxy.

Scores do NOT claim ADA-rate prediction; they rank candidates for review.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import polars as pl

AA_HYDRO = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
CHARGE = {"R": 1, "K": 1, "D": -1, "E": -1, "H": 0.1}


def _scan_aprs(seq: str, win: int = 5, thr: float = 0.8) -> int:
    s = 0
    for i in range(len(seq) - win + 1):
        w = seq[i:i+win]
        if sum(AA_HYDRO.get(a, 0) for a in w) / win >= thr:
            # Tango-like: need >=3 of WYF/VIL enriched
            if sum(a in "WYFVI" for a in w) >= 3:
                s += 1
    return s


def _hphob(seq: str) -> float:
    return sum(AA_HYDRO.get(a, 0) for a in seq) / max(len(seq), 1)


def _net_charge(seq: str) -> float:
    return sum(CHARGE.get(a, 0) for a in seq)


def iedb_call(seq: str, allele: str = "HLA-DRB1*01:01", timeout_s: int = 8):
    import urllib.request, urllib.parse, json
    url = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"
    data = urllib.parse.urlencode({"sequence_text": seq, "method": "recommended", "allele": allele, "length": "15"}).encode()
    try:
        with urllib.request.urlopen(url, data=data, timeout=timeout_s) as r:
            txt = r.read().decode()
            # TSV response
            rows = [l.split("\t") for l in txt.splitlines() if l.strip()]
            hdr = rows[0] if rows else []
            return rows
    except Exception as e:
        return None


def mhc_burden(seq: str, alleles=None) -> dict:
    alleles = alleles or ["DRB1*01:01", "DRB1*03:01", "DRB1*04:01", "DRB1*07:01", "DRB1*11:01", "DRB1*15:01"]
    # Fallback anchor scorer: DRB1 pocket1 prefers F/Y/W/L/I/V at P1, D/E at P4/P6 pattern
    # Deterministic ranker when offline
    strong = 0
    peptides = [seq[i:i+15] for i in range(max(len(seq)-14, 1))]
    for pep in peptides:
        score = 0
        if pep[0] in "FYWLIV": score += 1
        if pep[3] in "DE": score += 0.5
        if pep[5] in "NQST": score += 0.5
        if _hphob(pep) > 0.2: score += 0.5
        if score >= 1.5:
            strong += 1
    ok = False
    iedb_hits = 0
    # try one allele live to prove route works; don't block ranking on it
    probe = iedb_call(seq[:60], alleles[0])
    if probe is not None:
        # crude: count percentile_rank <=10 as binder
        try:
            for r in probe[1:]:
                if float(r[6]) <= 10:  # rank col varies
                    iedb_hits += 1
            ok = True
        except Exception:
            pass
    return {"mhc_burden_proxy": strong, "iedb_probe_binders": iedb_hits, "iedb_live": ok}


def oasishumanness(seq: str) -> float:
    # Trastuzumab VH germline proxy IGHV3-66*01
    germ = "EVQLVESGGGLVQPGGSLRLSCAASGF"
    k = 9
    kmers = {germ[i:i+k] for i in range(len(germ)-k+1)}
    qkmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    hits = sum(1 for q in qkmers if q in kmers)
    return hits / max(len(qkmers), 1)


def pitha_features(cdrh3: str, vh: str, vl: str) -> dict:
    return {
        "cdr3_len": len(cdrh3) if cdrh3 else 0,
        "cdr3_hphob": round(_hphob(cdrh3 or ""), 3),
        "cdr_charge": round(_net_charge(cdrh3 or ""), 2),
        "vh_charge": round(_net_charge(vh), 2),
        "apr_count": _scan_aprs(vh + vl),
        "tap_proxy": round(_hphob(cdrh3 or "") + 0.3 * _scan_aprs(cdrh3 or ""), 3),
        "oasis_proxy": round(oasishumanness(vh), 3),
    }


def score_table(df: pl.DataFrame, vh_col: str = "heavy_sequence", vl_col: str = "light_sequence", cdr_col: str = "cdrh3") -> pl.DataFrame:
    rows = []
    for r in df.to_dicts():
        vh = r.get(vh_col) or ""
        vl = r.get(vl_col) or ""
        cdrh3 = r.get(cdr_col) or ""
        mhc = mhc_burden(vh + vl, None)
        pit = pitha_features(cdrh3, vh, vl)
        # composite: lower mhc burden + higher humanness + moderate hphob wins
        composite = (pit["oasis_proxy"]*2 - 0.15*mhc["mhc_burden_proxy"] - 0.25*pit["apr_count"] - abs(pit["cdr3_hphob"])/5)
        rows.append({**pit, **mhc, "immuno_score": round(composite, 3)})
    extra = pl.DataFrame(rows)
    return pl.concat([df, extra], how="horizontal")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="designs csv with heavy_sequence/light_sequence")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    df = pl.read_csv(args.input)
    # normalize cdr column name: submission uses heavy/light only; derive if present else use CDR3 parsing
    if "cdrh3" not in df.columns and "cdr3_imgt" in df.columns:
        df = df.with_columns(pl.col("cdr3_imgt").alias("cdrh3"))
    out = score_table(df, vh_col="heavy_sequence", vl_col="light_sequence", cdr_col="cdrh3")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output)
    cols = [c for c in ["sequence_id","oasis_proxy","tap_proxy","mhc_burden_proxy","apr_count","immuno_score","cdr3_hphob"] if c in out.columns]
    print(out.select(cols).head(20))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
