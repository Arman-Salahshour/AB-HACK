#!/usr/bin/env python3
"""Build the final 30-design ranking and the hand-in CSV.

Policy (in order of authority):
  1. Optimisation parents: the best-50 designs that have a Boltz-2 refold, de-duplicated on the (VH, VL) pair,
     in composite order (binding models, interface PAE, pH-gated, liabilities).
  2. The top 4 parents are each followed by a framework twin (VH S85N+R87T, Ng et al. mAbs 2025) - the only
     literature-backed route to a thermostability gain. Twins share the parent's CDRH3 and binding evidence.
  3. De novo designs go last (AF2-IG judged none of them credible), in their existing evidence order.

Writes designs/final-30.csv (full detail, for the page) and designs/final-submission.csv (organisers' format:
sequence_id,kind,heavy_sequence,light_sequence,design_description,track). Does not touch designs/submission.csv.

    python scripts/final30.py
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

import importlib.util

from abhack.tickets import apply_framework_mutations

_spec = importlib.util.spec_from_file_location("results_radar", Path(__file__).with_name("results_radar.py"))
_rr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_rr)

N_PARENTS, N_TWINS, N_TOTAL = 18, 4, 30


def rd(p: str) -> pl.DataFrame:
    return pl.read_csv(p, infer_schema_length=10_000)


def mol1_lookup() -> dict[tuple[str, str], dict]:
    """Mol-1 (AF3-family, unguided) metrics keyed by (VH, VL) for every design that was checked."""
    out: dict[tuple[str, str], dict] = {}
    pairs = []
    a = Path("designs/af3-check.csv")
    if a.exists():
        chk = rd(str(a)).filter(pl.col("model") == "mol1")
        top12 = rd("designs/best50.csv").head(12)
        meta = dict(zip(top12["design_id"], zip(top12["vh_sequence"], top12["vl_sequence"], strict=True), strict=True))
        pairs.append((chk, meta))
    e = Path("designs/af3-extra.csv")
    if e.exists():
        chk = rd(str(e)).filter(pl.col("model") == "mol1")
        ex = rd("designs/_mol1_extra.csv")
        meta = dict(zip(ex["design_id"], zip(ex["vh_sequence"], ex["vl_sequence"], strict=True), strict=True))
        pairs.append((chk, meta))
    for chk, meta in pairs:
        for r in chk.to_dicts():
            if r["design_id"] in meta:
                out[meta[r["design_id"]]] = r
    return out


def why_parent(r: dict, rank: int, m: dict | None) -> str:
    bits = []
    kd = r["kd_pred_nM"]
    bits.append(f"binder probability {r['p_high']:.2f}, predicted K_D {kd:.0f} nM")
    bits.append(f"Boltz-2 interface PAE {r['min_ipae_boltz2']:.2f} A")
    if m:
        bits.append(f"Mol-1 (no template) places it in the trastuzumab pose, {m['pose_rmsd_A']:.1f} A")
    else:
        bits.append("Mol-1 check not run")
    sel = r.get("ph_log10_sel")
    if sel is not None:
        bits.append(f"pH selectivity {sel:+.2f} ({'slight preference for 6.5' if sel >= 0.1 else 'no clear preference'})")
    liab = int(r["liab_cdrh3"] or 0)
    bits.append("no CDRH3 liability" if liab == 0 else f"{liab} CDRH3 liabilit{'y' if liab == 1 else 'ies'}")
    return "; ".join(bits)


def main() -> None:
    best = rd("designs/best50.csv")
    parents = (best.filter(pl.col("min_ipae_boltz2").is_not_null())
                   .unique(subset=["vh_sequence", "vl_sequence"], keep="first", maintain_order=True)
                   .head(N_PARENTS))
    m1 = mol1_lookup()

    rows: list[dict] = []
    p_list = parents.to_dicts()
    seq_of: dict[tuple, int] = {}
    for i, r in enumerate(p_list):
        base = {
            "track": "optimization", "kind": "paired", "design_id": r["design_id"], "cdrh3": r["cdrh3"],
            "vh": r["vh_sequence"], "vl": r["vl_sequence"], "p_high": r["p_high"], "kd_pred_nM": r["kd_pred_nM"],
            "ipae_boltz2": r["min_ipae_boltz2"], "ph_sel": r.get("ph_log10_sel"), "liab_cdrh3": int(r["liab_cdrh3"] or 0),
            "composite": r["composite"], "framework": "trastuzumab",
        }
        _cif = _rr.refold_cif(r)
        base["loop_gap"] = _rr.loop_gap(_cif) if _cif else None
        mm = m1.get((r["vh_sequence"], r["vl_sequence"]))
        base["mol1_iptm"] = mm["iptm"] if mm else None
        base["mol1_pose_A"] = mm["pose_rmsd_A"] if mm else None
        base["mol1_epitope_overlap"] = mm["epitope_jaccard"] if mm else None
        base["why"] = why_parent(r, i + 1, mm)
        rows.append(base)
        seq_of[(r["vh_sequence"], r["vl_sequence"])] = len(rows) - 1
        if i < N_TWINS:
            tw = dict(base)
            tw.update(design_id=r["design_id"] + "-ft", vh=apply_framework_mutations(r["vh_sequence"]), framework="S85N+R87T")
            tw["why"] = ("Same CDRH3 and binding evidence as its parent, plus VH S85N+R87T: the literature "
                         "(Ng et al., mAbs 2025) reports higher Tm with HER2 binding retained. Creates an NLT "
                         "glycosylation sequon (flag; cell-free expression does not glycosylate).")
            rows.append(tw)

    # de novo, last: reuse the validated evidence order from the committed table
    fr = rd("designs/final-ranked.csv").filter(pl.col("track") == "de_novo").sort("submission_order")
    af = rd("designs/af2ig-denovo.csv")
    bg = rd("designs/boltzgen-scored.csv")
    for r in fr.head(N_TOTAL - len(rows)).to_dicts():
        did = r["design_id"]
        hit = af.filter(pl.col("design_id") == did)
        if hit.height:
            iptm, ipae, src = hit["af2ig_iptm"][0], hit["af2ig_ipae"][0], "RFD4 de novo VHH"
        else:
            b = bg.filter(pl.col("boltzgen_id") == did)
            iptm, ipae, src = (b["af2ig_iptm"][0], b["af2ig_ipae"][0], "BoltzGen de novo VHH") if b.height else (None, None, "de novo VHH")
        rows.append({
            "track": "denovo", "kind": "vhh", "design_id": did, "cdrh3": r.get("cdr3_imgt") or r.get("cdrh3") or "",
            "vh": r["vh_sequence"], "vl": "", "p_high": None, "kd_pred_nM": None, "ipae_boltz2": None,
            "ph_sel": None, "liab_cdrh3": None, "composite": None, "framework": "VHH",
            "mol1_iptm": None, "mol1_pose_A": None, "mol1_epitope_overlap": None, "loop_gap": None,
            "why": (f"{src}. AF2-IG from the designed coordinates: ipTM {iptm:.2f}, interface PAE {ipae:.0f} A, below the "
                    "credibility bar (ipTM >= 0.6, PAE <= 3 A). Included last, as a lottery ticket only."
                    if iptm is not None else f"{src}; no independent structure check passed. Included last."),
        })

    df = pl.DataFrame(rows, infer_schema_length=None).with_row_index("rank", offset=1)
    df = df.with_columns(rank=pl.col("rank").cast(pl.Int64))
    df.write_csv("designs/final-30.csv")

    sub = df.select(
        pl.format("d{}", pl.col("rank").cast(pl.String).str.zfill(2)).alias("sequence_id"),
        "kind", pl.col("vh").alias("heavy_sequence"), pl.col("vl").alias("light_sequence"),
        pl.when(pl.col("track") == "optimization")
          .then(pl.format("CDRH3 {} on trastuzumab framework{}", pl.col("cdrh3"),
                          pl.when(pl.col("framework") == "S85N+R87T").then(pl.lit(" + VH S85N,R87T")).otherwise(pl.lit(""))))
          .otherwise(pl.format("de novo VHH against HER2, rank {}", pl.col("rank"))).alias("design_description"),
        "track",
    )
    sub.write_csv("designs/final-submission.csv")
    print(f"{df.height} rows -> designs/final-30.csv, designs/final-submission.csv")
    print(df.group_by("track").len())
    with pl.Config(tbl_rows=40, tbl_cols=-1, fmt_str_lengths=22, float_precision=2):
        print(df.select("rank", "design_id", "track", "cdrh3", "p_high", "ipae_boltz2", "mol1_pose_A", "ph_sel", "liab_cdrh3"))


if __name__ == "__main__":
    main()
