#!/usr/bin/env python3
"""pH-selectivity (6.5 vs 7.4) for the designs already in the submission, both tracks.

Uses each design's own designed complex (optimisation: RFD4 cif; de novo: RFD4 or BoltzGen
cif), not a threaded trastuzumab model, so de novo scaffolds are scored on their own pose.
Writes designs/ph-submission.csv.
"""
import glob, tempfile
from pathlib import Path
import gemmi, polars as pl
from abhack.phscore import score_variant, PH_LOW, PH_HIGH

def cif_for(did, cif):
    if cif: return cif
    for f in glob.glob("designs/*.csv"):
        try: d = pl.read_csv(f, infer_schema_length=10000)
        except Exception: continue
        if {"design_id", "cif"} <= set(d.columns):
            hit = d.filter((pl.col("design_id") == did) & pl.col("cif").is_not_null())
            if hit.height: return hit["cif"][0]
    hit = glob.glob(f"runs/boltzgen/design/rank*_{did}.cif")
    return hit[0] if hit else None

def main():
    f = pl.read_csv("designs/final-ranked.csv", infer_schema_length=10000)
    rows = []
    for r in f.select("submission_order", "design_id", "format", "cif").to_dicts():
        cif = cif_for(r["design_id"], r["cif"])
        out = {**r, "cif": cif, "track": "denovo" if r["format"] == "vhh" else "optimisation"}
        if not cif or not Path(cif).exists():
            rows.append({**out, "ph_error": "no structure"}); continue
        st = gemmi.read_structure(cif); st.remove_ligands_and_waters()
        names = [c.name for c in st[0]]
        antigen = ("A",) if "shard" in r["design_id"] else ("C", "D")
        with tempfile.TemporaryDirectory() as d:
            pdb = Path(d) / "x.pdb"; st.write_pdb(str(pdb))
            try:
                s = score_variant("", "", pdb=pdb, antigen=antigen)
                rows.append({**out, "chains": ",".join(names), **{f"ph_{k}": v for k, v in s.items()}, "ph_error": ""})
            except Exception as e:  # noqa: BLE001
                rows.append({**out, "chains": ",".join(names), "ph_error": f"{type(e).__name__}: {e}"[:160]})
    df = pl.DataFrame(rows, infer_schema_length=None); df.write_csv("designs/ph-submission.csv")
    with pl.Config(tbl_rows=30, tbl_cols=-1, fmt_str_lengths=50):
        print(f"pH {PH_LOW} vs {PH_HIGH}")
        print(df.select("submission_order", "design_id", "track", "chains", "ph_log10_sel", "ph_his_summary", "ph_error"))
main()
