#!/usr/bin/env python
"""Export candidate structures as PDB files with their results beside them.

Two structures per candidate, which answer different questions:

* ``*_designed.pdb``  — what RFD4-Proteina drew: the Fv on HER2, the designed pose.
* ``*_boltz2.pdb``    — what Boltz-2 rebuilt from the sequence alone. Agreement
  between the two (backbone RMSD) is the evidence that the design is real.

Each file is named by submission rank, and ``results.csv`` plus ``README.md``
carry every metric for the pair.

    python scripts/export_structures.py designs/opt-s11-ranked-final.csv \
        --run runs/opt-novel-s11 --out designs/structures
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import gemmi
import polars as pl

#: Submission order: designs that rebuild their designed structure come first.
PREFERRED = ["her2_cdrh3_53", "her2_cdrh3_50", "her2_cdrh3_55",
             "her2_cdrh3_4", "her2_cdrh3_20", "her2_cdrh3_64", "her2_cdrh3_10"]

CHAIN_NOTE = {"A": "trastuzumab VL (unchanged)", "B": "VH, CDRH3 redesigned",
              "C": "HER2 domain IV", "D": "HER2 domain IV (second segment)"}


def to_pdb(src: Path, dest: Path) -> dict:
    st = gemmi.read_structure(str(src))
    st.setup_entities()
    st.assign_label_seq_id()
    dest.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(dest))
    return {ch.name: len(ch) for ch in st[0]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("table")
    ap.add_argument("--run", required=True)
    ap.add_argument("--out", default="designs/structures")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    run, out = Path(args.run), Path(args.out)
    df = pl.read_csv(args.table, infer_schema_length=10_000)

    rank = {d: i for i, d in enumerate(PREFERRED)}
    rows = sorted(df.to_dicts(), key=lambda r: (rank.get(r["design_id"], 99),
                                                r.get("min_ipae_boltz2") or 9))
    if args.limit:
        rows = rows[: args.limit]

    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    manifest = []
    for i, r in enumerate(rows, start=1):
        did = r["design_id"]
        short = did.replace("her2_cdrh3_", "D")
        stem = f"{i:02d}_{short}"

        designed = next(run.glob(f"**/inference_samples/{did}/{did}_generated.cif"), None)
        refold = next(run.glob(f"boltz2/{did}/**/*_model.cif"), None)

        chains = to_pdb(designed, out / f"{stem}_designed.pdb") if designed else {}
        if refold:
            to_pdb(refold, out / f"{stem}_boltz2.pdb")

        manifest.append({
            "rank": i,
            "design": short,
            "design_id": did,
            "cdrh3": r.get("cdrh3"),
            "designed_pdb": f"{stem}_designed.pdb" if designed else "",
            "refolded_pdb": f"{stem}_boltz2.pdb" if refold else "",
            "min_ipae": r.get("min_ipae_boltz2"),
            "iptm": r.get("iptm_boltz2"),
            "binder_plddt": r.get("binder_plddt_boltz2"),
            "backbone_rmsd": r.get("backbone_rmsd_boltz2"),
            "p_high": r.get("p_high"),
            "pred_kd_nM": r.get("kd_pred"),
            "mutations_from_trastuzumab": r.get("mutations_from_trastuzumab"),
            "min_mutations_to_provided": r.get("min_mutations_to_provided"),
            "chains": " ".join(f"{k}:{v}" for k, v in chains.items()),
        })

    m = pl.DataFrame(manifest)
    m.write_csv(out / "results.csv")

    lines = [
        "# HER2 candidate structures",
        "",
        "Trastuzumab variants with a redesigned CDRH3. Every design passes all",
        "`abguard` filters. Files are named by submission rank.",
        "",
        "## Two files per candidate",
        "",
        "| Suffix | What it is |",
        "| --- | --- |",
        "| `_designed.pdb` | what RFD4-Proteina drew — the designed pose on HER2 |",
        "| `_boltz2.pdb` | what Boltz-2 rebuilt **from the sequence alone** |",
        "",
        "Backbone RMSD between the two is the evidence that matters: under 1 Å means an",
        "independent model reproduced the design. A high RMSD with a good `min_ipae`",
        "means Boltz-2 is confident about a *different* interface than the one designed.",
        "",
        "## Chains",
        "",
        "| Chain | Contents |",
        "| --- | --- |",
    ]
    lines += [f"| {k} | {v} |" for k, v in CHAIN_NOTE.items()]
    lines += [
        "",
        "## Candidates",
        "",
        "| Rank | Design | CDRH3 | min_ipae | iptm | pLDDT | RMSD (Å) | P(high) | pred K_D (nM) | muts |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in manifest:
        def f(x, n=3):
            return "–" if x is None else (f"{x:.{n}f}" if isinstance(x, float) else str(x))
        lines.append(
            f"| {r['rank']} | {r['design']} | `{r['cdrh3']}` | {f(r['min_ipae'])} | "
            f"{f(r['iptm'])} | {f(r['binder_plddt'], 1)} | {f(r['backbone_rmsd'], 2)} | "
            f"{f(r['p_high'], 2)} | {f(r['pred_kd_nM'], 1)} | {r['mutations_from_trastuzumab']} |"
        )
    lines += [
        "",
        "## Reading the numbers",
        "",
        "- `min_ipae` — interface predicted aligned error, lower is better; under 0.5 is good.",
        "- `iptm` — interface confidence, higher is better.",
        "- `binder_plddt` — how well-formed the antibody itself is.",
        "- `RMSD` — backbone agreement between designed and refolded structure.",
        "- `P(high)` — FACS classifier (524k measured variants) on whether the loop binds.",
        "- `pred K_D` — SPR regressor (190 measured affinities). Trastuzumab measures",
        "  0.029–0.178 nM, and no design is predicted to beat it.",
        "",
        "Compare scores only within one folding model; cut-offs are calibrated per model.",
    ]
    (out / "README.md").write_text("\n".join(lines) + "\n")

    size = sum(p.stat().st_size for p in out.glob("*")) / 1e6
    print(f"{len(manifest)} candidates -> {out}  ({size:.1f} MB)")
    print(f"  {len(list(out.glob('*_designed.pdb')))} designed, "
          f"{len(list(out.glob('*_boltz2.pdb')))} refolded, + results.csv, README.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
