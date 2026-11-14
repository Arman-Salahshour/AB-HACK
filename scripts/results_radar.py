#!/usr/bin/env python3
"""Radar (spider) charts of the four judged-and-risk axes, as an HTML fragment.

Axes, each scaled 0-1 with 1 = better:

* Affinity        FACS binder model p_high (trained on 525k measured variants). Trastuzumab = 0.995.
* Thermostability MD thermal score (400 K) where run; otherwise a *proxy*: CDRH3 loop confidence from the
                  Boltz-2 refold (framework pLDDT minus loop pLDDT; a gap of 0 scores 1, 10 scores 0).
* pH selectivity  propka log10 K(6.5)/K(7.4), mapped -0.2 -> 0 and +0.6 -> 1 (trastuzumab control 0 -> 0.25).
* Safety proxy    sequence-only: CDRH3 liabilities, CDRH3 hydrophobicity, CDRH3 positive charge. It is NOT a
                  toxicity prediction; target-mediated cardiotoxicity is not captured by it.

Writes designs/radar-data.csv and prints the HTML fragment to the path given (default stdout is the file
docs/results-radar.html). Run after best50.py; re-run when MD results arrive.

    python scripts/results_radar.py -n 6
"""

from __future__ import annotations

import argparse
import glob
import html
import math
from pathlib import Path

import gemmi
import polars as pl
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from abhack.constants import CDRH3_WINDOW_END, CDRH3_WINDOW_START, TRASTUZUMAB_CDRH3_WINDOW, TRASTUZUMAB_VH
from abhack.model import MODELS_DIR, BinderClassifier

AXES = ["Affinity", "Thermostability", "pH selectivity", "Safety proxy"]


def clip(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def ph_score(sel: float | None) -> float | None:
    return None if sel is None else clip((sel + 0.2) / 0.8)


def safety(cdrh3: str, liab: float) -> float:
    gravy = ProteinAnalysis(cdrh3).gravy()
    pos = sum(cdrh3.count(a) for a in "KR") - sum(cdrh3.count(a) for a in "DE")
    parts = [clip(liab / 3.0), clip((gravy + 0.5) / 2.5), clip(max(0, pos) / 2.0)]
    return 1.0 - sum(parts) / len(parts)


def loop_gap(cif: str) -> float | None:
    """Framework pLDDT minus CDRH3-loop pLDDT on the VH chain of a Boltz-2 model (b-factor = pLDDT)."""
    try:
        st = gemmi.read_structure(cif)
        vh = st[0]["B"]
        loop, fw = [], []
        for i, r in enumerate(vh):
            b = sum(a.b_iso for a in r) / len(r)
            (loop if CDRH3_WINDOW_START <= i < CDRH3_WINDOW_END else fw).append(b)
        scale = 100.0 if max(fw) > 1.5 else 1.0
        return (sum(fw) / len(fw) - sum(loop) / len(loop)) * (100.0 / scale)
    except Exception:  # noqa: BLE001
        return None


def refold_cif(row: dict) -> str | None:
    did = row["design_id"].removesuffix("-ft")
    run = Path(row["cif"]).parts[1]
    hits = sorted(glob.glob(f"runs/{run}-gated/boltz2/{did}/seed_*/{did}_seed_*_sample_1_model.cif")
                  + glob.glob(f"runs/{run}/boltz2/{did}/seed_*/{did}_seed_*_sample_1_model.cif"))
    return hits[0] if hits else None


def collect(n: int) -> pl.DataFrame:
    best = pl.read_csv("designs/best50.csv", infer_schema_length=10_000).head(n)
    rows = []
    for r in best.to_dicts():
        cif = refold_cif(r)
        gap = loop_gap(cif) if cif else None
        rows.append({
            "design_id": r["design_id"], "cdrh3": r["cdrh3"], "p_high": r["p_high"], "kd_pred_nM": r["kd_pred_nM"],
            "ph_sel": r.get("ph_log10_sel"), "loop_gap": gap, "liab_cdrh3": r["liab_cdrh3"],
            "min_ipae": r.get("min_ipae_boltz2"),
            "affinity": r["p_high"],
            "thermostability": None if gap is None else clip(1 - gap / 10.0),
            "ph": ph_score(r.get("ph_log10_sel")),
            "safety": safety(r["cdrh3"], r["liab_cdrh3"] or 0),
        })
    return pl.DataFrame(rows, infer_schema_length=None)


# ---------------------------------------------------------------------------------------------- SVG
CX, CY, R = 160, 128, 62


def _pt(i: int, v: float) -> tuple[float, float]:
    ang = -math.pi / 2 + i * 2 * math.pi / len(AXES)
    return CX + R * v * math.cos(ang), CY + R * v * math.sin(ang)


def radar(series: list[dict], title: str, desc: str) -> str:
    """series: dicts with name, values (4 floats|None), css (stroke var), dash, labels (4 str)."""
    g = []
    for ring in (0.25, 0.5, 0.75, 1.0):
        pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in (_pt(i, ring) for i in range(len(AXES))))
        g.append(f'<polygon points="{pts}" class="rg"/>')
    for i, a in enumerate(AXES):
        x, y = _pt(i, 1.0)
        g.append(f'<line x1="{CX}" y1="{CY}" x2="{x:.1f}" y2="{y:.1f}" class="ax"/>')
        lx, ly = _pt(i, 1.28)
        anchor = "middle" if abs(lx - CX) < 8 else ("start" if lx > CX else "end")
        g.append(f'<text x="{lx:.1f}" y="{ly + 4:.1f}" text-anchor="{anchor}" class="al">{html.escape(a)}</text>')
    for s in series:
        vals = s["values"]
        if any(v is None for v in vals):
            continue
        pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in (_pt(i, v) for i, v in enumerate(vals)))
        dash = ' stroke-dasharray="5 4"' if s.get("dash") else ""
        fill = f'fill="{s["css"]}" fill-opacity="{s.get("fill", 0.14)}"' if not s.get("dash") else 'fill="none"'
        g.append(f'<polygon points="{pts}" {fill} stroke="{s["css"]}" stroke-width="2" stroke-linejoin="round"{dash}/>')
        for i, v in enumerate(vals):
            x, y = _pt(i, v)
            tip = html.escape(f'{s["name"]} - {AXES[i]}: {s["labels"][i]}')
            g.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{s["css"]}" stroke="var(--surface)" '
                     f'stroke-width="2"><title>{tip}</title></circle>')
    return (f'<svg viewBox="-12 0 388 256" role="img" aria-label="{html.escape(title)}: {html.escape(desc)}" '
            f'class="radar">{"".join(g)}</svg>')


def fmt(v: float | None, spec: str = ".2f") -> str:
    return "n/a" if v is None else format(v, spec)


def small_multiple(r: dict, rank: int) -> str:
    vals = [r["affinity"], r["thermostability"], r["ph"], r["safety"]]
    labels = [f'p(high) {fmt(r["p_high"])}, predicted K<sub>D</sub> {fmt(r["kd_pred_nM"], ".0f")} nM'.replace("<sub>", "").replace("</sub>", ""),
              f'loop gap {fmt(r["loop_gap"], ".1f")} pLDDT (proxy)',
              f'propka selectivity {fmt(r["ph_sel"], "+.2f")}',
              f'{int(r["liab_cdrh3"] or 0)} CDRH3 liabilities']
    svg = radar([{"name": r["design_id"], "values": vals, "css": "var(--s1)", "labels": labels}],
                f'Design {r["design_id"]}', ", ".join(f"{a} {fmt(v)}" for a, v in zip(AXES, vals, strict=True)))
    return (f'<figure class="rcard">{svg}<figcaption><strong>#{rank} {html.escape(r["design_id"])}</strong>'
            f'<span class="m">{html.escape(r["cdrh3"])}</span></figcaption></figure>')


def table(df: pl.DataFrame) -> str:
    head = ("<tr><th>#</th><th>Design</th><th>CDRH3</th><th class='m'>Affinity</th><th class='m'>Thermo (proxy)</th>"
            "<th class='m'>pH</th><th class='m'>Safety</th><th class='m'>Pred. K<sub>D</sub> nM</th>"
            "<th class='m'>pH sel.</th><th class='m'>Loop gap</th></tr>")
    body = "".join(
        f"<tr><td>{i}</td><td class='m'>{html.escape(r['design_id'])}</td><td class='m'>{html.escape(r['cdrh3'])}</td>"
        f"<td class='m num'>{fmt(r['affinity'])}</td><td class='m num'>{fmt(r['thermostability'])}</td>"
        f"<td class='m num'>{fmt(r['ph'])}</td><td class='m num'>{fmt(r['safety'])}</td>"
        f"<td class='m num'>{fmt(r['kd_pred_nM'], '.1f')}</td><td class='m num'>{fmt(r['ph_sel'], '+.2f')}</td>"
        f"<td class='m num'>{fmt(r['loop_gap'], '.1f')}</td></tr>"
        for i, r in enumerate(df.to_dicts(), 1)
    )
    return f'<div class="scroll"><table><thead>{head}</thead><tbody>{body}</tbody></table></div>'


# --------------------------------------------------------------------------- proof-of-concept panel
def poc_panel() -> tuple[str, str]:
    """Reference vs parent vs twin from real MD results; a status note while they are incomplete."""
    th = Path("runs/thermal-poc/thermal_scores.csv")
    if not th.exists():
        return "", "pending"
    t = pl.read_csv(th, infer_schema_length=10_000)
    ids = {"ref": "trastuzumab_ref", "par": None, "twin": None}
    did = Path("runs/thermal-poc/design_id.txt").read_text().strip()
    ids["par"], ids["twin"] = did, f"{did}-ft"
    got = {k: t.filter(pl.col("design_id") == v) for k, v in ids.items()}
    if any(g.is_empty() for g in got.values()):
        return "", "pending"
    ts = {k: float(g["thermal_score"][0]) for k, g in got.items()}
    best = pl.read_csv("designs/best50.csv", infer_schema_length=10_000).filter(pl.col("design_id") == did).row(0, named=True)
    clf = BinderClassifier.load(MODELS_DIR / "binder.joblib")
    p_ref = float(clf.p_high([TRASTUZUMAB_CDRH3_WINDOW])[0])
    sel = {"ref": 0.0, "par": best["ph_log10_sel"], "twin": best["ph_log10_sel"]}
    aff = {"ref": p_ref, "par": best["p_high"], "twin": best["p_high"]}
    saf = {"ref": safety(TRASTUZUMAB_CDRH3_WINDOW, 3.0), "par": safety(best["cdrh3"], best["liab_cdrh3"] or 0)}
    saf["twin"] = saf["par"]
    meta = {"ref": ("Trastuzumab (1N8Z)", "var(--ref)", True), "par": (did, "var(--s1)", False),
            "twin": (f"{did}-ft (S85N+R87T)", "var(--s2)", False)}
    series = []
    for k in ("ref", "par", "twin"):
        name, css, dash = meta[k]
        vals = [aff[k], ts[k], ph_score(sel[k]), saf[k]]
        labels = [f"p(high) {fmt(aff[k])}", f"thermal MD score {fmt(ts[k], '.3f')}",
                  f"propka selectivity {fmt(sel[k], '+.2f')}", "sequence proxy"]
        series.append({"name": name, "values": vals, "css": css, "dash": dash, "labels": labels})
    svg = radar(series, "Proof of concept", "trastuzumab reference, top design and its framework twin")
    legend = "".join(
        f'<li><span class="sw" style="background:{s["css"]};{"opacity:.55" if s["dash"] else ""}"></span>{html.escape(s["name"])}</li>'
        for s in series)
    return f'<figure class="rbig">{svg}<ul class="legend">{legend}</ul></figure>', "ready"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=6)
    ap.add_argument("-o", "--out", default="docs/results-radar.html")
    a = ap.parse_args()
    df = collect(a.n)
    df.write_csv("designs/radar-data.csv")
    poc, state = poc_panel()
    cards = "".join(small_multiple(r, i) for i, r in enumerate(df.to_dicts(), 1))
    frag = {"cards": cards, "table": table(df), "poc": poc, "poc_state": state}
    Path(a.out).write_text(
        f"<!--cards-->{frag['cards']}<!--table-->{frag['table']}<!--poc:{state}-->{frag['poc']}")
    print(f"wrote {a.out}: {df.height} designs, poc={state}")
    print(df.select("design_id", "affinity", "thermostability", "ph", "safety", "loop_gap"))


if __name__ == "__main__":
    main()
