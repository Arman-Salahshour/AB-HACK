#!/usr/bin/env python3
"""Per-frame trajectories of the MD runs, as a CSV and two SVG line charts (HTML fragment).

* pH runs:    fraction of the starting Fv-HER2 heavy-atom contacts (< 4.5 A) still within 5.5 A, per frame.
* thermal:    fraction of native intra-Fv contacts (same definition as analyze_thermal.q_native), per frame.

Molecules are re-imaged first (chains wrap across the box independently). Partial runs are drawn as far as
they got and labelled by their simulated time. Writes runs/md-timeseries.csv and docs/md-timeseries.html.

    python scripts/md_timeseries.py
"""

from __future__ import annotations

import glob
import html
import importlib.util
import json
import re
from pathlib import Path

import mdtraj as md
import numpy as np
import polars as pl

spec = importlib.util.spec_from_file_location("analyze_thermal", Path(__file__).with_name("analyze_thermal.py"))
_at = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_at)

FRAME_NS = 0.01


def load(stem: str) -> md.Trajectory:
    t = md.load(f"{stem}.dcd", top=f"{stem}_top.pdb")
    if t.unitcell_lengths is not None:
        t.image_molecules(inplace=True, anchor_molecules=[max(t.topology.find_molecules(), key=len)])
    return t


def interface_q(t: md.Trajectory, antigen=("C", "D")) -> np.ndarray:
    heavy = lambda chains: [a.index for a in t.topology.atoms  # noqa: E731
                            if a.element.symbol != "H" and (a.residue.chain.chain_id in chains) == chains_in]
    chains_in = False
    fv = heavy(antigen)
    chains_in = True
    ag = heavy(antigen)
    x0 = t.xyz[0]
    d0 = np.linalg.norm(x0[fv][:, None, :] - x0[ag][None, :, :], axis=2)
    ii, jj = np.where(d0 < 0.45)
    pairs = np.stack([np.array(fv)[ii], np.array(ag)[jj]], axis=1)
    d = md.compute_distances(t, pairs)
    return (d < 0.55).mean(axis=1)


def main() -> None:
    rows = []
    for stem in sorted(glob.glob("runs/phmd/*_pH*_r*.dcd")):
        stem = stem[:-4]
        m = re.search(r"([^/]+)_pH([\d.]+)_r(\d+)$", stem)
        meta_path = Path(f"{stem}.json")
        antigen = tuple(json.load(open(meta_path)).get("antigen", ["C", "D"])) if meta_path.exists() else ("C", "D")
        q = interface_q(load(stem), antigen)
        rows += [{"kind": "ph", "run": f"{m.group(1)}_pH{m.group(2)}_r{m.group(3)}", "design": m.group(1),
                  "cond": float(m.group(2)), "rep": int(m.group(3)), "t_ns": (i + 1) * FRAME_NS, "q": float(v)}
                 for i, v in enumerate(q)]
    for stem in sorted(glob.glob("runs/thermal-poc/*_T*_r*.dcd")):
        stem = stem[:-4]
        m = re.search(r"([^/]+)_T(\d+)_r(\d+)$", stem)
        q = _at.q_native(load(stem))
        rows += [{"kind": "thermal", "run": f"{m.group(1)}_T{m.group(2)}_r{m.group(3)}", "design": m.group(1),
                  "cond": float(m.group(2)), "rep": int(m.group(3)), "t_ns": (i + 1) * FRAME_NS, "q": float(v)}
                 for i, v in enumerate(q)]
    df = pl.DataFrame(rows)
    df.write_csv("runs/md-timeseries.csv")
    print(df.group_by("kind", "design", "cond").agg(pl.len().alias("frames"), pl.col("t_ns").max().alias("ns")).sort("kind", "design", "cond"))
    Path("docs/md-timeseries.html").write_text(charts(df))


# ----------------------------------------------------------------------------------------------- SVG
W, H, ML, MR, MT, MB = 560, 280, 46, 120, 14, 40


def line_chart(df: pl.DataFrame, series: list[dict], title: str, desc: str, ymin: float, ymax: float, xmax: float,
               ylabel: str) -> str:
    sx = lambda x: ML + (W - ML - MR) * x / xmax  # noqa: E731
    sy = lambda y: MT + (H - MT - MB) * (1 - (y - ymin) / (ymax - ymin))  # noqa: E731
    g = []
    n_y = 5
    for k in range(n_y + 1):
        v = ymin + (ymax - ymin) * k / n_y
        g.append(f'<line x1="{ML}" x2="{W - MR}" y1="{sy(v):.1f}" y2="{sy(v):.1f}" class="gl"/>')
        g.append(f'<text x="{ML - 8}" y="{sy(v) + 4:.1f}" text-anchor="end" class="tk">{v:.2f}</text>')
    for k in range(0, int(round(xmax * 10)) + 1, 2):
        v = k / 10
        if v <= xmax + 1e-9:
            g.append(f'<text x="{sx(v):.1f}" y="{H - MB + 18}" text-anchor="middle" class="tk">{v:.1f}</text>')
    g.append(f'<text x="{(ML + W - MR) / 2:.0f}" y="{H - 6}" text-anchor="middle" class="tk">simulated time after equilibration (ns)</text>')
    g.append(f'<text x="12" y="{(MT + H - MB) / 2:.0f}" transform="rotate(-90 12 {(MT + H - MB) / 2:.0f})" text-anchor="middle" class="tk">{html.escape(ylabel)}</text>')
    label_y = []
    for s in series:
        runs = df.filter(pl.col("run").is_in(s["runs"]))
        last = None
        for run in s["runs"]:
            r = runs.filter(pl.col("run") == run).sort("t_ns")
            if r.is_empty():
                continue
            pts = " ".join(f"{sx(t):.1f},{sy(min(max(v, ymin), ymax)):.1f}" for t, v in zip(r["t_ns"], r["q"], strict=True))
            tip = html.escape(f'{s["name"]}, {run}: {r["t_ns"].max():.2f} ns, final Q {r["q"][-1]:.3f}')
            g.append(f'<polyline points="{pts}" fill="none" stroke="{s["css"]}" stroke-width="2" stroke-linejoin="round" '
                     f'stroke-linecap="round" opacity=".9"><title>{tip}</title></polyline>')
            last = (r["t_ns"].max(), r["q"][-1])
        if last:
            label_y.append((sy(min(max(last[1], ymin), ymax)), s))
    # direct labels at the right edge, nudged apart
    label_y.sort(key=lambda a: a[0])
    prev = -99
    for y, s in label_y:
        y = max(y, prev + 14)
        prev = y
        g.append(f'<rect x="{W - MR + 10}" y="{y - 8:.1f}" width="10" height="10" rx="2" fill="{s["css"]}"/>')
        g.append(f'<text x="{W - MR + 25}" y="{y + 1:.1f}" class="dl">{html.escape(s["name"])}</text>')
    return (f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{html.escape(title)}: {html.escape(desc)}" class="ts">'
            f'{"".join(g)}</svg>')


def charts(df: pl.DataFrame) -> str:
    ph = df.filter(pl.col("kind") == "ph")
    runs = lambda design, cond: sorted(ph.filter((pl.col("design") == design) & (pl.col("cond") == cond))["run"].unique().to_list())  # noqa: E731
    ph_series = [
        {"name": "Reference, pH 6.5", "css": "var(--s1)", "runs": runs("ctrl", 6.5)},
        {"name": "Reference, pH 7.4", "css": "var(--s2)", "runs": runs("ctrl", 7.4)},
        {"name": "Top design, pH 6.5", "css": "var(--s3)", "runs": runs("her2_cdrh3_28", 6.5)},
        {"name": "Design 43, pH 6.5", "css": "var(--s4)", "runs": runs("opt43", 6.5)},
    ]
    th = df.filter(pl.col("kind") == "thermal")
    tr = lambda cond: sorted(th.filter(pl.col("cond") == cond)["run"].unique().to_list())  # noqa: E731
    th_series = [
        {"name": "Reference, 310 K", "css": "var(--s1)", "runs": tr(310.0)},
        {"name": "Reference, 400 K", "css": "var(--s2)", "runs": tr(400.0)},
    ]
    a = line_chart(ph, ph_series, "Interface contacts retained", "Fv to HER2 contacts kept over time at two pH values",
                   0.5, 1.0, 1.0, "interface contacts retained")
    b = line_chart(th, th_series, "Fv native contacts retained", "trastuzumab Fv at 310 K and 400 K",
                   0.8, 1.0, 0.7, "native contacts retained")
    return (f'<figure class="tsfig"><h4>Interface contacts over time (two-pH MD)</h4>{a}</figure>'
            f'<figure class="tsfig"><h4>Fv native contacts over time (thermal MD, reference only)</h4>{b}</figure>')


if __name__ == "__main__":
    main()
