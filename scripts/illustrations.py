#!/usr/bin/env python3
"""Antibody illustrations for the final page: (css, html).

The stylised IgG and VHH geometry and the loop colours are ported from the Oslo AI x Antibodies participant guide
(https://oslo-ai-antibodies-hackathon.vercel.app/participant_guide), where they are drawn by its own script. They are
re-drawn here as static SVG with this project's annotations added (the redesigned CDRH3 loop, the framework-twin
mutation site). The histidine curve is drawn from the Henderson-Hasselbalch equation.
"""

from __future__ import annotations

import math

CSS = """
  .figs { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 22px; }
  .abfig { margin: 0; min-width: 0; display: flex; flex-direction: column; gap: 8px;
    --ab-accent: #6b4eff; --ab-surface2: color-mix(in srgb, var(--fg) 7%, var(--surface)); --ab-ring: var(--bg); }
  @media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) .abfig { --ab-accent: #9482ff; } }
  :root[data-theme="dark"] .abfig { --ab-accent: #9482ff; }
  .abfig svg { width: 100%; height: auto; max-height: 360px; display: block; overflow: visible; }
  .abfig .chain { fill: var(--ab-accent); }
  .abfig .chain.light { opacity: .5; }
  .abfig .vdom { fill: none; stroke: var(--ab-accent); stroke-width: 1.5; stroke-dasharray: 3 4; opacity: .7; }
  .abfig .ag { fill: var(--ab-surface2); stroke: var(--muted); stroke-width: 1.4; }
  .abfig .lead { stroke: var(--muted); stroke-width: 1; fill: none; }
  .abfig .tie { stroke: var(--muted); stroke-width: 1.4; stroke-dasharray: 2 4; }
  .abfig .loop { stroke: var(--ab-ring); stroke-width: 2.5; }
  .abfig .loop.dim { opacity: .3; }
  .abfig .l-h1 { fill: #f2711c; } .abfig .l-h2 { fill: #e0457b; } .abfig .l-h3 { fill: #ffc53d; }
  .abfig .l-l1 { fill: #16a05b; } .abfig .l-l2 { fill: #0b86b8; } .abfig .l-l3 { fill: #a56bff; }
  .abfig .strand path { stroke: var(--ab-ring); stroke-width: 10; stroke-linecap: round; opacity: .3; fill: none; }
  .abfig .ssbar { stroke: var(--ab-ring); stroke-width: 3; stroke-linecap: round; fill: none; }
  .abfig .num text { font: 700 12px var(--font-mono); fill: #1a1a1a; }
  .abfig .lbl { font: 500 11px var(--font-mono); fill: var(--muted); letter-spacing: .04em; }
  .abfig .lbl.hi { fill: var(--fg); font-weight: 600; }
  .abfig .ring { fill: none; stroke: var(--fg); stroke-width: 1.6; stroke-dasharray: 3 3; }
  .abfig .twin { fill: #e5600b; stroke: var(--ab-ring); stroke-width: 2; }
  .abfig figcaption { font-size: 13.5px; color: var(--muted); text-align: center; max-width: none; }
  .abfig figcaption b { display: block; color: var(--fg); font-weight: 600; font-size: 14px; }
  .hcurve .ax { stroke: var(--muted); stroke-width: 1.2; fill: none; }
  .hcurve .gr { stroke: var(--line); stroke-width: 1; }
  .hcurve .tk { font: 400 11px var(--font-mono); fill: var(--muted); }
  .hcurve .cv-free { stroke: var(--muted); stroke-width: 2.5; fill: none; }
  .hcurve .cv-shift { stroke: var(--ab-accent); stroke-width: 2.5; fill: none; }
  .hcurve .vline { stroke: var(--fg); stroke-width: 1.2; stroke-dasharray: 4 4; opacity: .6; }
  .hcurve .dl { font: 600 11px var(--font-mono); }
"""


def _arm(side: str, highlight_h3: bool, twin: bool) -> str:
    L = side == "L"
    tx, rot = (192, -32) if L else (208, 32)
    light_x = -50 if L else 22
    l = {"l1": -43 if L else 43, "l2": -20 if L else 20, "l3": -32 if L else 32}
    h1 = -9 if L else 9
    h2 = 10 if L else -10
    dim = " dim" if highlight_h3 else ""
    ring = ('<circle class="ring" cx="0" cy="-219" r="17"/>' if highlight_h3 else "")
    mark = ('<path class="twin" d="M0 -106l7 7-7 7-7-7z"/>' if (twin and not L) else "")
    return f"""<g transform="translate({tx} 246) rotate({rot})">
      <rect class="chain" x="-16" y="-200" width="32" height="205" rx="15"/>
      <rect class="chain light" x="{light_x}" y="-200" width="28" height="125" rx="13"/>
      <rect class="vdom" x="{-54 if L else -20}" y="-206" width="74" height="64" rx="22"/>
      <circle class="loop l-l1{dim}" cx="{l['l1']}" cy="-207" r="7"/>
      <circle class="loop l-l2{dim}" cx="{l['l2']}" cy="-207" r="7"/>
      <circle class="loop l-l3{dim}" cx="{l['l3']}" cy="-217" r="8"/>
      <circle class="loop l-h1{dim}" cx="{h1}" cy="-207" r="7"/>
      <circle class="loop l-h2{dim}" cx="{h2}" cy="-207" r="7"/>
      <circle class="loop l-h3" cx="0" cy="-219" r="10"/>
      {ring}{mark}
    </g>"""


def igg_ours() -> str:
    svg = f"""<svg viewBox="0 0 400 430" role="img" aria-label="Stylised IgG antibody, a Y shape with six binding loops at each tip. The heavy-chain CDR3 loop, which our optimisation designs redesign, is circled; a diamond marks the framework twin mutation site.">
    <g class="chain"><rect x="170" y="240" width="28" height="150" rx="14"/><rect x="202" y="240" width="28" height="150" rx="14"/></g>
    <path class="ag" d="M30 40c4-22 30-30 46-16 18 3 26 22 12 36-10 14-34 14-46 4-8-6-14-14-12-24z"/>
    <path class="ag" d="M370 40c-4-22-30-30-46-16-18 3-26 22-12 36 10 14 34 14 46 4 8-6 14-14 12-24z"/>
    <path class="tie" d="M62 62L70 76M338 62L330 76"/>
    {_arm('L', True, False)}{_arm('R', True, True)}
    <text class="lbl" x="48" y="14">antigen</text><text class="lbl" x="318" y="14">antigen</text>
    <text class="lbl" x="200" y="414" text-anchor="middle">Fc (the stem)</text>
    <text class="lbl hi" x="200" y="70" text-anchor="middle">CDRH3: the loop we redesign</text>
    <path class="lead" d="M268 66L304 62"/>
    <text class="lbl hi" x="205" y="118" text-anchor="middle">twin: S85N + R87T</text>
    <path class="lead" d="M240 124L257 154"/>
  </svg>"""
    return (f'<figure class="abfig">{svg}<figcaption><b>Optimisation designs</b>Trastuzumab with one new CDRH3 loop. '
            f'The four twins also carry two framework mutations.</figcaption></figure>')


def vhh() -> str:
    svg = """<svg viewBox="0 0 400 430" role="img" aria-label="Stylised nanobody, one compact domain with three binding loops touching an antigen">
    <ellipse class="ag" cx="200" cy="72" rx="150" ry="52"/>
    <text class="lbl" x="200" y="76" text-anchor="middle">antigen</text>
    <path class="tie" d="M150 186L168 122M250 186L232 122M200 172L200 124"/>
    <rect class="chain" x="125" y="200" width="150" height="190" rx="46"/>
    <g class="strand"><path d="M158 236V356M186 236V356M214 236V356M242 236V356"/></g>
    <path class="ssbar" d="M186 302H214"/>
    <circle class="loop l-h1" cx="150" cy="197" r="12"/>
    <circle class="loop l-h2" cx="250" cy="197" r="12"/>
    <circle class="loop l-h3" cx="200" cy="186" r="15"/>
    <g class="num"><text x="150" y="201" text-anchor="middle">1</text><text x="250" y="201" text-anchor="middle">2</text><text x="200" y="191" text-anchor="middle">3</text></g>
    <text class="lbl" x="200" y="414" text-anchor="middle">VHH &middot; ~120 amino acids &middot; ~15 kDa</text>
  </svg>"""
    return ('<figure class="abfig">' + svg + '<figcaption><b>De novo nanobodies</b>All three loops are designed from '
            'scratch. No design passed an independent structure check.</figcaption></figure>')


def his_curve() -> str:
    W, H, ml, mr, mt, mb = 400, 330, 46, 14, 20, 56
    x = lambda ph: ml + (W - ml - mr) * (ph - 5.0) / 3.0  # noqa: E731
    y = lambda f: mt + (H - mt - mb) * (1 - f)  # noqa: E731
    frac = lambda ph, pka: 1 / (1 + 10 ** (ph - pka))  # noqa: E731
    g = []
    for f in (0, 0.25, 0.5, 0.75, 1.0):
        g.append(f'<line class="gr" x1="{ml}" x2="{W - mr}" y1="{y(f):.1f}" y2="{y(f):.1f}"/>'
                 f'<text class="tk" x="{ml - 8}" y="{y(f) + 4:.1f}" text-anchor="end">{f:.2f}</text>')
    for ph in (5, 5.5, 6, 6.5, 7, 7.4, 8):
        g.append(f'<text class="tk" x="{x(ph):.1f}" y="{H - mb + 18}" text-anchor="middle">{ph:g}</text>')
    g.append(f'<text class="tk" x="{(ml + W - mr) / 2:.0f}" y="{H - 8}" text-anchor="middle">pH</text>')
    g.append(f'<text class="tk" transform="rotate(-90 12 {(mt + H - mb) / 2:.0f})" x="12" y="{(mt + H - mb) / 2:.0f}" text-anchor="middle">fraction protonated</text>')
    for ph, lab in ((6.5, "tumour 6.5"), (7.4, "plasma 7.4")):
        g.append(f'<line class="vline" x1="{x(ph):.1f}" x2="{x(ph):.1f}" y1="{mt}" y2="{H - mb}"/>')
        g.append(f'<text class="tk" x="{x(ph) + (4 if ph < 7 else -4):.1f}" y="{mt + 12}" text-anchor="{"start" if ph < 7 else "end"}" style="fill:var(--fg)">{lab}</text>')
    pts = [5 + 3 * i / 60 for i in range(61)]
    for cls, pka in (("cv-free", 6.0), ("cv-shift", 7.0)):
        d = " ".join(f"{'M' if i == 0 else 'L'}{x(p):.1f},{y(frac(p, pka)):.1f}" for i, p in enumerate(pts))
        g.append(f'<path class="{cls}" d="{d}"/>')
    for pka, col, dy in ((6.0, "var(--muted)", 16), (7.0, "var(--ab-accent)", -8)):
        for ph in (6.5, 7.4):
            g.append(f'<circle cx="{x(ph):.1f}" cy="{y(frac(ph, pka)):.1f}" r="4.5" fill="{col}" stroke="var(--bg)" stroke-width="2"/>')
    d_free = frac(6.5, 6.0) - frac(7.4, 6.0)
    d_shift = frac(6.5, 7.0) - frac(7.4, 7.0)
    g.append(f'<text class="dl" x="{x(5.05):.1f}" y="{y(0.40):.1f}" style="fill:var(--muted)">free histidine (pKa 6.0)</text>')
    g.append(f'<text class="dl" x="{x(5.05):.1f}" y="{y(0.40) + 14:.1f}" style="fill:var(--muted)">gap {d_free:.2f}</text>')
    g.append(f'<text class="dl" x="{x(5.05):.1f}" y="{y(0.16):.1f}" style="fill:var(--ab-accent)">histidine lifted to pKa 7.0</text>')
    g.append(f'<text class="dl" x="{x(5.05):.1f}" y="{y(0.16) + 14:.1f}" style="fill:var(--ab-accent)">gap {d_shift:.2f}</text>')
    svg = (f'<svg class="hcurve" viewBox="0 0 {W} {H}" role="img" aria-label="Fraction of a histidine that is protonated against pH. '
           f'A free histidine changes by {d_free:.2f} between pH 6.5 and 7.4; one lifted to pKa 7.0 by a nearby acid changes by {d_shift:.2f}.">'
           f'{"".join(g)}</svg>')
    return ('<figure class="abfig">' + svg + '<figcaption><b>The pH idea</b>A histidine gains a charge as pH falls. '
            'A bigger gap between 6.5 and 7.4 needs its pKa pushed toward 7.</figcaption></figure>')


def build() -> tuple[str, str]:
    return CSS, f'<div class="figs">{igg_ours()}{vhh()}{his_curve()}</div>'


if __name__ == "__main__":
    css, html_ = build()
    print(len(css), len(html_))
    import re
    for tag in ("svg", "figure", "g"):
        print(tag, len(re.findall(f"<{tag}[ >]", html_)), len(re.findall(f"</{tag}>", html_)))
