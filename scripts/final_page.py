#!/usr/bin/env python3
"""Build the final ranking page (wrapper-less HTML for the Artifact tool) from designs/final-30.csv.

    python scripts/final_page.py            # writes docs/final-ranking.html
"""

from __future__ import annotations

import html
import re
from pathlib import Path

import polars as pl

OUT = Path("docs/final-ranking.html")
SRC_CSS = Path("/tmp/claude-1001/-home-hackathon/bfddba5e-cb2f-4f81-a0cb-929d3a8291b0/scratchpad/pipeline-share.html")


def css() -> str:
    t = SRC_CSS.read_text()
    m = re.search(r"<style>(.*?)</style>", t, flags=re.S)
    return m.group(1)


def f(v, spec=".2f", dash="–") -> str:
    return dash if v is None or (isinstance(v, float) and v != v) else format(v, spec)


def frag(path: str, start: str, end: str | None) -> str:
    t = Path(path).read_text()
    a = t.index(start) + len(start)
    return t[a:t.index(end, a)] if end else t[a:]


EXTRA_CSS = """
  .lede strong { color: var(--fg); font-weight: 600; }
  .policy { margin: 0; padding-left: 20px; }
  .rt td, .rt th { vertical-align: top; }
  .rt td.why { min-width: 340px; max-width: 460px; font-size: 13px; color: var(--muted); }
  .rt td.id { white-space: nowrap; }
  .rt tr.tw td:first-child { box-shadow: inset 3px 0 0 var(--s2); }
  .rt tr.dn td { color: var(--muted); }
  .rt .rk { font-weight: 600; font-variant-numeric: tabular-nums; }
  .rt .tag { font: 500 11px var(--font-mono); border: 1px solid currentColor; border-radius: 10px; padding: 0 6px; white-space: nowrap; }
  .t-opt { color: var(--accent); } .t-twin { color: var(--s2); } .t-dn { color: var(--muted); }
  .evid td:first-child { white-space: nowrap; font-weight: 500; }
  .nums { font-family: var(--font-mono); font-size: 12.5px; font-variant-numeric: tabular-nums; }
"""


def table(df: pl.DataFrame) -> str:
    head = ("<tr><th>#</th><th>Design</th><th>CDRH3</th><th>Why it is here</th>"
            "<th class='m'>Binder<br>prob.</th><th class='m'>Pred. K<sub>D</sub><br>nM</th><th class='m'>Boltz-2<br>iPAE &Aring;</th>"
            "<th class='m'>Mol-1 pose<br>&Aring;</th><th class='m'>pH<br>select.</th><th class='m'>Loop<br>gap</th></tr>")
    body = []
    for r in df.to_dicts():
        twin = r["design_id"].endswith("-ft")
        dn = r["track"] == "denovo"
        tag = ('<span class="tag t-dn">de novo</span>' if dn else
               ('<span class="tag t-twin">twin</span>' if twin else '<span class="tag t-opt">optimisation</span>'))
        cls = "tw" if twin else ("dn" if dn else "")
        sel = r["ph_sel"]
        body.append(
            f'<tr class="{cls}"><td class="rk">{r["rank"]}</td>'
            f'<td class="id"><span class="m">{html.escape(r["design_id"])}</span><br>{tag}</td>'
            f'<td class="m">{html.escape(r["cdrh3"] or "")}</td>'
            f'<td class="why">{html.escape(r["why"])}</td>'
            f'<td class="nums">{f(r["p_high"])}</td><td class="nums">{f(r["kd_pred_nM"], ".0f")}</td>'
            f'<td class="nums">{f(r["ipae_boltz2"])}</td><td class="nums">{f(r["mol1_pose_A"], ".1f")}</td>'
            f'<td class="nums">{f(sel, "+.2f")}</td><td class="nums">{f(r["loop_gap"], ".1f")}</td></tr>')
    return f'<div class="scroll"><table class="rt"><thead>{head}</thead><tbody>{"".join(body)}</tbody></table></div>'


import importlib.util

_s = importlib.util.spec_from_file_location("illustrations", Path(__file__).with_name("illustrations.py"))
_ill = importlib.util.module_from_spec(_s)
_s.loader.exec_module(_ill)


def main() -> None:
    df = pl.read_csv("designs/final-30.csv", infer_schema_length=10_000)
    n_opt = df.filter(pl.col("track") == "optimization").height
    n_par = df.filter((pl.col("track") == "optimization") & ~pl.col("design_id").str.ends_with("-ft")).height
    n_dn = df.filter(pl.col("track") == "denovo").height
    mol1_n = df.filter(pl.col("mol1_pose_A").is_not_null()).height
    cards = frag("docs/results-radar.html", "<!--cards-->", "<!--table-->")
    ts = frag("docs/md-timeseries.html", "", None) if False else Path("docs/md-timeseries.html").read_text()

    ill_css, ill_html = _ill.build()

    page = f"""<meta name="claude:viewport-fit-sides" content="cover">
<title>HER2 Final Ranking</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>{css()}{EXTRA_CSS}{ill_css}</style>
<main>
  <header>
    <div class="eyebrow">Oslo Antibodies &times; AI &middot; 14 Nov 2026</div>
    <h1>HER2 final ranking</h1>
    <p class="lede">
      <strong>{df.height} designs, ranked.</strong> Rank 1 is the one to make first. {n_opt} keep trastuzumab's framework
      with a new CDRH3 ({n_par} designs plus 4 framework twins); {n_dn} are de novo nanobodies, placed last because
      none passed an independent structure check. Nothing here is a measured result: affinity, thermostability and
      pH selectivity are predictions.
    </p>
  </header>

  <section>
    {ill_html}
  </section>

  <section>
    <h2>How the ranking was made</h2>
    <ol class="policy">
      <li><strong>Hard filters.</strong> The organisers' validator accepts all {df.height} (valid sequences, novel CDRH3, no unpaired
        cysteine). Of 225 generated optimisation designs, 139 were novel, 130 survived the pH gate, and 18 of the best 50 had a
        structure refold and went forward.</li>
      <li><strong>Binding evidence ranks first.</strong> The binder model, trained on 525,000 measured variants, and the
        Boltz-2 interface error decide the order. ipTM is not used: it barely separates binders from non-binders here.</li>
      <li><strong>Thermostability.</strong> The top 4 designs are each followed by a twin carrying two published
        framework mutations (VH S85N + R87T) reported to raise melting temperature while keeping HER2 binding. It is the
        only literature-backed route to beating trastuzumab on stability.</li>
      <li><strong>pH and safety break ties, they do not override.</strong> Designs that clearly bind better at plasma pH 7.4
        than at tumour pH 6.5 were removed. A sequence-liability count is a small penalty.</li>
      <li><strong>De novo last.</strong> Independent checks found no credible binder among them, so they fill the end of the
        list and are made only if extra allocation is won.</li>
    </ol>
  </section>

  <section>
    <h2>The ranked designs</h2>
    {table(df)}
    <p style="margin-top:12px">Binder probability and predicted K<sub>D</sub> come from models trained on measured data
      (trastuzumab itself scores 0.995 and about 0.2 nM). Mol-1 is an AlphaFold 3-family model run with no template; pose is how far
      its Fv sits from trastuzumab's crystal pose after aligning HER2 (shown for the {mol1_n} designs checked).
      pH select. is the predicted log<sub>10</sub> of binding at 6.5 over 7.4 (about &plusmn;0.1&ndash;0.3 noise).
      Loop gap is how much less confident the refold is about the new CDRH3 than the framework (lower is better).</p>
  </section>

  <section>
    <h2>Top six at a glance</h2>
    <p>Four axes, 0&ndash;1, higher is better. Thermostability here is the loop-confidence proxy, and safety is a
      sequence proxy (liabilities, hydrophobicity, charge), not a toxicity prediction.</p>
    <div class="rcards">{cards}</div>
  </section>

  <section>
    <h2>What is measured and what is not</h2>
    <div class="scroll"><table class="evid">
      <thead><tr><th>Question</th><th>What we have</th><th>How far to trust it</th></tr></thead>
      <tbody>
        <tr><td>Affinity</td><td>Binder and K<sub>D</sub> models trained on 525k measured variants; Boltz-2 interface error</td>
          <td>Predictions only. Our best design is predicted about 27&times; weaker than trastuzumab, so we do not expect to win on affinity.</td></tr>
        <tr><td>Thermostability</td><td>Framework twins from the literature; a loop-confidence proxy; thermal MD on trastuzumab only</td>
          <td>No design was simulated, because the MD was stopped early. The twins are the real bet.</td></tr>
        <tr><td>pH 6.5 versus 7.4</td><td>A pKa-based prediction for every design; two-pH MD on the reference</td>
          <td>Unvalidated and noisy; the MD gap on the reference is within replica noise.</td></tr>
        <tr><td>Safety</td><td>Sequence liabilities, hydrophobicity and charge</td>
          <td>A proxy. Trastuzumab's cardiac toxicity is target-mediated, so every design on this epitope inherits the risk.</td></tr>
        <tr><td>Structure</td><td>Boltz-2 (with the trastuzumab template); Mol-1 and OpenFold3 (no template)</td>
          <td>Mol-1 reproduces trastuzumab's pose and scores the designs the same, so it supports the binding mode but cannot tell good loops
            from bad. Hosted OpenFold3 failed on trastuzumab itself, so its result was discarded.</td></tr>
      </tbody>
    </table></div>
  </section>

  <section>
    <h2>The simulations</h2>
    <p>Each line is one replica, drawn as far as it got (10 ps per frame). 1 means every starting contact is still present. These runs were
      stopped early and cover the trastuzumab reference, not our designs, so they check the method and set a noise level; they do not rank anything.</p>
    <div class="tsgrid">{ts}</div>
  </section>

  <p class="foot">Hand-in file: <code>final-submission.csv</code> (30 rows, with the <code>track</code> column). Rankings come from
    <code>scripts/final30.py</code>; all design tables and trajectories are in the repo. Antibody schematics adapted from the
    <a href="https://oslo-ai-antibodies-hackathon.vercel.app/participant_guide">Oslo AI &times; Antibodies participant guide</a>.</p>
</main>
"""
    OUT.write_text(page)
    print(f"wrote {OUT} ({len(page):,} bytes)")


if __name__ == "__main__":
    main()
