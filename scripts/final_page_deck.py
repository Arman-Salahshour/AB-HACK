#!/usr/bin/env python3
"""Final ranking page in the look of the Odyssey presentation (wrapper-less HTML for the Artifact tool).

Structure and style follow the deck: navy illustrated 'night' slides with gold eyebrows and cream Cambria titles,
and cream content slides with navy text. Images are the deck's own, downscaled and embedded (the viewer only
allows data: images). Data tables come from designs/final-30.csv via final_page.py helpers.

    python scripts/final_page_deck.py        # writes docs/final-ranking.html
"""

from __future__ import annotations

import base64
import importlib.util
import re
from pathlib import Path

import polars as pl

HERE = Path(__file__).parent
WEB = Path("/tmp/claude-1001/-home-hackathon/bfddba5e-cb2f-4f81-a0cb-929d3a8291b0/scratchpad/deck/web")


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fp = load("final_page")
ill = load("illustrations")


def uri(name: str) -> str:
    mime = {"jpg": "image/jpeg", "webp": "image/webp", "png": "image/png"}[name.rsplit(".", 1)[1]]
    return f"data:{mime};base64," + base64.b64encode((WEB / name).read_bytes()).decode()


def base_css() -> str:
    c = fp.css()
    # the deck look is one committed palette: drop the base page's dark-mode token blocks
    i = c.index("@media (prefers-color-scheme: dark)")
    j = c.index(':root[data-theme="dark"]')
    k = c.index("}", j) + 1
    return c[:i] + c[k:]


ill_css = re.sub(r"\s*@media \(prefers-color-scheme: dark\)[^\n]*\n", "\n", ill.CSS)
ill_css = re.sub(r'\s*:root\[data-theme="dark"\] \.abfig \{[^}]*\}', "", ill_css)

DECK_CSS = """
  /* ---- Odyssey deck palette: navy 0E1F4A, deep blue 16306E, cream F4EFE4, light blue C4D6EE, orange F0872A, gold E9B44C ---- */
  :root {
    --bg: #f4efe4; --surface: #fbf8f1; --fg: #0e1f4a; --muted: #4a5879; --line: #d9cfba; --accent: #2463b5;
    --ok: #2f7a4d; --partial: #b5611a; --missing: #a33b2a;
    --s1: #2463b5; --s2: #f0872a; --s3: #1baf7a; --s4: #e9b44c; --ref: #4a5879;
    --font-serif: Cambria, Caladea, Georgia, "Times New Roman", serif;
    --font-sans: Calibri, Carlito, system-ui, -apple-system, "Segoe UI", sans-serif;
    --navy: #0e1f4a; --deep: #16306e; --cream: #f4efe4; --sky: #c4d6ee; --orange: #f0872a; --gold: #e9b44c; --brown: #b5611a;
    color-scheme: light;
  }
  html { color-scheme: light; }
  body { background: var(--cream); color: var(--navy); font-family: var(--font-sans); font-size: 17px; line-height: 1.55;
    padding: 0; }
  main { max-width: none; gap: 0; display: block; }
  .band { position: relative; overflow: hidden; }
  .wrap { max-width: 1040px; margin: 0 auto; padding-left: max(16px, env(safe-area-inset-left, 0px));
    padding-right: max(16px, env(safe-area-inset-right, 0px)); }
  .eyebrow { font: 700 12px var(--font-sans); letter-spacing: .32em; text-transform: uppercase; color: var(--gold); margin: 0 0 14px; }
  .band.paper .eyebrow { color: var(--brown); }
  .band h2.title { font: 400 clamp(30px, 5vw, 44px)/1.08 var(--font-serif); border: 0; padding: 0; margin: 0 0 18px; letter-spacing: 0;
    text-transform: none; color: inherit; text-wrap: balance; }
  .band h3.sub { font: 400 clamp(21px, 3vw, 26px)/1.2 var(--font-serif); margin: 38px 0 12px; color: inherit; }
  .band p { max-width: 62ch; }

  /* cover: the deck's Cover layout (storm) with the ship */
  .cover { background: #24498f url(@STORM@) center bottom / cover no-repeat; color: var(--cream); }
  .cover .wrap { min-height: clamp(470px, 82vh, 780px); display: flex; flex-direction: column; justify-content: center; position: relative; }
  .cover h1 { font: 400 clamp(64px, 13vw, 132px)/.95 var(--font-serif); margin: 0 0 20px; color: var(--cream); letter-spacing: -.01em; }
  .cover .sub { font: italic 400 clamp(19px, 2.6vw, 26px)/1.3 var(--font-serif); color: var(--sky); margin: 0 0 6px; }
  .cover .ship { position: absolute; right: max(16px, env(safe-area-inset-right, 0px)); bottom: 7%; width: clamp(170px, 31vw, 380px); height: auto; }
  @media (max-width: 640px) { .cover .wrap { justify-content: flex-start; padding-top: 56px; } .cover .ship { bottom: 4%; width: 46vw; } }

  /* night slides: the deck's Night layout (calm waves) */
  .night { background: #0e1f4d url(@CALM@) center bottom / 100% auto no-repeat; color: var(--cream);
    --fg: #f4efe4; --muted: #c4d6ee; --line: rgba(244,239,228,.24); --bg: #12285c; --surface: rgba(255,255,255,.07); }
  .night .wrap { padding-top: clamp(48px, 8vw, 88px); padding-bottom: calc(min(100vw, 1800px) * .15 + 30px); }
  .night p { color: var(--cream); }
  .night .say { font: italic 400 clamp(18px, 2.4vw, 22px)/1.4 var(--font-serif); color: var(--sky); max-width: 40ch; margin: 22px 0 0; }

  .split { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1.05fr); gap: clamp(20px, 4vw, 56px); align-items: center; }
  @media (max-width: 760px) { .split { grid-template-columns: minmax(0, 1fr); } }
  .cellheart { width: 100%; height: auto; display: block; }
  .twolabels { display: flex; justify-content: space-between; gap: 16px; margin-top: 6px; }
  .twolabels div { max-width: 48%; }
  .twolabels b { display: block; font: 700 15px var(--font-sans); color: var(--cream); }
  .twolabels i { font: italic 400 14px var(--font-sans); color: var(--gold); }

  /* paper slides: cream content slides */
  .paper { background: var(--cream); color: var(--navy); }
  .paper .wrap { padding-top: clamp(44px, 7vw, 80px); padding-bottom: clamp(44px, 7vw, 84px); }
  .paper + .paper .wrap { padding-top: 0; }
  .route { list-style: none; margin: 8px 0 0; padding: 0; display: grid; gap: 14px; counter-reset: stop; }
  .route li { counter-increment: stop; margin: 0; max-width: none; display: grid; grid-template-columns: 44px minmax(0, 1fr); gap: 4px 16px;
    padding: 14px 0; border-top: 1px solid var(--line); }
  .route li::before { content: counter(stop); grid-row: span 2; font: 400 34px/1 var(--font-serif); color: var(--brown); }
  .route b { font: 700 17px var(--font-sans); }
  .route span { color: var(--muted); font-size: 15.5px; }
  .route .n { font: 700 15px var(--font-sans); color: var(--navy); }
  .band.paper .scroll { background: var(--surface); border-color: var(--line); }
  .band.paper th { color: var(--muted); }
  .band.paper td, .band.paper th { border-bottom-color: var(--line); }
  .rt td.why { color: #3b4968; }
  .lede2 { font: italic 400 clamp(18px, 2.3vw, 21px)/1.4 var(--font-serif); color: var(--muted); max-width: 52ch; margin: 0 0 6px; }
  .band.paper .rcard, .band.paper .tsfig { color: var(--navy); }

  /* crew */
  .crewfig { width: 100%; height: auto; display: block; }
  .names { list-style: none; margin: 18px 0 0; padding: 0; display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 6px 24px; }
  .names li { margin: 0; max-width: none; font: 400 17px var(--font-serif); color: var(--cream); }
  .next { margin: 6px 0 0; padding-left: 20px; }
  .next li { margin-bottom: 9px; max-width: 60ch; color: var(--cream); }
  .night .foot { color: var(--sky); border-top-color: rgba(244,239,228,.24); margin-top: 34px; font-size: 13px; }
  .night .foot a { color: var(--gold); }
  .night .abfig svg { max-height: 300px; }
  .night .abfig { --ab-accent: #ab9dff; --ab-ring: #12285c; }
"""


def main() -> None:
    df = pl.read_csv("designs/final-30.csv", infer_schema_length=10_000)
    n_opt = df.filter(pl.col("track") == "optimization").height
    n_par = df.filter((pl.col("track") == "optimization") & ~pl.col("design_id").str.ends_with("-ft")).height
    n_dn = df.filter(pl.col("track") == "denovo").height
    cards = fp.frag("docs/results-radar.html", "<!--cards-->", "<!--table-->")
    ts = Path("docs/md-timeseries.html").read_text()

    css = (base_css() + fp.EXTRA_CSS + ill_css + DECK_CSS).replace("@STORM@", uri("storm.jpg")).replace("@CALM@", uri("calm.jpg"))

    page = f"""<meta name="claude:viewport-fit-sides" content="cover">
<title>Odyssey Final Ranking</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Caladea:ital,wght@0,400;0,700;1,400&family=Carlito:ital,wght@0,400;0,700;1,400&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>{css}</style>
<main>

  <section class="band cover">
    <div class="wrap">
      <div class="eyebrow">Antibody AI Hackathon &middot; Oslo 2026</div>
      <h1>Odyssey</h1>
      <p class="sub">Final ranking: {df.height} HER2 antibody designs</p>
      <p class="sub" style="font-size:clamp(16px,2vw,19px)">Rank 1 is the one to make first.</p>
      <img class="ship" src="{uri('ship.webp')}" alt="Orange-sailed ship">
    </div>
  </section>

  <section class="band night">
    <div class="wrap">
      <div class="eyebrow">I &middot; Unmet need</div>
      <h2 class="title">Strike the tumor,<br>spare the heart</h2>
      <div class="split">
        <div>
          <figure style="margin:0"><img class="cellheart" src="{uri('cell_heart.webp')}" alt="A purple tumor cell covered in HER2 receptors with antibodies attached, next to a heart that also carries HER2 receptors">
            <div class="twolabels"><div><b>Tumor microenvironment</b><i>Target HER2 here</i></div><div><b>Heart</b><i>Spare HER2 here</i></div></div>
          </figure>
          <p class="say">We target HER2 in the tumor microenvironment, not the HER2 the heart relies on.</p>
        </div>
        <div>
          {ill.his_curve()}
          <p style="margin-top:12px;font-size:15.5px">The tumor microenvironment is more acidic (around pH 6.5) than blood (7.4). A histidine gains a charge as pH falls, so a design can in
            principle hold on harder in the tumor. We predict this for every design; it is a prediction, not a measurement.</p>
        </div>
      </div>
    </div>
  </section>

  <section class="band paper">
    <div class="wrap">
      <div class="eyebrow">II &middot; Our voyage</div>
      <h2 class="title">Pipeline</h2>
      <div class="split" style="align-items:start">
        <ol class="route">
          <li><b>Generate</b><span>225 unique optimisation designs (a new CDRH3 loop on trastuzumab's framework) and 220 de novo nanobodies.</span></li>
          <li><b>Filter</b><span>139 pass the novelty rules and 130 pass the pH gate; only designs with a structure refold go forward ({n_par} parents).</span></li>
          <li><b>Cross-check</b><span>Boltz-2 interface error, then Mol-1 with no template: all 22 optimisation designs sit in trastuzumab's binding pose.</span></li>
          <li><b>Rank</b><span>Binding models decide the order. pH and safety break ties; they never override.</span></li>
          <li><b>Twin</b><span>The top 4 are followed by a framework twin (VH S85N + R87T) that the literature links to higher melting temperature.</span></li>
        </ol>
        <div class="figs" style="grid-template-columns:repeat(2,minmax(0,1fr));gap:14px">{ill.igg_ours()}{ill.vhh()}</div>
      </div>
    </div>
  </section>

  <section class="band paper">
    <div class="wrap">
      <div class="eyebrow">III &middot; Landfall</div>
      <h2 class="title">Results</h2>
      <p class="lede2">{df.height} designs, ranked: {n_opt} optimisation ({n_par} parents and 4 twins) and {n_dn} de novo, placed last because no independent check found them credible.
        Nothing here is a measured result: affinity, thermostability and pH selectivity are predictions.</p>

      <h3 class="sub">The ranked designs</h3>
      {fp.table(df)}
      <p style="margin-top:12px;font-size:15px;color:var(--muted)">Binder probability and predicted K<sub>D</sub> come from models trained on measured data (trastuzumab itself scores 0.995 and about 0.2 nM).
        Mol-1 is an AlphaFold 3-family model run with no template; pose is how far its antibody sits from trastuzumab's crystal pose after aligning HER2.
        pH select. is the predicted log<sub>10</sub> of binding at 6.5 over 7.4 (about &plusmn;0.1&ndash;0.3 noise). Loop gap is how much less confident the
        refold is about the new CDRH3 than about the framework (lower is better). Twin rows reuse their parent's structure numbers.</p>

      <h3 class="sub">Top six at a glance</h3>
      <p style="margin:0 0 14px;color:var(--muted)">Four axes, 0&ndash;1, higher is better. Thermostability is the loop-confidence proxy and safety is a sequence proxy, not a toxicity prediction.</p>
      <div class="rcards">{cards}</div>

      <h3 class="sub">What is measured and what is not</h3>
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
            <td>Mol-1 reproduces trastuzumab's pose and scores the designs the same, so it supports the binding mode but cannot tell good loops from bad.
              Hosted OpenFold3 failed on trastuzumab itself, so its result was discarded.</td></tr>
        </tbody>
      </table></div>

      <h3 class="sub">The simulations</h3>
      <p style="margin:0 0 14px;color:var(--muted)">Each line is one replica, drawn as far as it got (10 ps per frame). 1 means every starting contact is still present. These runs were stopped early
        and cover the trastuzumab reference, not our designs, so they check the method and set a noise level; they do not rank anything.</p>
      <div class="tsgrid">{ts}</div>
    </div>
  </section>

  <section class="band night">
    <div class="wrap">
      <div class="eyebrow">IV &middot; The crew and the journey ahead</div>
      <h2 class="title">Next steps and who we are</h2>
      <div class="split">
        <div>
          <ol class="next">
            <li><b>Measure the top designs.</b> Binding and thermostability in the lab; every number here is still a prediction.</li>
            <li><b>Finish the simulations on our designs.</b> Thermal and two-pH MD for the top design and its twin never completed.</li>
            <li><b>Test the pH idea directly.</b> Compare binding at 6.5 and 7.4 for the designs predicted to prefer 6.5.</li>
            <li><b>Add an energy estimate</b> (MM/GBSA on the MD frames) to the ranking; there is none yet.</li>
          </ol>
          <ul class="names">
            <li>Amir Amanzadi</li><li>Arman Salahshour</li><li>Joel Benjamin Heim</li>
            <li>Johanne Kloster Ellingsen</li><li>Paula Istvan</li><li>Tamjidmaa Khatanbaatar</li>
          </ul>
        </div>
        <img class="crewfig" src="{uri('crew.webp')}" alt="The six crew members on the ship">
      </div>
      <p class="foot">Hand-in file: <code>odyssey.csv</code> (30 rows, with the <code>track</code> column). Rankings come from <code>scripts/final30.py</code>; design tables and
        trajectories are in the repo. Antibody schematics adapted from the
        <a href="https://oslo-ai-antibodies-hackathon.vercel.app/participant_guide">Oslo AI &times; Antibodies participant guide</a>.</p>
    </div>
  </section>

</main>
"""
    page = page.replace("tumours", "tumors").replace("Tumours", "Tumors").replace("tumour", "tumor").replace("Tumour", "Tumor")
    Path("docs/final-ranking.html").write_text(page)
    print(f"wrote docs/final-ranking.html ({len(page) / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
