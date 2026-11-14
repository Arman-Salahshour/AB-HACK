# ab-hack — trastuzumab × HER2

Oslo Antibodies × AI hackathon, 14 November 2026. An agent-driven pipeline that designs
antibody variants against HER2, validates them with models that never saw the designs,
and ranks them on the only evidence that survives scrutiny.

**Submission:** [`odyssey.csv`](odyssey.csv) — 30 designs, 30/30 `abguard`-valid,
22 optimisation (18 parents + 4 framework twins) + 8 de novo, ranked. Rank 1 is the one to
make first. **Every number is a prediction; nothing here is a measured result.**
Results page: [`docs/final-ranking.html`](docs/final-ranking.html).

---

## The three findings

**1. Structure prediction cannot rank these designs.** On the organisers' own
cognate/non-cognate dataset (91 real pairs, ~8,100 shuffled), ipTM works as a binder
classifier at **AUROC 0.754 (AlphaFold3), 0.645 (Boltz-2), 0.568 (Chai-1)**. On our early 20
candidates the three predictors correlate **−0.15 to +0.13** with each other. Meanwhile a
classifier trained on the 524,346-variant FACS screen reaches **AUROC 0.989** on held-out
data, and the SPR regressor **0.930** at identifying sub-nanomolar binders. So structure
became a pass/fail gate and the measurement-trained models do the ranking.

**2. A generator cannot grade its own work.** 100 de novo designs scored **100/100** on
RFD4's own hotspot-obedience metric; **0 of 20** survived independent refolding — binder
pLDDT ~0.92 but 15–51 Å from the designed pose. AF2-IG agreed (ipTM 0.11–0.34 against a
0.5 threshold), and BoltzGen's own de novo pipeline agreed too. The de novo track was
dropped at 12:45 and the GPU moved to optimisation. This is why the pipeline has two
agents with opposed incentives.

**3. Using the measured epitope changes the designs.** Replacing the organisers' 4
hotspot atoms with the team's 17-residue mapped epitope, everything else identical:

| | 4 atoms | **17 residues** | 22-residue consensus |
| --- | --- | --- | --- |
| strong refolds | 13/30 | **30/30** | 24/30 |
| median `iptm` | 0.793 | **0.962** | 0.954 |
| median RMSD | 1.17 Å | **0.85 Å** | 0.94 Å |

---

## Results

```
225 unique optimisation designs → 139 pass novelty → 130 pass the pH gate
→ 18 with a Boltz-2 refold (parents) → +4 framework twins → 22
+ 8 de novo VHH, placed last                                   → 30 submitted
```

**Top of the final ranking** (full table: `designs/final-30.csv`; per-design reasoning is in
each `design_description` in `odyssey.csv`):

| rank | CDRH3 | P(high) | pred K_D | Boltz-2 iPAE | Mol-1 pose | pH sel. |
| --- | --- | --- | --- | --- | --- | --- |
| 1 (+twin 2) | `WPSGGTYGYD` | 0.87 | 10 nM | 0.32 | 1.5 Å | +0.03 |
| 3 (+twin 4) | `WYGRGWYVPD` | 0.88 | **9 nM** | 0.30 | 1.7 Å | +0.33 |
| 5 (+twin 6) | `YMSLYYYAFD` | 0.65 | 12 nM | 0.32 | 2.3 Å | −0.02 |
| 7 (+twin 8) | `TTSSVGYAHD` | 0.84 | 107 nM | 0.31 | 2.1 Å | −0.09 |
| — | *trastuzumab* | 0.995 | **0.22 nM** | — | — | — |

- **All 22 optimisation designs** sit in trastuzumab's binding pose when co-folded by
  Mol-1 with no template.
- **All 8 de novo designs** fail the credibility bar (AF2-IG ipTM 0.14–0.46 against ≥ 0.6),
  so they are included last as lottery tickets.

**The honest headline: no design is predicted to beat trastuzumab on affinity.** Best
predicted K_D is ~9 nM against the parent's 0.22 nM, roughly 40× weaker. The models score
the parent correctly (0.995, 0.22 nM), so this is a real gap, not a calibration artefact.
The realistic target is A-tier: epitope preserved, plus a thermostability and pH angle.

---

## The pipeline

```
 structures ─┐
 epitope     ├─► GENERATE ─► gate ─► pH gate ─► refold ─► cross-check ─► rank ─► twins ─► SUBMIT
 scaffold   ─┘   RFD4        │        │          Boltz-2    Mol-1 (no template)  │       S85N+R87T
                 BoltzGen    │        │          iPAE       AF2-IG (de novo)     │
                             │        └─ propka pH 6.5 vs 7.4 selectivity ≥ 0   └─ FACS classifier  AUROC 0.989
                             └─ abguard hard filters + novelty                       SPR regressor    AUROC 0.930
                                epitope contact ≥16 at ≥0.9 precision                liabilities
```

1. **Generate** — new CDRH3 loops on trastuzumab's framework (RFD4-Proteina, 17-residue
   mapped epitope), plus de novo VHH nanobodies (RFD4, BoltzGen).
2. **Gate** — `abguard` format and novelty rules; epitope engagement; interface PAE.
3. **pH gate** — tumour microenvironment (~pH 6.5) vs plasma (7.4): propka-based
   selectivity ≥ 0 keeps designs that do not lose binding when acidic. A prediction only.
4. **Refold and cross-check** — Boltz-2 refold of each design; Mol-1 with no template to
   confirm the binding pose; AF2-IG to judge de novo designs from their own coordinates.
5. **Rank** — binding models decide the order (FACS classifier, SPR regressor, interface
   PAE). pH and liabilities break ties; they never override.
6. **Twins** — each of the top 4 parents is followed by a framework twin with VH
   S85N + R87T, the one literature-backed route to higher Tm (Ng et al., mAbs 2025).
   Caveat: it creates an NLT glycosylation sequon (not an issue for cell-free expression).

```bash
uv sync && source env.sh
bash scripts/fetch_data.sh small

python -m abhack.pipeline generate --spec specs/opt-epitope/design.json --seed 1 --refold --compile true
python -m abhack.pipeline validate --run runs/opt-epitope-s1     # CPU, runs during generation
python scripts/best50.py -n 50          # funnel: gates → pH gate → composite → top 50
python scripts/final30.py               # final 30 + framework twins → designs/final-30.csv
python scripts/make_odyssey_csv.py      # hand-in CSV, cross-checked against the results page
bash scripts/validate.sh designs/odyssey.csv
python scripts/final_page_deck.py       # docs/final-ranking.html
```

**Always pass `--compile true` to generation:** measured **3.3×** (18.0 → 5.5 s/design).
The early 682 designs ran without it and cost ~2.5 h of GPU where ~45 min would have done.

### Two agents, opposed on purpose

`.claude/agents/` defines them; `docs/agents/` explains them.

- **`ab-generator`** proposes designs and is **forbidden from claiming they bind**.
- **`ab-evaluator`** judges with evidence the generator cannot produce, and is told **not
  to soften a null result**.

Reporting 0/40 at 12:45 freed three hours of GPU already committed to a configuration
that does not work.

---

## Documents

| | |
| --- | --- |
| [`docs/final-ranking.html`](docs/final-ranking.html) | the final ranking page — pipeline, 30 designs, simulations, next steps |
| [`docs/SUBMISSION.md`](docs/SUBMISSION.md) | earlier 20-design table with every metric (superseded by the 30) |
| [`docs/DEVELOPABILITY.md`](docs/DEVELOPABILITY.md) | liabilities, biophysics, toxicity, thermostability, head-to-head |
| [`docs/EFFICIENCY.md`](docs/EFFICIENCY.md) | why the filters ran in the wrong order, and the 3.3× |
| [`docs/TOOLS.md`](docs/TOOLS.md) | every tool used — and the seven that could not be run |
| [`docs/RUN-LOG.md`](docs/RUN-LOG.md) | every run, decision and failure, with times |
| [`docs/agents/`](docs/agents/) | generator and evaluator specifications |

---

## What is not measured

Stated plainly because it bounds every claim above:

- **Thermostability** — a judged axis. No Tm predictor was reachable; we report structural
  correlates and the literature-backed S85N+R87T twins. Thermal MD ran only on the
  trastuzumab reference, not on our designs.
- **pH selectivity** — predicted (propka). Two-pH MD was stopped early and covers the
  reference plus two designs; it checks the method and sets a noise level, it ranks nothing.
- **No energy estimate** — MM/GBSA on MD frames is not in the ranking yet.
- **Toxicity** — not predicted. We report the parent's *measured* FDA-label safety
  (anti-drug antibodies 1/903, CHF 0.4–3.2 %) and absence-of-red-flag checks.
- **Immunogenicity** — deliberately not predicted; tools in this class have limited
  validated power.
- **Solubility / aggregation** — none of Aggrescan3D, SoluProt, NetSolP or ProteinMPNN
  could be run here. See `docs/TOOLS.md`.
- **Affinity** — predicted, not measured. Calibrated against the parent, extrapolating
  beyond the training distribution by construction, since the novelty gate forces designs
  ≥3 mutations from every sequence the models learned from.

---

## Layout

```
odyssey.csv                   the submission
designs/                      tables, sequences, structures (PDB pairs per design)
specs/                        every design spec used, version controlled
scripts/                      design, refold, validate, structure prep, JapanFold clients
src/abhack/                   pipeline, novelty, epitope, models, submission
docs/                         the write-ups above
.claude/agents/               the runnable agents
```

`data/`, `runs/` and `models/` are gitignored — re-fetch with
`bash scripts/fetch_data.sh small`.

Data: Chinery et al. (FACS, SPR), Smorodina et al. (cognate/non-cognate ipTM),
PDB 1N8Z and 7MN8. Compute: 1× H100 (~5 h) and $0.44 of JapanFold credit.
