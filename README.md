# ab-hack — trastuzumab × HER2

Oslo Antibodies × AI hackathon, 14 November 2026. An agent-driven pipeline that designs
antibody variants against HER2, validates them with models that never saw the designs,
and ranks them on the only evidence that survives scrutiny.

**Submission:** [`odyssey.csv`](odyssey.csv) — 30 designs, 30/30 `abguard`-valid,
24 optimisation + 6 de novo, ranked.

---

## The three findings

**1. Structure prediction cannot rank these designs.** On the organisers' own
cognate/non-cognate dataset (91 real pairs, ~8,100 shuffled), ipTM works as a binder
classifier at **AUROC 0.754 (AlphaFold3), 0.645 (Boltz-2), 0.568 (Chai-1)**. On our 20
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
682 designs generated     348 abguard-eligible     78 structurally validated
```

**Top candidates** (full table in [`docs/SUBMISSION.md`](docs/SUBMISSION.md)):

| CDRH3 | P(high) | pred K_D | epitope | RMSD |
| --- | --- | --- | --- | --- |
| `WVYAGSYAHD` | 0.759 | **4.9 nM** | 18 | 1.04 Å |
| `YASLYVYAHD` | 0.761 | 17.5 nM | 18 | 0.91 Å |
| `YASLHDYYFD` | **0.883** | 87.7 nM | 17 | 0.76 Å |
| *trastuzumab* | **0.995** | **0.22 nM** | — | — |

**The honest headline: no design is predicted to beat trastuzumab on affinity.** Best
predicted K_D is 4.9 nM against a measured 0.029–0.178 nM. The models are calibrated —
they score the parent at 0.995 and 0.22 nM — so this is a real gap, not a modelling
artefact. S-tier via CDRH3 redesign alone looks out of reach; A-tier (epitope preserved
plus thermostability) is the realistic target.

---

## The pipeline

```
 structures ─┐
 epitope     ├──► GENERATE ──► gate ──► rank ──► SUBMIT
 scaffold   ─┘    RFD4          │        │
                                │        └─ FACS classifier  AUROC 0.989
                                │           SPR regressor    AUROC 0.930
                                └─ abguard (hard filters + novelty)
                                   Boltz-2 refold (pass/fail)
                                   epitope contact ≥16 at ≥0.9 precision
```

```bash
uv sync && source env.sh
bash scripts/fetch_data.sh small

python -m abhack.pipeline generate --spec specs/opt-epitope/design.json --seed 1 --refold
python -m abhack.pipeline validate --run runs/opt-epitope-s1     # CPU, runs during generation
python -m abhack.pipeline finalise --runs runs/... --limit 30
bash scripts/validate.sh designs/SUBMISSION.csv
```

**Always pass `--compile true` to generation:** measured **3.3×** (18.0 → 5.5 s/design).
We ran 682 designs without it and spent ~2.5 h of GPU where ~45 min would have done.

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
| [`docs/SUBMISSION.md`](docs/SUBMISSION.md) | the designs, every metric, full sequences |
| [`docs/DEVELOPABILITY.md`](docs/DEVELOPABILITY.md) | liabilities, biophysics, toxicity, thermostability, head-to-head |
| [`docs/EFFICIENCY.md`](docs/EFFICIENCY.md) | why the filters ran in the wrong order, and the 3.3× |
| [`docs/TOOLS.md`](docs/TOOLS.md) | every tool used — and the seven that could not be run |
| [`docs/RUN-LOG.md`](docs/RUN-LOG.md) | every run, decision and failure, with times |
| [`docs/agents/`](docs/agents/) | generator and evaluator specifications |

---

## What is not measured

Stated plainly because it bounds every claim above:

- **Thermostability** — a judged axis. No Tm predictor was reachable; we report structural
  correlates (VH–VL interface contacts, CDRH3 glycine content) and say so.
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
