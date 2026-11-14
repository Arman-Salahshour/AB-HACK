# Run log — Oslo Antibodies × AI, 14 November 2026

> **What is new here.** Most design write-ups report the batch that worked. This one
> reports the batch that did not, and when. The de novo track was killed at 12:45 by
> three independent methods agreeing on 0/40, and the three batches already queued
> behind it were cancelled. The timeline below is kept so that the *cost of finding out
> late* is visible: we queued 400 designs before the first validation returned. That
> ordering mistake, not the null result, is the transferable lesson.

What was built, what was run, what it produced, and why each decision was made.
Times are UTC on the hackathon VM (1× H100, 16 vCPU, 196 GB).

---

## 1. The pipeline

```
 INPUT                GENERATE              VALIDATE                 SUBMIT
 ─────                ────────              ────────                 ──────
 target .cif   ─┐
 scaffold .cif ─┼──►  RFD4-Proteina   ──►   Boltz-2 refold    ──►    abguard
 design spec   ─┘     module-run design     (sequence only)          rank + export
 (C_HOT epitope)      ~18 s / design        ~14 s / design           ms
                           │                     │
                           │                     ├─ abhack.epitope   (contacts)
                           │                     ├─ FACS classifier  (P(high))
                           │                     ├─ SPR regressor    (pred K_D)
                           │                     └─ OpenDDE-abag     (JapanFold)
                           ▼                     ▼
                    self-assessed          independent evidence
                    (NOT evidence)
```

One command per stage:

```bash
python -m abhack.pipeline generate --spec <spec> --seed N --refold
python -m abhack.pipeline validate --run runs/<name>       # CPU, runs during generation
python -m abhack.pipeline finalise --runs ... --limit 20
python -m abhack.pipeline status
```

Two agents in `.claude/agents/`: **ab-generator** proposes and is forbidden from
claiming its output binds; **ab-evaluator** judges with evidence the generator cannot
produce and is instructed not to soften a null result.

### Why the separation exists

The generator's own metrics (`hotspot_obedience`, `clash_metrics`, `pre_oracle_pass`)
describe the geometry of a proposal, not whether it binds. Measured: **100 de novo
designs scored 100/100 on hotspot obedience; 0 of 20 survived independent refolding.**

---

## 2. Timeline

| Time | Event |
| --- | --- |
| 10:00 | Repo, pipeline, env, data fetched. Smoke test: design → refold → score → abguard, 3/3 valid |
| 10:40 | Probed `abguard` directly; found the novelty gate is stricter than the PDF implies |
| 11:20 | De novo VHH track built (HER2 + caplacizumab framework + hotspots), seeds 1–3 |
| 12:10 | **Seed 1 refold: 0/20.** Seed 2: 0/20. Binder pLDDT ~0.92 but 15–51 Å from the designed pose |
| 12:30 | AF2-IG on JapanFold confirms: ipTM 0.11–0.34 vs 0.5 threshold. BoltzGen's own de novo: 0.16–0.29 |
| 12:45 | **Dropped de novo. Moved GPU to the optimisation track** (seeds 11, 12) |
| 13:20 | Seed 11: 61/64 abguard-valid, 13/30 strong refold. OpenDDE-abag ipTM 0.63–0.96 |
| 13:35 | User supplied the visualised epitope (17 HER2 residues) and 7MN8 |
| 14:10 | **Seed 21, steered with the 17 residues: 66 valid, 30/30 strong refold** |
| 14:30 | Consensus run (22 residues) generating |

---

## 3. Runs

| Run | Weight | Spec | Designs | abguard-valid | Refolded | Strong |
| --- | --- | --- | --- | --- | --- | --- |
| `denovo-vhh-s1` | antibody | de novo VHH | 100 | – | 20 | **0** |
| `denovo-vhh-s2` | antibody | de novo VHH | 100 | – | 20 | **0** |
| `denovo-vhh-s3` | antibody | de novo VHH | 27 | – | 0 | – |
| `opt-novel-s11` | her2-novel | organisers' 4-atom | 100 | 61 | 30 | 13 |
| `opt-novel-s12` | her2-novel | organisers' 4-atom | 100 | 64 | 0 | – |
| `opt-explore-s13` | her2 tsr-k 1 | organisers' 4-atom | 55 | – | 0 | – |
| **`opt-epitope-s21`** | her2-novel | **17-residue epitope** | 100 | **66** | 30 | **30** |
| `opt-consensus-s22` | her2-novel | 22-residue consensus | running | – | – | – |
| `obedience-test` | her2-novel | control experiment | 10 | – | 10 | 3 |
| `inverse-bench` | her2-novel | inverse folding | 20 | – | 0 | – |

`opt-novel-s12-gated` and `opt-explore-s13-gated` were run by the team, not by this
session; their parameters are in their own `run.log`.

---

## 4. Results

### The decisive comparison

| | de novo VHH | optimisation (4-atom) | **optimisation (steered)** |
| --- | --- | --- | --- |
| median `min_ipae` | 2.48 – 15+ | 0.439 | **0.316** |
| median `iptm` | 0.43 – 0.82 | 0.793 | **0.962** |
| median backbone RMSD | 15 – 51 Å | 1.17 Å | **0.85 Å** |
| strong refolds | **0 / 40** | 13 / 30 | **30 / 30** |
| AF2-IG ipTM | 0.11 – 0.34 | – | – |
| OpenDDE-abag ipTM | – | 0.63 – 0.96 | – |

### Best candidates

| Design | CDRH3 | epitope | prec | min_ipae | iptm | RMSD | P(high) | pred K_D |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `s21_94` | WASSYYYYFD | 20 | 1.00 | 0.314 | 0.967 | 0.70 Å | 0.45 | **5.9 nM** |
| `s21_85` | WASLYEYYFD | 20 | 1.00 | 0.301 | 0.971 | 0.99 Å | 0.58 | 68.9 nM |
| `s21_51` | WSFSVDYAFD | 20 | 1.00 | 0.337 | 0.956 | 0.84 Å | – | – |
| `s11_53` | WTSSVGYAHD | 15/17 | 0.88 | 0.475 | 0.777 | 0.82 Å | **0.78** | 22.2 nM |

`WASSYYYYFD` appeared independently in seed 11 (`_50`) and seed 21 (`_94`) under
different steering — the strongest reproducibility signal in the dataset.

### Pool

```
482 designs generated   →   191 abguard-eligible   →   47 refolded   →   43 strong
```

---

## 5. Decisions and their evidence

**Dropped the de novo track (12:45).** 0/40 across three independent methods (Boltz-2
twice, AF2-IG, plus BoltzGen's internal scoring of its own designs). Adaptyv's published
EGFR competition agrees on direction: optimised binders 18–28 % hit rate, best de novo
method 9 %, RFdiffusion one binder in total. Their recommendation is explicit: start
from an existing binder and redesign its non-interface regions.

**Steered with the measured epitope (13:40).** The organisers' spec used 4 atoms in one
corner of the interface. The team's visualisation gave 17 residues; 7MN8 independently
confirmed 14 of them and added 5 more. Result: strong refolds went from 13/30 to 30/30
and median `iptm` from 0.79 to 0.96.

**Ranked structure-reproduction above interface confidence.** Two seed-11 designs
(`_64`, `_10`) had the best `min_ipae` but 9–14 Å RMSD — Boltz-2 confident about an
interface that is not the designed one. The epitope data, which knows nothing about
RMSD, independently ranked those two last (8/17 and 5/17 coverage). Two unrelated lines
of evidence, same conclusion.

**Did not pursue millisecond decoding.** Inverse folding (`--inverse-folding true`)
measured 21.1 s/sequence versus 18.0 s for full generation — RFD4 runs its full trunk
either way, so there is no lightweight decoder to isolate. Millisecond sequence design
needs a different model (ProteinMPNN, 1.7M params, ~20 ms). Inverse folding remains
useful for refining a validated backbone: 20 sequences from `_53` gave 12 unique loops,
all 6 mutations from trastuzumab, 9 passing novelty.

---

## 6. Open limitations

- **No design is predicted to beat trastuzumab on affinity.** Best predicted 5.9 nM
  against a measured 0.029–0.178 nM. Verified this is not model clipping: the regressor
  predicts as low as 0.05 nM on its own training data, and 25 of its 190 training points
  are sub-nanomolar. S-tier is unlikely through CDRH3 redesign alone.
- **`pass` is empty, not `True`.** `refold_metrics.csv` could not compute target-aligned
  RMSD ("fewer than 3 target atoms fixed in absolute coordinates"), so the scorer
  declined to judge. "Strong refold" is our threshold on the metrics the docs recommend,
  not the tool's verdict.
- **Thermostability is unmeasured.** No data in the provided tables and no predictor
  built. The free lever is that the novelty gate only inspects CDRH3, so framework
  mutations cost nothing against it — identified but not executed.
- **Only domain IV was modelled.** 109 of HER2's 607 residues. 7MN8 shows the 56 HER2
  residues contacting HER3 do not overlap the trastuzumab epitope, so dimerisation
  blocking is not a concern, but no design was folded in full-receptor context.
- **The obedience control was underpowered.** 30 % strong at obedience 0.75 vs 43 % at
  1.00, n=10 vs n=30 — inconclusive.

---

## 7. Failures hit, and the fixes

| Symptom | Cause | Fix |
| --- | --- | --- |
| `Invalid RFD4-Proteina design input` | PDB has no `label_seq_id`/entity records | `scripts/prepare_structure.py` |
| `Failed to compute leaving groups for '12-21'` | scaffold CDRs encoded as residues *named* for their length | strip non-standard residues |
| Same error, with a spec | a `"_comment"` key — every top-level key is read as a design id | no comments in the JSON |
| `Provide at least one finite hotspot seed, all on the same target chain` | `--resample-hotspots` with hotspots on chains C and D | drop the flag |
| `NaN values detected in input or generated coords` | inverse folding driven by a spec | pass a structure directly |
| 0 epitope contacts on every design | generated structures renumber chains from 1 | shift to reference numbering |
| novelty table all null | CDRH3 found by framework anchors a VHH lacks | number with `immunum` |
| refold metrics would not join | `refold_metrics.csv` keys on `query_id`, not `design_id` | join on `query_id` |
| `Input is not a design output` | hand-built run subset missing `result.json` | copy and filter the manifest |
| designs rejected that our filter passed | abguard uses a 610k reference incl. patents, aligns with indels | abguard is the authority |

---

## 8. Cost

| | |
| --- | --- |
| GPU | ~4 h on one H100 |
| JapanFold | **$0.20** of $1,100 (AF2-IG, BoltzGen, OpenDDE-abag) |
| Designs generated | 482 |

---

## 9. Where everything is

```
designs/best/                    the single best candidate, designed + refolded PDB
designs/structures/              seed 11 top 17, two PDBs each + results.csv
designs/structures-s21/          seed 21 top 12, two PDBs each
designs/her2-candidates.fasta    sequences with metrics in the headers
designs/opt-*.csv                per-run tables: metrics, novelty, epitope coverage
specs/                           every design spec used, version controlled
reports/her2-candidates.html     shareable candidate page
docs/agents/                     generator and evaluator specifications
.claude/agents/                  the runnable agents
src/abhack/                      pipeline, novelty, epitope, models, submission
```

Raw runs (`runs/`) and trained models (`models/`) are gitignored and die with the VM.
Everything above is in the repository.
