# CLAUDE.md — working notes for this repo

Oslo Antibodies × AI hackathon. Goal: designs that beat trastuzumab on **affinity**
(SPR), **epitope preservation** (competitive binding) and **thermostability**.
Submission: ≤ 20 ranked designs, free split between the de novo and optimisation tracks.
Deadline: the VM is shut down at the end of **Saturday 14 November 2026** — copy out
anything worth keeping (`designs/`, `reports/`, `src/`, this repo) before then.

## Environment

- VM: 1× H100 80 GB, 16 vCPU, 196 GB RAM, ~358 GB free disk. No internet-facing services.
- Python: `uv sync` → `.venv`. Run things as `.venv/bin/python -m abhack.<module>`.
- Models run as **docker images**, not via `apheris-foundry` (there is no Foundry Hub here).
  `source env.sh` defines every image tag plus `abrun`, which mounts the repo root at `/work`.
- `abguard` (uv tool, on PATH) is the organisers' validator and the final authority on filters.
- `foundry-toolkit module list` lists the images; `foundry-toolkit module describe --name … --version …`
  shows every flag of one.

## Pipeline

```
scripts/design.sh    → runs/<name>/design       RFD4-Proteina design
scripts/refold.sh    → runs/<name>/prepared, /<fold>, /scored-<fold>
abhack.designs       → designs/<name>.csv       sequences + metrics + novelty verdict
abhack.evaluate      → ranked table             composite z-score
abhack.submission    → designs/submission.csv   abguard schema, ≤ 20 rows
scripts/validate.sh  → abguard verdict          run this before any submission
abhack.calibrate     → reports/calibration.md   what a confidence score is worth
```

Two things about this flow that cost us designs:

- **`refold.sh` selects before any gate.** It passes `--top-n`, so RFD4's own
  `pre_oracle_pass`/H-bond ordering decides what gets refolded — novelty, epitope and
  confidence are never consulted. In the smoke run that refolded the two *worse*-scoring
  designs and skipped the two that passed novelty most comfortably. **Gate first, then
  refold the survivors.** Refolding costs ~25 s/design/model + 20 s apo, so the order is
  free to fix and expensive to get wrong.
- **`abhack.evaluate` is orphaned.** Its `composite` never reaches `abhack.submission`,
  which does its own ranking. Two ranking paths exist; only `submission` writes the file.
  Don't tune `evaluate` expecting it to change a submission.

## Facts worth not re-deriving

- Trastuzumab CDRH3 window (what the data calls `cdrh3`): `WGGDGFYAMD`, 10 residues,
  between the fixed anchors `YYCSR` and `YWGQG` in VH. Sequences in `abhack.constants`.
- Baseline SPR KD of trastuzumab in the provided table: 0.029–0.178 nM (median 0.134, n=8).
  Beating that on affinity is hard; thermostability and epitope preservation are the
  cheaper axes to win.
- Provided data for the novelty gate = FACS ∪ SPR CDRH3s = **525,063** unique sequences.
  A design matching any of them fails. Known FACS binders therefore always fail.
- **How `abguard` really decides** (verified on a probe table, 2026-11-14 — see
  `src/abhack/novelty.py`):
  - it compares the **IMGT CDR3**, `SRWGGDGFYAMDY` for trastuzumab — the 10-residue
    window plus `SR`/`Y` flanks (`immunum.Annotator(chains=["H","K","L"], scheme="imgt")`);
  - `identity_to_trastuzumab` is over the **whole domain**, not the loop, so any CDRH3-only
    redesign scores ~0.97 and lands in the `novel` (= optimisation) class;
  - that class needs **≥ 3 mutations from every known CDRH3**, and abguard's *bundled*
    reference already contains the 524k FACS table, the SPR table and trastuzumab —
    no `--reference-db` needed;
  - the `de_novo` class (domain identity < 70 %) needs only 2.
  - So: 4 mutations from trastuzumab is **not** sufficient; the design must be 3 mutations
    from all 525,063 provided loops. Expect many `her2-2026-11-07` designs to fail here.
  - The task description also asks for CDRH3 identity to trastuzumab < 70 % (≥ 4 of the 10
    window positions changed). The shipped abguard does not check it; `classify()` does,
    unless `enforce_cdrh3_identity=False`.
  - **The 0.70 line is inclusive**: `novelty.py:104` `IDENTITY_THRESHOLD = 0.70` with
    `identity >= threshold` → `novel`. So identity *exactly* 0.70 is an optimisation design
    at the 3-mutation bar, not de novo. Measured: 2 of the 100 seed-1 VHHs land on 0.70
    exactly, and `abguard novelty-report` confirms `novelty_class = novel` for them.
    `abhack.designs.add_novelty` mirrors this; it used to call them `"unknown"`.
  - **`MUTATIONS_REQUIRED` is not configurable** (`novelty.py:146-149`, `{novel: 3,
    de_novo: 2}`) and distances are **Levenshtein** via `polars-distance`, so an indel
    counts as one mutation rather than invalidating the comparison.
  - The bundled reference is `abguard/reference_data/abguard_cdrh3.parquet`, **610,109**
    IMGT CDRH3s (`manifest.json`: scheme imgt, min_confidence 0.1), resolved by
    `reference.bundled_reference()` at `reference.py:673` and used automatically by
    `team-validate`. That is larger than our 525,063 provided loops — it carries prior art
    beyond the hackathon tables, so a design can pass our local check and still fail here.
- **Only `abguard team-validate` is trustworthy.** `novelty-report` prints a
  `--no-novelty: prior art was not checked` warning and leaves `is_novel_enough` null *even
  when `--reference-db` is passed*, and on `designs/_probe.csv` it disagreed outright:
  it said "2 of 3 novel enough" while `team-validate` rejected all three on novelty.
  Read `novelty-report` only for its distance and nearest-neighbour columns.
- **`designs/_probe.csv`'s annotations are wrong.** T03 is labelled "expect pass" (4
  mutations from trastuzumab) but **fails**: it is 2 mutations from `SRWGGDGFVAMVY`, a
  *different* FACS loop. Useful as a regression test — it demonstrates that distance from
  trastuzumab says nothing about distance from the other 525k loops — but don't trust the
  `design_description` text as the expected verdict.
- `abguard liability-report` additionally flags LAP sequence liabilities (isomerisation,
  deamidation, Trp/Met oxidation, N-glycosylation …). Trastuzumab itself carries several,
  so they are not disqualifying — but avoiding new ones in CDRH3 is free developability.
- Expected yield per weight (organisers' numbers, 300 designs each):
  `antibody-2026-10-30` 100 % pass novelty but 0 % known binders, 8.2 Å interface error;
  `her2-2026-11-07` balanced (`--sequence-tsr-k 2`) 6 % pass, exploratory (`1`) 74 % pass;
  `her2-novel-2026-11-09` 86 % pass, 6.8 Å — best default.
- `refold-score` leaves `pass` empty on the HER2 example (HER2 and the framework are fixed
  relative to each other, so RMSD is not meaningful). Rank on `min_ipae` (lower better),
  `iptm`, `binder_plddt`, and only within one folding model.
- Do not rank on CDRH3 refold agreement: folding models do not reproduce trastuzumab's own
  CDRH3 loop, so a mismatch there means little.
- **De novo needs a different molecule, not a different loop.** `abguard` splits the
  tracks on identity to trastuzumab alone (`novelty.py`: `identity = 1 - (d_H + d_L) /
  (max(|VH|,|ref VH|) + max(|VL|,|ref VL|))`, VHH uses heavy terms only). Measured:
  caplacizumab VHH framework 0.60, our generated VHHs 0.65, adalimumab VH 0.63 / VL 0.72
  — so a VHH clears the 0.70 line but a paired Fv on adalimumab is borderline.
- Two traps when writing a design spec, both cost a failed run:
  - the tutorial's VHH file encodes CDRs as residues *named* for their length
    (`13-13`, `9-10`, `12-21` at 23/38/76), which works as a direct `--input` but
    collides with a contig that states lengths itself — strip them, as
    `scripts/build_denovo_input.py` does, or the module dies with "Failed to compute
    leaving groups";
  - generated structures **renumber every chain from 1**, so HER2 C489→1 and D591→1.
    Anything comparing against 1N8Z numbering must shift (`abhack.epitope` does).
- File-format gotchas, found while wiring this up: `refold_metrics.csv` keys designs by
  **`query_id`** (its `design` column is a path); `per_example_metrics.csv` keys by
  `example_id`; the FACS zip also contains a `__MACOSX/` entry. All handled in `abhack`.
- Measured on this VM: smoke run (4 designs, refold 2 with Boltz-2) ≈ 5 min end to end;
  `min_ipae` ≈ 0.31–0.33, `iptm` ≈ 0.95, `binder_plddt` ≈ 96–99 for her2-novel designs.
- Baselines trained here: FACS 3-class classifier 0.91 held-out accuracy (135k rows,
  one-hot + HistGB); SPR pKD regressor CV R² 0.56 ± 0.13 on 190 points.

## Confidence calibration — done, 2026-11-14

`abhack.calibrate` (→ `reports/calibration.md`) fits the filter on
`data/raw/zenodo-18390239-ab-ag-complexes/iptm_confidence_scores.zip`: one 91×91 matrix per
predictor, **diagonal = real cognate pair, off-diagonal = shuffled pairing** (non-binder by
construction). ~91 positives and ~8,100 negatives per predictor.

| predictor | AUROC | enrichment @ ipTM 0.8 | cutoff for 5 % FPR | real kept there |
| --- | --- | --- | --- | --- |
| AF3 | 0.754 | 28.8× | 0.560 | 49.5 % |
| **Boltz-2** (ours) | **0.645** | **1.10×** | **0.951** | **17.6 %** |
| Chai-1 | 0.568 | 89× | 0.393 | 16.9 % |

- **Do not rank on `iptm_boltz2`.** At the conventional 0.8 cutoff Boltz-2 admits 81.7 % of
  known non-binders against 90.1 % of real pairs: 1.10× enrichment, effectively no
  information. Its usable cutoff is **0.951**, which keeps under a fifth of true binders.
  Our refolds all run on Boltz-2, so this is the column we had been ranking on.
- Chai-1's 89× is a trap: it keeps 1.1 % of real pairs, i.e. one. Enrichment alone never
  picks a filter — read it with the recall beside it.
- Rank instead on **interface PAE** against the organisers' vendored antibody criterion,
  `min_ipae ≤ 2.0 Å` (`data/docs/06-evaluate.md`; protein binder is pTM ≥ 0.8, iPAE ≤ 1.5 Å,
  RMSD ≤ 2.5 Å). The tutorial gives **no ipTM threshold anywhere**, consistent with the above.
  This benchmark contains ipTM only, so the PAE threshold is vendored, not calibrated here.
- Independent corroboration: Adaptyv's EGFR analysis puts these metrics at ~0.2 correlation
  with measured outcome. Two unrelated datasets, same conclusion — **a coarse filter, never
  a fine ranking.**

Honesty constraints that must survive into any write-up: the negatives are *computational*
shuffled pairings, not measured non-binders; the "454 VHH–antigen SPR true negatives" named
in the task PDF **do not exist in the released data** (what exists is 12 non-cognate rows
marked `no reportable fit`, which the organisers' README explicitly forbids relabelling as
negatives, plus **4** cognate fitted KDs, not 15); and a shuffled pair is an *easier*
negative than an optimised design, so every FPR above is a lower bound on ours.

AUROC is hand-rolled (no sklearn in this env) and pinned in `tests/test_calibrate.py`
against five hand-computed cases including tie handling — a wrong AUROC would make a useless
predictor look usable, which is the error the module exists to prevent.

## Submission-path invariants — do not re-break these

`abhack.submission` lost designs silently four different ways. Each fix has a test in
`tests/test_submission.py`; if a test there fails, read this list before "fixing" it.

- **`passes = None` means "undecided, abguard decides" — never "failed".** A de novo VHH has
  no trastuzumab CDRH3 window, so only abguard can grade its CDR3. `df.filter(pl.col("passes"))`
  drops nulls, which deleted the **entire de novo track** before ranking.
- **Dedupe on sequences, not on `cdrh3`.** `cdrh3` is null for every VHH and polars treats
  all nulls as one value, so the whole VHH set collapsed to a single row. Key is
  `(kind, format, vh_sequence, vl_sequence)`.
- **`kind` is per row**, derived from the `format` column `parse_design` emits. One global
  stamp mislabels a track in a mixed submission and fails abguard's per-partition schema.
- **A missing `rank_by` raises.** It used to skip the sort and return concatenation order,
  which looks exactly like a ranking and is not one. This is the failure mode when a VHH
  table (no `*_boltz2` columns) is concatenated with a refolded one.
- **The 250 aa cap and the unpaired-cysteine rule are enforced, not warned.** Both were
  already computed and discarded; the cysteine check now covers VL too, which it never did.

## Conventions

- Raw data stays in `data/raw/` (gitignored); parquet caches in `data/processed/`.
- One run folder per design config: `runs/<weight>-s<seed>[-<setting>]`, with `run.log`.
- Every design table keeps `design_id`, `cdrh3`, `vh_sequence`, `vl_sequence` so it can be
  concatenated with any other table.
- Ranking decisions and their reasons go in `reports/` — 40 % of the pitch score is judgment,
  so write down *why* a design was picked, not just the number.

## Decision log

**2026-11-14, ~12:45 — dropped the de novo VHH track, moved GPU to optimisation.**

Evidence, internal:

| Batch | min_ipae | iptm | binder_pLDDT | RMSD to design | pass |
| --- | --- | --- | --- | --- | --- |
| CDRH3 redesign, trastuzumab framework | 0.31–0.33 | 0.95 | 96–99 | — | n/a |
| De novo VHH, top 20 of 100 (seed 1) | 2.48–15+ | 0.43–0.82 | 89–93 | **15–51 Å** | **0/20** |

The de novo VHHs fold well (pLDDT ~90) but Boltz-2 docks them 15–51 Å from the
designed pose — they are well-formed nanobodies that bind somewhere else. All
100 had scored 100/100 on hotspot obedience, which is RFD4 grading its own work.

Evidence, external — Adaptyv's EGFR competition (601 designs characterised):

- optimised existing binders **18 %** hit rate, diversified **28 %**, best de
  novo method (BindCraft) **9 %**; RFdiffusion produced one binder in total;
- their explicit recommendation: *start from an existing binder and redesign its
  non-interface regions*;
- top affinity came from 10 framework mutations to cetuximab, not de novo;
- ranking metrics (ipTM, iPAE, ESM PLL) correlate ~0.2 with outcome and not at
  all within the top 200 — a coarse filter, never a fine ranking;
- expression ~95 %, "essentially solved" in cell-free, so developability
  scoring is a marginal gain here.

**The tactic that follows:** the novelty gate only examines CDRH3, so framework
mutations are free. A design can carry a CDRH3 novel enough to pass the filter,
the inherited paratope for epitope preservation, and stabilising framework
mutations for thermostability — one design improving all three axes, which is
what S-tier requires.

For ranking, prefer the models trained on this system (FACS classifier 0.91
accuracy, SPR regressor R² 0.56) over generic ESM pseudo-likelihood (AUC 0.72
in Adaptyv's analysis).

Sources: adaptyvbio.com/blog/po104, blog.booleanbiotech.com/what-we-learned-adaptyv,
biorxiv 2025.04.17.648362.

**Post-mortem on the same batch (independent analysis, ~13:00) — two explanations ruled
out by measurement, one still open.** Read the open one first: it limits how strongly the
0/20 can be blamed on the designs.

- **NOT ruled out: the refold ran without an MSA.** `runs/denovo-vhh-s1/prepared/holo/queries.json`
  has `msa_path: None` on **every** chain (C, D and H), so Boltz-2 predicted single-sequence
  despite `missing_msas: fail` in its config. Boltz-2 is an MSA-dependent model, so its pose
  prediction was handicapped and 0/20 overstates how bad the designs are. Two mitigations:
  the target still reproduced to 0.45 Å because the template carried it (below), and a *de
  novo* binder has no natural MSA to find anyway — but the docking judgement is still made
  on a weakened model. **The right control is AF2-IG** (`scripts/japanfold_af2ig.py`), which
  starts from the designed coordinates instead of rediscovering the pose from sequence.
  Until that runs, treat "the de novo designs do not bind where designed" as probable, not
  established.

Exact columns, seed-1 top 20, `runs/denovo-vhh-s1/scored-boltz2/refold_metrics.csv`:
`binder_ptm` 0.952–0.985, `binder_plddt` 89.2–94.4, `min_ipae` 2.48–20.84,
`target_aligned_binder_rmsd` 10.9–55.2 Å, `backbone_rmsd` 6.2–17.7 Å, `pass` False ×20.
Note `pass` carries real reasons here (`refolded_target_aligned_binder_rmsd.rmsd=27.5`),
unlike the smoke run — **the de novo spec does fix enough group-0 target atoms**, so there
is no spec fix to make on that account.

- **Not the split target.** The hotspots (`C572 C573 D591 D593`) straddle the artificial
  chain break at the unresolved 581–590 loop, so the obvious suspect was C and D drifting
  apart during refolding. Measured, superposing each refold on `prepared/holo/template.cif`:
  chain C self-RMSD **0.45 Å** (0.38–0.57), chain D displaced **2.0 Å** (1.3–3.3), C572–D591
  gap **15.8 Å** against 16.25 Å designed. The target reproduces almost exactly. Splitting
  HER2 domain IV across two chains did no harm — don't spend time "fixing" it.
- **Not foldability.** `binder_ptm` 0.95–0.99 and `binder_plddt` ~90: these are well-formed
  nanobodies. They simply dock 10.9–55.2 Å from the designed pose, i.e. they bind elsewhere.
- **The metrics contradict each other, which is the decisive tell.** `_15` has the best
  `min_ipae` (2.48) *and* best `iptm` (0.82) but sits 39 Å from its designed pose; `_70` has
  the best RMSD (10.9 Å) with `min_ipae` 18.2. No design is good on both, so none is a real
  binder — Boltz-2 is confidently predicting *an* interface, just not the designed one.
- **Root cause at the level of the weights.** De novo ran on `antibody-2026-10-30`, the
  *general* antibody checkpoint, which knows nothing about HER2. Both HER2-finetuned
  checkpoints only redesign CDRH3 on a fixed framework (weights README: `--model`
  her2-conservative/balanced/exploratory/her2-novel; `12-her2.md`: "designs CDR loops on a
  fixed antibody framework"). So the de novo track structurally **could not use** the one
  4.2 GB asset actually trained on this target. That is the headline reason it underperformed.
- Weak predictors of docking accuracy on those 20 (`target_aligned_binder_rmsd`, Spearman,
  **n=20 so suggestive only**): `n_contacts` ρ = −0.58, `hbond_atom_pair_count_coi` ρ = −0.44,
  `vh_length` ρ = **+0.50** (longer docks worse — the `12-21` CDR3 range is too permissive).
  `hotspot_obedience` ρ = −0.24 and `epitope_fraction` ρ = +0.16 are both uninformative here.
- Caveat on the comparison in the table above: CDRH3 redesign refolds against a **templated
  known crystal complex**, so `min_ipae` 0.31 there partly reflects the answer being supplied.
  The optimisation track's advantage is real, but those particular numbers overstate it.

## Open work

- ~~Calibrate the confidence filter on the Smorodina ipTM data~~ — **done**, see
  *Confidence calibration* above and `reports/calibration.md`. What remains is **wiring the
  result into ranking**: `abhack.submission` should rank on `min_ipae ≤ 2.0 Å` and stop
  consulting `iptm_boltz2`.
- **Join `epitope_fraction` into the design table.** `abhack.epitope` is imported by no other
  module; its output sits in `designs/epitope-s1.csv`, keyed by `design_id`, ready to join in
  `designs.collect`. It discriminates as a *measurement* (37 of 100 designs below 1.0, range
  0.464–1.0) and it is judged axis 2 — but on the refolded 20 it did **not** predict docking
  success, so treat it as a measure of the judged property, not as a binding predictor.
  `contacts_hotspot` is True for all 100 and useless as a filter.
- **Collapse the two ranking paths** (`evaluate.py` vs `submission.py`). One must go.
- Thermostability has no provided data — needs a proxy (e.g. framework germlining,
  Aggrescan3D/SoluProt-style developability scores, or an ESM pseudo-likelihood). The **130
  liability occurrences across 20 designs** that `abguard liability-report` already prints
  (0 hard errors) are free input and are currently parsed by nothing.
- Expression/solubility scoring (Aggrescan 3D, SoluProt, NetSolP, ProteinMPNN soluble
  weights) is recommended by the organisers but not installed yet. Note `08-adaptyv.md` does
  **not** document these — they come from the task PDF; `08-adaptyv.md` is a miniprotein
  funnel case study with no developability method in it.
- **A second-opinion refold from a different architecture.** Every design we have passed
  through one design model and one folding model, and the calibration says that folding model
  is the weakest measured. Nothing currently separates "the designs are bad" from "Boltz-2 is
  wrong about them". Untouched and already on disk: `of3ob` (OpenBind), `of3p2`, `protenix_v2`,
  `mol1` — all template-capable. Not calibratable from the benchmark (it covers only
  AF3/Boltz-2/Chai-1), so treat as a consistency check. `esmfold2` and `boltz1` take no
  template and fold the target unguided — don't use them for binders.
- **JapanFold is installed but never called** (`~/.claude/skills/japanfold/SKILL.md`, added
  12:27, verbatim upstream including an `allowed-tools` grant of `Bash(curl *)` and
  `Bash(python3 *)`; needs `JAPANFOLD_API_KEY` in the shell). Its **AF2-IG** scores a binder
  you already designed (structure + binder sequence → ipTM and interface PAE) and is
  architecturally independent of everything in Foundry, which is exactly the orthogonal check
  the calibration argues for, since AF3-family ipTM is the only kind that separated.

## Housekeeping

- `designs/denovo-s1.csv` and `designs/denovo-vhh-s1.csv` are the **same 100 designs with
  identical column sets** — collected twice in parallel. Delete one.
- A diagram of this pipeline, its three gates and the evidence above lives at
  **`docs/pipeline.html`** — standalone, open it in a browser, no server or build step.
  Mirrored at https://claude.ai/artifact/7EFQa8tE1fePRTgGhWH18Y (anyone with the link). The
  repo copy is the one that survives the VM; keep the two in step if either changes.
- Tests: 35 passing (`.venv/bin/python -m pytest tests/ -q`) across `test_novelty.py`,
  `test_submission.py`, `test_calibrate.py`.
