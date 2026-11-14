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
```

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

## Open work

- Calibrate the confidence filter on the Smorodina cognate/non-cognate ipTM data
  (`data/raw/zenodo-18390239-ab-ag-complexes/iptm_confidence_scores.zip`) before trusting
  ipTM/ipAE to rank designs.
- Thermostability has no provided data — needs a proxy (e.g. framework germlining,
  Aggrescan3D/SoluProt-style developability scores, or an ESM pseudo-likelihood).
- Expression/solubility scoring (Aggrescan 3D, SoluProt, NetSolP, ProteinMPNN soluble
  weights) is recommended by the organisers but not installed yet.
