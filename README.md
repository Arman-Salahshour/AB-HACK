# ab-hack — Oslo Antibodies × AI hackathon (HER2 / trastuzumab)

Pipeline for improving trastuzumab against HER2, judged on three axes measured by Adaptyv:

1. **Affinity** — SPR binding to HER2, vs. the trastuzumab baseline.
2. **Epitope preservation** — competitive binding against trastuzumab.
3. **Thermostability** — direct measurement.

A design that beats the baseline on all three is S-tier; two → A-tier; one → B-tier.
Up to 20 ranked designs are submitted, split freely between the **de novo** and
**optimisation** tracks.

## Hard filters (enforced by `abguard`)

| Filter | Rule |
| --- | --- |
| Format | VHH, scFv or Fab (paired) |
| Antibody-ness | IMGT numbering confidence ≥ 0.1 with `immunum` |
| Length | ≤ 250 aa total (single chain, or VH+VL combined) |
| Cysteines | no unpaired cysteines |
| Novelty (de novo) | domain identity < 70 % to any known antibody, ≥ 2 CDRH3 mutations |
| Novelty (optimisation) | ≥ 3 CDRH3 mutations from **every** provided sequence; task description also asks for < 70 % CDRH3 identity to trastuzumab |

`abguard` decides the track from identity over the **whole domain**, so a CDRH3-only
redesign on trastuzumab's framework (~97 % identical) is always an optimisation design.
It compares the IMGT CDR3 — `SRWGGDGFYAMDY`, i.e. the 10-residue window `WGGDGFYAMD`
plus its `SR`/`Y` flanks — against a bundled reference that already contains the FACS
and SPR tables. Being 4 mutations from trastuzumab is not enough: the loop must be
3 mutations away from all 525,063 provided CDRH3s. A design is never moved between
tracks; if it misses its track's criteria it is rejected.

## Layout

```
data/raw/          hackathon data from s3://apheris-oslo-hackathon/data/ (not in git)
data/docs/         RFD4-Proteina tutorial: Data/, Settings/, Scripts/, Tutorial.html
data/processed/    parquet caches built by src/abhack/data.py
runs/              docker run outputs: design / prepared / folded / scored (not in git)
designs/           candidate tables and the ranked submission
scripts/           thin shell wrappers around the Foundry docker images
src/abhack/        python package: data, novelty, design parsing, scoring, models
tests/             unit tests for the filter logic
reports/           figures and notes for the pitch
```

## Quickstart

```bash
uv sync                                  # create .venv from pyproject.toml
source env.sh                            # image tags + `abrun` docker helper
nvidia-smi                               # H100 should be visible

bash scripts/fetch_data.sh small         # HER2 tables + tutorial docs (~12 MB)
bash scripts/smoke_test.sh               # 4 designs → refold → score (~5 min)

# full design run (100 designs, ~30 s each)
WEIGHT=her2-novel N=100 SEED=1 OUT=runs/novel-s1 bash scripts/design.sh
python -m abhack.designs collect runs/novel-s1 -o designs/novel-s1.csv
# ONE ranking path: gates (novelty, ipae) BEFORE refolding, then a weighted
# z-composite; see src/abhack/ranking.py for the components and why iptm is absent
python -m abhack.ranking rank designs/novel-s1.csv -o designs/opt-ranked.csv
# gate-first refold: only survivors, selected by the ranking — not RFD4's top-n
python scripts/refold_gated.py runs/novel-s1 designs/opt-ranked.csv --out runs/novel-s1-gated
OUT=runs/novel-s1-gated TOPN=<n> bash scripts/refold.sh
# thermostability tickets: paper-validated framework twins of the top designs
python -m abhack.tickets designs/opt-ranked.csv -n 4 -o designs/opt-tickets.csv
# quota-merge the tracks (12 optimisation / 8 de novo default) and convert
python -m abhack.ranking merge designs/opt-ranked.csv designs/denovo-ranked.csv -o designs/final-ranked.csv
python -m abhack.submission build designs/final-ranked.csv --rank-by submission_order --descending -o designs/submission.csv
bash scripts/validate.sh designs/submission.csv
```

## De novo track: VHH against the trastuzumab epitope

The HER2-finetuned weights only redesign CDRH3 on trastuzumab's framework, so they
produce ~97 %-identical designs that are always optimisation. For de novo you need a
different molecule: `antibody-2026-10-30` designing all three CDRs of a VHH, docked at
trastuzumab's epitope. A VHH sits at ~0.65 identity to trastuzumab's VH — under the
0.70 line, so `abguard` classes it `de_novo`, where only 2 CDRH3 mutations are needed.

**Verdict, 2026-11-14: do not rerun this.** Both seeds (200 designs) refolded 0/20
at the designed pose under Boltz-2, and the AF2-IG control — which starts from the
designed coordinates and cannot be fooled by a missing MSA — agrees: ipTM 0.10–0.34,
interface PAE 19–27 Å on all 40 tested (`designs/af2ig-denovo.csv`). Two independent
architectures, same answer: well-folded nanobodies that bind somewhere else. The de
novo track instead draws on BoltzGen (`scripts/japanfold_boltzgen.py`, off-box) with
the same AF2-IG check (`scripts/japanfold_af2ig_boltzgen.py`); best candidate ipTM
0.46 / PAE 16 Å — still below the credibility bar, so de novo entries are ranked
lottery tickets for a thin track, not evidence-backed binders
(`scripts/build_denovo_table.py` → `designs/denovo-ranked.csv`).

```bash
python scripts/build_denovo_input.py          # writes data/denovo-vhh/
WEIGHT=antibody N=100 SEED=1 INPUT=data/denovo-vhh/design.json \
  OUT=runs/denovo-vhh-s1 bash scripts/design.sh
OUT=runs/denovo-vhh-s1 TOPN=20 FOLD=boltz2 WEIGHT=antibody bash scripts/refold.sh
python -m abhack.epitope runs/denovo-vhh-s1   # is it binding where trastuzumab binds?
python -m abhack.designs collect runs/denovo-vhh-s1 -o designs/denovo-s1.csv
```

The input folder (`scripts/build_denovo_input.py`) holds three files, per
`data/docs/15-antibody-design.md`:

| File | What |
| --- | --- |
| `target_her2.cif` | HER2 domain IV, chains C (489–580) + D (591–607) from the 1N8Z crop, trastuzumab's own chains dropped |
| `framework_vhh.cif` | caplacizumab VHH framework, 83 residues with CDR gaps at 23/38/76 |
| `design.json` | tracked in `specs/denovo-vhh/`; contig designs CDR1 13, CDR2 9–10, CDR3 12–21 |

`C_HOT` pins the four HER2 atoms nearest trastuzumab's CDRH3 in 1N8Z, which is what
makes the design compete with trastuzumab in the competitive-binding assay — the
epitope-preservation axis. `abhack.epitope` then measures whether a design actually
landed there: contacts within 4.5 Å, scored against the 44 HER2 residues within 10 Å
of those hotspots.

Unlike the HER2 CDRH3 exercise, `refold-score`'s `pass` column is meaningful here: a de
novo VHH has a real designed pose, so RMSD-to-design is measurable.

## Models on this VM

Design (`apheris-rfd4-proteina`, tasks `design`, `refold-prepare`, `refold-score`, `metrics`, `finetune`):

| Name in `env.sh` | Weight | Use |
| --- | --- | --- |
| `RFD4_ANTIBODY` | `antibody-2026-10-30` | any CDR / VHH / Fv on a fixed framework |
| `RFD4_HER2` | `her2-2026-11-07` | trastuzumab CDRH3; `--sequence-tsr-k 2` balanced, `1` exploratory |
| `RFD4_HER2_NOVEL` | `her2-novel-2026-11-09` | mostly novel CDRH3 loops, lowest interface error (6.8 Å) |

Folding: `FOLD_BOLTZ2`, `FOLD_OF3OB`, `FOLD_OF3P2`, `FOLD_PROTENIX_V2`, `FOLD_PROTENIX_V1`,
`FOLD_MOL1`, `FOLD_BOLTZ1` (no template), `FOLD_ESMFOLD2` (no template).

For the HER2 example `refold-score` leaves `pass` empty (the example fixes HER2 and the
framework relative to each other), so rank on `min_ipae` (lower better), `iptm` and
`binder_plddt`, and only within one folding model.

## Data

| Path | Contents |
| --- | --- |
| `data/raw/her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip` | 524,346 CDRH3 variants, FACS class high/medium/low |
| `data/raw/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv` | 768 designs, 190 with KD (nM) |
| `data/raw/zenodo-18390239-ab-ag-complexes/iptm_confidence_scores.zip` | ipTM matrices, cognate vs non-cognate — calibrate the confidence filter |
| `data/raw/sigma2-ab-ag-champloo/spr-supplementary/` | SPR cross-reactivity tables |

Both trastuzumab tables together are the "provided data" the novelty gate compares against
(525,063 unique CDRH3s). Cite Chinery et al. when using them.

## References

- Task description: `ai-antibodies-hackathon-apheris.pdf`
- Tutorial: `data/docs/Tutorial.html`, VM guide `data/docs/18-vm.md`, HER2 models `data/docs/12-her2.md`
