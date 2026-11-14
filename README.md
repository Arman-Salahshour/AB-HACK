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
OUT=runs/novel-s1 TOPN=20 bash scripts/refold.sh
python -m abhack.designs collect runs/novel-s1 -o designs/novel-s1.csv
python -m abhack.submission build designs/*.csv -o designs/submission.csv
bash scripts/validate.sh designs/submission.csv
```

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
