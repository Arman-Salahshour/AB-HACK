#!/usr/bin/env bash
# Build the final ranked submission from the collected design tables.
#
# The whole ordering policy in one place (see src/abhack/ranking.py for the
# reasoning behind each component):
#
#   0b. pH     — propka selectivity (6.5 vs 7.4) for every design -> designs/ph-scores.csv
#   1. rank    — gates (novelty, ipae, pH >= 0) before scoring; weighted z-composite.
#                pH credit needs two-pH MD agreement; thermal MD joins if its results exist
#   2. tickets — paper-validated thermostability twins of the top 4 designs
#   3. merge   — 12 optimisation / 8 de novo slots, spill on evidence
#   4. build   — abguard schema, submission_order ascending (rank 1 = row 1)
#   5. validate— abguard team-validate is the authority; run it or do not submit
#
#   bash scripts/final_submission.sh            # everything, writes submission.csv
#   SINCE_GATED=0 bash scripts/final_submission.sh   # skip joining gated refolds
#   PH_GATE=0 ...      # drop the pH gate;  RUN_PH_MD=1 SYSTEMS="ctrl:C,D opt43:C,D" ...  # run two-pH MD first
#   RUN_THERMAL=1 THERMAL_PDBS="designs/structures/*_boltz2.pdb" ...  # run thermal MD on a FINAL shortlist
#   REQUIRE_THERMAL=1 ...  # refuse to rank if any design lacks thermal results or fails the 310 K control
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source env.sh

OPT_TABLES=(designs/opt-s11.csv designs/opt-s12.csv designs/opt-s13.csv)
DENOVO=designs/denovo-ranked.csv
GATED_RUNS=(runs/opt-novel-s12-gated runs/opt-explore-s13-gated)
[ "${SINCE_GATED:-1}" = 1 ] || GATED_RUNS=()

# 0. fold gated-refold metrics back into their design tables
for run in "${GATED_RUNS[@]}"; do
  metrics="$run/scored-boltz2/refold_metrics.csv"
  [ -f "$metrics" ] || { echo "skip $run (no refold metrics yet)"; continue; }
  case "$run" in
    *novel-s12*) src_table=designs/opt-s12.csv ;;
    *explore-s13*) src_table=designs/opt-s13.csv ;;
    *) echo "no table mapping for $run; skip"; continue ;;
  esac
  echo "joining $metrics -> $src_table"
  .venv/bin/python - "$metrics" "$src_table" <<'PY'
import sys
import polars as pl

metrics, table = sys.argv[1], sys.argv[2]
s = pl.read_csv(metrics)
keep = {c: f"{c}_boltz2" for c in ("min_ipae", "iptm", "binder_plddt") if c in s.columns}
sdf = s.select(["query_id", *keep]).rename({"query_id": "design_id", **keep})
t = pl.read_csv(table, infer_schema_length=10_000)
t = t.drop([c for c in keep.values() if c in t.columns])
t.join(sdf, on="design_id", how="left").write_csv(table)
PY
done

# 0b. pH selectivity for every design (propka, own designed complex); unscored = undecided
.venv/bin/python scripts/ph_score_pool.py "${OPT_TABLES[@]}" "$DENOVO" -o designs/ph-scores.csv

# optional MD stages (GPU, slow): two-pH interface MD and thermal unfolding MD
if [ "${RUN_PH_MD:-0}" = 1 ]; then bash scripts/run_ph_batch.sh; fi
if [ "${RUN_THERMAL:-0}" = 1 ]; then
  : "${THERMAL_PDBS:?set THERMAL_PDBS to the final shortlist structures}"
  # shellcheck disable=SC2086
  bash scripts/run_all_thermal.sh $THERMAL_PDBS
fi

RANK_EXTRA=(--ph-scores designs/ph-scores.csv)
[ -f runs/phmd/ph_md_summary.csv ] && RANK_EXTRA+=(--ph-md runs/phmd/ph_md_summary.csv)
[ "${PH_GATE:-1}" = 1 ] || RANK_EXTRA+=(--no-ph-gate)
if [ -f runs/thermal/thermal_scores.csv ]; then
  RANK_EXTRA+=(--thermal-results runs/thermal/thermal_scores.csv)
  [ -f runs/thermal/thermal_control.csv ] && RANK_EXTRA+=(--thermal-control runs/thermal/thermal_control.csv)
  [ "${REQUIRE_THERMAL:-0}" = 1 ] && RANK_EXTRA+=(--require-thermal)
fi

# 1-2. rank + tickets (twins inherit their parent's evidence by construction)
.venv/bin/python -m abhack.ranking rank "${OPT_TABLES[@]}" "${RANK_EXTRA[@]}" -o designs/opt-ranked.csv
.venv/bin/python -m abhack.tickets designs/opt-ranked.csv -n 4 -o designs/opt-tickets.csv
.venv/bin/python - <<'PY'
import polars as pl
ranked = pl.read_csv("designs/opt-ranked.csv", infer_schema_length=10_000)
tickets = pl.read_csv("designs/opt-tickets.csv", infer_schema_length=10_000)
pl.concat([ranked, tickets], how="diagonal_relaxed").sort(
    "composite", descending=True, nulls_last=True
).write_csv("designs/opt-ranked.csv")
print(f"opt pool: {ranked.height} + {tickets.height} tickets")
PY

# 2b. apply the same pH gate/credit to the de novo table so both tracks face one policy
.venv/bin/python -m abhack.ranking rank "$DENOVO" "${RANK_EXTRA[@]}" -o designs/denovo-ph-ranked.csv

# 3. quota-merge the tracks
.venv/bin/python -m abhack.ranking merge designs/opt-ranked.csv designs/denovo-ph-ranked.csv \
  --opt-slots 12 --denovo-slots 8 -o designs/final-ranked.csv

# 4. convert to the abguard schema — ascending: submission_order 1 is row 1
.venv/bin/python -m abhack.submission build designs/final-ranked.csv \
  --rank-by submission_order -o designs/submission.csv

# 5. the authority
bash scripts/validate.sh designs/submission.csv
