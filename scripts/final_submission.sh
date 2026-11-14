#!/usr/bin/env bash
# Build the final ranked submission from the collected design tables.
#
# The whole ordering policy in one place (see src/abhack/ranking.py for the
# reasoning behind each component):
#
#   1. rank    — gates (novelty, ipae) before scoring; weighted z-composite
#   2. tickets — paper-validated thermostability twins of the top 4 designs
#   3. merge   — 12 optimisation / 8 de novo slots, spill on evidence
#   4. build   — abguard schema, submission_order ascending (rank 1 = row 1)
#   5. validate— abguard team-validate is the authority; run it or do not submit
#
#   bash scripts/final_submission.sh            # everything, writes submission.csv
#   SINCE_GATED=0 bash scripts/final_submission.sh   # skip joining gated refolds
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
t = pl.read_csv(table, infer_schema_length=10_000).drop(
    [c for c in keep.values() if c in t.columns]
)
t.join(sdf, on="design_id", how="left").write_csv(table)
PY
done

# 1-2. rank + tickets (twins inherit their parent's evidence by construction)
.venv/bin/python -m abhack.ranking rank "${OPT_TABLES[@]}" -o designs/opt-ranked.csv
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

# 3. quota-merge the tracks
.venv/bin/python -m abhack.ranking merge designs/opt-ranked.csv "$DENOVO" \
  --opt-slots 12 --denovo-slots 8 -o designs/final-ranked.csv

# 4. convert to the abguard schema — ascending: submission_order 1 is row 1
.venv/bin/python -m abhack.submission build designs/final-ranked.csv \
  --rank-by submission_order -o designs/submission.csv

# 5. the authority
bash scripts/validate.sh designs/submission.csv
