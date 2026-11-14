#!/usr/bin/env bash
# Thermal screen for paired Fv candidates: keep A (VL) + B (VH), strip HER2 C/D.
set -euo pipefail
OUTDIR=${OUTDIR:-runs/thermal}; NS_HOT=${NS_HOT:-5}; NS_CTRL=${NS_CTRL:-2}; T_HOT=${T_HOT:-400}
REPLICAS=${REPLICAS:-3}; CONCURRENCY=${CONCURRENCY:-4}; BINDER_CHAINS=${BINDER_CHAINS:-A,B}
(( $# )) || { echo "Usage: $0 designs/structures/*_boltz2.pdb" >&2; exit 2; }
mkdir -p "$OUTDIR" jobs; : > jobs/thermal_joblist.txt
for pdb in "$@"; do
  stem=$(basename "${pdb%.pdb}"); design_id=$(sed -E 's/^[0-9]+_//; s/_(boltz2|designed)$//' <<<"$stem")
  printf '%q ' python scripts/run_thermal.py --pdb "$pdb" --out "$OUTDIR/${design_id}_T310_r1" --design-id "$design_id" --binder-chains "$BINDER_CHAINS" --format fv --temp 310 --ns "$NS_CTRL" --seed 1 >> jobs/thermal_joblist.txt; printf '\n' >> jobs/thermal_joblist.txt
  for seed in $(seq 1 "$REPLICAS"); do
    printf '%q ' python scripts/run_thermal.py --pdb "$pdb" --out "$OUTDIR/${design_id}_T${T_HOT}_r${seed}" --design-id "$design_id" --binder-chains "$BINDER_CHAINS" --format fv --temp "$T_HOT" --ns "$NS_HOT" --seed "$seed" >> jobs/thermal_joblist.txt; printf '\n' >> jobs/thermal_joblist.txt
  done
done
xargs -P "$CONCURRENCY" -I{} bash -c '{}' < jobs/thermal_joblist.txt
python scripts/analyze_thermal.py --runs "$OUTDIR/*_T${T_HOT}_*" --out "$OUTDIR/thermal_scores.csv"
python scripts/analyze_thermal.py --runs "$OUTDIR/*_T310_*" --out "$OUTDIR/thermal_control.csv"
