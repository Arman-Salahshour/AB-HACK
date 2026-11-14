#!/usr/bin/env bash
# Two-pH complex MD batch. Usage: bash scripts/run_ph_batch.sh   (env: NS, REPS, PHS, CONCURRENCY)
set -euo pipefail
cd "$(dirname "$0")/.."
NS=${NS:-1}; REPS=${REPS:-2}; PHS=${PHS:-"6.5 7.4"}; CONCURRENCY=${CONCURRENCY:-3}; D=runs/phmd
: > "$D/joblist.txt"
for spec in ctrl:C,D opt43:C,D opt28:C,D phsw0037:C,D dn8shard:A; do
  name=${spec%%:*}; ag=${spec##*:}
  for ph in $PHS; do for r in $(seq 1 "$REPS"); do
    echo ".venv/bin/python scripts/run_ph_md.py --pdb $D/inputs/$name.pdb --out $D/${name}_pH${ph}_r${r} --ph $ph --ns $NS --seed $r --antigen $ag" >> "$D/joblist.txt"
  done; done
done
xargs -P "$CONCURRENCY" -I{} bash -c '{} || echo "FAILED: {}"' < "$D/joblist.txt"
.venv/bin/python scripts/analyze_ph_md.py "$D/*_pH*_r*" -o "$D/ph_md_summary.csv"
