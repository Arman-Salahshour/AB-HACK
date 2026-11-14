#!/usr/bin/env bash
# Two-pH (6.5 / 7.4) complex MD batch.
#   SYSTEMS="ctrl:C,D opt43:C,D"  name:antigen-chains, each needs runs/phmd/inputs/<name>.pdb
#   env: NS (ns per run), REPS, PHS, CONCURRENCY. Default SYSTEMS is the trastuzumab control only:
#   the control is what defines the noise band the ranking uses, so it must always be in the batch.
set -euo pipefail
cd "$(dirname "$0")/.."
NS=${NS:-1}; REPS=${REPS:-2}; PHS=${PHS:-"6.5 7.4"}; CONCURRENCY=${CONCURRENCY:-3}; D=runs/phmd
: > "$D/joblist.txt"
SYSTEMS=${SYSTEMS:-"ctrl:C,D"}
case " $SYSTEMS " in *" ctrl:"*) ;; *) SYSTEMS="ctrl:C,D $SYSTEMS" ;; esac
for spec in $SYSTEMS; do
  name=${spec%%:*}; ag=${spec##*:}
  for ph in $PHS; do for r in $(seq 1 "$REPS"); do
    [ "${SKIP_EXISTING:-1}" = 1 ] && [ -f "$D/${name}_pH${ph}_r${r}.json" ] && continue
    echo ".venv/bin/python scripts/run_ph_md.py --pdb $D/inputs/$name.pdb --out $D/${name}_pH${ph}_r${r} --ph $ph --ns $NS --seed $r --antigen $ag" >> "$D/joblist.txt"
  done; done
done
xargs -P "$CONCURRENCY" -I{} bash -c '{} || echo "FAILED: {}"' < "$D/joblist.txt"
.venv/bin/python scripts/analyze_ph_md.py "$D/*_pH*_r*" -o "$D/ph_md_summary.csv"
