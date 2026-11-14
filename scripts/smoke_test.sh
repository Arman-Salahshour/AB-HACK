#!/usr/bin/env bash
# End-to-end check that the toolchain works: 4 designs -> refold 2 -> score -> collect.
# About 5 minutes on the H100.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

OUT=${OUT:-runs/smoke}
WEIGHT=${WEIGHT:-her2-novel} N=${N:-4} SEED=${SEED:-1} OUT="$OUT" bash scripts/design.sh
OUT="$OUT" TOPN=${TOPN:-2} FOLD=${FOLD:-boltz2} WEIGHT="${WEIGHT:-her2-novel}" bash scripts/refold.sh

if [[ -x .venv/bin/python ]]; then
  .venv/bin/python -m abhack.designs collect "$OUT" -o "designs/smoke.csv"
  .venv/bin/python -m abhack.evaluate "designs/smoke.csv" --model "${FOLD:-boltz2}"
fi
echo "smoke test OK"
