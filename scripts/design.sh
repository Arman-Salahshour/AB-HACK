#!/usr/bin/env bash
# Generate designs with RFD4-Proteina.
#
#   WEIGHT=her2-novel N=100 SEED=1 OUT=runs/novel-s1 bash scripts/design.sh
#
#   WEIGHT   antibody | her2 | her2-novel          (default her2-novel)
#   N        designs per input, max 200            (default 20)
#   SEED     random seed                           (default 1)
#   TSR_K    sequence randomness: 1 varied, 4 default; only for WEIGHT=her2
#   INPUT    design spec                           (default the HER2 CDRH3 example)
#   OUT      output folder                         (default runs/<weight>-s<seed>)
#   EXTRA    extra flags passed to module-run design
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source env.sh

WEIGHT="${WEIGHT:-her2-novel}"
N="${N:-20}"
SEED="${SEED:-1}"
INPUT="${INPUT:-data/docs/Data/12-her2/her2_cdrh3.yaml}"
OUT="${OUT:-runs/${WEIGHT}-s${SEED}}"
IMG="$(rfd4_image "$WEIGHT")"

args=(--input "$INPUT" --output "$OUT/design"
      --num-replicates "$N" --seed "$SEED"
      --recycle-strategy high/front_loaded)
[[ -n "${TSR_K:-}" ]] && args+=(--sequence-tsr-k "$TSR_K")
[[ -n "${EXTRA:-}" ]] && read -r -a extra <<< "$EXTRA" && args+=("${extra[@]}")

echo "design: weight=$WEIGHT n=$N seed=$SEED tsr_k=${TSR_K:-default} -> $OUT/design"
mkdir -p "$OUT"
rm -rf "$OUT/design"
abrun "$IMG" module-run design "${args[@]}"

printf '%s\n' "weight=$WEIGHT n=$N seed=$SEED tsr_k=${TSR_K:-} input=$INPUT date=$(date -Is)" \
  >> "$OUT/run.log"
echo "done. Next: OUT=$OUT bash scripts/refold.sh"
