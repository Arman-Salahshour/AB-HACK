#!/usr/bin/env bash
# Refold designs and score them: refold-prepare -> fold -> refold-score.
#
#   OUT=runs/novel-s1 TOPN=20 FOLD=boltz2 bash scripts/refold.sh
#
#   OUT     run folder holding design/           (required)
#   TOPN    designs to refold                    (default 20)
#   FOLD    boltz2|of3ob|of3p2|protenix-v2|protenix-v1|mol1|boltz1|esmfold2 (default boltz2)
#   WEIGHT  design image used for prepare/score  (default her2-novel)
#   SEEDS   folding seeds                        (default 42)
#   SAMPLES diffusion samples per query          (default 1)
#
# The HER2 example leaves `pass` empty; rank on min_ipae (lower better), iptm
# and binder_plddt, and only within one folding model.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source env.sh

OUT="${OUT:?set OUT to the run folder, e.g. OUT=runs/novel-s1}"
TOPN="${TOPN:-20}"
FOLD="${FOLD:-boltz2}"
WEIGHT="${WEIGHT:-her2-novel}"
SEEDS="${SEEDS:-42}"
SAMPLES="${SAMPLES:-1}"

RFD="$(rfd4_image "$WEIGHT")"
IFS='|' read -r FOLD_IMG FOLD_MODEL <<< "$(fold_image "$FOLD")"

echo "[1/3] prepare top $TOPN from $OUT/design"
rm -rf "$OUT/prepared"
abrun "$RFD" module-run refold-prepare \
  --input "$OUT/design" --output "$OUT/prepared" --top-n "$TOPN"

echo "[2/3] fold with $FOLD ($FOLD_MODEL)"
rm -rf "$OUT/$FOLD"
predict=(predict --input "$OUT/prepared/holo/queries.json" --output "$OUT/$FOLD"
         --write-detailed-scores true --num-diffusion-samples "$SAMPLES"
         --seeds "$SEEDS" --query-failures skip)
if fold_has_template "$FOLD"; then
  predict+=(--template-paths "$OUT/prepared/holo/template.cif")
fi
CMD=$(docker inspect -f '{{join .Config.Cmd " "}}' "$FOLD_IMG")
# shellcheck disable=SC2086
abrun "$FOLD_IMG" $CMD "${predict[@]}"

echo "[3/3] score"
rm -rf "$OUT/scored-$FOLD"
abrun "$RFD" module-run refold-score \
  --input "$OUT/design" --output "$OUT/scored-$FOLD" \
  --predictions "$OUT/$FOLD" --model "$FOLD_MODEL"

column -s, -t < "$OUT/scored-$FOLD/refold_metrics.csv" | cut -c1-200
echo "done. Next: python -m abhack.designs collect $OUT -o designs/$(basename "$OUT").csv"
