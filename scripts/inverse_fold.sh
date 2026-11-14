#!/usr/bin/env bash
# Sequence decode only — no diffusion over structure.
#
# RFD4-Proteina keeps every backbone coordinate and runs just its sequence head,
# so you get new sequences for a structure you already have. You decide what is
# decoded by what you leave out of C_SEQ in the spec; see
# specs/inverse-fold/design.json.
#
#   bash scripts/inverse_fold.sh                       # the HER2 CDRH3 example
#   SPEC=specs/mine/design.json N=50 bash scripts/inverse_fold.sh
#   STRUCT=runs/opt-novel-s11/design/.../x_generated.cif bash scripts/inverse_fold.sh
#
#   SPEC    design spec (default specs/inverse-fold/design.json)
#   STRUCT  decode a single structure directly, no spec — every chain kept,
#           the model decides what to decode from the file's own annotation
#   N       sequences to generate      (default 20)
#   SEED    random seed                (default 1)
#   OUT     output folder              (default runs/inverse-<seed>)
#   WEIGHT  antibody | her2 | her2-novel  (default her2-novel)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source env.sh

N="${N:-20}"
SEED="${SEED:-1}"
WEIGHT="${WEIGHT:-her2-novel}"
OUT="${OUT:-runs/inverse-s${SEED}}"
IMG="$(rfd4_image "$WEIGHT")"

if [[ -n "${STRUCT:-}" ]]; then
  INPUT="$STRUCT"
else
  SPEC="${SPEC:-specs/inverse-fold/design.json}"
  WORK="data/$(basename "$(dirname "$SPEC")")"
  mkdir -p "$WORK"
  cp "$SPEC" "$WORK/design.json"
  # the spec names its structure relative to its own folder
  STRUCT_NAME=$(.venv/bin/python -c "import json,sys; d=json.load(open('$SPEC')); print(next(v['input'] for k,v in d.items() if isinstance(v,dict) and 'input' in v))")
  [[ -f "$WORK/$STRUCT_NAME" ]] || cp "data/docs/Data/12-her2/$STRUCT_NAME" "$WORK/" 2>/dev/null || true
  INPUT="$WORK/design.json"
fi

echo "inverse folding: input=$INPUT n=$N seed=$SEED weight=$WEIGHT -> $OUT"
mkdir -p "$OUT"
rm -rf "$OUT/design"
abrun "$IMG" module-run design \
  --input "$INPUT" --output "$OUT/design" \
  --inverse-folding true \
  --num-replicates "$N" --seed "$SEED" \
  --recycle-strategy high/front_loaded

printf '%s\n' "inverse-folding n=$N seed=$SEED input=$INPUT date=$(date -Is)" >> "$OUT/run.log"
echo "done. Sequences:"
echo "  .venv/bin/python -m abhack.designs collect $OUT"
