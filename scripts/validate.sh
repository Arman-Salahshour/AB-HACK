#!/usr/bin/env bash
# Check a submission table with abguard — the organisers' authority on the hard filters.
#
#   bash scripts/validate.sh designs/submission.csv [--drop]
#
# Columns required: sequence_id, kind, heavy_sequence, light_sequence, design_description.
# Uses abguard's bundled reference database; pass REFERENCE_DB=<path> to override it.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

TABLE="${1:?usage: scripts/validate.sh <submission.csv> [extra abguard flags]}"
shift || true

ref=()
[[ -n "${REFERENCE_DB:-}" ]] && ref=(--reference-db "$REFERENCE_DB")

echo "== novelty report =="
abguard novelty-report "${ref[@]}" "$TABLE" || true

echo
echo "== liabilities =="
abguard liability-report "$TABLE" || true

echo
echo "== validate (hard filters + novelty gate) =="
abguard team-validate "${ref[@]}" "$TABLE" "$@"
