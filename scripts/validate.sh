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

# abguard ships a reference bundle (610,109 IMGT CDRH3s) inside its own package, and
# `team-validate` falls back to it with no flag. Pass it explicitly anyway so the
# reports agree about what they compared against. abguard is a uv tool, so its package
# is not importable from this venv — resolve it on disk.
if [[ -z "${REFERENCE_DB:-}" ]]; then
  REFERENCE_DB=$(ls -d "$HOME"/.local/share/uv/tools/abguard/lib/python*/site-packages/abguard/reference_data 2>/dev/null | head -1 || true)
fi

ref=()
if [[ -n "${REFERENCE_DB:-}" && -d "${REFERENCE_DB}" ]]; then
  ref=(--reference-db "$REFERENCE_DB")
else
  echo "WARNING: no reference database resolved; team-validate will use its own bundle" >&2
fi

# NOTE: `novelty-report` prints a "--no-novelty: prior art was not checked" warning and
# leaves is_novel_enough null even WITH --reference-db, and on a probe table it disagreed
# with team-validate (said 2 of 3 novel enough; team-validate rejected all 3 on novelty).
# Trust `team-validate` — it is the gate the organisers run — and read novelty-report only
# for its per-design distances and nearest-neighbour columns.

echo "== novelty report =="
abguard novelty-report "${ref[@]}" "$TABLE" || true

echo
echo "== liabilities =="
abguard liability-report "$TABLE" || true

echo
echo "== validate (hard filters + novelty gate) =="
abguard team-validate "${ref[@]}" "$TABLE" "$@"
