#!/usr/bin/env bash
# Copy hackathon data out of the read-only bucket.
#
#   bash scripts/fetch_data.sh small   # HER2 tables + tutorial docs + ipTM scores (~12 MB)
#   bash scripts/fetch_data.sh ngs     # ready-made FACS train/val/test splits (2.5 GB)
#   bash scripts/fetch_data.sh champloo# per-replica confidence archives (1.8 GB)
#   bash scripts/fetch_data.sh all
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

BUCKET=s3://apheris-oslo-hackathon
RAW=data/raw
mkdir -p "$RAW"

small() {
  aws s3 cp "$BUCKET/data/README.md" "$RAW/README.md"
  aws s3 sync "$BUCKET/data/her2-trastuzumab/" "$RAW/her2-trastuzumab/"
  aws s3 sync "$BUCKET/rfd4-proteina/docs/" data/docs/
  aws s3 cp "$BUCKET/data/zenodo-18390239-ab-ag-complexes/iptm_confidence_scores.zip" \
    "$RAW/zenodo-18390239-ab-ag-complexes/iptm_confidence_scores.zip"
  aws s3 sync "$BUCKET/data/sigma2-ab-ag-champloo/spr-supplementary/" \
    "$RAW/sigma2-ab-ag-champloo/spr-supplementary/"
}

ngs() { aws s3 sync "$BUCKET/data/zenodo-10831512-trastuzumab-ngs/" "$RAW/zenodo-10831512-trastuzumab-ngs/"; }

champloo() { aws s3 sync "$BUCKET/data/sigma2-ab-ag-champloo/" "$RAW/sigma2-ab-ag-champloo/"; }

accuracy() {
  aws s3 cp "$BUCKET/data/zenodo-18390239-ab-ag-complexes/topmodel_structure_accuracy.zip" \
    "$RAW/zenodo-18390239-ab-ag-complexes/topmodel_structure_accuracy.zip"
}

case "${1:-small}" in
  small) small ;;
  ngs) ngs ;;
  champloo) champloo ;;
  accuracy) accuracy ;;
  all) small; accuracy; ngs; champloo ;;
  *) sed -n '2,9p' "$0"; exit 1 ;;
esac
du -sh "$RAW"/* 2>/dev/null || true
