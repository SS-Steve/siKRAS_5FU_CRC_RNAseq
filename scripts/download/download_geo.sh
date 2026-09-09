#!/usr/bin/env bash
# ============================================================================
# Download the 5-FU resistance dataset from GEO.
#   GSE153412: HCT116 5-FU resistant (116RU) vs sensitive (116SU), 3v3.
# Used by notebook 10 (fivefu_mechanism.py) for the resistance-reversal test.
# Outputs are saved to Data/5FU_resistance/.
# ============================================================================
set -euo pipefail

OUT_DIR="Data/5FU_resistance"
mkdir -p "$OUT_DIR"

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153412/suppl"

echo "Downloading GSE153412 estimated counts..."
curl -L -o "$OUT_DIR/GSE153412_est_counts.tsv.gz" \
  "$BASE/GSE153412_est_counts.tsv.gz"

echo "Downloading GSE153412 TPM..."
curl -L -o "$OUT_DIR/GSE153412_tpm.tsv.gz" \
  "$BASE/GSE153412_tpm.tsv.gz"

echo "Done. Files saved to $OUT_DIR/"
