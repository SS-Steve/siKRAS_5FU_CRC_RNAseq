#!/usr/bin/env bash
# ============================================================================
# Download TCGA CRC (COAD + READ) data from the Zenodo mirror of UCSC Xena,
# then merge COAD + READ into combined CRC files used by notebook 08.
# Outputs are saved to Data/TCGA/.
# ============================================================================
set -euo pipefail

OUT_DIR="Data/TCGA"
mkdir -p "$OUT_DIR"

ZENODO_BASE="https://zenodo.org/api/records/15122700/files"

echo "Downloading COAD raw counts..."
curl -L -o "$OUT_DIR/COAD_counts.tsv.gz" \
  "$ZENODO_BASE/Colon.Adenocarcinoma_TCGA_rsem.genes_expected_count_clean_raw_counts.tabular.gz/content"

echo "Downloading COAD survival..."
curl -L -o "$OUT_DIR/COAD_survival.tsv.gz" \
  "$ZENODO_BASE/Colon.Adenocarcinoma_TCGA_survival_clean.tabular.gz/content"

echo "Downloading READ raw counts..."
curl -L -o "$OUT_DIR/READ_counts.tsv.gz" \
  "$ZENODO_BASE/Rectum.Adenocarcinoma_TCGA_rsem.genes_expected_count_clean_raw_counts.tabular.gz/content"

echo "Downloading READ survival..."
curl -L -o "$OUT_DIR/READ_survival.tsv.gz" \
  "$ZENODO_BASE/Rectum.Adenocarcinoma_TCGA_survival_clean.tabular.gz/content"

echo "Merging COAD + READ into combined CRC files..."
python - <<'PY'
import pandas as pd

coad_c = pd.read_csv("Data/TCGA/COAD_counts.tsv.gz", sep="\t", index_col=0)
read_c = pd.read_csv("Data/TCGA/READ_counts.tsv.gz", sep="\t", index_col=0)
crc_c = pd.concat([coad_c, read_c], axis=1)
crc_c.to_csv("Data/TCGA/CRC_raw_counts.tsv.gz", sep="\t")

coad_s = pd.read_csv("Data/TCGA/COAD_survival.tsv.gz", sep="\t", index_col=0)
read_s = pd.read_csv("Data/TCGA/READ_survival.tsv.gz", sep="\t", index_col=0)
crc_s = pd.concat([coad_s, read_s], axis=0)
crc_s.to_csv("Data/TCGA/CRC_survival.tsv.gz", sep="\t")

print(f"  CRC_raw_counts.tsv.gz: {crc_c.shape[0]} genes x {crc_c.shape[1]} samples")
print(f"  CRC_survival.tsv.gz:   {crc_s.shape[0]} samples")
PY

echo "Done. Files saved to $OUT_DIR/"
