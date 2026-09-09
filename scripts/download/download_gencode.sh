#!/usr/bin/env bash
# ============================================================================
# Download GENCODE v47 GTF (used for TCGA TPM gene-length calculation).
# ============================================================================
set -euo pipefail

OUT_DIR="Data/GENCODE"
mkdir -p "$OUT_DIR"

echo "Downloading GENCODE v47 annotation GTF..."
curl -L -o "$OUT_DIR/gencode.v47.annotation.gtf.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz"

echo "Done. Saved to $OUT_DIR/gencode.v47.annotation.gtf.gz"
