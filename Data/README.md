# Data directory

This directory holds the input and output data for the pipeline.

## Tracked in git

- `gene_sets/` — local GMT caches (Reactome 2022, LINCS L1000 Chemical
  Perturbation Consensus) used by the 5-FU mechanistic analyses (notebook
  `10`). These are tracked so that those analyses run fully offline.

## Generated at runtime (git-ignored)

The remaining subdirectories are produced by running the notebooks and are
excluded from version control:

| Path | Created by | Contents |
|------|-----------|----------|
| `samplesheet.csv` | notebook `00` | sample metadata |
| `fastqc/`, `multiqc/` | notebook `01` | quality control |
| `AlignBAM/` | notebook `02` | STAR alignments |
| `count_matrix_raw.tsv`, `count_matrix_processed.tsv` | notebooks `03`–`04` | count matrices |
| `tables/` | notebooks `05`–`10` | DESeq2 / GSEA / Cox result tables |
| `TCGA/` | `scripts/download/download_tcga.sh` | TCGA COAD/READ data |
| `5FU_resistance/` | `scripts/download/download_geo.sh` | GSE153412 data |
| `kras_allele_counts.tsv` | notebook `09` | KRAS allele counts |
| `logs/` | background jobs | log files |
