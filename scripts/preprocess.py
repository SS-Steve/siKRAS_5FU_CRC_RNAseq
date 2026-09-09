"""Step 04: Preprocess count matrix (filter zero-count genes, transpose)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.helpers import ensure_dir


def preprocess(input_file: str, output: str) -> pd.DataFrame:
    """
    Filter zero-count genes and transpose to samples x genes.

    Only genes with zero total counts across ALL samples are removed
    (no aggressive low-expression filtering; DESeq2 handles that internally).
    """
    counts = pd.read_csv(input_file, sep="\t", index_col=0)  # genes x samples
    n_before = counts.shape[0]
    counts = counts[counts.sum(axis=1) > 0]
    n_after = counts.shape[0]
    print(f"Removed {n_before - n_after} zero-count genes "
          f"({n_before} -> {n_after})")

    # Transpose to samples x genes for DESeq2
    counts = counts.T
    ensure_dir(str(Path(output).parent))
    counts.to_csv(output, sep="\t")
    print(f"Processed matrix: {counts.shape[0]} samples x {counts.shape[1]} genes -> {output}")
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess count matrix.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="Data/count_matrix_processed.tsv")
    args = parser.parse_args()
    preprocess(args.input, args.output)


if __name__ == "__main__":
    main()
