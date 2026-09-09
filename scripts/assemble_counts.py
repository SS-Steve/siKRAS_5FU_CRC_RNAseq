"""Step 03: Assemble a genes x samples count matrix from STAR output."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.helpers import ensure_dir


def assemble_counts(align_dir: str, output: str,
                    strandedness: str = "unstranded") -> pd.DataFrame:
    """
    Assemble STAR ReadsPerGene.out.tab files into a genes x samples matrix.

    Parameters
    ----------
    align_dir : str
        Directory with per-sample STAR output subdirectories.
    output : str
        Output TSV path (genes x samples).
    strandedness : str
        One of "unstranded", "forward", "reverse".
        STAR column mapping:
          col2 = unstranded, col3 = 1st-strand, col4 = 2nd-strand.
    """
    col_map = {"unstranded": 1, "forward": 2, "reverse": 3}
    if strandedness not in col_map:
        raise ValueError(f"strandedness must be one of {list(col_map)}")
    count_col = col_map[strandedness]

    align_path = Path(align_dir)
    matrices = []
    for sample_dir in sorted(align_path.iterdir()):
        if not sample_dir.is_dir():
            continue
        count_file = sample_dir / "ReadsPerGene.out.tab"
        if not count_file.exists():
            print(f"[Skip] {sample_dir.name}: no ReadsPerGene.out.tab")
            continue
        df = pd.read_csv(count_file, sep="\t", header=None, index_col=0,
                         usecols=[0, count_col], names=["gene", sample_dir.name])
        matrices.append(df)

    if not matrices:
        raise FileNotFoundError("No STAR count files found.")

    matrix = pd.concat(matrices, axis=1)
    # Remove STAR summary rows (N_unmapped, N_multimapping, etc.)
    matrix = matrix[~matrix.index.str.startswith("N_")]
    ensure_dir(str(Path(output).parent))
    matrix.to_csv(output, sep="\t")
    print(f"Count matrix: {matrix.shape[0]} genes x {matrix.shape[1]} samples -> {output}")
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser(description="Assemble STAR count matrix.")
    parser.add_argument("--align-dir", required=True)
    parser.add_argument("--output", default="Data/count_matrix_raw.tsv")
    parser.add_argument("--strandedness", default="unstranded",
                        choices=["unstranded", "forward", "reverse"])
    args = parser.parse_args()
    assemble_counts(args.align_dir, args.output, args.strandedness)


if __name__ == "__main__":
    main()
