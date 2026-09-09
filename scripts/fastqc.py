"""Step 01: Run FastQC quality control on all FASTQ files."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.helpers import ensure_dir


def run_fastqc(samplesheet: str, output_dir: str, fastqc_bin: str = "fastqc",
               n_jobs: int = 6) -> None:
    """
    Run FastQC on every FASTQ file listed in the samplesheet.

    Parameters
    ----------
    samplesheet : str
        Path to samplesheet CSV (from step 00).
    output_dir : str
        Directory for FastQC reports.
    fastqc_bin : str
        Path or name of the FastQC executable.
    n_jobs : int
        Number of parallel FastQC jobs.
    """
    sheet = pd.read_csv(samplesheet)
    fastq_files = sheet["fastq_1"].tolist()
    if "fastq_2" in sheet.columns:
        fastq_files += sheet["fastq_2"].dropna().tolist()

    ensure_dir(output_dir)
    processes = []
    for f in fastq_files:
        p = subprocess.Popen(
            [fastqc_bin, "--extract", "-q", "-o", output_dir, f],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        processes.append(p)
        if len(processes) >= n_jobs:
            for proc in processes:
                proc.wait()
            processes = []
    for proc in processes:
        proc.wait()
    print(f"FastQC complete: {len(fastq_files)} files -> {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run FastQC on all FASTQ files.")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--output-dir", default="Data/fastqc")
    parser.add_argument("--fastqc-bin", default="fastqc")
    parser.add_argument("--jobs", type=int, default=6)
    args = parser.parse_args()
    run_fastqc(args.samplesheet, args.output_dir, args.fastqc_bin, args.jobs)


if __name__ == "__main__":
    main()
