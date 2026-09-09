"""Step 02: STAR alignment with --quantMode GeneCounts."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.helpers import ensure_dir


def run_star(samplesheet: str, genome_dir: str, out_dir: str,
             threads: int = 30, quant_mode: str = "GeneCounts",
             star_bin: str = "STAR", force: bool = False) -> None:
    """
    Align each sample with STAR and produce per-gene counts.

    Outputs are written to out_dir/<sample>/Aligned.sortedByCoord.out.bam
    and ReadsPerGene.out.tab.
    """
    sheet = pd.read_csv(samplesheet)
    ensure_dir(out_dir)

    for _, row in sheet.iterrows():
        sample = row["sample"]
        sample_out = Path(out_dir) / sample
        expected = sample_out / "ReadsPerGene.out.tab"
        if expected.exists() and not force:
            print(f"[Skip] {sample}: output exists")
            continue

        sample_out.mkdir(parents=True, exist_ok=True)
        read_files = [row["fastq_1"]]
        if row.get("layout") == "paired" and row.get("fastq_2"):
            read_files.append(row["fastq_2"])

        cmd = [
            star_bin,
            "--runThreadN", str(threads),
            "--genomeDir", genome_dir,
            "--readFilesIn", *read_files,
            "--outFileNamePrefix", f"{sample_out}/",
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--quantMode", quant_mode,
        ]
        # Handle gzipped FASTQ files
        if any(str(f).endswith(".gz") for f in read_files):
            cmd.extend(["--readFilesCommand", "zcat"])

        print(f"[Run] {sample}: STAR ...")
        subprocess.run(cmd, check=True)

    print(f"STAR alignment complete -> {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="STAR alignment.")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--genome-dir", required=True)
    parser.add_argument("--out-dir", default="Data/AlignBAM")
    parser.add_argument("--threads", type=int, default=30)
    parser.add_argument("--quant-mode", default="GeneCounts")
    parser.add_argument("--star-bin", default="STAR")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    run_star(args.samplesheet, args.genome_dir, args.out_dir,
             args.threads, args.quant_mode, args.star_bin, args.force)


if __name__ == "__main__":
    main()
