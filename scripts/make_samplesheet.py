"""Step 00: Generate a samplesheet from raw FASTQ directories."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.helpers import ensure_dir


def _find_read_files(fq: list, read_number: int) -> list:
    """Find R1/R2 files by matching the read number at the END of the filename.

    The read number must appear immediately before the extension (e.g.
    ``Sample_1.fq.gz``, ``Sample_R2.fastq.gz``). This avoids false matches when
    the sample name itself contains the other read number (e.g. ``Control_2_1.fq.gz``
    must match R1, not R2).
    """
    pattern = re.compile(rf"[_.]R?{read_number}\.(fq|fastq)(\.gz)?$")
    return [f for f in fq if pattern.search(f.name)]


def detect_layout(subdir: Path) -> str:
    """Detect paired-end vs single-end from FASTQ files in a directory."""
    files = sorted(subdir.glob("*.fq*")) + sorted(subdir.glob("*.fastq*"))
    has_r1 = len(_find_read_files(files, 1)) > 0
    has_r2 = len(_find_read_files(files, 2)) > 0
    return "paired" if (has_r1 and has_r2) else "single"


def make_samplesheet(raw_dir: str, output: str) -> pd.DataFrame:
    """
    Scan raw_dir for sample subdirectories and build a samplesheet.

    The 'group' column is auto-inferred by stripping a trailing replicate
    suffix (e.g. Control_1 -> Control). Verify the result before DESeq2.
    """
    raw_path = Path(raw_dir)
    rows = []
    for subdir in sorted(raw_path.iterdir()):
        if not subdir.is_dir():
            continue
        sample = subdir.name
        layout = detect_layout(subdir)
        fq = sorted(subdir.glob("*.fq*")) + sorted(subdir.glob("*.fastq*"))
        # Auto-infer group by stripping a trailing replicate suffix (_1, _2, ...)
        group = re.sub(r"[_\-.]\d+$", "", sample)
        if layout == "paired":
            r1 = _find_read_files(fq, 1)[0]
            r2 = _find_read_files(fq, 2)[0]
            rows.append({
                "sample": sample,
                "layout": layout,
                "fastq_1": str(r1),
                "fastq_2": str(r2),
                "group": group,
            })
        else:
            rows.append({
                "sample": sample,
                "layout": layout,
                "fastq_1": str(fq[0]),
                "fastq_2": "",
                "group": group,
            })
    df = pd.DataFrame(rows)
    ensure_dir(str(Path(output).parent))
    df.to_csv(output, index=False)
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a samplesheet from raw FASTQ.")
    parser.add_argument("--raw-dir", required=True, help="Directory with sample subdirectories.")
    parser.add_argument("--output", default="Data/samplesheet.csv")
    args = parser.parse_args()

    df = make_samplesheet(args.raw_dir, args.output)
    print(f"Generated samplesheet with {len(df)} samples -> {args.output}")
    print("NOTE: verify the auto-inferred 'group' column before running DESeq2.")


if __name__ == "__main__":
    main()
