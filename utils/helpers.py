"""Common helper functions for the siKRAS RNA-seq pipeline."""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import pandas as pd
import yaml


def load_config(config_path: str = "config.yaml") -> dict:
    """Load the YAML configuration file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def read_samplesheet(
    samplesheet: str, group_col: str = "group"
) -> pd.DataFrame:
    """
    Read a samplesheet CSV and return sample metadata indexed by sample name.

    Parameters
    ----------
    samplesheet : str
        Path to the samplesheet CSV.
    group_col : str
        Column name defining experimental groups.

    Returns
    -------
    pd.DataFrame
        Metadata with a 'condition' column derived from group_col.
    """
    sheet = pd.read_csv(samplesheet)
    if "sample" not in sheet.columns:
        raise ValueError("Samplesheet must contain a 'sample' column.")
    meta = sheet.set_index("sample")
    if group_col not in meta.columns:
        raise ValueError(f"Group column '{group_col}' not found in samplesheet.")
    meta = meta.rename(columns={group_col: "condition"})
    return meta


def read_table(path: str, **kwargs) -> pd.DataFrame:
    """Read a TSV/CSV table, transparently handling gzip compression."""
    if str(path).endswith(".gz"):
        return pd.read_csv(path, compression="gzip", **kwargs)
    return pd.read_csv(path, **kwargs)


def ensure_dir(path: str) -> None:
    """Create a directory if it does not exist."""
    Path(path).mkdir(parents=True, exist_ok=True)


def launch_background(cmd: list, log_file: str) -> None:
    """
    Launch a command in the background (detached, like nohup).

    The command survives the notebook kernel, so it is suitable for
    long-running steps (FastQC, STAR alignment). Progress is written to
    `log_file`, which the user can monitor with `tail -f`.

    Parameters
    ----------
    cmd : list
        Command and arguments (e.g. ["python", "scripts/fastqc.py", "..."]).
    log_file : str
        Path to the log file.
    """
    import subprocess
    ensure_dir(str(Path(log_file).parent))
    log = open(log_file, "w")
    proc = subprocess.Popen(
        cmd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True
    )
    print(f"Launched in background (PID {proc.pid})")
    print(f"Monitor progress: tail -f {log_file}")


def strip_tcga_barcode(sample_id: str, length: int = 12) -> str:
    """Strip a TCGA barcode to its patient-level identifier (e.g., TCGA-XX-XXXX)."""
    return str(sample_id)[:length]


def list_sample_dirs(raw_dir: str) -> list[str]:
    """List sample subdirectories under a raw data directory."""
    return sorted(
        d.name for d in Path(raw_dir).iterdir() if d.is_dir()
    )


def read_gtf_gene_lengths(gtf_path: str) -> dict[str, int]:
    """
    Parse a GENCODE GTF and return gene_name -> union exon length.

    Parameters
    ----------
    gtf_path : str
        Path to a gzipped GTF file.

    Returns
    -------
    dict
        Mapping of gene symbol to total exon length (bp).
    """
    from collections import defaultdict

    gene_exons = defaultdict(set)
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9 or parts[2] != "exon":
                continue
            attrs = parts[8]
            gene_name = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("gene_name "):
                    gene_name = kv.split(" ")[1].strip('"')
                    break
            if gene_name:
                gene_exons[gene_name].add((int(parts[3]), int(parts[4])))

    gene_lengths = {}
    for gene, exons in gene_exons.items():
        total = 0
        cur_start = cur_end = None
        for start, end in sorted(exons):
            if cur_start is None:
                cur_start, cur_end = start, end
            elif start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                total += cur_end - cur_start + 1
                cur_start, cur_end = start, end
        if cur_start is not None:
            total += cur_end - cur_start + 1
        gene_lengths[gene] = total
    return gene_lengths
