"""Step 09: KRAS G13D allele-specific analysis via pysam pileup."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils import plotting
from utils.helpers import ensure_dir


def count_kras_alleles(bam_dir: str, samplesheet: str,
                       chrom: str = "NC_000012.12",
                       position: int = 25245347,
                       out_file: str = "Data/kras_allele_counts.tsv") -> pd.DataFrame:
    """
    Count WT vs G13D reads at the KRAS c.38 position in each BAM.

    The G13D mutation is c.38G>A. On the minus strand, the reference genome
    has C (WT) and the mutant reads carry T (G13D).

    Parameters
    ----------
    bam_dir : str
        Directory with per-sample STAR output (Aligned.sortedByCoord.out.bam).
    samplesheet : str
        Sample metadata CSV.
    chrom, position : str, int
        Genomic coordinate of the variant (1-based).
    out_file : str
        Output TSV with WT/G13D/total counts per sample.
    """
    sheet = pd.read_csv(samplesheet)
    rows = []
    for sample in sheet["sample"]:
        bam_path = Path(bam_dir) / sample / "Aligned.sortedByCoord.out.bam"
        if not bam_path.exists():
            print(f"[Skip] {sample}: BAM not found")
            continue

        bam = pysam.AlignmentFile(str(bam_path), "rb")
        n_wt = n_g13d = 0
        for col in bam.pileup(chrom, position - 1, position, truncate=True,
                              min_base_quality=0, min_mapping_quality=0):
            if col.pos != position - 1:
                continue
            for read in col.pileups:
                if read.is_del or read.is_refskip:
                    continue
                base = read.alignment.query_sequence[read.query_position].upper()
                if base == "C":
                    n_wt += 1
                elif base == "T":
                    n_g13d += 1
        bam.close()
        rows.append({"sample": sample, "WT": n_wt, "G13D": n_g13d,
                     "total": n_wt + n_g13d})

    df = pd.DataFrame(rows)
    ensure_dir(str(Path(out_file).parent))
    df.to_csv(out_file, sep="\t", index=False)
    print(f"KRAS allele counts -> {out_file}")
    return df


def plot_allele_counts(allele_df: pd.DataFrame, out_path: str) -> None:
    """Plot WT vs G13D read counts and G13D% by group."""
    plotting.setup_style()
    df = allele_df.copy()
    df["group"] = df["sample"].str.replace(r"[_\-.]\d+$", "", regex=True)
    df["G13D_pct"] = df["G13D"] / df["total"] * 100

    group_order = ["Control", "Empty", "siScr", "siKRAS"]
    group_order = [g for g in group_order if g in df["group"].unique()]
    # Reorder samples by group_order so Panel A matches Panel B
    df["_group_rank"] = df["group"].map({g: i for i, g in enumerate(group_order)})
    df = df.sort_values(["_group_rank", "sample"]).reset_index(drop=True)
    df = df.drop(columns=["_group_rank"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: stacked bars of WT vs G13D
    ax = axes[0]
    x = np.arange(len(df))
    w = 0.35
    ax.bar(x - w/2, df["WT"], w, color=plotting.COLOR_DOWN, label="KRAS WT", alpha=0.9)
    ax.bar(x + w/2, df["G13D"], w, color=plotting.COLOR_UP, label="KRAS G13D", alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(df["sample"], rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Read count")
    ax.set_title("KRAS allele counts")
    plotting.style_legend(ax)
    plotting.add_grid(ax, axis="y")

    # Panel B: G13D% by group
    ax = axes[1]
    means = df.groupby("group")["G13D_pct"].mean().reindex(group_order)
    stds = df.groupby("group")["G13D_pct"].std().reindex(group_order)
    colors = {"Control": "#000000", "Empty": "#A9A9A9",
              "siScr": "#88ADFF", "siKRAS": "#FF7E79"}
    x = np.arange(len(group_order))
    ax.bar(x, means, color=[colors[g] for g in group_order], width=0.5, alpha=0.9)
    ax.errorbar(x, means, yerr=stds, fmt="none", color="black", capsize=5)
    for i, g in enumerate(group_order):
        sub = df[df["group"] == g]["G13D_pct"]
        ax.scatter([i]*len(sub), sub, color="black", s=40, zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels(group_order)
    ax.set_ylabel("G13D allele (%)")
    ax.set_title("Allele balance")
    plotting.add_grid(ax, axis="y")

    plotting.save_figure(fig, out_path)
    print(f"KRAS allele plot -> {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="KRAS allele counting.")
    parser.add_argument("--bam-dir", required=True)
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--chrom", default="NC_000012.12")
    parser.add_argument("--position", type=int, default=25245347)
    parser.add_argument("--output", default="Data/kras_allele_counts.tsv")
    args = parser.parse_args()
    count_kras_alleles(args.bam_dir, args.samplesheet,
                       args.chrom, args.position, args.output)


if __name__ == "__main__":
    main()
