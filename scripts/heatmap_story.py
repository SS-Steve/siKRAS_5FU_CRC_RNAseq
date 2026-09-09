#!/usr/bin/env python3
"""Generate the Figure 5F heatmap: representative KRAS-pathway and
proliferation/repair genes (row Z-scored) across treatment groups."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from utils import plotting
from utils.helpers import ensure_dir

# Genes selected from the GSEA leading edges and classic KRAS feedback markers.
UP_GENES = [
    "SPRY2", "SPRY4", "DUSP6", "DUSP4", "FOSL1", "HBEGF",
    "EPHA2", "ITGA2", "PLAU", "PLAUR",
]
DOWN_GENES = [
    # G2/M checkpoint
    "AURKA", "CCNB2", "PTTG1", "CDKN3", "NEK2", "KATNA1", "HMMR", "LMNB1",
    # DNA repair
    "POLR2C", "TAF13", "DGUOK", "TP53", "POLR2A", "DDB2", "GTF2A2",
    # mTORC1 / metabolism
    "MTHFD2", "HSPA5", "TXNRD1", "G6PD", "GCLC",
]


def main() -> None:
    plotting.setup_style()

    counts = pd.read_csv("Data/count_matrix_processed.tsv", sep="\t", index_col=0)  # samples x genes
    sheet = pd.read_csv("Data/samplesheet.csv").set_index("sample")

    genes = [g for g in UP_GENES + DOWN_GENES if g in counts.columns]
    missing = [g for g in UP_GENES + DOWN_GENES if g not in counts.columns]
    if missing:
        print(f"  [heatmap] genes not found: {missing}")

    # Sample order: Control -> Empty -> siScr -> siKRAS
    group_order = ["Control", "Empty", "siScr", "siKRAS"]
    ordered_samples = []
    for g in group_order:
        ordered_samples += [s for s in counts.index if sheet.loc[s, "group"] == g]

    expr = np.log1p(counts.loc[ordered_samples, genes])
    # Row Z-score (per gene across samples)
    expr_z = expr.subtract(expr.mean(axis=0), axis=1).div(expr.std(axis=0), axis=1)
    expr_z = expr_z.T  # genes x samples

    fig, ax = plt.subplots(figsize=(8, max(6, len(genes) * 0.32)))
    sns.heatmap(expr_z, cmap=plotting.BWR_CMAP, center=0,
                cbar_kws={"label": "Row Z-score", "shrink": 0.6},
                linewidths=0.4, linecolor="white", ax=ax)
    ax.set_xticklabels(ordered_samples, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(genes, fontsize=9)

    # Separate up/down gene blocks with a horizontal line
    n_up = len([g for g in UP_GENES if g in counts.columns])
    ax.axhline(n_up, color="black", linewidth=1.2)
    ax.text(-0.5, n_up / 2, "KRAS feedback", rotation=90, va="center",
            ha="right", fontsize=9, fontweight="bold")
    ax.text(-0.5, n_up + (len(genes) - n_up) / 2, "Proliferation / repair",
            rotation=90, va="center", ha="right", fontsize=9, fontweight="bold")

    ax.set_title("KRAS silencing suppresses proliferation and repair programmes")
    plotting.save_figure(fig, "Figures/06_visualization/heatmap_story.png")
    print(f"  Heatmap -> Figures/06_visualization/heatmap_story.png ({len(genes)} genes, {len(ordered_samples)} samples)")


if __name__ == "__main__":
    main()
