"""Step 06: Visualization — MA, volcano, PCA, heatmaps."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils import plotting
from utils.helpers import ensure_dir

# Housekeeping genes for MA/volcano annotation
HK_GENES = {
    "ACTB", "GAPDH", "B2M", "HPRT1", "PPIA", "RPLP0", "TBP", "UBC",
    "YWHAZ", "PGK1", "RPL13A", "GUSB", "SDHA", "HMBS",
    "TFRC", "RPL19", "IPO8", "CYC1", "EIF4A2",
}


def ma_plot(res: pd.DataFrame, out_path: str, shrunken: bool = False) -> None:
    """MA plot: log2FC vs log2 baseMean.

    Parameters
    ----------
    res : DataFrame
        DESeq2 results (must contain baseMean, log2FoldChange, padj).
    out_path : str
        Output PNG path.
    shrunken : bool
        If True, label the y-axis as apeglm-shrunken log2 fold change.
    """
    plotting.setup_style()
    res = res.copy()
    res["A"] = np.log2(res["baseMean"] + 1)
    res["M"] = res["log2FoldChange"]

    sig_up = (res["padj"] < 0.05) & (res["log2FoldChange"] > 0.5)
    sig_dn = (res["padj"] < 0.05) & (res["log2FoldChange"] < -0.5)
    hk = res.index.intersection(HK_GENES)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(res["A"], res["M"], alpha=0.3, s=5, color=plotting.COLOR_NS, label="NS")
    ax.scatter(res.loc[sig_up, "A"], res.loc[sig_up, "M"], alpha=0.9, s=8,
               color=plotting.COLOR_UP, label=f"Up: {sig_up.sum()}")
    ax.scatter(res.loc[sig_dn, "A"], res.loc[sig_dn, "M"], alpha=0.9, s=8,
               color=plotting.COLOR_DOWN, label=f"Down: {sig_dn.sum()}")
    if len(hk) > 0:
        ax.scatter(res.loc[hk, "A"], res.loc[hk, "M"], s=20,
                   facecolors="white", edgecolors="black", linewidth=0.8,
                   label="Housekeeping", zorder=5)

    ax.axhline(0, color="black", linewidth=1)
    ax.axhline(0.5, color="black", linestyle="--")
    ax.axhline(-0.5, color="black", linestyle="--")
    y_lim = np.ceil(max(abs(res["M"].min()), abs(res["M"].max())))
    ax.set_ylim(-y_lim, y_lim)
    ax.set_xlabel("Log\u2082 Average Expression (A)")
    if shrunken:
        ax.set_ylabel("apeglm-shrunken Log\u2082 Fold Change (M)")
    else:
        ax.set_ylabel("Log\u2082 Fold Change (M)")
    plotting.add_grid(ax)
    plotting.style_legend(ax)
    plotting.save_figure(fig, out_path)
    print(f"  MA plot -> {out_path}")


def volcano_plot(res: pd.DataFrame, out_path: str) -> None:
    """Volcano plot: -log10(padj) vs log2FC."""
    plotting.setup_style()
    res = res.copy()
    res["neglog10p"] = -np.log10(res["padj"].clip(lower=1e-300))

    sig_up = (res["padj"] < 0.05) & (res["log2FoldChange"] > 0.5)
    sig_dn = (res["padj"] < 0.05) & (res["log2FoldChange"] < -0.5)
    hk = res.index.intersection(HK_GENES)

    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.scatter(res["log2FoldChange"], res["neglog10p"], alpha=0.3, s=5,
               color=plotting.COLOR_NS, label="NS")
    ax.scatter(res.loc[sig_up, "log2FoldChange"], res.loc[sig_up, "neglog10p"],
               alpha=0.9, s=8, color=plotting.COLOR_UP, label=f"Up: {sig_up.sum()}")
    ax.scatter(res.loc[sig_dn, "log2FoldChange"], res.loc[sig_dn, "neglog10p"],
               alpha=0.9, s=8, color=plotting.COLOR_DOWN, label=f"Down: {sig_dn.sum()}")
    if len(hk) > 0:
        ax.scatter(res.loc[hk, "log2FoldChange"], res.loc[hk, "neglog10p"], s=20,
                   facecolors="white", edgecolors="black", linewidth=0.8,
                   label="Housekeeping", zorder=5)

    ax.axvline(0, color="black", linewidth=1)
    ax.axvline(0.5, color="black", linestyle="--")
    ax.axvline(-0.5, color="black", linestyle="--")
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", alpha=0.5)
    x_lim = np.ceil(max(abs(res["log2FoldChange"].min()), abs(res["log2FoldChange"].max())))
    ax.set_xlim(-x_lim, x_lim)
    ax.set_xlabel("Log\u2082 Fold Change")
    ax.set_ylabel("\u2212Log\u2081\u2080 Adjusted P-value")
    plotting.add_grid(ax)
    plotting.style_legend(ax, loc="upper right")
    plotting.save_figure(fig, out_path)
    print(f"  Volcano -> {out_path}")


def pca_plot(counts: pd.DataFrame, meta: pd.DataFrame, out_path: str,
             group_order: list = None) -> None:
    """PCA plot coloured by group."""
    plotting.setup_style()
    dds = DeseqDataSet(counts=counts, metadata=meta[["condition"]],
                       design_factors="condition")
    dds.deseq2()
    sc.tl.pca(dds)
    variance = dds.uns["pca"]["variance_ratio"][:2] * 100

    groups = group_order or sorted(meta["condition"].unique())
    color_map = {
        "Control": "#000000", "Empty": "#A9A9A9",
        "siScr": "#88ADFF", "siKRAS": "#FF7E79",
    }
    palette = {g: color_map.get(g, "#999999") for g in groups}

    sc.pl.pca(dds, color="condition", size=200, palette=palette, show=False)
    ax = plt.gca()
    ax.set_xlabel(f"PCA1 ({variance[0]:.2f}% variance)")
    ax.set_ylabel(f"PCA2 ({variance[1]:.2f}% variance)")
    ax.set_title("")
    # Restore axis ticks (scanpy hides them by default) and add grid
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.tick_params(labelsize=10)
    ax.ticklabel_format(style='scientific', axis='both', scilimits=(0, 0))
    plotting.add_grid(ax)
    # Style the legend (white background + border)
    leg = ax.get_legend()
    if leg is not None:
        leg.set_frame_on(True)
        leg.get_frame().set_facecolor('white')
        leg.get_frame().set_edgecolor('black')
    plotting.save_figure(ax.figure, out_path)
    print(f"  PCA -> {out_path}")


def heatmap(counts: pd.DataFrame, genes: list, out_path: str,
            meta: pd.DataFrame = None, cmap=None, title: str = "") -> None:
    """Z-scored heatmap of selected genes."""
    plotting.setup_style()
    genes = [g for g in genes if g in counts.columns]
    expr = np.log1p(counts[genes].T)
    expr_z = expr.subtract(expr.mean(axis=1), axis=0).div(expr.std(axis=1), axis=0)

    cmap = cmap or plotting.BWR_CMAP
    fig, ax = plt.subplots(figsize=(10, max(6, len(genes) * 0.3)))
    sns.heatmap(expr_z, cmap=cmap, ax=ax,
                cbar_kws={"label": "Z-score", "shrink": 0.5})
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_title(title)
    plotting.save_figure(fig, out_path)
    print(f"  Heatmap -> {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize DESeq2 results.")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--results", required=True, help="DESeq2 full results TSV")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--out-dir", default="Figures")
    parser.add_argument("--prefix", default="")
    args = parser.parse_args()

    ensure_dir(args.out_dir)
    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    res = pd.read_csv(args.results, sep="\t", index_col=0)
    sheet = pd.read_csv(args.samplesheet).set_index("sample")
    common = counts.index.intersection(sheet.index)
    meta = sheet.loc[common].rename(columns={"group": "condition"})

    plotting.setup_style()
    ma_plot(res, f"{args.out_dir}/{args.prefix}MA_plot.png")
    volcano_plot(res, f"{args.out_dir}/{args.prefix}volcano_plot.png")
    pca_plot(counts.loc[common], meta, f"{args.out_dir}/{args.prefix}PCA_plot.png")


if __name__ == "__main__":
    main()
