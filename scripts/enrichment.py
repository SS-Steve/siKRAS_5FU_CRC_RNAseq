"""Step 07: Enrichment analysis — GSEA (prerank) and ORA (Enrichr).

Generates result tables AND figures:
  - GSEA running enrichment plots (top N terms per gene set)
  - Enrichr barplot
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gseapy.plot import gseaplot

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils import plotting
from utils.helpers import ensure_dir


def run_prerank_gsea(results_file: str, gene_sets: list,
                     permutation_num: int = 1000, seed: int = 6,
                     out_dir: str = "Data/tables",
                     fig_dir: str = "Figures/07_enrichment",
                     top_n: int = 10,
                     story_pathways: list = None) -> pd.DataFrame:
    """
    Prerank GSEA. Genes are ranked by the DESeq2 Wald statistic.

    Generates two sets of GSEA running plots:
      - ``fig_dir/top_10/``: the top ``top_n`` terms by FDR (unbiased)
      - ``fig_dir/siKRAS_5FU/``: story-relevant pathways (if provided)
    """
    res = pd.read_csv(results_file, sep="\t", index_col=0)
    ranking = res[["stat"]].dropna().sort_values("stat", ascending=False)
    ranking.index.name = "Geneid"
    ranking = ranking[~ranking.index.duplicated(keep="first")]

    ensure_dir(out_dir)
    top5_dir = f"{fig_dir}/top_10"
    story_dir = f"{fig_dir}/siKRAS_5FU"
    ensure_dir(top5_dir)
    ensure_dir(story_dir)
    plotting.setup_style()

    all_rows = []
    for gs in gene_sets:
        print(f"[GSEA] {gs} ...")
        pre = gp.prerank(rnk=ranking, gene_sets=gs, seed=seed,
                         permutation_num=permutation_num,
                         no_plot=True, verbose=False)
        for term, d in pre.results.items():
            all_rows.append({
                "gene_set": gs, "Term": term,
                "NES": d["nes"], "fdr": d["fdr"], "es": d["es"], "pval": d["pval"],
            })

        # --- Top N by FDR (unbiased) ---
        df_gs = pd.DataFrame([
            {"Term": t, "fdr": pre.results[t]["fdr"]}
            for t in pre.results
        ]).sort_values("fdr")
        for term in df_gs.head(top_n)["Term"]:
            safe = term[:50].replace("/", "_").replace("(", "").replace(")", "").replace(" ", "_")
            d = pre.results[term]
            gseaplot(rank_metric=pre.ranking, term=term,
                     ofname=f"{top5_dir}/gsea_{safe}.png", **d)
            plt.close("all")

        # --- Story-relevant pathways ---
        if story_pathways:
            for story_name in story_pathways:
                # Match the story pathway name to an actual term
                match = None
                for term in pre.results:
                    if story_name.lower() == term.lower() or story_name.lower() in term.lower():
                        match = term
                        break
                if match is None:
                    print(f"    [story] '{story_name}' not found in {gs}")
                    continue
                safe = story_name[:50].replace("/", "_").replace("(", "").replace(")", "").replace(" ", "_")
                d = pre.results[match]
                gseaplot(rank_metric=pre.ranking, term=match,
                         ofname=f"{story_dir}/gsea_{safe}.png", **d)
                plt.close("all")

    df = pd.DataFrame(all_rows).sort_values("fdr")
    df.to_csv(f"{out_dir}/prerank_gsea_results.tsv", sep="\t", index=False)
    print(f"GSEA complete: {len(df)} terms -> {out_dir}/prerank_gsea_results.tsv")
    return df


def run_enrichr(significant_file: str, gene_sets: list,
                out_dir: str = "Data/tables",
                fig_dir: str = "Figures/07_enrichment") -> pd.DataFrame:
    """Enrichr over-representation analysis on significant DEGs (with barplot)."""
    sig = pd.read_csv(significant_file, sep="\t", index_col=0)
    gene_list = sig.index.tolist()

    enr = gp.enrichr(gene_list=gene_list, gene_sets=gene_sets,
                     organism="human", outdir=None)
    ensure_dir(out_dir)
    ensure_dir(fig_dir)
    plotting.setup_style()
    enr.results.to_csv(f"{out_dir}/enrichr_results.tsv", sep="\t", index=False)

    # Barplot
    unique_sets = enr.results["Gene_set"].unique()
    cmap = plt.get_cmap("Dark2", len(unique_sets))
    color_map = {gs: cmap(i) for i, gs in enumerate(unique_sets)}
    ax = gp.barplot(
        enr.results, column="Adjusted P-value", group="Gene_set",
        size=10, top_term=10, figsize=(8, max(6, len(gene_sets) * 2.5)),
        color=color_map,
    )
    ax.set_xlabel("\u2212Log\u2081\u2080 Adjusted P-value", fontsize=13, fontweight="bold")
    ax.grid(color="gray", alpha=0.3, axis="x")
    ax.set_axisbelow(True)
    leg = ax.get_legend()
    if leg is not None:
        leg.set_frame_on(True)
        leg.get_frame().set_facecolor('white')
        leg.get_frame().set_edgecolor('black')
    plotting.save_figure(ax.figure, f"{fig_dir}/enrichr_barplot.png")
    print(f"Enrichr complete: {len(enr.results)} terms -> {out_dir}/enrichr_results.tsv")
    return enr.results


def main() -> None:
    parser = argparse.ArgumentParser(description="GSEA + Enrichr enrichment.")
    parser.add_argument("--results", required=True, help="DESeq2 full results TSV")
    parser.add_argument("--significant", required=True, help="Significant DEGs TSV")
    parser.add_argument("--gene-sets", nargs="+",
                        default=["MSigDB_Hallmark_2020"])
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=6)
    parser.add_argument("--out-dir", default="Data/tables")
    parser.add_argument("--fig-dir", default="Figures/07_enrichment")
    parser.add_argument("--top-n", type=int, default=5)
    args = parser.parse_args()

    run_prerank_gsea(args.results, args.gene_sets, args.permutations,
                     args.seed, args.out_dir, args.fig_dir, args.top_n)
    run_enrichr(args.significant, args.gene_sets, args.out_dir, args.fig_dir)


if __name__ == "__main__":
    main()
