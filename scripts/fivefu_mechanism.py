"""Step 10: Supplementary analyses — DNA repair sub-pathways, 5-FU response."""

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
from scipy import stats
from gseapy.plot import gseaplot

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils import plotting
from utils.helpers import ensure_dir

def _load_gmt(gmt_file: str) -> dict:
    """Load a GMT file into a {name: [genes]} dict."""
    gmt = {}
    with open(gmt_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                gmt[parts[0]] = parts[2:]
    return gmt

def dna_repair_subpathways(results_file: str, out_file: str,
                           permutation_num: int = 1000,
                           fig_dir: str = "Figures/10_5FU_mechanism/gsea_repair",
                           gmt_file: str = "Data/gene_sets/Reactome_2022.gmt") -> pd.DataFrame:
    """
    GSEA of DNA repair sub-pathways (Reactome) using DESeq2 Wald ranking.

    Also generates GSEA running plots for the key repair sub-pathways.
    """
    res = pd.read_csv(results_file, sep="\t", index_col=0)
    ranking = res[["stat"]].dropna().sort_values("stat", ascending=False)
    ranking.index.name = "Geneid"
    ranking = ranking[~ranking.index.duplicated(keep="first")]

    pre = gp.prerank(rnk=ranking, gene_sets=gmt_file, seed=6,
                     permutation_num=permutation_num, no_plot=True, verbose=False)

    # Extract DNA-repair-related terms
    keywords = ["repair", "excision", "recombination", "nonhomologous",
                "homologous", "replication", "checkpoint", "DNA damage"]
    rows = []
    for term, d in pre.results.items():
        if any(k in term.lower() for k in keywords):
            rows.append({"Term": term, "NES": d["nes"], "fdr": d["fdr"]})
    df = pd.DataFrame(rows).sort_values("fdr")
    ensure_dir(str(Path(out_file).parent))
    df.to_csv(out_file, sep="\t", index=False)

    # GSEA running plots for the key repair sub-pathways.
    # Exact Reactome term names (incl. R-HSA IDs) to guarantee figure-text
    # consistency with the reported NES values.
    key_terms = {
        "NER": "Nucleotide Excision Repair R-HSA-5696398",
        "BER": "Base-Excision Repair, AP Site Formation R-HSA-73929",
        "NHEJ": "Nonhomologous End-Joining (NHEJ) R-HSA-5693571",
        "HR": "Defective HDR Thru Homologous Recombination (HRR) Due To BRCA1 Loss-Of-Function R-HSA-9701192",
        "DNA_Replication": "DNA Replication R-HSA-69306",
        "DNA_Repair": "DNA Repair R-HSA-73894",
    }
    plotting.setup_style()
    ensure_dir(fig_dir)
    for label, term in key_terms.items():
        if term not in pre.results:
            print(f"    [repair plot] '{label}' ({term}) not found")
            continue
        d = pre.results[term]
        gseaplot(rank_metric=pre.ranking, term=term,
                 ofname=f"{fig_dir}/gsea_{label}.png", **d)
        plt.close("all")

    print(f"DNA repair sub-pathways: {len(df)} terms -> {out_file}")
    return df

def fivefu_response_overlap(results_file: str, out_file: str,
                           fig_dir: str = "Figures/10_5FU_mechanism/gsea_5fu_response",
                           gmt_file: str = "Data/gene_sets/LINCS_L1000_Chem_Pert_Consensus_Sigs.gmt") -> pd.DataFrame:
    """
    Compare siKRAS signature with LINCS L1000 5-FU response signatures.

    Generates GSEA running plots for Fluorouracil Up and Down.
    """
    res = pd.read_csv(results_file, sep="\t", index_col=0)
    ranking = res[["stat"]].dropna().sort_values("stat", ascending=False)
    ranking.index.name = "Geneid"
    ranking = ranking[~ranking.index.duplicated(keep="first")]

    gs = _load_gmt(gmt_file)
    fu_sets = {g: genes for g, genes in gs.items() if "Fluorouracil" in g}

    plotting.setup_style()
    ensure_dir(fig_dir)
    rows = []
    for name, genes in fu_sets.items():
        pre = gp.prerank(rnk=ranking, gene_sets={name: genes}, seed=6,
                         permutation_num=1000, no_plot=True, verbose=False)
        if name in pre.results:
            d = pre.results[name]
            rows.append({"geneset": name, "NES": d["nes"], "fdr": d["fdr"]})
            safe = name.replace(" ", "_")
            gseaplot(rank_metric=pre.ranking, term=name,
                     ofname=f"{fig_dir}/gsea_{safe}.png", **d)
            plt.close("all")
    df = pd.DataFrame(rows)
    ensure_dir(str(Path(out_file).parent))
    df.to_csv(out_file, sep="\t", index=False)
    print(f"5-FU response overlap -> {out_file}")
    return df

def fivefu_resistance_reversal(results_file: str, out_file: str,
                               counts_file: str = "Data/5FU_resistance/GSE153412_est_counts.tsv.gz") -> tuple:
    """
    Check if siKRAS reverses the 5-FU resistance signature.

    Resistance signature is computed from GSE153412 HCT116 (resistant vs
    sensitive, 3v3) using the SAME method as our own data:
    DESeq2 (Wald test) + lfcShrink (apeglm).

    Returns
    -------
    tuple
        (summary_df, resistance_results, siKRAS_results)
    """
    import gzip
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # --- 1. Load GSE153412 HCT116 (resistant vs sensitive, 3v3) ---
    with gzip.open(counts_file, 'rt') as f:
        counts = pd.read_csv(f, sep='\t', index_col=0)
    res_cols = ['116RU1', '116RU2', '116RU3']   # HCT116 resistant untreated
    sens_cols = ['116SU1', '116SU2', '116SU3']  # HCT116 sensitive untreated
    sub = counts[res_cols + sens_cols].round().astype(int)
    meta = pd.DataFrame({'condition': ['resistant']*3 + ['sensitive']*3},
                        index=sub.columns)

    dds = DeseqDataSet(counts=sub.T, metadata=meta, design_factors='condition')
    dds.deseq2()
    res = DeseqStats(dds, contrast=['condition', 'resistant', 'sensitive'])
    res.summary()
    coeff = [c for c in res.LFC.columns if 'resistant' in c and 'T.' in c]
    if coeff:
        res.lfc_shrink(coeff=coeff[0])
    resistance_df = res.results_df  # resistance log2FC (shrunk) + padj

    # --- 2. siKRAS results ---
    siKRAS_df = pd.read_csv(results_file, sep='\t', index_col=0)

    # --- 3. Reversal analysis ---
    res_up = resistance_df[
        (resistance_df.padj < 0.05) & (resistance_df.log2FoldChange > 0.5)]
    res_dn = resistance_df[
        (resistance_df.padj < 0.05) & (resistance_df.log2FoldChange < -0.5)]

    up_in = siKRAS_df.reindex(res_up.index).dropna(subset=['log2FoldChange'])
    dn_in = siKRAS_df.reindex(res_dn.index).dropna(subset=['log2FoldChange'])

    rows = []
    if len(up_in) > 10:
        t, p = stats.ttest_1samp(up_in['log2FoldChange'], 0)
        rows.append({'signature': 'resistance_up',
                     'mean_log2FC_siKRAS': up_in['log2FoldChange'].mean(),
                     'reversal_pct': (up_in['log2FoldChange'] < 0).mean() * 100,
                     'p_value': p})
    if len(dn_in) > 10:
        t, p = stats.ttest_1samp(dn_in['log2FoldChange'], 0)
        rows.append({'signature': 'resistance_down',
                     'mean_log2FC_siKRAS': dn_in['log2FoldChange'].mean(),
                     'reversal_pct': (dn_in['log2FoldChange'] > 0).mean() * 100,
                     'p_value': p})

    df = pd.DataFrame(rows)
    ensure_dir(str(Path(out_file).parent))
    df.to_csv(out_file, sep='\t', index=False)
    print(f"5-FU resistance reversal -> {out_file}")
    return df, resistance_df, siKRAS_df

def plot_resistance_scatter(resistance_df: pd.DataFrame, siKRAS_df: pd.DataFrame,
                           out_path: str) -> None:
    """
    Scatter: 5-FU resistance log2FC vs siKRAS log2FC.
    Both axes are DESeq2 + lfcShrink results (same method).
    Coloured by concordance between the two signatures.
    """
    plotting.setup_style()
    common = resistance_df.index.intersection(siKRAS_df.index)
    x = resistance_df.loc[common, "log2FoldChange"]
    y = siKRAS_df.loc[common, "log2FoldChange"]
    res_padj = resistance_df.loc[common, "padj"]
    sik_padj = siKRAS_df.loc[common, "padj"]

    # Significance masks (both sides: padj < 0.05, |log2FC| > 0.5)
    res_up = (res_padj < 0.05) & (x > 0.5)
    res_dn = (res_padj < 0.05) & (x < -0.5)
    sik_up = (sik_padj < 0.05) & (y > 0.5)
    sik_dn = (sik_padj < 0.05) & (y < -0.5)

    both_up = res_up & sik_up
    both_dn = res_dn & sik_dn
    discord = (res_up & sik_dn) | (res_dn & sik_up)  # opposite directions

    # Housekeeping genes (marked like the MA plot)
    HK_GENES = {
        "ACTB", "GAPDH", "B2M", "HPRT1", "PPIA", "RPLP0", "TBP", "UBC",
        "YWHAZ", "PGK1", "RPL13A", "GUSB", "SDHA", "HMBS",
        "TFRC", "RPL19", "IPO8", "CYC1", "EIF4A2",
    }
    hk = common.intersection(HK_GENES)

    fig, ax = plt.subplots(figsize=(7, 7))
    # Background: ALL genes (grey)
    ax.scatter(x, y, s=5, alpha=0.3, color=plotting.COLOR_NS, label="NS")
    # Overlay significant categories
    ax.scatter(x[both_up], y[both_up], s=8, alpha=0.9,
               color=plotting.COLOR_UP, label="Both up")
    ax.scatter(x[both_dn], y[both_dn], s=8, alpha=0.9,
               color=plotting.COLOR_DOWN, label="Both down")
    ax.scatter(x[discord], y[discord], s=8, alpha=0.9,
               color="#bb99ff", label="Reversal candidates")
    # Housekeeping genes (white fill, black edge)
    if len(hk) > 0:
        ax.scatter(x[hk], y[hk], s=20, facecolors="white",
                   edgecolors="black", linewidth=0.8,
                   label="Housekeeping", zorder=5)

    ax.axhline(0, color="black", linewidth=1)
    ax.axvline(0, color="black", linewidth=1)
    # Dashed lines for BOTH significances (black, like the MA plot)
    for v in [0.5, -0.5]:
        ax.axhline(v, color="black", linestyle="--")
        ax.axvline(v, color="black", linestyle="--")
    ax.set_xlabel("Log2FC (5-FU resistant vs sensitive, HCT116)")
    ax.set_ylabel("Log2FC (siKRAS vs siScr)")
    ax.set_title("5-FU resistance reversal")
    plotting.style_legend(ax, loc="upper left")
    plotting.add_grid(ax)
    plotting.save_figure(fig, out_path)
    print(f"Resistance scatter -> {out_path}")

def plot_pearson_correlation(resistance_df: pd.DataFrame, siKRAS_df: pd.DataFrame,
                             out_path: str) -> None:
    """
    Pearson correlation between resistance and siKRAS log2FC, computed for
    three groups: all genes, significant genes (union), and double-significant
    genes (significant in both datasets). Three fitted lines are drawn.
    """
    from scipy.stats import pearsonr, linregress
    plotting.setup_style()

    common = resistance_df.index.intersection(siKRAS_df.index)
    x = resistance_df.loc[common, "log2FoldChange"]
    y = siKRAS_df.loc[common, "log2FoldChange"]
    res_padj = resistance_df.loc[common, "padj"]
    sik_padj = siKRAS_df.loc[common, "padj"]

    # Drop genes with missing log2FC on either side
    mask = x.notna() & y.notna()
    x, y = x[mask], y[mask]
    res_padj, sik_padj = res_padj[mask], sik_padj[mask]

    res_sig = (res_padj < 0.05) & (abs(x) > 0.5)
    sik_sig = (sik_padj < 0.05) & (abs(y) > 0.5)
    union = res_sig | sik_sig
    both = res_sig & sik_sig

    # Pearson for the three groups
    rho_all, p_all = pearsonr(x, y)
    rho_union, p_union = pearsonr(x[union], y[union])
    rho_both, p_both = pearsonr(x[both], y[both])

    fig, ax = plt.subplots(figsize=(7, 7))
    # Scatter: non-significant -> union-only -> double-significant
    ax.scatter(x[~union], y[~union], s=6, alpha=0.12, color=plotting.COLOR_NS,
               label=f"Non-significant (n={(~union).sum()})")
    union_only = union & ~both
    ax.scatter(x[union_only], y[union_only], s=12, alpha=0.5, color="#888888",
               label=f"Significant (n={int(union.sum())})")
    ax.scatter(x[both], y[both], s=22, alpha=0.85, color="#333333",
               label=f"Double-significant (n={int(both.sum())})")

    # Three fitted lines
    for m, color, rho, lbl in [
        (None, plotting.COLOR_UP, rho_all, f"All genes (r={rho_all:.3f})"),
        (union, plotting.COLOR_DOWN, rho_union, f"Significant (r={rho_union:.3f})"),
        (both, "#1B7837", rho_both, f"Double-significant (r={rho_both:.3f})"),
    ]:
        xx = x if m is None else x[m]
        yy = y if m is None else y[m]
        slope, intercept, _, _, _ = linregress(xx, yy)
        xs = np.linspace(xx.min(), xx.max(), 100)
        ax.plot(xs, slope * xs + intercept, color=color, linewidth=3, label=lbl)

    ax.axhline(0, color="black", linewidth=1, linestyle="--")
    ax.axvline(0, color="black", linewidth=1, linestyle="--")
    ax.set_xlabel("Log2FC (5-FU resistant vs sensitive, HCT116)")
    ax.set_ylabel("Log2FC (siKRAS vs siScr)")
    ax.set_title(f"Pearson r: all={rho_all:.3f}, sig={rho_union:.3f}, double={rho_both:.3f}")
    plotting.style_legend(ax)
    plotting.add_grid(ax)
    plotting.save_figure(fig, out_path)
    print(f"Pearson correlation -> {out_path}")

def main() -> None:
    parser = argparse.ArgumentParser(description="5-FU mechanism analyses.")
    parser.add_argument("--results", required=True, help="siKRAS vs siScr DESeq2 results")
    parser.add_argument("--resistance-counts",
                        default="Data/5FU_resistance/GSE153412_est_counts.tsv.gz",
                        help="GSE153412 estimated counts (optional)")
    parser.add_argument("--out-dir", default="Data/tables")
    args = parser.parse_args()

    dna_repair_subpathways(args.results, f"{args.out_dir}/dna_repair_subpathways.tsv")
    fivefu_response_overlap(args.results, f"{args.out_dir}/5fu_response_overlap.tsv")
    fivefu_resistance_reversal(args.results,
                               f"{args.out_dir}/5fu_resistance_reversal.tsv",
                               counts_file=args.resistance_counts)

if __name__ == "__main__":
    main()
