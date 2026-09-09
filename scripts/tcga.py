"""Step 08: TCGA validation — TPM, KM survival, Cox regression."""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils import plotting
from utils.helpers import ensure_dir, read_gtf_gene_lengths, strip_tcga_barcode


def counts_to_tpm(counts_file: str, gtf_file: str, out_file: str) -> pd.DataFrame:
    """
    Convert raw counts to log2(TPM+1) using GENCODE gene lengths.
    """
    counts = pd.read_csv(counts_file, sep="\t", index_col=0, compression="gzip")
    counts.columns = [strip_tcga_barcode(c) for c in counts.columns]
    counts = counts.loc[:, ~counts.columns.duplicated()]

    gene_lengths = read_gtf_gene_lengths(gtf_file)
    lengths = np.array([gene_lengths.get(g, np.nan) for g in counts.index])
    valid = ~np.isnan(lengths)
    counts = counts.loc[valid]
    lengths = lengths[valid]

    rpk = counts.div(lengths, axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    log2tpm = np.log2(tpm + 1)

    ensure_dir(str(Path(out_file).parent))
    log2tpm.to_csv(out_file, sep="\t", compression="gzip")
    print(f"TPM: {log2tpm.shape[0]} genes x {log2tpm.shape[1]} samples -> {out_file}")
    return log2tpm


def km_plot(gene: str, tpm: pd.DataFrame, survival: pd.DataFrame,
            out_path: str) -> float:
    """Kaplan-Meier plot with median split and log-rank test."""
    e = tpm.loc[gene].dropna()
    common = e.index.intersection(survival.index)
    e, s = e[common], survival.loc[common]
    mask = s["OS.time"].notna() & s["OS"].notna()
    e, s = e[mask], s[mask]
    events = s["OS"].astype(float).astype(int)
    median = e.median()
    high = e >= median

    plotting.setup_style()
    fig, ax = plt.subplots(figsize=(7, 6))
    kmf = KaplanMeierFitter()
    kmf.fit(s.loc[high, "OS.time"], events.loc[high], label=f"High (n={high.sum()})")
    kmf.plot(ax=ax, color=plotting.COLOR_UP, linewidth=2)
    kmf.fit(s.loc[~high, "OS.time"], events.loc[~high], label=f"Low (n={(~high).sum()})")
    kmf.plot(ax=ax, color=plotting.COLOR_DOWN, linewidth=2)

    lr = logrank_test(s.loc[high, "OS.time"], s.loc[~high, "OS.time"],
                      events.loc[high], events.loc[~high])
    ax.set_title(f"{gene}  \u2014  Log-rank adj. p={lr.p_value:.4f}")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Overall Survival")
    plotting.add_grid(ax)
    plotting.style_legend(ax)
    plotting.save_figure(fig, out_path)
    return lr.p_value


def cox_regression(gene: str, tpm: pd.DataFrame, survival: pd.DataFrame) -> dict:
    """Univariate and multivariate Cox regression (adjust age/gender/stage)."""
    e = tpm.loc[gene].dropna()
    common = e.index.intersection(survival.index)
    s = survival.loc[common]
    gz = (e[common] - e[common].mean()) / e[common].std()

    df = pd.DataFrame({
        "gene": gz, "OS.time": s["OS.time"].astype(float),
        "OS": s["OS"].astype(float).astype(int),
        "age": s["age_at_initial_pathologic_diagnosis"].astype(float),
        "gender": (s["gender"] == "MALE").astype(int),
    })

    def collapse_stage(x):
        if isinstance(x, str) and x != "[Discrepancy]":
            for r in ["IV", "III", "II", "I"]:
                if f"Stage {r}" in x:
                    return r
        return np.nan

    df["stage"] = s["ajcc_pathologic_tumor_stage"].apply(collapse_stage)
    stage_dummies = pd.get_dummies(df["stage"], prefix="stage", drop_first=True)
    df = pd.concat([df, stage_dummies], axis=1)

    # Univariate
    sub_uni = df[["gene", "OS.time", "OS"]].dropna()
    cph_uni = CoxPHFitter().fit(sub_uni, "OS.time", "OS")
    # Multivariate
    cov_cols = ["gene", "age", "gender"] + [c for c in df.columns if c.startswith("stage_")]
    sub_mv = df[cov_cols + ["OS.time", "OS"]].dropna()
    cph_mv = CoxPHFitter().fit(sub_mv, "OS.time", "OS")

    return {
        "gene": gene,
        "HR_uni": np.exp(cph_uni.params_["gene"]),
        "p_uni": cph_uni.summary.loc["gene", "p"],
        "HR_mv": np.exp(cph_mv.params_["gene"]),
        "p_mv": cph_mv.summary.loc["gene", "p"],
    }


def screen_degs(deg_file: str, tpm: pd.DataFrame, survival: pd.DataFrame,
                out_file: str) -> pd.DataFrame:
    """Systematic univariate Cox screening of all DEGs with FDR."""
    degs = pd.read_csv(deg_file, sep="\t", index_col=0)
    base = pd.DataFrame({
        "OS.time": survival["OS.time"].astype(float),
        "OS": survival["OS"].astype(float).astype(int),
    })
    rows = []
    for g in degs.index:
        if g not in tpm.index:
            continue
        e = tpm.loc[g].dropna()
        common = e.index.intersection(base.index)
        if len(common) < 50:
            continue
        gz = (e[common] - e[common].mean()) / e[common].std()
        sub = pd.DataFrame({
            "gene": gz, "OS.time": base.loc[common, "OS.time"],
            "OS": base.loc[common, "OS"]}).dropna()
        if sub["OS"].sum() < 10:
            continue
        cph = CoxPHFitter().fit(sub, "OS.time", "OS")
        rows.append({
            "gene": g, "log2FC": degs.loc[g, "log2FoldChange"],
            "HR": np.exp(cph.params_["gene"]), "cox_p": cph.summary.loc["gene", "p"],
        })
    df = pd.DataFrame(rows).sort_values("cox_p")
    df["cox_fdr"] = multipletests(df["cox_p"], method="fdr_bh")[1]
    ensure_dir(str(Path(out_file).parent))
    df.to_csv(out_file, sep="\t", index=False)
    print(f"Screening: {len(df)} genes -> {out_file}")
    return df


def plot_cox_volcano(screening_df: pd.DataFrame, out_path: str) -> None:
    """
    Volcano-style overview of the TCGA prognostic screening.
    x = log2(HR), y = -log10(cox_p). Highlights FDR-significant genes.
    """
    import matplotlib.pyplot as plt
    from utils import plotting

    plotting.setup_style()
    df = screening_df.copy()
    df["log2HR"] = np.log2(df["HR"])
    df["neglog10p"] = -np.log10(df["cox_p"].clip(lower=1e-300))

    sig = df["cox_fdr"] < 0.05
    risk = sig & (df["HR"] > 1)
    prot = sig & (df["HR"] < 1)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df["log2HR"], df["neglog10p"], s=5, alpha=0.3,
               color=plotting.COLOR_NS, label="NS")
    ax.scatter(df.loc[risk, "log2HR"], df.loc[risk, "neglog10p"],
               s=20, alpha=0.9, color=plotting.COLOR_UP, label=f"Risk (n={risk.sum()})")
    ax.scatter(df.loc[prot, "log2HR"], df.loc[prot, "neglog10p"],
               s=20, alpha=0.9, color=plotting.COLOR_DOWN,
               label=f"Protective (n={prot.sum()})")

    # Annotate FDR-significant genes
    for gene in df[sig].index:
        g = df.loc[gene]
        ax.annotate(gene, (g["log2HR"], g["neglog10p"]),
                    xytext=(5, 5), textcoords="offset points", fontsize=9,
                    fontweight="bold")

    ax.axvline(0, color="black", linewidth=1)
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", alpha=0.5)
    x_lim = np.ceil(max(abs(df["log2HR"].min()), abs(df["log2HR"].max())))
    ax.set_xlim(-x_lim, x_lim)
    ax.set_xlabel("Log₂ Hazard Ratio (protective \u2190 \u2192 risk)")
    ax.set_ylabel("\u2212Log₁\u2080 P-value")
    ax.set_title("TCGA prognostic screening (siKRAS DEGs)")
    plotting.style_legend(ax)
    plotting.add_grid(ax)
    plotting.save_figure(fig, out_path)
    print(f"Cox volcano -> {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="TCGA validation.")
    parser.add_argument("--counts", required=True, help="TCGA raw counts (gzip)")
    parser.add_argument("--survival", required=True, help="TCGA survival (gzip)")
    parser.add_argument("--gtf", required=True, help="GENCODE GTF for gene lengths")
    parser.add_argument("--deg-file", required=True, help="siKRAS vs siScr DEGs")
    parser.add_argument("--genes", nargs="+", default=["PCOLCE2", "GDF15", "HMMR"])
    parser.add_argument("--out-dir", default="Data/tables")
    parser.add_argument("--fig-dir", default="Figures")
    args = parser.parse_args()

    tpm = counts_to_tpm(args.counts, args.gtf, f"{args.out_dir}/CRC_log2tpm.tsv.gz")
    with gzip.open(args.survival, "rt") as f:
        survival = pd.read_csv(f, sep="\t", index_col=0)
    survival.index = [strip_tcga_barcode(i) for i in survival.index]

    # Cox for candidate genes
    cox_rows = []
    for g in args.genes:
        if g in tpm.index:
            r = cox_regression(g, tpm, survival)
            cox_rows.append(r)
            km_plot(g, tpm, survival, f"{args.fig_dir}/KM_{g}.png")
            print(f"  {g}: HR_mv={r['HR_mv']:.3f}, p_mv={r['p_mv']:.4f}")
    if cox_rows:
        pd.DataFrame(cox_rows).to_csv(f"{args.out_dir}/cox_results.tsv", sep="\t", index=False)

    screen_degs(args.deg_file, tpm, survival, f"{args.out_dir}/cox_screening.tsv")


if __name__ == "__main__":
    main()
