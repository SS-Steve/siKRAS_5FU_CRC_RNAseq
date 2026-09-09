"""Step 05: DESeq2 differential expression with lfcShrink (apeglm)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.helpers import ensure_dir


def run_deseq2(counts_file: str, samplesheet: str, contrasts: list,
               significance: tuple = (10, 0.05, 0.5),
               out_dir: str = "Data/tables",
               dispersion_plot: str | None = None) -> dict:
    """
    Run DESeq2 with lfcShrink for each contrast.

    Parameters
    ----------
    counts_file : str
        Processed count matrix (samples x genes).
    samplesheet : str
        Sample metadata CSV.
    contrasts : list
        List of [factor, numerator, denominator] contrasts.
    significance : tuple
        (baseMean, padj, |log2FC|) thresholds.
    out_dir : str
        Output directory for results tables.
    dispersion_plot : str, optional
        If provided, save a DESeq2 dispersion plot to this path.

    Returns
    -------
    dict
        Mapping of contrast label -> (dds, results_df, significant_df).
    """
    base_mean_thr, padj_thr, lfc_thr = significance

    counts = pd.read_csv(counts_file, sep="\t", index_col=0)
    sheet = pd.read_csv(samplesheet).set_index("sample")
    common = counts.index.intersection(sheet.index)
    counts = counts.loc[common]
    meta = sheet.loc[common]

    # Reference level: set the FIRST denominator as reference so that
    # lfcShrink coefficient naming is predictable.
    all_levels = sorted(meta["group"].unique())
    meta = meta.copy()
    meta["group"] = pd.Categorical(meta["group"], categories=all_levels)

    dds = DeseqDataSet(counts=counts, metadata=meta[["group"]],
                       design_factors="group")
    dds.deseq2()

    # Optional: save the dispersion plot (QC of mean-variance modelling).
    if dispersion_plot:
        Path(dispersion_plot).parent.mkdir(parents=True, exist_ok=True)
        dds.plot_dispersions(save_path=dispersion_plot)
        print(f"Dispersion plot saved: {dispersion_plot}")

    results = {}
    for contrast in contrasts:
        factor, num, denom = contrast
        label = f"{num}_vs_{denom}"
        print(f"\n=== Contrast: {label} ===")

        # Relevel so the denominator is reference (enables lfcShrink coeff).
        meta_rl = meta.copy()
        meta_rl["group"] = pd.Categorical(
            meta_rl["group"], categories=[denom] + [g for g in all_levels if g != denom]
        )
        dds_rl = DeseqDataSet(counts=counts, metadata=meta_rl[["group"]],
                              design_factors="group")
        dds_rl.deseq2()

        res = DeseqStats(dds_rl, contrast=["group", num, denom])
        res.summary()

        # Save unshrunk results (before lfcShrink)
        unshrunk_df = res.results_df.copy()
        ensure_dir(out_dir)
        unshrunk_df.to_csv(f"{out_dir}/deseq2_results_unshrunk_{label}.tsv", sep="\t")

        # lfcShrink (apeglm) on the numerator coefficient
        coeff = [c for c in res.LFC.columns if num in c and "T." in c]
        if coeff:
            res.lfc_shrink(coeff=coeff[0])

        res_df = res.results_df
        sig = res_df[
            (res_df.baseMean >= base_mean_thr)
            & (res_df.padj < padj_thr)
            & (abs(res_df.log2FoldChange) > lfc_thr)
        ]

        res_df.to_csv(f"{out_dir}/deseq2_results_{label}.tsv", sep="\t")
        sig.to_csv(f"{out_dir}/deseq2_significant_{label}.tsv", sep="\t")

        n_up = (sig.log2FoldChange > 0).sum()
        n_down = (sig.log2FoldChange < 0).sum()
        print(f"  DEGs: {len(sig)} ({n_up} up, {n_down} down)")
        results[label] = (dds_rl, res_df, sig)

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="DESeq2 + lfcShrink.")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--contrasts", nargs="+", required=True,
                        help="Each contrast as 'num,denom' (factor is 'group').")
    parser.add_argument("--baseMean", type=float, default=10)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--log2FC", type=float, default=0.5)
    parser.add_argument("--out-dir", default="Data/tables")
    args = parser.parse_args()

    contrasts = [["group"] + c.split(",") for c in args.contrasts]
    run_deseq2(args.counts, args.samplesheet, contrasts,
               (args.baseMean, args.padj, args.log2FC), args.out_dir)


if __name__ == "__main__":
    main()
