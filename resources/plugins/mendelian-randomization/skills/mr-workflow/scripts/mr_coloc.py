#!/usr/bin/env python3
"""mr_coloc: Exact single-causal-variant Bayesian colocalization (Wakefield 2009
Approximate Bayes Factor / Giambartolomei et al. 2014 coloc.abf) across a shared
genomic locus, outputting PP.H0-PP.H4, conditional SNP.PP.H4, 95% credible set,
and a regional mirror plot.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, read_summary_stats, write_artifact


def wakefield_log_abf(beta: np.ndarray, se: np.ndarray, prior_sd: float) -> np.ndarray:
    """Wakefield (2009) log Approximate Bayes Factor (natural log) for H1 vs H0."""
    v = se * se
    w = prior_sd * prior_sd
    r = w / (v + w)
    z = beta / se
    return 0.5 * np.log(1.0 - r) + 0.5 * r * (z * z)


def filter_region(df: pd.DataFrame, chrom: str | None, start: int | None, end: int | None) -> pd.DataFrame:
    out = df
    if chrom is not None and out["CHR"].notna().any():
        ch_clean = str(chrom).replace("chr", "").replace("CHR", "")
        out = out[out["CHR"].astype(str) == ch_clean]
    if start is not None and out["POS"].notna().any():
        out = out[out["POS"] >= start]
    if end is not None and out["POS"].notna().any():
        out = out[out["POS"] <= end]
    return out.reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trait1", required=True)
    parser.add_argument("--trait2", required=True)
    parser.add_argument("--chrom", default=None)
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)
    parser.add_argument("--trait1_type", choices=["quant", "cc"], default="quant")
    parser.add_argument("--trait2_type", choices=["quant", "cc"], default="cc")
    parser.add_argument("--p1", type=float, default=1e-4)
    parser.add_argument("--p2", type=float, default=1e-4)
    parser.add_argument("--p12", type=float, default=1e-5)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    df1, _ = read_summary_stats(args.trait1)
    df2, _ = read_summary_stats(args.trait2)
    df1 = filter_region(df1, args.chrom, args.start, args.end)
    df2 = filter_region(df2, args.chrom, args.start, args.end)

    shared_snps = [s for s in df1["SNP"] if s in set(df2["SNP"])]
    if len(shared_snps) < 10:
        fail(f"mr_coloc requires >= 10 shared SNPs in the locus; found {len(shared_snps)}")

    s1 = df1.set_index("SNP").loc[shared_snps]
    s2 = df2.set_index("SNP").loc[shared_snps]

    # Standard coloc.abf prior SDs: 0.15 for SDY=1 standardized quantitative traits, 0.20 for binary log-OR
    sd1 = 0.15 if args.trait1_type == "quant" else 0.20
    sd2 = 0.15 if args.trait2_type == "quant" else 0.20

    l1 = wakefield_log_abf(s1["BETA"].to_numpy(dtype=float), s1["SE"].to_numpy(dtype=float), sd1)
    l2 = wakefield_log_abf(s2["BETA"].to_numpy(dtype=float), s2["SE"].to_numpy(dtype=float), sd2)
    lsum = l1 + l2

    # Log unnormalized posterior scores for H0..H4 (Giambartolomei et al., PLoS Genet 2014 Eq. 2)
    lh0 = 0.0
    lh1 = np.log(args.p1) + float(logsumexp(l1))
    lh2 = np.log(args.p2) + float(logsumexp(l2))
    lh4 = np.log(args.p12) + float(logsumexp(lsum))

    # H3: sum_{i != j} ABF1_i * ABF2_j = (sum ABF1)*(sum ABF2) - sum(ABF1_i * ABF2_i)
    log_prod_all = float(logsumexp(l1) + logsumexp(l2))
    log_diag = float(logsumexp(lsum))
    if log_prod_all > log_diag + 1e-12:
        log_offdiag = log_prod_all + np.log1p(-np.exp(log_diag - log_prod_all))
        lh3 = np.log(args.p1) + np.log(args.p2) + float(log_offdiag)
    else:
        lh3 = -1e30

    all_lh = np.array([lh0, lh1, lh2, lh3, lh4], dtype=float)
    denom = float(logsumexp(all_lh))
    pps = np.exp(all_lh - denom)
    posteriors = {f"PP.H{i}": round(float(pps[i]), 5) for i in range(5)}

    # Conditional per-SNP posterior probability under H4 + 95% credible set
    snp_pp_h4 = np.exp(lsum - float(logsumexp(lsum)))
    snp_df = pd.DataFrame(
        {
            "SNP": shared_snps,
            "CHR": s1["CHR"].to_numpy(),
            "POS": s1["POS"].to_numpy(),
            "beta_trait1": np.round(s1["BETA"].to_numpy(dtype=float), 5),
            "se_trait1": np.round(s1["SE"].to_numpy(dtype=float), 5),
            "pval_trait1": s1["P"].to_numpy(dtype=float),
            "beta_trait2": np.round(s2["BETA"].to_numpy(dtype=float), 5),
            "se_trait2": np.round(s2["SE"].to_numpy(dtype=float), 5),
            "pval_trait2": s2["P"].to_numpy(dtype=float),
            "log_abf1": np.round(l1, 4),
            "log_abf2": np.round(l2, 4),
            "SNP.PP.H4": np.round(snp_pp_h4, 6),
        }
    ).sort_values("SNP.PP.H4", ascending=False).reset_index(drop=True)

    cum_pp = np.cumsum(snp_df["SNP.PP.H4"].to_numpy(dtype=float))
    n_cs = int(np.searchsorted(cum_pp, 0.95, side="left")) + 1
    snp_df["in_credible_set_95"] = np.arange(len(snp_df)) < n_cs

    pph4 = float(pps[4])
    pph3 = float(pps[3])
    if pph4 >= 0.80:
        call = "colocalized_H4_strong"
    elif pph4 >= 0.50:
        call = "colocalized_H4_moderate"
    elif pph3 >= 0.80:
        call = "distinct_causal_variants_H3"
    else:
        call = "inconclusive_or_single_trait"

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    snps_path = outdir / "coloc_snps.tsv"
    snp_df.to_csv(snps_path, sep="\t", index=False)

    # Regional mirror plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8.0, 7.2), sharex=True)
    pos_arr = snp_df["POS"].to_numpy(dtype=float)
    x_coord = pos_arr / 1e6 if np.isfinite(pos_arr).all() else np.arange(len(snp_df), dtype=float)
    x_lbl = f"Position on chr{snp_df['CHR'].iloc[0]} (Mb)" if np.isfinite(pos_arr).all() else "Variant Rank"

    nlp1 = -np.log10(np.clip(snp_df["pval_trait1"].to_numpy(dtype=float), 1e-300, 1.0))
    nlp2 = -np.log10(np.clip(snp_df["pval_trait2"].to_numpy(dtype=float), 1e-300, 1.0))
    ax1.scatter(x_coord, nlp1, s=20, color="#2563EB", alpha=0.82)
    ax1.scatter([x_coord[0]], [nlp1[0]], s=55, color="#E11D48", zorder=5, label=f"Lead H4: {snp_df['SNP'].iloc[0]}")
    ax1.set_ylabel("Trait 1 −log10(P)")
    ax1.set_title(f"Bayesian Colocalization: PP.H4 = {pph4:.3f}, PP.H3 = {pph3:.3f} ({call})", fontsize=10)
    ax1.legend(frameon=False, fontsize=8, loc="upper right")

    ax2.scatter(x_coord, nlp2, s=20, color="#10B981", alpha=0.82)
    ax2.scatter([x_coord[0]], [nlp2[0]], s=55, color="#E11D48", zorder=5)
    ax2.set_ylabel("Trait 2 −log10(P)")

    cs_mask = snp_df["in_credible_set_95"].to_numpy(dtype=bool)
    ax3.scatter(x_coord[~cs_mask], snp_df["SNP.PP.H4"].to_numpy()[~cs_mask], s=18, color="#94A3B8", alpha=0.7)
    ax3.scatter(x_coord[cs_mask], snp_df["SNP.PP.H4"].to_numpy()[cs_mask], s=32, color="#E11D48", alpha=0.9, label=f"95% Credible Set (n={n_cs})")
    ax3.set_ylabel("SNP.PP.H4")
    ax3.set_xlabel(x_lbl)
    ax3.legend(frameon=False, fontsize=8, loc="upper right")
    for ax in (ax1, ax2, ax3):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig_path = outdir / "coloc-locus.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {"trait1_type": args.trait1_type, "trait2_type": args.trait2_type, "p1": args.p1, "p2": args.p2, "p12": args.p12}
    artifacts = [
        write_artifact(fig_path, "mr_coloc", f"Regional colocalization plot (PP.H4={pph4:.3f})", inputs=(args.trait1, args.trait2), parameters=params_meta),
        write_artifact(snps_path, "mr_coloc", "Per-SNP Wakefield ABF and SNP.PP.H4 table", inputs=(args.trait1, args.trait2), parameters=params_meta, shape=snp_df.shape),
    ]

    warnings: list[str] = []
    if len(shared_snps) < 50:
        warnings.append(f"locus contains only {len(shared_snps)} shared SNPs; verify dense unfiltered summary stats were passed rather than pre-clumped lead SNPs")

    emit(
        {
            "tool": "mr_coloc",
            "n_shared_snps": len(shared_snps),
            "priors": {"p1": args.p1, "p2": args.p2, "p12": args.p12, "sd1": sd1, "sd2": sd2},
            "posteriors": posteriors,
            "coloc_call": call,
            "top_snp": {
                "SNP": str(snp_df["SNP"].iloc[0]),
                "SNP.PP.H4": float(snp_df["SNP.PP.H4"].iloc[0]),
            },
            "credible_set_95": snp_df.head(min(n_cs, 25))[["SNP", "SNP.PP.H4"]].to_dict(orient="records"),
            "outputs": {"snps_table": str(snps_path)},
            "figures": [str(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
