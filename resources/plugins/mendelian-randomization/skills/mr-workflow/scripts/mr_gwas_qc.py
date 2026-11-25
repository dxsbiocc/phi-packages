#!/usr/bin/env python3
"""mr_gwas_qc: Standardize GWAS/QTL summary statistics, compute genomic inflation
factor lambda_GC, extract physical-window-clumped lead variants, and render
Manhattan + QQ plots.
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import clump_by_distance, emit, fail, read_summary_stats, write_artifact


def chrom_sort_key(ch: str) -> tuple[int, str]:
    s = str(ch).upper()
    if s.isdigit():
        return (int(s), "")
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (order.get(s, 99), s)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--p_cutoff", type=float, default=5e-8)
    parser.add_argument("--clump_kb", type=int, default=10000)
    parser.add_argument("--min_maf", type=float, default=0.01)
    parser.add_argument("--sample_size", type=int, default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    # Count raw lines quickly
    raw_df_peek = pd.read_csv(args.gwas, sep="," if args.gwas.lower().endswith((".csv", ".csv.gz")) else "\t")
    n_input = int(len(raw_df_peek))

    clean_df, scale = read_summary_stats(args.gwas, default_n=args.sample_size)
    if args.min_maf > 0 and clean_df["EAF"].notna().any():
        maf = np.minimum(clean_df["EAF"], 1.0 - clean_df["EAF"])
        clean_df = clean_df[clean_df["EAF"].isna() | (maf >= args.min_maf)].reset_index(drop=True)

    if clean_df.empty:
        fail("0 valid biallelic SNP records remained after GWAS QC")

    z = clean_df["BETA"].to_numpy(dtype=float) / clean_df["SE"].to_numpy(dtype=float)
    chi2 = z * z
    lambda_gc = float(np.median(chi2) / 0.454936423119572)

    gw_sig = int((clean_df["P"] <= args.p_cutoff).sum())
    leads = clump_by_distance(clean_df, args.p_cutoff, args.clump_kb)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    clean_path = outdir / "gwas_clean.tsv"
    clean_df.to_csv(clean_path, sep="\t", index=False)

    leads_path = outdir / "lead_instruments.tsv"
    leads.to_csv(leads_path, sep="\t", index=False)

    # Render 2-panel Manhattan + QQ plot
    fig, (ax_m, ax_q) = plt.subplots(1, 2, figsize=(12.0, 4.5), gridspec_kw={"width_ratios": [2.2, 1.0]})

    has_coords = clean_df["CHR"].notna().any() and clean_df["POS"].notna().any()
    nlp_all = -np.log10(np.clip(clean_df["P"].to_numpy(dtype=float), 1e-300, 1.0))
    if has_coords:
        plot_df = clean_df.dropna(subset=["CHR", "POS"]).copy()
        # Downsample non-significant background points if > 50k for fast crisp plotting
        sig_part = plot_df[plot_df["P"] <= 1e-3]
        bg_part = plot_df[plot_df["P"] > 1e-3]
        if len(bg_part) > 40000:
            bg_part = bg_part.sample(n=40000, random_state=42)
        plot_df = pd.concat([sig_part, bg_part], ignore_index=True)
        chroms = sorted(plot_df["CHR"].astype(str).unique(), key=chrom_sort_key)
        offset = 0.0
        ticks = []
        tick_lbls = []
        colors = ("#1E3A8A", "#38BDF8")
        for idx_c, ch in enumerate(chroms):
            sub = plot_df[plot_df["CHR"].astype(str) == ch]
            pos = sub["POS"].to_numpy(dtype=float)
            p_nlp = -np.log10(np.clip(sub["P"].to_numpy(dtype=float), 1e-300, 1.0))
            x_cum = pos + offset
            ax_m.scatter(x_cum, p_nlp, s=9, color=colors[idx_c % 2], alpha=0.75, linewidths=0)
            ticks.append(offset + 0.5 * (float(pos.min()) + float(pos.max())))
            tick_lbls.append(ch)
            offset += float(pos.max()) + 2e6
        ax_m.set_xticks(ticks[:24])
        ax_m.set_xticklabels(tick_lbls[:24], fontsize=7.5)
        ax_m.set_xlabel("Chromosome")
    else:
        idx_x = np.arange(len(nlp_all))
        ax_m.scatter(idx_x, nlp_all, s=10, color="#2563EB", alpha=0.7)
        ax_m.set_xlabel("Variant Index")

    ax_m.axhline(-np.log10(args.p_cutoff), color="#E11D48", linestyle="--", linewidth=1.0, label=f"P = {args.p_cutoff:.0e}")
    ax_m.set_ylabel("−log10(P)")
    ax_m.set_title(f"Manhattan Plot ({len(leads)} clumped leads, λ_GC = {lambda_gc:.3f})", fontsize=10)
    ax_m.legend(frameon=False, fontsize=8, loc="upper right")
    ax_m.spines["top"].set_visible(False)
    ax_m.spines["right"].set_visible(False)

    # QQ plot
    obs_p = np.sort(np.clip(clean_df["P"].to_numpy(dtype=float), 1e-300, 1.0))
    n_p = len(obs_p)
    if n_p > 50000:
        # Keep top 5000 smallest P + uniform quantile grid
        keep_idx = np.unique(np.concatenate([np.arange(5000), np.linspace(5000, n_p - 1, 10000, dtype=int)]))
        exp_nlp = -np.log10((keep_idx + 0.5) / n_p)
        obs_nlp = -np.log10(obs_p[keep_idx])
    else:
        exp_nlp = -np.log10((np.arange(n_p) + 0.5) / n_p)
        obs_nlp = -np.log10(obs_p)

    lim_max = max(float(exp_nlp.max()), float(obs_nlp.max())) * 1.05
    ax_q.plot([0, float(exp_nlp.max())], [0, float(exp_nlp.max())], color="#94A3B8", linestyle="--", linewidth=1.0)
    ax_q.scatter(exp_nlp, obs_nlp, s=12, color="#0F172A", alpha=0.75)
    ax_q.set_xlabel("Expected −log10(P)")
    ax_q.set_ylabel("Observed −log10(P)")
    ax_q.set_title(f"QQ Plot (λ_GC = {lambda_gc:.3f})", fontsize=10)
    ax_q.spines["top"].set_visible(False)
    ax_q.spines["right"].set_visible(False)
    fig.tight_layout()
    fig_path = outdir / "gwas-manhattan-qq.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {"p_cutoff": args.p_cutoff, "clump_kb": args.clump_kb, "min_maf": args.min_maf}
    artifacts = [
        write_artifact(fig_path, "mr_gwas_qc", f"GWAS Manhattan & QQ plot (λ_GC={lambda_gc:.3f})", inputs=(args.gwas,), parameters=params_meta),
        write_artifact(leads_path, "mr_gwas_qc", f"Distance-clumped lead variants (n={len(leads)})", inputs=(args.gwas,), parameters=params_meta, shape=leads.shape),
        write_artifact(clean_path, "mr_gwas_qc", "Standardized GWAS summary statistics", inputs=(args.gwas,), parameters=params_meta, shape=clean_df.shape),
    ]

    warnings: list[str] = [
        f"lead loci clumped via physical window ({args.clump_kb} kb) without an external 1000G PLINK LD r^2 reference panel"
    ]
    if lambda_gc > 1.20:
        warnings.append(f"elevated genomic inflation factor lambda_GC={lambda_gc:.3f} (> 1.20); check population stratification vs polygenicity")
    if len(leads) == 0:
        warnings.append(f"0 variants passed p_cutoff={args.p_cutoff}; check statistical power")

    emit(
        {
            "tool": "mr_gwas_qc",
            "n_variants_input": n_input,
            "n_variants_clean": int(len(clean_df)),
            "effect_scale_input": scale,
            "lambda_gc": round(lambda_gc, 4),
            "n_genome_wide_sig": gw_sig,
            "n_lead_loci": int(len(leads)),
            "top_leads": leads.head(15).to_dict(orient="records"),
            "outputs": {
                "clean_gwas": str(clean_path),
                "lead_instruments": str(leads_path),
            },
            "figures": [str(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
