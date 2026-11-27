#!/usr/bin/env python3
"""gwas_sumstats_harmonize: Standardize PLINK2, REGENIE, SAIGE, BOLT-LMM, or GWAS
Catalog summary statistics, enforce IUPAC alleles and SE > 0, filter MAF/INFO,
compute lambda_GC and lambda_1000, extract independent lead loci, and render a
2-panel Manhattan + QQ figure.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    chrom_sort_key,
    clump_by_distance,
    emit,
    fail,
    is_valid_iupac_allele,
    parse_sumstats,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sumstats", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--n_cases", type=int, default=None)
    ap.add_argument("--n_controls", type=int, default=None)
    ap.add_argument("--min_maf", type=float, default=0.01)
    ap.add_argument("--min_info", type=float, default=0.80)
    ap.add_argument("--lead_window_kb", type=int, default=500)
    ap.add_argument("--p_threshold", type=float, default=5e-8)
    args = ap.parse_args()

    default_n: float | None = None
    n_eff_cc: float | None = None
    if args.n_cases is not None and args.n_controls is not None:
        if args.n_cases <= 0 or args.n_controls <= 0:
            fail("--n_cases and --n_controls must be positive integers")
        default_n = float(args.n_cases + args.n_controls)
        n_eff_cc = float(4.0 * args.n_cases * args.n_controls / (args.n_cases + args.n_controls))

    df, meta = parse_sumstats(
        args.sumstats,
        require_alleles=True,
        require_beta_se=True,
        require_n=False,
        default_n=default_n,
    )
    n_input = int(meta["n_input"])

    if n_eff_cc is not None:
        df["N"] = df["N"].where(df["N"].notna() & (df["N"] > 0), default_n)
        df["N_EFF"] = n_eff_cc

    # 1) Allele & summary-statistic validity
    ea_ok = df["EA"].map(is_valid_iupac_allele)
    nea_ok = df["NEA"].map(is_valid_iupac_allele)
    beta_arr = df["BETA"].to_numpy(dtype=float, copy=True)
    se_arr = df["SE"].to_numpy(dtype=float, copy=True)
    z_arr = df["Z"].to_numpy(dtype=float, copy=True)
    p_arr = df["P"].to_numpy(dtype=float, copy=True)
    snp_ok = df["SNP"].notna() & (df["SNP"].astype(str) != "")

    core_valid = (
        snp_ok
        & ea_ok
        & nea_ok
        & (df["EA"] != df["NEA"])
        & np.isfinite(beta_arr)
        & np.isfinite(se_arr)
        & (se_arr > 0.0)
        & np.isfinite(z_arr)
        & np.isfinite(p_arr)
        & (p_arr > 0.0)
        & (p_arr <= 1.0)
    )
    n_invalid_core = int((~core_valid).sum())
    df_cur = df.loc[core_valid].copy()

    # 2) MAF filter when EAF is available
    n_maf_filtered = 0
    if args.min_maf > 0 and df_cur["EAF"].notna().any():
        eaf_arr = df_cur["EAF"].to_numpy(dtype=float, copy=True)
        maf = np.minimum(eaf_arr, 1.0 - eaf_arr)
        maf_ok = np.isnan(eaf_arr) | ((eaf_arr > 0.0) & (eaf_arr < 1.0) & (maf >= args.min_maf))
        n_maf_filtered = int((~maf_ok).sum())
        df_cur = df_cur.loc[maf_ok].copy()

    # 3) INFO filter when INFO is available
    n_info_filtered = 0
    if args.min_info > 0 and df_cur["INFO"].notna().any():
        info_arr = df_cur["INFO"].to_numpy(dtype=float, copy=True)
        info_ok = np.isnan(info_arr) | (info_arr >= args.min_info)
        n_info_filtered = int((~info_ok).sum())
        df_cur = df_cur.loc[info_ok].copy()

    if df_cur.empty:
        fail("0 valid variants remained after GWAS harmonization and QC filters")

    # 4) Deduplicate keeping lowest P
    before_dedup = int(len(df_cur))
    df_cur = df_cur.sort_values("P", ascending=True).drop_duplicates(subset=["SNP"], keep="first")
    n_dups = before_dedup - int(len(df_cur))

    # Sort canonically by CHR, POS
    if meta["has_coords"]:
        df_cur["_chr_key"] = df_cur["CHR"].map(chrom_sort_key)
        df_cur = df_cur.sort_values(["_chr_key", "POS", "P"], ascending=[True, True, True]).drop(columns=["_chr_key"])
    df_cur = df_cur.reset_index(drop=True)

    # 5) Genomic inflation metrics
    z_clean = df_cur["Z"].to_numpy(dtype=float, copy=True)
    chi2 = z_clean * z_clean
    lambda_gc = float(np.median(chi2) / 0.454936423119572)

    lambda_1000: float | None = None
    if args.n_cases is not None and args.n_controls is not None:
        lambda_1000 = float(1.0 + 500.0 * (lambda_gc - 1.0) * (1.0 / args.n_cases + 1.0 / args.n_controls))
    elif df_cur["N"].notna().any():
        n_med = float(np.nanmedian(df_cur["N"].to_numpy(dtype=float, copy=True)))
        if n_med > 0:
            lambda_1000 = float(1.0 + (lambda_gc - 1.0) * (1000.0 / n_med))

    # 6) Extract independent genome-wide significant lead loci
    gw_sig = int((df_cur["P"] <= args.p_threshold).sum())
    leads = clump_by_distance(df_cur, args.p_threshold, args.lead_window_kb)

    out_cols = ["SNP", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "Z", "P", "N", "N_EFF"]
    if meta["has_info"]:
        out_cols.append("INFO")
    df_out = df_cur[out_cols].copy()
    leads_out = leads[out_cols].copy() if not leads.empty else pd.DataFrame(columns=out_cols)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    params_meta = {
        "n_cases": args.n_cases,
        "n_controls": args.n_controls,
        "min_maf": args.min_maf,
        "min_info": args.min_info,
        "lead_window_kb": args.lead_window_kb,
        "p_threshold": args.p_threshold,
    }

    harm_path = outdir / "harmonized_sumstats.tsv"
    leads_path = outdir / "lead_loci.tsv"
    fig_path = outdir / "gwas_manhattan_qq.png"

    # 7) 2-panel Manhattan + QQ figure
    plt = setup_matplotlib()
    fig, (ax_m, ax_q) = plt.subplots(1, 2, figsize=(12.2, 4.6), gridspec_kw={"width_ratios": [2.2, 1.0]})

    has_coords = bool((df_out["CHR"] != "").any() and df_out["POS"].notna().any())
    if has_coords:
        plot_df = df_out[(df_out["CHR"] != "") & df_out["POS"].notna()].copy()
        sig_part = plot_df[plot_df["P"] <= 1e-3]
        bg_part = plot_df[plot_df["P"] > 1e-3]
        if len(bg_part) > 40000:
            bg_part = bg_part.sample(n=40000, random_state=42)
        plot_df = pd.concat([sig_part, bg_part], ignore_index=True)

        chroms = sorted(plot_df["CHR"].astype(str).unique(), key=chrom_sort_key)
        offset = 0.0
        ticks: list[float] = []
        tick_lbls: list[str] = []
        colors = ("#1E3A8A", "#38BDF8")
        chr_offsets: dict[str, float] = {}
        for idx_c, ch in enumerate(chroms):
            sub = plot_df[plot_df["CHR"].astype(str) == ch]
            pos = sub["POS"].to_numpy(dtype=float, copy=True)
            p_nlp = -np.log10(np.clip(sub["P"].to_numpy(dtype=float, copy=True), 1e-300, 1.0))
            x_cum = pos + offset
            ax_m.scatter(x_cum, p_nlp, s=9, color=colors[idx_c % 2], alpha=0.78, linewidths=0)
            chr_offsets[ch] = offset
            ticks.append(offset + 0.5 * (float(pos.min()) + float(pos.max())))
            tick_lbls.append(ch)
            offset += float(pos.max()) + 2e6
        ax_m.set_xticks(ticks[:24])
        ax_m.set_xticklabels(tick_lbls[:24], fontsize=7.5)
        ax_m.set_xlabel("Chromosome")

        if not leads_out.empty:
            lead_x = [
                float(r["POS"]) + chr_offsets.get(str(r["CHR"]), 0.0)
                for _, r in leads_out.head(8).iterrows()
                if pd.notna(r["POS"]) and str(r["CHR"]) in chr_offsets
            ]
            lead_y = [
                -math.log10(max(float(r["P"]), 1e-300))
                for _, r in leads_out.head(8).iterrows()
                if pd.notna(r["POS"]) and str(r["CHR"]) in chr_offsets
            ]
            lead_names = [
                str(r["SNP"])
                for _, r in leads_out.head(8).iterrows()
                if pd.notna(r["POS"]) and str(r["CHR"]) in chr_offsets
            ]
            if lead_x:
                ax_m.scatter(lead_x, lead_y, s=34, color="#E11D48", zorder=5, edgecolors="white", linewidths=0.6)
                for lx, ly, nm in zip(lead_x[:5], lead_y[:5], lead_names[:5]):
                    ax_m.annotate(nm, (lx, ly), xytext=(0, 5), textcoords="offset points", ha="center", fontsize=7.2, color="#0F172A")
    else:
        p_nlp = -np.log10(np.clip(df_out["P"].to_numpy(dtype=float, copy=True), 1e-300, 1.0))
        ax_m.scatter(np.arange(len(p_nlp)), p_nlp, s=10, color="#2563EB", alpha=0.72)
        ax_m.set_xlabel("Variant Index")

    ax_m.axhline(-math.log10(args.p_threshold), color="#E11D48", linestyle="--", linewidth=1.0, label=f"P = {args.p_threshold:.0e}")
    ax_m.set_ylabel("−log10(P)")
    l1000_txt = f", λ_1000 = {lambda_1000:.3f}" if lambda_1000 is not None else ""
    ax_m.set_title(f"Manhattan Plot ({len(leads_out)} lead loci, λ_GC = {lambda_gc:.3f}{l1000_txt})", fontsize=10)
    ax_m.legend(frameon=False, fontsize=8, loc="upper right")

    # QQ panel
    obs_p = np.sort(np.clip(df_out["P"].to_numpy(dtype=float, copy=True), 1e-300, 1.0))
    n_p = len(obs_p)
    if n_p > 50000:
        keep_idx = np.unique(np.concatenate([np.arange(5000), np.linspace(5000, n_p - 1, 10000, dtype=int)]))
        exp_nlp = -np.log10((keep_idx + 0.5) / n_p)
        obs_nlp = -np.log10(obs_p[keep_idx])
    else:
        exp_nlp = -np.log10((np.arange(n_p) + 0.5) / n_p)
        obs_nlp = -np.log10(obs_p)

    ax_q.plot([0, float(exp_nlp.max())], [0, float(exp_nlp.max())], color="#94A3B8", linestyle="--", linewidth=1.0)
    ax_q.scatter(exp_nlp, obs_nlp, s=12, color="#0F172A", alpha=0.75)
    ax_q.set_xlabel("Expected −log10(P)")
    ax_q.set_ylabel("Observed −log10(P)")
    ax_q.set_title(f"QQ Plot (λ_GC = {lambda_gc:.3f}{l1000_txt})", fontsize=10)

    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    artifacts = [
        write_tsv(df_out, harm_path, "gwas_sumstats_harmonize", "Harmonized GWAS summary statistics", inputs=(args.sumstats,), parameters=params_meta),
        write_tsv(leads_out, leads_path, "gwas_sumstats_harmonize", f"Independent lead loci (n={len(leads_out)})", inputs=(args.sumstats,), parameters=params_meta),
        write_artifact(fig_path, "gwas_sumstats_harmonize", f"GWAS Manhattan & QQ plot (λ_GC={lambda_gc:.3f})", inputs=(args.sumstats,), parameters=params_meta),
    ]

    warnings: list[str] = []
    if (args.n_cases is None) != (args.n_controls is None):
        warnings.append("Only one of --n_cases or --n_controls was provided; case-control N_EFF and lambda_1000 require both.")
    if not meta["has_eaf"]:
        warnings.append("EAF column absent; MAF filtering was skipped.")
    if not meta["has_info"]:
        warnings.append("Imputation INFO column absent; INFO filtering was skipped.")
    if lambda_gc > 1.15 and (lambda_1000 is None or lambda_1000 > 1.05):
        warnings.append(f"Elevated genomic inflation (lambda_GC={lambda_gc:.3f}); verify polygenicity vs confounding via gwas_ldsc_interpret.")
    if len(leads_out) == 0:
        warnings.append(f"No variants reached P <= {args.p_threshold:.1e}.")

    n_dropped = n_input - int(len(df_out))
    emit(
        {
            "tool": "gwas_sumstats_harmonize",
            "n_variants_input": n_input,
            "n_variants_harmonized": int(len(df_out)),
            "n_dropped_qc": int(max(0, n_dropped)),
            "effect_scale_input": meta["effect_scale_input"],
            "lambda_gc": round(lambda_gc, 4),
            "lambda_1000": round(lambda_1000, 4) if lambda_1000 is not None else None,
            "n_cases": args.n_cases,
            "n_controls": args.n_controls,
            "n_eff": round(n_eff_cc, 2) if n_eff_cc is not None else None,
            "n_genome_wide_sig": gw_sig,
            "n_lead_loci": int(len(leads_out)),
            "qc_breakdown": {
                "invalid_alleles_or_stats": n_invalid_core,
                "maf_filtered": n_maf_filtered,
                "info_filtered": n_info_filtered,
                "duplicates_removed": n_dups,
            },
            "top_leads": leads_out.head(15).to_dict(orient="records"),
            "outputs": {
                "harmonized_sumstats": str(harm_path),
                "lead_loci": str(leads_path),
            },
            "figures": [str(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
