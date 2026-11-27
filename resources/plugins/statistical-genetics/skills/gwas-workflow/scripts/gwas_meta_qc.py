#!/usr/bin/env python3
"""gwas_meta_qc: Align effect/non-effect alleles across >=2 GWAS cohorts, flag
cross-cohort EAF discordance, run Fixed-Effect IVW and DerSimonian-Laird
Random-Effects meta-analysis, compute Cochran's Q, I^2, and direction strings,
and export heterogeneity and EAF QC diagnostics.
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
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    chrom_sort_key,
    emit,
    fail,
    is_palindromic,
    is_valid_iupac_allele,
    parse_str_list,
    parse_sumstats,
    setup_matplotlib,
    strand_flip,
    write_artifact,
    write_tsv,
)


def _align_to_ref(
    ea: str,
    nea: str,
    beta: float,
    eaf: float,
    ref_ea: str,
    ref_nea: str,
    ref_eaf: float,
) -> tuple[str, float, float]:
    """Align a cohort record to (ref_ea, ref_nea). Returns (status, aligned_beta, aligned_eaf)."""
    pal = is_palindromic(ref_ea, ref_nea)
    if pal:
        if {ea, nea} != {ref_ea, ref_nea}:
            return ("incompatible", np.nan, np.nan)
        # Palindromic SNPs require unambiguous EAF away from 0.5 in both cohorts
        if not (np.isfinite(eaf) and np.isfinite(ref_eaf)):
            return ("palindromic_ambiguous", np.nan, np.nan)
        if abs(eaf - 0.5) < 0.08 or abs(ref_eaf - 0.5) < 0.08:
            return ("palindromic_ambiguous", np.nan, np.nan)
        same_eaf_diff = abs(eaf - ref_eaf)
        flip_eaf_diff = abs((1.0 - eaf) - ref_eaf)
        if same_eaf_diff <= flip_eaf_diff:
            return ("aligned", beta, eaf)
        return ("flipped", -beta, 1.0 - eaf)

    if ea == ref_ea and nea == ref_nea:
        return ("aligned", beta, eaf)
    if ea == ref_nea and nea == ref_ea:
        return ("flipped", -beta, 1.0 - eaf if np.isfinite(eaf) else np.nan)
    f_ea, f_nea = strand_flip(ea), strand_flip(nea)
    if f_ea == ref_ea and f_nea == ref_nea:
        return ("strand_aligned", beta, eaf)
    if f_ea == ref_nea and f_nea == ref_ea:
        return ("flipped", -beta, 1.0 - eaf if np.isfinite(eaf) else np.nan)
    return ("incompatible", np.nan, np.nan)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohorts", nargs="+", action="extend", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--cohort_names", nargs="+", action="extend", default=None)
    ap.add_argument("--eaf_diff_max", type=float, default=0.20)
    args = ap.parse_args()

    cohort_paths = parse_str_list(args.cohorts)
    if len(cohort_paths) < 2:
        fail("gwas_meta_qc requires at least 2 cohort files via --cohorts")

    names_arg = parse_str_list(args.cohort_names)
    if names_arg and len(names_arg) != len(cohort_paths):
        fail(f"Length of --cohort_names ({len(names_arg)}) must match --cohorts ({len(cohort_paths)})")
    cohort_names = names_arg if names_arg else [Path(p).stem.split(".")[0] for p in cohort_paths]

    cohort_dfs: list[pd.DataFrame] = []
    per_cohort_counts: dict[str, int] = {}
    for cname, cpath in zip(cohort_names, cohort_paths):
        cdf, _ = parse_sumstats(cpath, require_alleles=True, require_beta_se=True, require_n=False)
        valid = (
            cdf["SNP"].notna()
            & (cdf["SNP"].astype(str) != "")
            & cdf["EA"].map(is_valid_iupac_allele)
            & cdf["NEA"].map(is_valid_iupac_allele)
            & (cdf["EA"] != cdf["NEA"])
            & np.isfinite(cdf["BETA"].to_numpy(dtype=float, copy=True))
            & np.isfinite(cdf["SE"].to_numpy(dtype=float, copy=True))
            & (cdf["SE"].to_numpy(dtype=float, copy=True) > 0.0)
        )
        cdf = cdf.loc[valid].sort_values("P", ascending=True).drop_duplicates(subset=["SNP"], keep="first").reset_index(drop=True)
        if cdf.empty:
            fail(f"Cohort '{cname}' ({cpath}) has 0 valid SNP records")
        per_cohort_counts[cname] = int(len(cdf))
        cohort_dfs.append(cdf)

    # Index SNPs across cohorts using first appearance as reference allele frame
    snp_records: dict[str, dict[str, Any]] = {}
    n_flipped = 0
    n_incompatible = 0
    n_pal_ambig = 0

    for c_idx, cdf in enumerate(cohort_dfs):
        snps = cdf["SNP"].to_numpy()
        chrs = cdf["CHR"].to_numpy()
        poss = cdf["POS"].to_numpy(dtype=float, copy=True)
        eas = cdf["EA"].to_numpy()
        neas = cdf["NEA"].to_numpy()
        betas = cdf["BETA"].to_numpy(dtype=float, copy=True)
        ses = cdf["SE"].to_numpy(dtype=float, copy=True)
        eafs = cdf["EAF"].to_numpy(dtype=float, copy=True)
        ns = cdf["N"].to_numpy(dtype=float, copy=True)

        for i in range(len(cdf)):
            snp = str(snps[i])
            ea = str(eas[i])
            nea = str(neas[i])
            b_val = float(betas[i])
            se_val = float(ses[i])
            eaf_val = float(eafs[i])
            n_val = float(ns[i])

            if snp not in snp_records:
                # If palindromic with ambiguous EAF in reference itself, track ref cleanly
                snp_records[snp] = {
                    "SNP": snp,
                    "CHR": str(chrs[i]),
                    "POS": poss[i],
                    "EA": ea,
                    "NEA": nea,
                    "ref_eaf": eaf_val,
                    "entries": {c_idx: (b_val, se_val, eaf_val, n_val)},
                }
            else:
                rec = snp_records[snp]
                status, abeta, aeaf = _align_to_ref(ea, nea, b_val, eaf_val, rec["EA"], rec["NEA"], rec["ref_eaf"])
                if status == "incompatible":
                    n_incompatible += 1
                    continue
                if status == "palindromic_ambiguous":
                    n_pal_ambig += 1
                    continue
                if status == "flipped":
                    n_flipped += 1
                rec["entries"][c_idx] = (abeta, se_val, aeaf, n_val)
                if not rec["CHR"] and str(chrs[i]):
                    rec["CHR"] = str(chrs[i])
                if pd.isna(rec["POS"]) and np.isfinite(poss[i]):
                    rec["POS"] = poss[i]

    k_tot = len(cohort_dfs)
    meta_rows: list[dict[str, Any]] = []
    eaf_pairs_x: list[float] = []
    eaf_pairs_y: list[float] = []
    eaf_pairs_out: list[bool] = []

    for snp, rec in snp_records.items():
        entries = rec["entries"]
        k = len(entries)
        betas_k = np.array([entries[ci][0] for ci in sorted(entries)], dtype=float)
        ses_k = np.array([entries[ci][1] for ci in sorted(entries)], dtype=float)
        eafs_k = np.array([entries[ci][2] for ci in sorted(entries)], dtype=float)
        ns_k = np.array([entries[ci][3] for ci in sorted(entries)], dtype=float)

        w_fe = 1.0 / (ses_k * ses_k)
        sum_w = float(np.sum(w_fe))
        beta_fe = float(np.sum(w_fe * betas_k) / sum_w)
        se_fe = float(1.0 / math.sqrt(sum_w))
        z_fe = beta_fe / se_fe
        p_fe = float(max(2.0 * stats.norm.sf(abs(z_fe)), 1e-300))

        if k >= 2:
            q_stat = float(np.sum(w_fe * (betas_k - beta_fe) ** 2))
            q_p = float(stats.chi2.sf(q_stat, df=k - 1))
            i2 = float(max(0.0, (q_stat - (k - 1)) / q_stat) * 100.0) if q_stat > 0.0 else 0.0
            c_denom = sum_w - float(np.sum(w_fe * w_fe)) / sum_w
            tau2 = float(max(0.0, (q_stat - (k - 1)) / c_denom)) if c_denom > 0.0 else 0.0
            w_re = 1.0 / (ses_k * ses_k + tau2)
            sum_wre = float(np.sum(w_re))
            beta_re = float(np.sum(w_re * betas_k) / sum_wre)
            se_re = float(1.0 / math.sqrt(sum_wre))
            z_re = beta_re / se_re
            p_re = float(max(2.0 * stats.norm.sf(abs(z_re)), 1e-300))
        else:
            q_stat, q_p, i2, tau2 = 0.0, 1.0, 0.0, 0.0
            beta_re, se_re, p_re = beta_fe, se_fe, p_fe

        dir_chars = []
        for ci in range(k_tot):
            if ci not in entries:
                dir_chars.append("?")
            elif entries[ci][0] > 0:
                dir_chars.append("+")
            elif entries[ci][0] < 0:
                dir_chars.append("-")
            else:
                dir_chars.append("0")
        direction = "".join(dir_chars)

        fin_eaf = eafs_k[np.isfinite(eafs_k)]
        if fin_eaf.size >= 2:
            eaf_min = float(np.min(fin_eaf))
            eaf_max = float(np.max(fin_eaf))
            eaf_diff = eaf_max - eaf_min
            eaf_mean = float(np.mean(fin_eaf))
            eaf_outlier = bool(eaf_diff > args.eaf_diff_max)
            if 0 in entries and 1 in entries and np.isfinite(entries[0][2]) and np.isfinite(entries[1][2]):
                eaf_pairs_x.append(float(entries[0][2]))
                eaf_pairs_y.append(float(entries[1][2]))
                eaf_pairs_out.append(eaf_outlier)
        elif fin_eaf.size == 1:
            eaf_min = eaf_max = eaf_mean = float(fin_eaf[0])
            eaf_diff = 0.0
            eaf_outlier = False
        else:
            eaf_min = eaf_max = eaf_mean = eaf_diff = np.nan
            eaf_outlier = False

        n_tot = float(np.nansum(ns_k)) if np.isfinite(ns_k).any() else np.nan
        meta_rows.append(
            {
                "SNP": snp,
                "CHR": rec["CHR"],
                "POS": rec["POS"],
                "EA": rec["EA"],
                "NEA": rec["NEA"],
                "EAF_mean": eaf_mean,
                "EAF_diff": eaf_diff,
                "EAF_outlier": eaf_outlier,
                "n_cohorts": k,
                "N_total": n_tot,
                "Direction": direction,
                "BETA_FE": beta_fe,
                "SE_FE": se_fe,
                "P_FE": p_fe,
                "BETA_RE": beta_re,
                "SE_RE": se_re,
                "P_RE": p_re,
                "Cochran_Q": q_stat,
                "Q_P": q_p,
                "I2": i2,
                "tau2_DL": tau2,
            }
        )

    meta_df = pd.DataFrame(meta_rows)
    if meta_df.empty:
        fail("0 SNPs remained after cross-cohort allele alignment")

    if (meta_df["CHR"] != "").any() and meta_df["POS"].notna().any():
        meta_df["_ck"] = meta_df["CHR"].map(chrom_sort_key)
        meta_df = meta_df.sort_values(["_ck", "POS", "P_FE"], ascending=[True, True, True]).drop(columns=["_ck"]).reset_index(drop=True)

    meta_multi = meta_df[meta_df["n_cohorts"] >= 2]
    n_snps_meta = int(len(meta_multi))
    n_eaf_outliers = int(meta_df["EAF_outlier"].sum())
    n_high_het = int((meta_multi["I2"] >= 50.0).sum()) if n_snps_meta > 0 else 0
    n_gws_fe = int((meta_df["P_FE"] <= 5e-8).sum())
    n_gws_re = int((meta_df["P_RE"] <= 5e-8).sum())
    median_i2 = float(np.median(meta_multi["I2"].to_numpy(dtype=float, copy=True))) if n_snps_meta > 0 else 0.0

    top_hits = meta_df.sort_values("P_FE", ascending=True)
    gws_hits = top_hits[top_hits["P_FE"] <= 5e-8]
    if len(gws_hits) >= 10:
        top_hits_df = gws_hits.reset_index(drop=True)
    else:
        top_hits_df = top_hits.head(25).reset_index(drop=True)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    meta_path = outdir / "gwas_meta_sumstats.tsv"
    top_path = outdir / "gwas_meta_top_hits.tsv"
    fig_path = outdir / "gwas_meta_heterogeneity.png"

    # 2-panel diagnostic figure: heterogeneity I2 vs -log10(P_FE), and Cohort 1 vs Cohort 2 EAF QC
    plt = setup_matplotlib()
    fig, (ax_het, ax_eaf) = plt.subplots(1, 2, figsize=(11.6, 4.6))

    if n_snps_meta > 0:
        nlp_fe = -np.log10(np.clip(meta_multi["P_FE"].to_numpy(dtype=float, copy=True), 1e-300, 1.0))
        i2_vals = meta_multi["I2"].to_numpy(dtype=float, copy=True)
        out_mask = meta_multi["EAF_outlier"].to_numpy(dtype=bool, copy=True)
        ax_het.scatter(nlp_fe[~out_mask], i2_vals[~out_mask], s=14, color="#2563EB", alpha=0.72, label="Consistent EAF")
        if out_mask.any():
            ax_het.scatter(nlp_fe[out_mask], i2_vals[out_mask], s=24, color="#E11D48", alpha=0.88, label="EAF outlier")
        ax_het.axhline(50.0, color="#F59E0B", linestyle="--", linewidth=1.0, label="I² = 50%")
        ax_het.axvline(-math.log10(5e-8), color="#64748B", linestyle=":", linewidth=1.0, label="P = 5e-8")
        ax_het.set_xlabel("−log10(P_FE)")
        ax_het.set_ylabel("Heterogeneity I² (%)")
        ax_het.set_ylim(-3.0, 105.0)
        ax_het.set_title(f"Meta-Analysis Heterogeneity (median I² = {median_i2:.1f}%)", fontsize=10)
        ax_het.legend(frameon=False, loc="upper right")
    else:
        ax_het.text(0.5, 0.5, "No multi-cohort overlapping SNPs", ha="center", va="center", transform=ax_het.transAxes)
        ax_het.set_title("Meta-Analysis Heterogeneity", fontsize=10)

    if eaf_pairs_x:
        ex = np.asarray(eaf_pairs_x, dtype=float)
        ey = np.asarray(eaf_pairs_y, dtype=float)
        eo = np.asarray(eaf_pairs_out, dtype=bool)
        ax_eaf.plot([0, 1], [0, 1], color="#64748B", linestyle="--", linewidth=1.0)
        ax_eaf.fill_between([0, 1], [max(0, 0 - args.eaf_diff_max), max(0, 1 - args.eaf_diff_max)], [min(1, 0 + args.eaf_diff_max), min(1, 1 + args.eaf_diff_max)], color="#E2E8F0", alpha=0.45)
        ax_eaf.scatter(ex[~eo], ey[~eo], s=14, color="#0F766E", alpha=0.72, label="Concordant")
        if eo.any():
            ax_eaf.scatter(ex[eo], ey[eo], s=26, color="#E11D48", alpha=0.90, label=f"|ΔEAF| > {args.eaf_diff_max:.2f}")
        ax_eaf.set_xlabel(f"{cohort_names[0]} EAF")
        ax_eaf.set_ylabel(f"{cohort_names[1]} aligned EAF")
        ax_eaf.set_xlim(-0.02, 1.02)
        ax_eaf.set_ylim(-0.02, 1.02)
        ax_eaf.set_title(f"Cross-Cohort EAF Concordance ({n_eaf_outliers} outliers)", fontsize=10)
        ax_eaf.legend(frameon=False, loc="upper left")
    else:
        ax_eaf.text(0.5, 0.5, "EAF column unavailable across cohorts", ha="center", va="center", transform=ax_eaf.transAxes)
        ax_eaf.set_title("Cross-Cohort EAF Concordance", fontsize=10)

    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    params_meta = {"cohort_names": cohort_names, "eaf_diff_max": args.eaf_diff_max}
    artifacts = [
        write_tsv(meta_df, meta_path, "gwas_meta_qc", f"GWAS meta-analysis summary statistics (K={k_tot})", inputs=cohort_paths, parameters=params_meta),
        write_tsv(top_hits_df, top_path, "gwas_meta_qc", f"Top GWAS meta-analysis hits (n={len(top_hits_df)})", inputs=cohort_paths, parameters=params_meta),
        write_artifact(fig_path, "gwas_meta_qc", "Meta-analysis heterogeneity and EAF concordance plot", inputs=cohort_paths, parameters=params_meta),
    ]

    warnings: list[str] = []
    if n_eaf_outliers > 0:
        warnings.append(f"Flagged {n_eaf_outliers} SNPs with cross-cohort EAF discrepancy > {args.eaf_diff_max:.2f}.")
    if n_incompatible > 0:
        warnings.append(f"Dropped {n_incompatible} strand/allele-incompatible cohort records.")
    if n_pal_ambig > 0:
        warnings.append(f"Dropped {n_pal_ambig} ambiguous palindromic records (missing EAF or |EAF - 0.5| < 0.08).")
    if n_high_het > 0:
        warnings.append(f"{n_high_het} meta-analyzed SNPs show substantial heterogeneity (I^2 >= 50%); evaluate RE_DL alongside FE_IVW.")

    emit(
        {
            "tool": "gwas_meta_qc",
            "n_cohorts": k_tot,
            "cohort_names": cohort_names,
            "n_snps_total": int(len(meta_df)),
            "n_snps_meta": n_snps_meta,
            "n_alleles_flipped": n_flipped,
            "n_incompatible_dropped": n_incompatible,
            "n_palindromic_ambiguous_dropped": n_pal_ambig,
            "n_eaf_outliers": n_eaf_outliers,
            "n_high_heterogeneity": n_high_het,
            "n_genome_wide_sig_fe": n_gws_fe,
            "n_genome_wide_sig_re": n_gws_re,
            "median_i2": round(median_i2, 3),
            "per_cohort_counts": per_cohort_counts,
            "top_hits": top_hits_df.head(15).to_dict(orient="records"),
            "outputs": {
                "meta_sumstats": str(meta_path),
                "meta_top_hits": str(top_path),
            },
            "figures": [str(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
