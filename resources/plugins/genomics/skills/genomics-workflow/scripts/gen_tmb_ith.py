#!/usr/bin/env python3
"""gen_tmb_ith: FoCR-compliant TMB and Mroz/Rocco MATH intratumor heterogeneity."""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    NON_SYNONYMOUS_CLASSES,
    emit,
    extract_vaf,
    fail,
    load_clinical,
    read_maf,
    read_seg,
)


def variant_in_cnv(sample: str, chrom: str, pos: float, seg_by_sample: dict[str, pd.DataFrame]) -> bool:
    sub = seg_by_sample.get(sample)
    if sub is None or sub.empty or not np.isfinite(pos):
        return False
    c = str(chrom).replace("chr", "")
    hit = sub[(sub["Chromosome"] == c) & (sub["Start"] <= pos) & (sub["End"] >= pos)]
    return not hit.empty


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--capture_size_mb", type=float, required=True)
    parser.add_argument("--min_vaf", type=float, default=0.05)
    parser.add_argument("--vaf_col", default=None)
    parser.add_argument("--seg", default=None)
    parser.add_argument("--cnv_seg_cutoff", type=float, default=0.2)
    parser.add_argument("--clinical", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    if args.capture_size_mb <= 0:
        fail("capture_size_mb must be > 0")

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_maf(args.maf)
    all_samples = sorted(df["Tumor_Sample_Barcode"].unique())
    vaf, vaf_source = extract_vaf(df, args.vaf_col)
    df["_vaf"] = vaf if vaf is not None else np.nan

    warnings: list[str] = []
    if vaf is None:
        warnings.append("no VAF column detected; TMB computed without min_vaf filter and MATH score is NA")

    seg_by_sample: dict[str, pd.DataFrame] = {}
    cnv_filtered = False
    if args.seg:
        seg_df = read_seg(args.seg)
        alt_seg = seg_df[seg_df["Segment_Mean"].abs() >= args.cnv_seg_cutoff]
        for s, grp in alt_seg.groupby("Sample"):
            seg_by_sample[str(s)] = grp
        cnv_filtered = True

    nonsyn = df[df["Variant_Classification"].isin(NON_SYNONYMOUS_CLASSES)].copy()
    if vaf is not None:
        tmb_df = nonsyn[(nonsyn["_vaf"].isna()) | (nonsyn["_vaf"] >= args.min_vaf)]
    else:
        tmb_df = nonsyn

    rows = []
    n_cnv_excluded_total = 0
    for s in all_samples:
        n_mut = int((tmb_df["Tumor_Sample_Barcode"] == s).sum())
        tmb_val = n_mut / args.capture_size_mb

        # MATH score on somatic VAFs >= min_vaf (excluding CNV-altered segments if --seg provided)
        math_val = np.nan
        n_vaf_used = 0
        if vaf is not None:
            sub_s = df[(df["Tumor_Sample_Barcode"] == s) & (df["_vaf"].notna()) & (df["_vaf"] >= args.min_vaf)]
            if cnv_filtered and not sub_s.empty:
                keep_mask = []
                for _, r in sub_s.iterrows():
                    in_cnv = variant_in_cnv(s, r["Chromosome"], float(r["Start_Position"]), seg_by_sample)
                    if in_cnv:
                        n_cnv_excluded_total += 1
                    keep_mask.append(not in_cnv)
                sub_s = sub_s.loc[keep_mask]
            v_arr = sub_s["_vaf"].to_numpy(dtype=float)
            n_vaf_used = int(len(v_arr))
            if n_vaf_used >= 3:
                med = float(np.median(v_arr))
                if med > 0:
                    mad = float(np.median(np.abs(v_arr - med)))
                    math_val = 100.0 * mad / med

        rows.append(
            {
                "Tumor_Sample_Barcode": s,
                "n_coding_nonsyn": n_mut,
                "tmb_per_mb": round(tmb_val, 4),
                "n_vaf_loci": n_vaf_used,
                "math_score": round(math_val, 4) if np.isfinite(math_val) else np.nan,
            }
        )

    out_df = pd.DataFrame(rows)

    group_tests = None
    if args.clinical and args.group_col:
        clin, scol = load_clinical(args.clinical, args.sample_col)
        if args.group_col not in clin.columns:
            fail(f"group_col '{args.group_col}' not found in clinical table")
        gmap = dict(zip(clin[scol], clin[args.group_col].astype(str)))
        out_df[args.group_col] = out_df["Tumor_Sample_Barcode"].map(gmap)
        valid_g = out_df.dropna(subset=[args.group_col])
        groups = sorted(valid_g[args.group_col].unique())
        if len(groups) >= 2:
            group_tests = {"group_col": args.group_col, "groups": groups}
            for metric in ("tmb_per_mb", "math_score"):
                arrays = [valid_g.loc[valid_g[args.group_col] == g, metric].dropna().to_numpy(dtype=float) for g in groups]
                arrays_nonempty = [a for a in arrays if len(a) >= 2]
                if len(arrays_nonempty) == 2:
                    u_stat, p_val = stats.mannwhitneyu(arrays_nonempty[0], arrays_nonempty[1], alternative="two-sided")
                    group_tests[metric] = {"test": "mannwhitneyu", "stat": round(float(u_stat), 4), "pvalue": float(p_val)}
                elif len(arrays_nonempty) > 2:
                    h_stat, p_val = stats.kruskal(*arrays_nonempty)
                    group_tests[metric] = {"test": "kruskal", "stat": round(float(h_stat), 4), "pvalue": float(p_val)}

    out_path = os.path.join(args.output_dir, "tmb_math_per_sample.tsv")
    out_df.to_csv(out_path, sep="\t", index=False)

    # Plot TMB and MATH distributions
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    tmb_vals = out_df["tmb_per_mb"].dropna().to_numpy()
    axes[0].hist(tmb_vals, bins=min(20, max(5, len(tmb_vals) // 2)), color="#38BDF8", edgecolor="#0F172A", alpha=0.85)
    axes[0].axvline(10.0, color="#EF4444", linestyle="--", linewidth=1.1, label="10 mut/Mb ref")
    axes[0].set_xlabel("TMB (mutations / Mb)")
    axes[0].set_ylabel("Samples")
    axes[0].set_title(f"TMB distribution (capture = {args.capture_size_mb} Mb)")
    axes[0].legend(frameon=False)

    math_vals = out_df["math_score"].dropna().to_numpy()
    if len(math_vals) > 0:
        axes[1].hist(math_vals, bins=min(20, max(5, len(math_vals) // 2)), color="#A78BFA", edgecolor="#0F172A", alpha=0.85)
        axes[1].set_xlabel("MATH score (100 * MAD / median VAF)")
        axes[1].set_ylabel("Samples")
        axes[1].set_title("Intratumor heterogeneity (MATH)")
    else:
        axes[1].text(0.5, 0.5, "MATH unavailable (no VAF)", ha="center", va="center")
        axes[1].set_axis_off()
    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "tmb-math-distribution.png")
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)

    emit(
        {
            "maf": os.path.abspath(args.maf),
            "table": os.path.abspath(out_path),
            "capture_size_mb": args.capture_size_mb,
            "min_vaf": args.min_vaf,
            "vaf_source": vaf_source,
            "cnv_filtered_for_math": cnv_filtered,
            "n_variants_excluded_by_cnv": n_cnv_excluded_total,
            "n_samples": len(all_samples),
            "tmb_median": round(float(np.median(tmb_vals)), 4) if len(tmb_vals) else None,
            "tmb_high_ge10_samples": int(np.sum(tmb_vals >= 10.0)),
            "math_median": round(float(np.median(math_vals)), 4) if len(math_vals) else None,
            "group_comparison": group_tests,
            "figures": [os.path.abspath(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
