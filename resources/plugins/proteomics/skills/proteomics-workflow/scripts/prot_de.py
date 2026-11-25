#!/usr/bin/env python3
"""prot_de: empirical-Bayes moderated t-test / Welch / Mann-Whitney differential abundance."""

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
    bh_fdr,
    clean_sample_name,
    detect_format_and_columns,
    emit,
    fail,
    infer_scale,
    load_metadata,
    moderated_t_test,
    read_table,
    write_artifact,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--id_col", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default="condition")
    parser.add_argument("--case", required=True)
    parser.add_argument("--control", required=True)
    parser.add_argument("--method", choices=["moderated_t", "welch", "mannwhitney"], default="moderated_t")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--log2fc", type=float, default=1.0)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.input)
    det = detect_format_and_columns(df, args.id_col)
    id_col = det["id_col"]
    sample_cols = det["sample_cols"]

    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    mat.columns = [clean_sample_name(c) for c in sample_cols]
    scale = infer_scale(mat)
    warnings: list[str] = []
    if scale["scale"] == "raw_intensity":
        warnings.append("input appears to be linear raw intensities; run prot_qc_impute to log2-transform first")

    meta, scol = load_metadata(args.metadata, args.sample_col)
    if args.group_col not in meta.columns:
        fail(f"group_col '{args.group_col}' not found in metadata")
    meta_map = dict(zip(meta[scol], meta[args.group_col].astype(str)))

    case_cols = [c for c in mat.columns if meta_map.get(c) == args.case]
    ctrl_cols = [c for c in mat.columns if meta_map.get(c) == args.control]
    if len(case_cols) < 2 or len(ctrl_cols) < 2:
        fail(
            f"need >= 2 samples per group; found {len(case_cols)} for '{args.case}' and {len(ctrl_cols)} for '{args.control}'"
        )
    if len(case_cols) < 3 or len(ctrl_cols) < 3:
        warnings.append(
            f"small replicate count ({args.case}: n={len(case_cols)}, {args.control}: n={len(ctrl_cols)}); interpret FDR cautiously"
        )

    x_case = mat[case_cols].to_numpy(dtype=float)
    x_ctrl = mat[ctrl_cols].to_numpy(dtype=float)
    ids = df[id_col].astype(str).str.strip().to_numpy()

    eb_meta = None
    if args.method == "moderated_t":
        res = moderated_t_test(x_case, x_ctrl)
        log2fc = res["log2fc"]
        mean_case = res["mean_case"]
        mean_ctrl = res["mean_control"]
        ave_expr = res["ave_expr"]
        stat_vals = res["stat"]
        pvals = res["pvalue"]
        padj = res["padj"]
        eb_meta = {"d0": res["d0"], "s0": res["s0"], "shrunk": res["shrunk"]}
    else:
        n_feat = len(ids)
        mean_case = np.nanmean(x_case, axis=1)
        mean_ctrl = np.nanmean(x_ctrl, axis=1)
        log2fc = mean_case - mean_ctrl
        ave_expr = 0.5 * (mean_case + mean_ctrl)
        stat_vals = np.full(n_feat, np.nan)
        pvals = np.full(n_feat, np.nan)
        for i in range(n_feat):
            a = x_case[i][np.isfinite(x_case[i])]
            b = x_ctrl[i][np.isfinite(x_ctrl[i])]
            if len(a) < 2 or len(b) < 2:
                continue
            if args.method == "welch":
                t_val, p_val = stats.ttest_ind(a, b, equal_var=False)
            else:
                t_val, p_val = stats.mannwhitneyu(a, b, alternative="two-sided")
            stat_vals[i] = float(t_val)
            pvals[i] = float(p_val)
        padj = bh_fdr(pvals)

    direction = np.where(
        (padj <= args.fdr) & (log2fc >= args.log2fc),
        "up",
        np.where((padj <= args.fdr) & (log2fc <= -args.log2fc), "down", "ns"),
    )

    out_df = pd.DataFrame(
        {
            "protein_id": ids,
            "log2fc": np.round(log2fc, 5),
            "ave_expr": np.round(ave_expr, 5),
            "mean_case": np.round(mean_case, 5),
            "mean_control": np.round(mean_ctrl, 5),
            "stat": np.round(stat_vals, 5),
            "pvalue": pvals,
            "padj": padj,
            "regulation": direction,
        }
    ).sort_values(["padj", "pvalue"], na_position="last")

    de_path = os.path.join(args.output_dir, "de_results.tsv")
    out_df.to_csv(de_path, sep="\t", index=False)

    n_up = int((out_df["regulation"] == "up").sum())
    n_down = int((out_df["regulation"] == "down").sum())

    # Volcano plot
    fig, ax = plt.subplots(figsize=(6.0, 4.8))
    neg_log_p = -np.log10(np.clip(out_df["pvalue"].to_numpy(dtype=float), 1e-300, 1.0))
    fc = out_df["log2fc"].to_numpy(dtype=float)
    reg = out_df["regulation"].to_numpy()
    ax.scatter(fc[reg == "ns"], neg_log_p[reg == "ns"], s=14, color="#94A3B8", alpha=0.55, label="ns")
    ax.scatter(fc[reg == "up"], neg_log_p[reg == "up"], s=20, color="#EF4444", alpha=0.85, label=f"up ({n_up})")
    ax.scatter(fc[reg == "down"], neg_log_p[reg == "down"], s=20, color="#3B82F6", alpha=0.85, label=f"down ({n_down})")
    ax.axvline(args.log2fc, color="#64748B", linestyle="--", linewidth=0.9)
    ax.axvline(-args.log2fc, color="#64748B", linestyle="--", linewidth=0.9)
    for _, row in out_df[out_df["regulation"] != "ns"].head(10).iterrows():
        yval = -np.log10(max(float(row["pvalue"]), 1e-300))
        ax.annotate(str(row["protein_id"]), (float(row["log2fc"]), yval), fontsize=7, alpha=0.85)
    ax.set_xlabel(f"log2FC ({args.case} vs {args.control})")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"Differential protein abundance ({args.method})")
    ax.legend(frameon=False)
    fig.tight_layout()
    volcano_path = os.path.join(args.output_dir, "de-volcano.png")
    fig.savefig(volcano_path, dpi=150)
    plt.close(fig)

    # MA plot
    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    ave = out_df["ave_expr"].to_numpy(dtype=float)
    ax.scatter(ave[reg == "ns"], fc[reg == "ns"], s=14, color="#94A3B8", alpha=0.55)
    ax.scatter(ave[reg == "up"], fc[reg == "up"], s=20, color="#EF4444", alpha=0.85)
    ax.scatter(ave[reg == "down"], fc[reg == "down"], s=20, color="#3B82F6", alpha=0.85)
    ax.axhline(0.0, color="#334155", linewidth=1.0)
    ax.set_xlabel("Average log2 intensity")
    ax.set_ylabel(f"log2FC ({args.case} vs {args.control})")
    ax.set_title("Proteomics MA plot")
    fig.tight_layout()
    ma_path = os.path.join(args.output_dir, "de-ma.png")
    fig.savefig(ma_path, dpi=150)
    plt.close(fig)

    top_hits = [
        {
            "protein_id": str(r["protein_id"]),
            "log2fc": round(float(r["log2fc"]), 4),
            "padj": float(r["padj"]) if pd.notna(r["padj"]) else None,
            "regulation": str(r["regulation"]),
        }
        for _, r in out_df[out_df["regulation"] != "ns"].head(15).iterrows()
    ]

    params_meta = {
        "contrast": f"{args.case}_vs_{args.control}",
        "method": args.method,
        "fdr": args.fdr,
        "log2fc": args.log2fc,
    }
    in_files = (args.input, args.metadata)
    artifacts = [
        write_artifact(
            volcano_path,
            "prot_de",
            f"Proteomics volcano plot ({args.case} vs {args.control}, {args.method})",
            inputs=in_files,
            parameters=params_meta,
        ),
        write_artifact(
            ma_path,
            "prot_de",
            f"Proteomics MA plot ({args.case} vs {args.control})",
            inputs=in_files,
            parameters=params_meta,
        ),
        write_artifact(
            de_path,
            "prot_de",
            f"Differential protein abundance ({args.case} vs {args.control})",
            inputs=in_files,
            parameters=params_meta,
            shape=out_df.shape,
        ),
    ]
    de_abs = os.path.abspath(de_path)

    emit(
        {
            "input": os.path.abspath(args.input),
            "de_table": de_abs,
            "contrast": f"{args.case}_vs_{args.control}",
            "method": args.method,
            "n_case": len(case_cols),
            "n_control": len(ctrl_cols),
            "n_tested": int(out_df["pvalue"].notna().sum()),
            "fdr_cutoff": args.fdr,
            "log2fc_cutoff": args.log2fc,
            "n_up": n_up,
            "n_down": n_down,
            "empirical_bayes": eb_meta,
            "top_regulated": top_hits,
            "figures": [os.path.abspath(volcano_path), os.path.abspath(ma_path)],
            "scientificContext": {
                de_abs: {
                    "valueSemantics": "de_stats",
                }
            },
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
