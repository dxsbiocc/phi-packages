#!/usr/bin/env python3
"""meth_qc: global DNA methylation bimodality, coverage summary, logit M-value PCA,
and sample outlier detection.
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
from _common import beta_to_m, emit, fail, m_to_beta, read_methylation_matrix, read_table, write_artifact

PALETTE = ["#2563EB", "#E11D48", "#10B981", "#F59E0B", "#8B5CF6", "#06B6D4"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--coverage_matrix", default=None)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    mat_df, scale = read_methylation_matrix(args.matrix)
    samples = list(mat_df.columns)
    n_cpgs, n_samples = mat_df.shape
    if n_samples < 2:
        fail("meth_qc requires at least 2 sample columns")

    raw_arr = mat_df.to_numpy(dtype=float)
    beta_arr = raw_arr if scale == "beta" else m_to_beta(raw_arr)
    m_arr = beta_to_m(beta_arr) if scale == "beta" else raw_arr

    cov_df = None
    if args.coverage_matrix:
        cov_df = read_table(args.coverage_matrix).set_index(read_table(args.coverage_matrix).columns[0])

    sample_rows = []
    for j, sname in enumerate(samples):
        col = beta_arr[:, j]
        fin = col[np.isfinite(col)]
        n_fin = max(1, fin.size)
        row = {
            "sample": sname,
            "valid_cpgs": int(fin.size),
            "missing_fraction": round(float(1.0 - fin.size / max(1, n_cpgs)), 4),
            "mean_beta": round(float(np.mean(fin)), 4) if fin.size > 0 else None,
            "median_beta": round(float(np.median(fin)), 4) if fin.size > 0 else None,
            "frac_hypo_lt02": round(float(np.sum(fin < 0.2) / n_fin), 4),
            "frac_intermediate_02_08": round(float(np.sum((fin >= 0.2) & (fin <= 0.8)) / n_fin), 4),
            "frac_hyper_gt08": round(float(np.sum(fin > 0.8) / n_fin), 4),
        }
        if cov_df is not None and sname in cov_df.columns:
            cvals = pd.to_numeric(cov_df[sname], errors="coerce").to_numpy(dtype=float)
            pos_c = cvals[np.isfinite(cvals) & (cvals > 0)]
            row["median_coverage"] = round(float(np.median(pos_c)), 2) if pos_c.size > 0 else 0.0
        sample_rows.append(row)

    # M-value PCA on top 5000 most variable CpGs (median-imputed)
    row_med = np.nanmedian(m_arr, axis=1)
    row_med = np.where(np.isfinite(row_med), row_med, 0.0)
    m_imp = m_arr.copy()
    na_i = np.where(~np.isfinite(m_imp))
    if na_i[0].size > 0:
        m_imp[na_i] = row_med[na_i[0]]

    rv = np.var(m_imp, axis=1)
    top_idx = np.argsort(rv)[::-1][: min(5000, n_cpgs)]
    X = m_imp[top_idx].T  # n_samples x top_cpgs
    X = X - np.mean(X, axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    scores = U * S
    tot_var = float(np.sum(S * S))
    var_pct = [round(100.0 * float(s * s) / tot_var, 2) if tot_var > 1e-12 else 0.0 for s in S[: min(5, len(S))]]

    # Outlier detection on PC1/PC2 robust distance (> 3.5 MADs from median)
    pc12 = scores[:, : min(2, scores.shape[1])]
    med_pc = np.median(pc12, axis=0)
    dists = np.linalg.norm(pc12 - med_pc, axis=1)
    mad_d = float(np.median(np.abs(dists - np.median(dists)))) * 1.4826
    outliers = []
    for j, sname in enumerate(samples):
        is_out = bool(mad_d > 1e-9 and (dists[j] - np.median(dists)) / mad_d > 3.5)
        sample_rows[j]["pc1"] = round(float(scores[j, 0]), 4)
        sample_rows[j]["pc2"] = round(float(scores[j, 1]), 4) if scores.shape[1] >= 2 else 0.0
        sample_rows[j]["pca_outlier"] = is_out
        if is_out:
            outliers.append(sname)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    qc_df = pd.DataFrame(sample_rows)
    qc_path = outdir / "sample_qc.tsv"
    qc_df.to_csv(qc_path, sep="\t", index=False)

    # Group mapping if metadata provided
    grp_map: dict[str, str] = {}
    if args.metadata and args.group_col:
        meta = read_table(args.metadata)
        scol = args.sample_col or str(meta.columns[0])
        if scol in meta.columns and args.group_col in meta.columns:
            grp_map = {
                str(k).strip(): str(v).strip()
                for k, v in zip(meta[scol], meta[args.group_col])
            }

    # Plot: (A) Beta density curves + (B) M-value PCA
    fig, (ax_b, ax_p) = plt.subplots(1, 2, figsize=(11.0, 4.5))
    grid = np.linspace(0.0, 1.0, 60)
    uniq_grps = sorted(set(grp_map.values())) if grp_map else []
    grp_colors = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(uniq_grps)}

    for j, sname in enumerate(samples):
        fin = beta_arr[:, j]
        fin = fin[np.isfinite(fin)]
        if fin.size < 5:
            continue
        hist, edges = np.histogram(fin, bins=30, range=(0.0, 1.0), density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        col = grp_colors.get(grp_map.get(sname, ""), "#2563EB")
        ax_b.plot(centers, hist, color=col, alpha=0.65, linewidth=1.3)
    ax_b.set_xlabel("Methylation β-value")
    ax_b.set_ylabel("Density")
    ax_b.set_title("Per-Sample β-Value Distribution")
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    pc2_vals = scores[:, 1] if scores.shape[1] >= 2 else np.zeros(n_samples)
    if uniq_grps:
        for g in uniq_grps:
            idx_g = [i for i, s in enumerate(samples) if grp_map.get(s) == g]
            if idx_g:
                ax_p.scatter(scores[idx_g, 0], pc2_vals[idx_g], s=46, color=grp_colors[g], label=g, alpha=0.88)
        ax_p.legend(frameon=False, fontsize=8.5, loc="best")
    else:
        ax_p.scatter(scores[:, 0], pc2_vals, s=46, color="#2563EB", alpha=0.88)
    for i, sname in enumerate(samples):
        if sname in outliers:
            ax_p.annotate(sname, (scores[i, 0], pc2_vals[i]), fontsize=7.5, color="#E11D48")
    p1_lbl = f"PC1 ({var_pct[0]:.1f}%)" if var_pct else "PC1"
    p2_lbl = f"PC2 ({var_pct[1]:.1f}%)" if len(var_pct) >= 2 else "PC2"
    ax_p.set_xlabel(p1_lbl)
    ax_p.set_ylabel(p2_lbl)
    ax_p.set_title("Logit M-Value PCA")
    ax_p.spines["top"].set_visible(False)
    ax_p.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = outdir / "meth-qc.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    in_files = (args.matrix, args.coverage_matrix, args.metadata)
    artifacts = [
        write_artifact(fig_path, "meth_qc", "DNA methylation β-density & M-value PCA QC", inputs=in_files),
        write_artifact(qc_path, "meth_qc", "Per-sample DNA methylation QC summary", inputs=in_files, shape=qc_df.shape),
    ]

    warnings: list[str] = []
    if outliers:
        warnings.append(f"PCA outlier sample(s) flagged (>3.5 MAD): {outliers}")

    emit(
        {
            "tool": "meth_qc",
            "input_scale": scale,
            "n_samples": n_samples,
            "n_cpgs": n_cpgs,
            "pca_variance_explained_pct": var_pct,
            "per_sample": sample_rows,
            "outliers": outliers,
            "outputs": {"sample_qc": str(qc_path)},
            "figures": [str(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
