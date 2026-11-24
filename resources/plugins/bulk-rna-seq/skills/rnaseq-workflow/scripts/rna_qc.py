#!/usr/bin/env python3
"""rna_qc: sample-level QC, library size & detected gene summary, log2(CPM+1)
PCA ordination, pairwise Pearson correlation heatmap, and outlier flagging.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import emit, fail, read_expression_matrix, read_metadata, setup_matplotlib, write_artifact


def compute_log2_cpm(df):
    x = df.to_numpy(dtype=float)
    if (x < 0).any():
        # Already log/VST-centered
        return x, "input_as_is"
    col_sums = x.sum(axis=0)
    if np.median(col_sums) <= 50.0:
        return x, "input_log_scale"
    safe_sums = np.where(col_sums > 0, col_sums, 1.0)
    cpm = x / safe_sums * 1e6
    return np.log2(cpm + 1.0), "log2_cpm_plus_1"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--metadata")
    ap.add_argument("--sample_col", default=None)
    ap.add_argument("--group_col", default="condition")
    ap.add_argument("--top_var", type=int, default=500)
    ap.add_argument("--min_corr_outlier", type=float, default=0.80)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    df = read_expression_matrix(args.input)
    samples = df.columns.tolist()
    if len(samples) < 3:
        fail(f"rna_qc requires at least 3 samples (got {len(samples)})")

    meta = read_metadata(args.metadata, samples=samples, sample_col=args.sample_col) if args.metadata else None
    labels = None
    if meta is not None and args.group_col in meta.columns:
        labels = meta[args.group_col].astype(str).to_numpy()

    x_log, transform_used = compute_log2_cpm(df)
    lib_sizes = df.sum(axis=0).to_numpy(dtype=float)
    detected = (df.to_numpy(dtype=float) > 0).sum(axis=0)

    # Sample-to-sample Pearson correlation
    corr = np.corrcoef(x_log.T)
    mean_corr = []
    for i in range(len(samples)):
        others = [corr[i, j] for j in range(len(samples)) if j != i]
        mean_corr.append(float(np.mean(others)))

    # PCA on top variable genes
    gene_var = x_log.var(axis=1)
    k_var = min(args.top_var, x_log.shape[0])
    top_idx = np.argsort(gene_var)[::-1][:k_var]
    x_sub = x_log[top_idx].T  # samples x k_var
    x_c = x_sub - x_sub.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(x_c, full_matrices=False)
    coords = u * s
    eig = (s * s)
    var_exp = (eig / eig.sum()).tolist() if eig.sum() > 0 else []

    # Outlier heuristic: mean correlation below floor OR PC1/PC2 distance > 3 SD
    pc_dist = np.linalg.norm(coords[:, : min(2, coords.shape[1])], axis=1)
    pc_cut = float(pc_dist.mean() + 3.0 * pc_dist.std()) if pc_dist.std() > 0 else np.inf
    qc_rows = []
    flagged = []
    for i, samp in enumerate(samples):
        is_out = bool(mean_corr[i] < args.min_corr_outlier or pc_dist[i] > pc_cut)
        if is_out:
            flagged.append(samp)
        row = {
            "sample": samp,
            "library_size": round(float(lib_sizes[i]), 2),
            "detected_genes": int(detected[i]),
            "mean_pearson_r": round(mean_corr[i], 4),
            "PC1": round(float(coords[i, 0]), 4) if coords.shape[1] >= 1 else 0.0,
            "PC2": round(float(coords[i, 1]), 4) if coords.shape[1] >= 2 else 0.0,
            "outlier_flag": is_out,
        }
        if labels is not None:
            row["group"] = str(labels[i])
        qc_rows.append(row)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = [args.input] + ([args.metadata] if args.metadata else [])
    params = {
        "transform": transform_used,
        "top_var": k_var,
        "min_corr_outlier": args.min_corr_outlier,
    }

    qc_df = pd.DataFrame(qc_rows)
    qc_path = outdir / "sample-qc.tsv"
    qc_df.to_csv(qc_path, sep="\t", index=False)

    corr_df = pd.DataFrame(np.round(corr, 4), index=samples, columns=samples)
    corr_df.index.name = "sample"
    corr_path = outdir / "sample-correlation.tsv"
    corr_df.to_csv(corr_path, sep="\t")

    plt = setup_matplotlib()
    # PCA figure
    fig, ax = plt.subplots(figsize=(6, 5))
    if labels is not None:
        for g in np.unique(labels):
            m = labels == g
            ax.scatter(coords[m, 0], coords[m, 1], label=str(g), s=55, alpha=0.85)
        ax.legend(fontsize=8)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=55, alpha=0.85)
    for i, s_name in enumerate(samples):
        ax.annotate(s_name, (coords[i, 0], coords[i, 1]), fontsize=7)
    ax.set_xlabel(f"PC1 ({var_exp[0] * 100:.1f}%)" if len(var_exp) > 0 else "PC1")
    ax.set_ylabel(f"PC2 ({var_exp[1] * 100:.1f}%)" if len(var_exp) > 1 else "PC2")
    ax.set_title(f"Bulk RNA-seq PCA (top {k_var} variable genes)")
    fig.tight_layout()
    pca_path = outdir / "qc-pca.png"
    fig.savefig(pca_path, dpi=180)
    plt.close(fig)

    # Correlation heatmap
    fig, ax = plt.subplots(figsize=(max(5, 0.45 * len(samples) + 1.5), max(4.5, 0.45 * len(samples) + 1.2)))
    im = ax.imshow(corr, vmin=max(0.0, float(np.min(corr)) - 0.05), vmax=1.0, cmap="YlGnBu")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=7)
    ax.set_title("Sample-to-sample Pearson correlation")
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    hm_path = outdir / "qc-correlation-heatmap.png"
    fig.savefig(hm_path, dpi=180)
    plt.close(fig)

    artifacts = [
        write_artifact(qc_path, "rna_qc", "Bulk RNA-seq sample QC summary", inputs=inputs, parameters=params, shape=qc_df.shape),
        write_artifact(corr_path, "rna_qc", "Sample-to-sample Pearson correlation matrix", inputs=inputs, parameters=params, shape=corr_df.shape),
        write_artifact(pca_path, "rna_qc", "Bulk RNA-seq PCA plot", inputs=inputs, parameters=params),
        write_artifact(hm_path, "rna_qc", "Sample correlation heatmap", inputs=inputs, parameters=params),
    ]

    emit(
        {
            "tool": "rna_qc",
            "n_samples": int(len(samples)),
            "n_genes": int(df.shape[0]),
            "transform": transform_used,
            "pca_variance_explained": [round(float(v), 4) for v in var_exp[:5]],
            "flagged_outliers": flagged,
            "outputs": {
                "sample_qc": str(qc_path),
                "correlation_matrix": str(corr_path),
                "pca_figure": str(pca_path),
                "correlation_figure": str(hm_path),
            },
            "artifacts": artifacts,
        }
    )


if __name__ == "__main__":
    main()
