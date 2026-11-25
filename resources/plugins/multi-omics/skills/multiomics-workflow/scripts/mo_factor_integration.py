#!/usr/bin/env python3
"""mo_factor_integration: Block-balanced Multi-Omics Joint Latent Factor
Decomposition (Multiple Factor Analysis / Consensus SVD with 1/sigma_1(X_k)
block-inertia weighting) across 2+ matched-sample omics matrices.
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
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    bh_fdr,
    emit,
    fail,
    parse_labeled_paths,
    read_matrix,
    read_metadata_aligned,
    transform_block,
    write_artifact,
)

PALETTE = ["#2563EB", "#E11D48", "#10B981", "#F59E0B", "#8B5CF6", "#06B6D4"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrices", action="extend", nargs="+", default=[], required=True)
    parser.add_argument("--transforms", action="extend", nargs="+", default=[])
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default=None)
    parser.add_argument("--n_factors", type=int, default=5)
    parser.add_argument("--max_features_per_block", type=int, default=2000)
    parser.add_argument("--top_loadings", type=int, default=25)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    labeled = parse_labeled_paths(args.matrices, prefix="modality")
    if len(labeled) < 2:
        fail("mo_factor_integration requires at least 2 omics matrices")

    tf_list = [t.strip() for item in args.transforms for t in str(item).split(",") if t.strip()]
    if tf_list and len(tf_list) != len(labeled):
        fail(f"--transforms length ({len(tf_list)}) must match --matrices length ({len(labeled)})")

    raw_blocks = [(lbl, path, read_matrix(path)) for lbl, path in labeled]
    shared_set = set.intersection(*(set(df.columns) for _, _, df in raw_blocks))
    shared = [s for s in raw_blocks[0][2].columns if s in shared_set]
    if len(shared) < 6:
        fail(f"mo_factor_integration requires >= 6 shared samples across all blocks; found {len(shared)}")

    k_factors = min(args.n_factors, len(shared) - 1)
    if k_factors < 1:
        fail("unable to extract >= 1 latent factor from shared samples")

    scaled_blocks: list[np.ndarray] = []
    block_meta: list[dict] = []
    feature_names_per_block: list[list[str]] = []

    for idx, (lbl, path, df) in enumerate(raw_blocks):
        tf_req = tf_list[idx] if tf_list else "auto"
        t_df, tf_used = transform_block(df[shared], tf_req)
        variances = t_df.var(axis=1)
        non_const = variances[variances > 1e-12].sort_values(ascending=False)
        if non_const.empty:
            fail(f"modality '{lbl}' has 0 non-constant features across shared samples")
        keep_feats = non_const.head(args.max_features_per_block).index.tolist()

        # X_k: n_samples x p_k, feature-standardized
        xk = t_df.loc[keep_feats].to_numpy(dtype=float).T
        xk = xk - np.mean(xk, axis=0, keepdims=True)
        sd = np.std(xk, axis=0, ddof=1, keepdims=True)
        xk = xk / np.where(sd > 1e-12, sd, 1.0)

        # MFA Escofier-Pagès block weight: divide by leading singular value sigma_1(X_k)
        s_vals = np.linalg.svd(xk, compute_uv=False)
        sigma1 = float(s_vals[0]) if s_vals.size > 0 and s_vals[0] > 1e-12 else 1.0
        xk_weighted = xk / sigma1

        scaled_blocks.append(xk_weighted)
        feature_names_per_block.append(keep_feats)
        block_meta.append(
            {
                "modality": lbl,
                "path": str(path.resolve()),
                "transform": tf_used,
                "n_features_retained": len(keep_feats),
                "sigma1_weight": round(1.0 / sigma1, 6),
            }
        )

    # Joint SVD on concatenated weighted blocks [X_1, ..., X_K]
    X_concat = np.hstack(scaled_blocks)
    U, S, Vt = np.linalg.svd(X_concat, full_matrices=False)
    U_k = U[:, :k_factors]
    S_k = S[:k_factors]
    scores = U_k * S_k  # n_samples x k_factors
    factor_names = [f"Factor{f + 1}" for f in range(k_factors)]

    # Per-block R^2 variance explained by each orthogonal shared factor u_f:
    # Since columns of U_k are orthonormal, projection of X_b onto u_f is u_f (u_f^T X_b),
    # with Frobenius squared norm ||u_f^T X_b||_2^2!
    r2_mat = np.zeros((len(scaled_blocks), k_factors), dtype=float)
    loading_rows: list[dict] = []

    col_offset = 0
    for b_idx, xk_w in enumerate(scaled_blocks):
        mod_name = block_meta[b_idx]["modality"]
        feats = feature_names_per_block[b_idx]
        p_b = xk_w.shape[1]
        total_ss = float(np.sum(xk_w * xk_w))
        v_block = Vt[:k_factors, col_offset : col_offset + p_b]  # k_factors x p_b
        col_offset += p_b

        for f_idx in range(k_factors):
            uf = U_k[:, f_idx]
            proj = uf @ xk_w  # 1 x p_b
            ss_f = float(np.sum(proj * proj))
            r2_mat[b_idx, f_idx] = (100.0 * ss_f / total_ss) if total_ss > 1e-12 else 0.0

            w_vec = v_block[f_idx]
            order = np.argsort(np.abs(w_vec))[::-1][: min(args.top_loadings, p_b)]
            for rank_i, feat_i in enumerate(order, start=1):
                loading_rows.append(
                    {
                        "factor": factor_names[f_idx],
                        "modality": mod_name,
                        "rank": rank_i,
                        "feature": feats[feat_i],
                        "loading": round(float(w_vec[feat_i]), 5),
                        "abs_loading": round(float(abs(w_vec[feat_i])), 5),
                    }
                )

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    scores_df = pd.DataFrame(np.round(scores, 5), index=shared, columns=factor_names)
    scores_df.index.name = "sample"
    scores_path = outdir / "shared_factor_scores.tsv"
    scores_df.to_csv(scores_path, sep="\t")

    r2_df = pd.DataFrame(
        np.round(r2_mat, 3),
        index=[m["modality"] for m in block_meta],
        columns=factor_names,
    )
    r2_df["Total_Pct"] = np.round(r2_df.sum(axis=1), 3)
    r2_df.index.name = "modality"
    r2_path = outdir / "block_variance_explained.tsv"
    r2_df.to_csv(r2_path, sep="\t")

    load_df = pd.DataFrame(loading_rows)
    load_path = outdir / "top_feature_loadings.tsv"
    load_df.to_csv(load_path, sep="\t", index=False)

    # Factor characterization (shared vs modality-specific)
    factor_summaries = []
    for f_idx, fname in enumerate(factor_names):
        per_mod = {
            block_meta[b]["modality"]: round(float(r2_mat[b, f_idx]), 3)
            for b in range(len(block_meta))
        }
        active_mods = [m for m, pct in per_mod.items() if pct >= 5.0]
        factor_summaries.append(
            {
                "factor": fname,
                "singular_value": round(float(S_k[f_idx]), 4),
                "variance_pct_by_modality": per_mod,
                "kind": "shared_multiomics" if len(active_mods) >= 2 else "modality_dominant",
            }
        )

    meta_aligned = read_metadata_aligned(args.metadata, shared, args.sample_col) if args.metadata else None
    grp_assoc = None
    if args.group_col and meta_aligned is not None:
        if args.group_col not in meta_aligned.columns:
            fail(f"group_col '{args.group_col}' not found in metadata")
        grps = meta_aligned[args.group_col].astype(str)
        uniq = sorted(grps.unique())
        if len(uniq) >= 2:
            g_rows = []
            for f_idx, fname in enumerate(factor_names):
                fvals = scores[:, f_idx]
                arrays = [fvals[(grps == g).to_numpy()] for g in uniq if (grps == g).sum() >= 2]
                if len(arrays) >= 2:
                    stat_v, p_v = stats.kruskal(*arrays)
                    g_rows.append({"factor": fname, "test": "Kruskal-Wallis", "statistic": round(float(stat_v), 4), "pvalue": float(p_v)})
            if g_rows:
                p_arr = bh_fdr(np.array([r["pvalue"] for r in g_rows], dtype=float))
                for r, q in zip(g_rows, p_arr):
                    r["fdr"] = float(q)
                grp_assoc = g_rows

    # Two-panel diagnostic figure: (A) Per-block R^2 heatmap + (B) Factor1 vs Factor2 sample scatter
    fig, (ax_r2, ax_sc) = plt.subplots(1, 2, figsize=(11.2, 4.6))
    im = ax_r2.imshow(r2_mat, cmap="YlGnBu", aspect="auto")
    ax_r2.set_xticks(np.arange(k_factors))
    ax_r2.set_xticklabels(factor_names, rotation=30, ha="right")
    ax_r2.set_yticks(np.arange(len(block_meta)))
    ax_r2.set_yticklabels([m["modality"] for m in block_meta])
    for bi in range(len(block_meta)):
        for fi in range(k_factors):
            ax_r2.text(fi, bi, f"{r2_mat[bi, fi]:.1f}%", ha="center", va="center", fontsize=8, color="#0F172A")
    ax_r2.set_title("Variance Explained (%) per Modality")
    fig.colorbar(im, ax=ax_r2, fraction=0.04, pad=0.03)

    f2_y = scores[:, 1] if k_factors >= 2 else np.zeros(len(shared))
    if args.group_col and meta_aligned is not None and args.group_col in meta_aligned.columns:
        grps = meta_aligned[args.group_col].astype(str)
        for g_idx, g in enumerate(sorted(grps.unique())):
            m = (grps == g).to_numpy()
            ax_sc.scatter(scores[m, 0], f2_y[m], s=42, color=PALETTE[g_idx % len(PALETTE)], alpha=0.88, label=g)
        ax_sc.legend(frameon=False, fontsize=8.5, loc="best")
    else:
        ax_sc.scatter(scores[:, 0], f2_y, s=42, color="#2563EB", alpha=0.88)
    ax_sc.axhline(0.0, color="#CBD5E1", linestyle=":", linewidth=0.9)
    ax_sc.axvline(0.0, color="#CBD5E1", linestyle=":", linewidth=0.9)
    ax_sc.set_xlabel("Factor 1 Score")
    ax_sc.set_ylabel("Factor 2 Score" if k_factors >= 2 else "Zero")
    ax_sc.set_title(f"Joint Multi-Omics Sample Space (n={len(shared)})")
    ax_sc.spines["top"].set_visible(False)
    ax_sc.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = outdir / "multiomics-factors.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    in_files = [p for _, p in labeled] + ([args.metadata] if args.metadata else [])
    params_meta = {"n_factors": k_factors, "max_features_per_block": args.max_features_per_block}
    artifacts = [
        write_artifact(fig_path, "mo_factor_integration", "Joint multi-omics factor R² & sample space", inputs=in_files, parameters=params_meta),
        write_artifact(r2_path, "mo_factor_integration", "Per-modality variance explained (%) table", inputs=in_files, parameters=params_meta, shape=r2_df.shape),
        write_artifact(scores_path, "mo_factor_integration", "Shared multi-omics sample factor scores", inputs=in_files, parameters=params_meta, shape=scores_df.shape),
        write_artifact(load_path, "mo_factor_integration", "Top multi-omics feature loadings per block", inputs=in_files, parameters=params_meta, shape=load_df.shape),
    ]

    emit(
        {
            "tool": "mo_factor_integration",
            "n_shared_samples": len(shared),
            "n_factors": k_factors,
            "modalities": block_meta,
            "factor_variance_explained": factor_summaries,
            "group_associations": grp_assoc,
            "outputs": {
                "factor_scores": str(scores_path),
                "variance_explained": str(r2_path),
                "top_loadings": str(load_path),
            },
            "figures": [str(fig_path)],
            "scientificContext": {
                str(scores_path): {
                    "valueSemantics": "normalized",
                    "sampleIds": shared,
                }
            },
            "artifacts": artifacts,
            "warnings": [],
        }
    )


if __name__ == "__main__":
    main()
