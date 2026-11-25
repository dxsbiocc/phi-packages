#!/usr/bin/env python3
"""mo_pair_assoc: Two-block global coupling (Symmetric Procrustes & Mantel
permutation tests) and covariate-adjusted Partial Spearman cross-block network.
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
from scipy.spatial.distance import pdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import bh_fdr, emit, fail, read_matrix, read_metadata_aligned, transform_block, write_artifact

PALETTE = ["#2563EB", "#E11D48", "#10B981", "#F59E0B", "#8B5CF6", "#06B6D4"]


def pca_coords(sample_by_feat: np.ndarray, k: int) -> np.ndarray:
    xc = sample_by_feat - np.mean(sample_by_feat, axis=0, keepdims=True)
    sd = np.std(xc, axis=0, ddof=1, keepdims=True)
    xc = xc / np.where(sd > 1e-12, sd, 1.0)
    u, s, _ = np.linalg.svd(xc, full_matrices=False)
    k_eff = min(k, u.shape[1])
    return u[:, :k_eff] * s[:k_eff]


def symmetric_procrustes(X: np.ndarray, Y: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Symmetric Procrustes M^2 statistic and aligned coordinates (X_std, Y_rot)."""
    Xc = X - np.mean(X, axis=0, keepdims=True)
    Yc = Y - np.mean(Y, axis=0, keepdims=True)
    nx = np.linalg.norm(Xc)
    ny = np.linalg.norm(Yc)
    if nx < 1e-12 or ny < 1e-12:
        return 1.0, Xc, Yc
    Xs = Xc / nx
    Ys = Yc / ny
    # Match column width
    d = max(Xs.shape[1], Ys.shape[1])
    if Xs.shape[1] < d:
        Xs = np.pad(Xs, ((0, 0), (0, d - Xs.shape[1])))
    if Ys.shape[1] < d:
        Ys = np.pad(Ys, ((0, 0), (0, d - Ys.shape[1])))
    u, s, vt = np.linalg.svd(Xs.T @ Ys, full_matrices=False)
    R = vt.T @ u.T
    scale = float(np.sum(s))
    m2 = float(np.clip(1.0 - scale * scale, 0.0, 1.0))
    Y_aligned = scale * (Ys @ R)
    return m2, Xs, Y_aligned


def build_covariate_design(meta: pd.DataFrame, cov_cols: list[str]) -> tuple[np.ndarray, int]:
    cols = [np.ones(len(meta), dtype=float)]
    for c in cov_cols:
        if c not in meta.columns:
            fail(f"covariate '{c}' not found in metadata columns: {list(meta.columns)}")
        s = meta[c]
        num = pd.to_numeric(s, errors="coerce")
        if num.notna().all() and num.nunique() > 2:
            v = num.to_numpy(dtype=float)
            cols.append((v - v.mean()) / max(v.std(ddof=1), 1e-12))
        else:
            cats = s.fillna("MISSING").astype(str)
            lvls = sorted(cats.unique())
            for lvl in lvls[1:]:
                cols.append((cats == lvl).to_numpy(dtype=float))
    Z = np.column_stack(cols)
    df_used = int(np.linalg.matrix_rank(Z) - 1)
    return Z, max(0, df_used)


def rank_residualize(mat_feat_by_sample: np.ndarray, Z: np.ndarray | None) -> np.ndarray:
    """Rank-transform each feature across samples, then residualize against design Z if given."""
    ranked = np.apply_along_axis(stats.rankdata, 1, mat_feat_by_sample).astype(float)
    if Z is not None:
        # R = Y - Y Z (Z^+)^T => transpose: S = ranked.T (n x p)
        S = ranked.T
        beta, _, _, _ = np.linalg.lstsq(Z, S, rcond=None)
        res = (S - Z @ beta).T
    else:
        res = ranked - np.mean(ranked, axis=1, keepdims=True)
    norms = np.linalg.norm(res, axis=1, keepdims=True)
    return res / np.where(norms > 1e-12, norms, 1.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix_x", required=True)
    parser.add_argument("--matrix_y", required=True)
    parser.add_argument("--name_x", default="block_X")
    parser.add_argument("--name_y", default="block_Y")
    parser.add_argument("--transform_x", choices=["auto", "none", "log1p", "clr", "zscore"], default="auto")
    parser.add_argument("--transform_y", choices=["auto", "none", "log1p", "clr", "zscore"], default="auto")
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default=None)
    parser.add_argument("--covariates", action="extend", nargs="+", default=[])
    parser.add_argument("--top_x", type=int, default=100)
    parser.add_argument("--top_y", type=int, default=100)
    parser.add_argument("--fdr_cutoff", type=float, default=0.05)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    covs = [c.strip() for item in args.covariates for c in str(item).split(",") if c.strip()]
    if covs and not args.metadata:
        fail("--covariates requires --metadata")

    raw_x = read_matrix(args.matrix_x)
    raw_y = read_matrix(args.matrix_y)
    shared = [s for s in raw_x.columns if s in set(raw_y.columns)]
    if len(shared) < 6:
        fail(f"mo_pair_assoc requires >= 6 shared samples between matrix_x and matrix_y; found {len(shared)}")

    tx_df, tx_used = transform_block(raw_x[shared], args.transform_x)
    ty_df, ty_used = transform_block(raw_y[shared], args.transform_y)

    meta_aligned = read_metadata_aligned(args.metadata, shared, args.sample_col) if args.metadata else None
    Z_cov, n_cov_df = (build_covariate_design(meta_aligned, covs) if (covs and meta_aligned is not None) else (None, 0))
    if len(shared) - 2 - n_cov_df < 3:
        fail(f"too few residual degrees of freedom ({len(shared) - 2 - n_cov_df}) after adjusting for {covs}")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # 1) Procrustes on top PCA axes (up to min(5, n-1))
    n_axes = min(5, len(shared) - 1, tx_df.shape[0], ty_df.shape[0])
    px = pca_coords(tx_df.to_numpy(dtype=float).T, n_axes)
    py = pca_coords(ty_df.to_numpy(dtype=float).T, n_axes)
    m2_obs, x_ord, y_ord = symmetric_procrustes(px, py)
    r_proc = float(np.sqrt(max(0.0, 1.0 - m2_obs)))
    perm_le = 0
    for _ in range(args.permutations):
        perm_idx = rng.permutation(len(shared))
        m2_p, _, _ = symmetric_procrustes(px, py[perm_idx])
        if m2_p <= m2_obs + 1e-12:
            perm_le += 1
    p_proc = float((perm_le + 1) / (args.permutations + 1))

    # 2) Mantel test on Euclidean pairwise sample distances in transformed space
    dx_vec = pdist(tx_df.to_numpy(dtype=float).T, metric="euclidean")
    dy_vec = pdist(ty_df.to_numpy(dtype=float).T, metric="euclidean")
    r_mantel_raw, _ = stats.spearmanr(dx_vec, dy_vec)
    r_mantel = float(r_mantel_raw) if np.isfinite(r_mantel_raw) else 0.0
    # Prebuild square distance matrix for fast row/col permutation
    from scipy.spatial.distance import squareform
    dy_sq = squareform(dy_vec)
    tri_i, tri_j = np.triu_indices(len(shared), k=1)
    mantel_ge = 0
    for _ in range(args.permutations):
        p_idx = rng.permutation(len(shared))
        dy_perm = dy_sq[np.ix_(p_idx, p_idx)][tri_i, tri_j]
        rp, _ = stats.spearmanr(dx_vec, dy_perm)
        if np.isfinite(rp) and rp >= r_mantel - 1e-12:
            mantel_ge += 1
    p_mantel = float((mantel_ge + 1) / (args.permutations + 1))

    # 3) Pairwise Spearman / Partial Spearman network on top-variance features
    var_x = tx_df.var(axis=1)
    var_y = ty_df.var(axis=1)
    keep_x = var_x[var_x > 1e-12].sort_values(ascending=False).head(args.top_x).index.tolist()
    keep_y = var_y[var_y > 1e-12].sort_values(ascending=False).head(args.top_y).index.tolist()
    if not keep_x or not keep_y:
        fail("one of the blocks has 0 non-constant features")

    ux = rank_residualize(tx_df.loc[keep_x].to_numpy(dtype=float), Z_cov)
    uy = rank_residualize(ty_df.loc[keep_y].to_numpy(dtype=float), Z_cov)
    rmat = np.clip(ux @ uy.T, -0.999999, 0.999999)
    df_t = max(1, len(shared) - 2 - n_cov_df)
    tmat = rmat * np.sqrt(df_t / np.maximum(1e-12, 1.0 - rmat * rmat))
    pmat = 2.0 * stats.t.sf(np.abs(tmat), df=df_t)
    qflat = bh_fdr(pmat.ravel()).reshape(pmat.shape)

    rows = []
    for i, fx in enumerate(keep_x):
        for j, fy in enumerate(keep_y):
            rows.append(
                {
                    "feature_x": fx,
                    "feature_y": fy,
                    "cor": round(float(rmat[i, j]), 4),
                    "pvalue": float(pmat[i, j]),
                    "fdr": float(qflat[i, j]),
                }
            )
    assoc_df = pd.DataFrame(rows).sort_values("pvalue", ascending=True)
    assoc_path = outdir / "cross_block_associations.tsv"
    assoc_df.to_csv(assoc_path, sep="\t", index=False)

    sig_df = assoc_df[assoc_df["fdr"] < args.fdr_cutoff]
    edges_path = outdir / "significant_network_edges.tsv"
    sig_df.to_csv(edges_path, sep="\t", index=False)

    # Procrustes residual arrow plot (2D)
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    x2 = x_ord[:, :2] if x_ord.shape[1] >= 2 else np.column_stack([x_ord[:, 0], np.zeros(len(shared))])
    y2 = y_ord[:, :2] if y_ord.shape[1] >= 2 else np.column_stack([y_ord[:, 0], np.zeros(len(shared))])
    for i in range(len(shared)):
        ax.annotate(
            "",
            xy=(y2[i, 0], y2[i, 1]),
            xytext=(x2[i, 0], x2[i, 1]),
            arrowprops={"arrowstyle": "->", "color": "#94A3B8", "lw": 0.9, "alpha": 0.75},
        )
    if args.group_col and meta_aligned is not None and args.group_col in meta_aligned.columns:
        grps = meta_aligned[args.group_col].astype(str)
        for g_idx, g in enumerate(sorted(grps.unique())):
            m = (grps == g).to_numpy()
            col = PALETTE[g_idx % len(PALETTE)]
            ax.scatter(x2[m, 0], x2[m, 1], color=col, marker="o", s=36, label=f"{g} ({args.name_x})")
            ax.scatter(y2[m, 0], y2[m, 1], color=col, marker="^", s=38, alpha=0.85)
    else:
        ax.scatter(x2[:, 0], x2[:, 1], color="#2563EB", marker="o", s=36, label=args.name_x)
        ax.scatter(y2[:, 0], y2[:, 1], color="#E11D48", marker="^", s=38, label=args.name_y)
    ax.set_xlabel("Procrustes Dimension 1")
    ax.set_ylabel("Procrustes Dimension 2")
    ax.set_title(
        f"Procrustes {args.name_x} → {args.name_y} (M²={m2_obs:.3f}, r={r_proc:.3f}, p={p_proc:.3g})",
        fontsize=9.8,
    )
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    proc_fig = outdir / "procrustes-arrows.png"
    fig.savefig(proc_fig, dpi=160)
    plt.close(fig)

    # Top cross-block correlation heatmap
    top_fx = keep_x[: min(20, len(keep_x))]
    top_fy = keep_y[: min(20, len(keep_y))]
    sub_r = pd.DataFrame(rmat, index=keep_x, columns=keep_y).loc[top_fx, top_fy]
    fig2, ax2 = plt.subplots(figsize=(max(5.5, 0.34 * len(top_fy) + 2.0), max(4.5, 0.28 * len(top_fx) + 1.8)))
    im = ax2.imshow(sub_r.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="auto")
    ax2.set_xticks(np.arange(len(top_fy)))
    ax2.set_xticklabels([s[:20] for s in top_fy], rotation=45, ha="right", fontsize=7.5)
    ax2.set_yticks(np.arange(len(top_fx)))
    ax2.set_yticklabels([s[:22] for s in top_fx], fontsize=7.5)
    kind_lbl = f"Partial Spearman (adj: {','.join(covs)})" if covs else "Spearman"
    ax2.set_title(f"{args.name_x} × {args.name_y} {kind_lbl} Heatmap", fontsize=9.8)
    fig2.colorbar(im, ax=ax2, fraction=0.035, pad=0.03, label="r")
    fig2.tight_layout()
    hm_fig = outdir / "cross-assoc-heatmap.png"
    fig2.savefig(hm_fig, dpi=160)
    plt.close(fig2)

    in_files = (args.matrix_x, args.matrix_y, args.metadata)
    params_meta = {
        "transform_x": tx_used,
        "transform_y": ty_used,
        "covariates": covs,
        "fdr_cutoff": args.fdr_cutoff,
        "permutations": args.permutations,
    }
    artifacts = [
        write_artifact(proc_fig, "mo_pair_assoc", f"Procrustes superimposition ({args.name_x} vs {args.name_y})", inputs=in_files, parameters=params_meta),
        write_artifact(hm_fig, "mo_pair_assoc", f"Cross-omics {kind_lbl} heatmap", inputs=in_files, parameters=params_meta),
        write_artifact(assoc_path, "mo_pair_assoc", "All-pairs cross-block correlation table", inputs=in_files, parameters=params_meta, shape=assoc_df.shape),
        write_artifact(edges_path, "mo_pair_assoc", f"FDR-significant cross-block network edges (n={len(sig_df)})", inputs=in_files, parameters=params_meta, shape=sig_df.shape),
    ]

    warnings: list[str] = []
    if len(shared) < 15:
        warnings.append(f"sample size n={len(shared)} (< 15): cross-block correlation estimates have wide confidence intervals")
    if not covs and args.group_col:
        warnings.append(
            f"group_col '{args.group_col}' was supplied for plotting without --covariates; pass --covariates {args.group_col} if you want within-group partial correlations"
        )

    emit(
        {
            "tool": "mo_pair_assoc",
            "n_shared_samples": len(shared),
            "block_x": {"name": args.name_x, "transform": tx_used, "n_features_tested": len(keep_x)},
            "block_y": {"name": args.name_y, "transform": ty_used, "n_features_tested": len(keep_y)},
            "covariates_adjusted": covs,
            "procrustes": {
                "m2": round(m2_obs, 4),
                "r": round(r_proc, 4),
                "pvalue": round(p_proc, 5),
                "n_axes": n_axes,
            },
            "mantel": {
                "r": round(r_mantel, 4),
                "pvalue": round(p_mantel, 5),
                "method": "spearman_euclidean",
            },
            "network": {
                "method": "partial_spearman" if covs else "spearman",
                "n_pairs_tested": int(len(assoc_df)),
                "fdr_cutoff": args.fdr_cutoff,
                "n_sig_edges": int(len(sig_df)),
                "top_edges": sig_df.head(10).to_dict(orient="records"),
            },
            "outputs": {
                "associations_table": str(assoc_path),
                "significant_edges_table": str(edges_path),
            },
            "figures": [str(proc_fig), str(hm_fig)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
