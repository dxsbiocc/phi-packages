#!/usr/bin/env python3
"""mo_subtype: Similarity Network Fusion (Wang et al., Nature Methods 2014) for multi-omics molecular subtyping."""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    parse_labeled_paths,
    read_matrix,
    read_metadata_aligned,
    transform_block,
    write_artifact,
)


def _scaled_exponential_affinity(X_samples_by_feat: np.ndarray, k_nn: int, mu: float = 0.5) -> np.ndarray:
    """Compute Wang et al. (2014) scaled exponential similarity matrix W."""
    n = X_samples_by_feat.shape[0]
    diff = X_samples_by_feat[:, None, :] - X_samples_by_feat[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    np.fill_diagonal(dist, 0.0)

    k_eff = max(1, min(k_nn, n - 1))
    d_sorted = np.sort(dist, axis=1)
    # Exclude self-distance at index 0; average over k_eff nearest neighbors
    eps_i = np.mean(d_sorted[:, 1 : 1 + k_eff], axis=1)
    eps_ij = (eps_i[:, None] + eps_i[None, :] + dist) / 3.0
    eps_ij = np.maximum(eps_ij, 1e-8)

    exponent = -(dist * dist) / (max(mu, 1e-4) * eps_ij)
    W = np.exp(np.clip(exponent, -60.0, 0.0))
    W = 0.5 * (W + W.T)
    np.fill_diagonal(W, 1.0)
    return W


def _normalize_full_kernel(W: np.ndarray) -> np.ndarray:
    """Build full SNF status matrix P with 0.5 self-weight and 0.5 row-normalized off-diagonal."""
    n = W.shape[0]
    P = W.copy()
    np.fill_diagonal(P, 0.0)
    row_sums = P.sum(axis=1, keepdims=True)
    off = np.where(row_sums > 1e-12, P / np.maximum(row_sums, 1e-12), 1.0 / max(1, n - 1))
    np.fill_diagonal(off, 0.0)
    out = 0.5 * off
    np.fill_diagonal(out, 0.5)
    return 0.5 * (out + out.T)


def _knn_sparse_kernel(W: np.ndarray, k_nn: int) -> np.ndarray:
    """Build local K-NN normalized affinity matrix S."""
    n = W.shape[0]
    k_eff = max(1, min(k_nn, n - 1))
    W_nodiag = W.copy()
    np.fill_diagonal(W_nodiag, 0.0)
    S = np.zeros_like(W_nodiag)
    for i in range(n):
        idx = np.argsort(W_nodiag[i])[::-1][:k_eff]
        vals = W_nodiag[i, idx]
        tot = float(np.sum(vals))
        if tot > 1e-12:
            S[i, idx] = vals / tot
        else:
            S[i, idx] = 1.0 / k_eff
    return S


def _spectral_embedding(W: np.ndarray, k_max: int) -> tuple[np.ndarray, np.ndarray]:
    """Compute symmetric normalized Laplacian eigenvectors & ascending eigenvalues."""
    n = W.shape[0]
    W_sym = 0.5 * (W + W.T)
    np.fill_diagonal(W_sym, 0.0)
    deg = np.sum(W_sym, axis=1)
    d_inv_sqrt = np.where(deg > 1e-12, 1.0 / np.sqrt(np.maximum(deg, 1e-12)), 0.0)
    L_norm = np.eye(n) - (d_inv_sqrt[:, None] * W_sym) * d_inv_sqrt[None, :]
    L_norm = 0.5 * (L_norm + L_norm.T)
    evals, evecs = np.linalg.eigh(L_norm)
    k_take = min(k_max, n)
    return evals[:k_take], evecs[:, :k_take]


def _spectral_labels(W: np.ndarray, k: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    _, evecs = _spectral_embedding(W, k)
    norms = np.linalg.norm(evecs, axis=1, keepdims=True)
    U = np.where(norms > 1e-12, evecs / np.maximum(norms, 1e-12), 0.0)
    km = KMeans(n_clusters=k, n_init=20, random_state=seed)
    labels = km.fit_predict(U) + 1
    return labels, U


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrices", action="extend", nargs="+", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--transforms", action="extend", nargs="+", default=None)
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--sample_col", default=None)
    ap.add_argument("--group_col", default=None)
    ap.add_argument("--k_clusters", type=int, default=None)
    ap.add_argument("--n_neighbors", type=int, default=15)
    ap.add_argument("--n_iters", type=int, default=15)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    labeled = parse_labeled_paths(args.matrices, prefix="modality")
    if len(labeled) < 2:
        fail("mo_subtype requires >= 2 omics feature x sample matrices.")
    if args.k_clusters is not None and not (2 <= args.k_clusters <= 8):
        fail(f"k_clusters must be in [2, 8] (got {args.k_clusters}).")
    if args.n_iters < 1:
        fail("n_iters must be >= 1.")

    raw_blocks = [(lbl, read_matrix(p), str(p)) for lbl, p in labeled]
    shared_samples = list(raw_blocks[0][1].columns)
    for _, df_b, _ in raw_blocks[1:]:
        col_set = set(df_b.columns)
        shared_samples = [s for s in shared_samples if s in col_set]

    n_shared = len(shared_samples)
    if n_shared < 6:
        fail(f"mo_subtype requires >= 6 shared samples across all modalities (found {n_shared}).")

    tr_specs: list[str] = []
    if args.transforms:
        for item in args.transforms:
            for tok in str(item).split(","):
                tok = tok.strip()
                if tok:
                    tr_specs.append(tok)

    meta_df = None
    if args.metadata:
        meta_df = read_metadata_aligned(args.metadata, shared_samples, sample_col=args.sample_col)
        if args.group_col and args.group_col not in meta_df.columns:
            fail(f"group_col '{args.group_col}' not found in metadata columns: {list(meta_df.columns)}")

    k_nn = max(3, min(int(args.n_neighbors), max(3, n_shared - 2)))
    warnings: list[str] = []
    if k_nn < int(args.n_neighbors):
        warnings.append(f"Clamped n_neighbors from {args.n_neighbors} to {k_nn} for N={n_shared} shared samples.")

    W_list: list[np.ndarray] = []
    P_list: list[np.ndarray] = []
    S_list: list[np.ndarray] = []
    mod_meta: list[dict[str, Any]] = []

    for idx, (lbl, df_b, raw_p) in enumerate(raw_blocks):
        sub = df_b[shared_samples]
        t_req = tr_specs[idx] if idx < len(tr_specs) else "auto"
        t_df, t_used = transform_block(sub, t_req)
        x = t_df.to_numpy(dtype=float, copy=True).T  # samples x features
        var_f = np.var(x, axis=0)
        keep_f = var_f > 1e-12
        if int(keep_f.sum()) < 2:
            fail(f"Modality '{lbl}' has < 2 non-constant features across shared samples.")
        x = x[:, keep_f]
        mu = np.mean(x, axis=0, keepdims=True)
        sd = np.std(x, axis=0, ddof=1, keepdims=True)
        x_z = (x - mu) / np.maximum(sd, 1e-8)

        W_m = _scaled_exponential_affinity(x_z, k_nn=k_nn)
        W_list.append(W_m)
        P_list.append(_normalize_full_kernel(W_m))
        S_list.append(_knn_sparse_kernel(W_m, k_nn=k_nn))
        mod_meta.append(
            {
                "modality": lbl,
                "path": raw_p,
                "n_features_used": int(x_z.shape[1]),
                "transform": t_used,
            }
        )

    m_blocks = len(W_list)
    P_cur = [p.copy() for p in P_list]
    for _ in range(int(args.n_iters)):
        P_next = []
        for m in range(m_blocks):
            other_avg = np.mean([P_cur[v] for v in range(m_blocks) if v != m], axis=0)
            P_m = S_list[m] @ other_avg @ S_list[m].T
            P_m = _normalize_full_kernel(P_m)
            P_next.append(P_m)
        P_cur = P_next

    W_fused = np.mean(P_cur, axis=0)
    W_fused = 0.5 * (W_fused + W_fused.T)
    np.fill_diagonal(W_fused, 1.0)

    # Distance matrix for silhouette score
    w_off_max = float(np.max(W_fused - np.diag(np.diag(W_fused)))) if n_shared > 1 else 1.0
    w_scale = max(w_off_max, 1e-8)
    dist_fused = np.clip(1.0 - (W_fused / w_scale), 0.0, 1.0)
    np.fill_diagonal(dist_fused, 0.0)

    max_k_cand = min(6, max(2, n_shared - 2))
    evals, _ = _spectral_embedding(W_fused, max_k_cand + 2)
    k_candidates: list[dict[str, Any]] = []
    if args.k_clusters is not None:
        k_opt = min(int(args.k_clusters), max(2, n_shared - 1))
        labels, emb_u = _spectral_labels(W_fused, k_opt, int(args.seed))
        sil_opt = float(silhouette_score(dist_fused, labels, metric="precomputed"))
        k_candidates.append({"k": k_opt, "silhouette": round(sil_opt, 4), "eigengap": None})
        selection_mode = "user_specified"
    else:
        best_score = -1e9
        k_opt = 2
        labels = np.ones(n_shared, dtype=int)
        emb_u = np.zeros((n_shared, 2), dtype=float)
        sil_opt = 0.0
        for kc in range(2, max_k_cand + 1):
            lbl_c, u_c = _spectral_labels(W_fused, kc, int(args.seed))
            if len(np.unique(lbl_c)) < 2 or len(np.unique(lbl_c)) >= n_shared:
                continue
            sil_c = float(silhouette_score(dist_fused, lbl_c, metric="precomputed"))
            egap = float(evals[kc] - evals[kc - 1]) if kc < len(evals) else 0.0
            comb = sil_c + 0.15 * max(0.0, egap)
            k_candidates.append({"k": kc, "silhouette": round(sil_c, 4), "eigengap": round(egap, 5)})
            if comb > best_score:
                best_score = comb
                k_opt = kc
                labels = lbl_c
                emb_u = u_c
                sil_opt = sil_c
        selection_mode = "auto_eigengap_silhouette"

    # Compute per-modality NMI concordance against the fused subtypes
    for m_idx, m_info in enumerate(mod_meta):
        lbl_m, _ = _spectral_labels(W_list[m_idx], k_opt, int(args.seed))
        m_info["modality_concordance_nmi"] = round(float(normalized_mutual_info_score(labels, lbl_m)), 4)

    # Group / external label validation if provided
    group_assoc = None
    if meta_df is not None and args.group_col:
        ext_vals = meta_df[args.group_col].astype(str).to_numpy(copy=True)
        ari = float(adjusted_rand_score(ext_vals, labels))
        nmi = float(normalized_mutual_info_score(ext_vals, labels))
        ctab = pd.crosstab(pd.Series(labels, name="snf_subtype"), pd.Series(ext_vals, name=args.group_col))
        chi2_stat, chi2_p, _, _ = stats.chi2_contingency(ctab)
        group_assoc = {
            "group_col": args.group_col,
            "ari": round(ari, 4),
            "nmi": round(nmi, 4),
            "chi2_stat": round(float(chi2_stat), 4),
            "chi2_pvalue": float(chi2_p),
        }

    # Build subtype summary counts
    subtype_names = [f"Subtype_{int(lb)}" for lb in labels]
    counts_ser = pd.Series(subtype_names).value_counts().sort_index()
    subtype_counts = {str(k): int(v) for k, v in counts_ser.items()}

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    assign_df = pd.DataFrame(
        {
            "sample": shared_samples,
            "snf_cluster": labels.astype(int),
            "snf_subtype": subtype_names,
            "spectral_dim1": np.round(emb_u[:, 0], 6),
            "spectral_dim2": np.round(emb_u[:, 1] if emb_u.shape[1] > 1 else np.zeros(n_shared), 6),
        }
    )
    if meta_df is not None and args.group_col:
        assign_df[args.group_col] = meta_df[args.group_col].astype(str).to_numpy(copy=True)

    assign_path = out_dir / "snf_cluster_assignments.tsv"
    assign_df.to_csv(assign_path, sep="\t", index=False)

    fused_df = pd.DataFrame(np.round(W_fused, 6), index=shared_samples, columns=shared_samples)
    fused_df.index.name = "sample"
    fused_path = out_dir / "snf_fused_similarity.tsv"
    fused_df.to_csv(fused_path, sep="\t")

    # 2-panel diagnostic figure: ordered SNF similarity matrix + spectral embedding
    order_idx = np.argsort(labels, kind="mergesort")
    W_ord = W_fused[np.ix_(order_idx, order_idx)].copy()
    np.fill_diagonal(W_ord, np.nan)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.8), dpi=160)
    im = ax1.imshow(W_ord, cmap="YlGnBu", aspect="auto")
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04, label="Fused SNF Similarity")
    ax1.set_title(f"Fused SNF Similarity Matrix (K={k_opt}, Sil={sil_opt:.3f})")
    ax1.set_xlabel("Samples (ordered by SNF subtype)")
    ax1.set_ylabel("Samples (ordered by SNF subtype)")

    palette = ["#2563eb", "#dc2626", "#059669", "#d97706", "#7c3aed", "#0891b2", "#be123c", "#475569"]
    u2 = emb_u[:, 1] if emb_u.shape[1] > 1 else np.zeros(n_shared)
    for c_id in sorted(np.unique(labels)):
        m_c = labels == c_id
        ax2.scatter(
            emb_u[m_c, 0],
            u2[m_c],
            color=palette[(int(c_id) - 1) % len(palette)],
            s=45,
            alpha=0.9,
            label=f"Subtype_{int(c_id)} (n={int(m_c.sum())})",
        )
    ax2.set_xlabel("Normalized Laplacian Eigenvector 1")
    ax2.set_ylabel("Normalized Laplacian Eigenvector 2")
    ax2.set_title("SNF Spectral Embedding")
    ax2.legend(frameon=False, fontsize=8.2)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    fig.tight_layout()
    fig_path = out_dir / "snf_fused_heatmap_embedding.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    inputs_all = [m["path"] for m in mod_meta] + ([args.metadata] if args.metadata else [])
    params_meta = {
        "k_clusters": k_opt,
        "selection_mode": selection_mode,
        "n_neighbors": k_nn,
        "n_iters": int(args.n_iters),
    }
    artifacts = [
        write_artifact(assign_path, "mo_subtype", "SNF Multi-Omics Subtype Assignments", inputs=inputs_all, parameters=params_meta, shape=assign_df.shape),
        write_artifact(fused_path, "mo_subtype", "SNF Fused Sample Similarity Matrix", inputs=inputs_all, parameters=params_meta, shape=fused_df.shape),
        write_artifact(fig_path, "mo_subtype", "SNF Fused Similarity Heatmap & Spectral Embedding", inputs=inputs_all, parameters=params_meta),
    ]

    emit(
        {
            "status": "ok",
            "tool": "mo_subtype",
            "n_shared_samples": n_shared,
            "n_modalities": m_blocks,
            "k_clusters": int(k_opt),
            "k_selection_mode": selection_mode,
            "silhouette_score": round(sil_opt, 4),
            "subtype_counts": subtype_counts,
            "k_candidates": k_candidates,
            "modalities": mod_meta,
            "group_association": group_assoc,
            "outputs": {
                "cluster_assignments": str(assign_path),
                "fused_similarity_matrix": str(fused_path),
                "snf_plot": str(fig_path),
            },
            "figures": [str(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
