#!/usr/bin/env python3
"""prot_network: protein co-abundance correlation network, modules, and hub scoring."""

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
from scipy.cluster.hierarchy import leaves_list, linkage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    bh_fdr,
    clean_sample_name,
    detect_format_and_columns,
    emit,
    fail,
    read_table,
)


def connected_components(nodes: list[str], edges: list[tuple[str, str]]) -> dict[str, int]:
    adj: dict[str, set[str]] = {n: set() for n in nodes}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    visited: set[str] = set()
    comp_map: dict[str, int] = {}
    cid = 0
    for n in nodes:
        if n in visited:
            continue
        cid += 1
        stack = [n]
        visited.add(n)
        while stack:
            cur = stack.pop()
            comp_map[cur] = cid
            for nxt in adj[cur]:
                if nxt not in visited:
                    visited.add(nxt)
                    stack.append(nxt)
    return comp_map


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--de_table", default=None)
    parser.add_argument("--id_col", default=None)
    parser.add_argument("--top_n", type=int, default=80)
    parser.add_argument("--method", choices=["spearman", "pearson"], default="spearman")
    parser.add_argument("--min_abs_corr", type=float, default=0.7)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.input)
    det = detect_format_and_columns(df, args.id_col)
    id_col = det["id_col"]
    sample_cols = det["sample_cols"]

    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    mat.columns = [clean_sample_name(c) for c in sample_cols]
    mat.index = df[id_col].astype(str).str.strip()
    mat = mat[~mat.index.duplicated(keep="first")].dropna(axis=0)
    if mat.shape[0] < 4 or mat.shape[1] < 4:
        fail("co-abundance network requires >= 4 complete proteins and >= 4 samples")

    # Prioritize regulated proteins if de_table is provided, else top variance
    if args.de_table:
        de = read_table(args.de_table)
        did = "protein_id" if "protein_id" in de.columns else de.columns[0]
        de_ids = [str(x).strip() for x in de.sort_values("pvalue")[did] if str(x).strip() in mat.index]
        chosen = list(dict.fromkeys(de_ids))[: args.top_n]
        if len(chosen) < 4:
            chosen = list(mat.var(axis=1).sort_values(ascending=False).head(args.top_n).index)
    else:
        chosen = list(mat.var(axis=1).sort_values(ascending=False).head(args.top_n).index)

    sub = mat.loc[chosen]
    X = sub.to_numpy(dtype=float)
    n_p = X.shape[0]

    if args.method == "spearman":
        corr_mat, p_mat = stats.spearmanr(X, axis=1)
        if np.ndim(corr_mat) == 0:
            corr_mat = np.array([[1.0, float(corr_mat)], [float(corr_mat), 1.0]])
            p_mat = np.array([[0.0, float(p_mat)], [float(p_mat), 0.0]])
    else:
        corr_mat = np.corrcoef(X)
        n_s = X.shape[1]
        r_clip = np.clip(corr_mat, -0.999999, 0.999999)
        t_val = r_clip * np.sqrt((n_s - 2) / np.maximum(1.0 - r_clip ** 2, 1e-12))
        p_mat = 2.0 * stats.t.sf(np.abs(t_val), df=n_s - 2)
        np.fill_diagonal(p_mat, 0.0)

    iu, ju = np.triu_indices(n_p, k=1)
    r_vals = corr_mat[iu, ju]
    p_vals = p_mat[iu, ju]
    padj_vals = bh_fdr(p_vals)

    keep = np.isfinite(r_vals) & (np.abs(r_vals) >= args.min_abs_corr) & (padj_vals <= args.fdr)
    edge_rows = []
    edge_pairs = []
    strength = {p: 0.0 for p in chosen}
    degree = {p: 0 for p in chosen}
    for idx in np.flatnonzero(keep):
        u = chosen[int(iu[idx])]
        v = chosen[int(ju[idx])]
        rv = float(r_vals[idx])
        edge_pairs.append((u, v))
        strength[u] += abs(rv)
        strength[v] += abs(rv)
        degree[u] += 1
        degree[v] += 1
        edge_rows.append(
            {
                "protein_a": u,
                "protein_b": v,
                "corr": round(rv, 4),
                "pvalue": float(p_vals[idx]),
                "padj": float(padj_vals[idx]),
            }
        )

    edges_df = pd.DataFrame(edge_rows)
    if len(edges_df):
        edges_df = edges_df.sort_values("corr", key=lambda s: s.abs(), ascending=False)
    else:
        edges_df = pd.DataFrame(columns=["protein_a", "protein_b", "corr", "pvalue", "padj"])
    edges_path = os.path.join(args.output_dir, "coabundance_edges.tsv")
    edges_df.to_csv(edges_path, sep="\t", index=False)

    comp = connected_components(chosen, edge_pairs)
    hub_df = pd.DataFrame(
        [
            {
                "protein_id": p,
                "module": int(comp[p]),
                "degree": int(degree[p]),
                "strength": round(float(strength[p]), 4),
            }
            for p in chosen
        ]
    ).sort_values(["strength", "degree"], ascending=False)
    hubs_path = os.path.join(args.output_dir, "protein_hubs.tsv")
    hub_df.to_csv(hubs_path, sep="\t", index=False)

    # Clustered correlation heatmap
    corr_clean = np.nan_to_num(corr_mat, nan=0.0)
    np.fill_diagonal(corr_clean, 1.0)
    order = np.arange(n_p)
    if n_p >= 3:
        try:
            Z = linkage(1.0 - np.clip(corr_clean, -1.0, 1.0), method="average")
            order = leaves_list(Z)
        except Exception:
            pass
    ordered_corr = corr_clean[np.ix_(order, order)]
    fig, ax = plt.subplots(figsize=(6.0, 5.5))
    im = ax.imshow(ordered_corr, cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="auto")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=f"{args.method} r")
    ax.set_title(f"Protein co-abundance ({n_p} proteins, {len(edges_df)} edges)")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    hm_path = os.path.join(args.output_dir, "coabundance-heatmap.png")
    fig.savefig(hm_path, dpi=150)
    plt.close(fig)

    n_nontrivial_modules = int(sum(1 for _, grp in hub_df.groupby("module") if len(grp) >= 2))
    emit(
        {
            "input": os.path.abspath(args.input),
            "edges_table": os.path.abspath(edges_path),
            "hubs_table": os.path.abspath(hubs_path),
            "method": args.method,
            "n_proteins": n_p,
            "n_edges": int(len(edges_df)),
            "n_modules_ge2": n_nontrivial_modules,
            "top_hubs": [
                {
                    "protein_id": str(r["protein_id"]),
                    "module": int(r["module"]),
                    "degree": int(r["degree"]),
                    "strength": float(r["strength"]),
                }
                for _, r in hub_df.head(10).iterrows()
            ],
            "figures": [os.path.abspath(hm_path)],
        }
    )


if __name__ == "__main__":
    main()
