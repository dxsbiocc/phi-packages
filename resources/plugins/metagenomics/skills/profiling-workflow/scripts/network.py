#!/usr/bin/env python3
"""mg_network: Spearman co-occurrence network analysis across taxa or
functional features with BH FDR control, hub detection, and network figure.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from _common import emit, fail, read_matrix, setup_matplotlib, write_artifact


def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    if m == 0:
        return np.array([], dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out = np.empty(m, dtype=float)
    out[order] = q
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--min_prevalence", type=float, default=0.2)
    ap.add_argument("--top", type=int, default=100)
    ap.add_argument("--min_abs_rho", type=float, default=0.4)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--features", nargs="+", help="Optional explicit subset of feature names to analyze.")
    args = ap.parse_args()

    matrix = read_matrix(args.input)  # features x samples
    n_samples = matrix.shape[1]
    if n_samples < 5:
        fail(f"co-occurrence network needs at least 5 samples (got {n_samples})")

    df = matrix.copy()
    if args.features:
        wanted = []
        for item in args.features:
            for part in str(item).split(","):
                part = part.strip()
                if part:
                    wanted.append(part)
        present = [f for f in wanted if f in df.index]
        if len(present) < 2:
            fail(f"fewer than 2 requested features found in matrix ({len(present)} matched)")
        df = df.loc[present]
    else:
        prev = (df > 0).sum(axis=1) / float(n_samples)
        df = df.loc[prev >= args.min_prevalence]
        if args.top and len(df) > args.top:
            order = df.mean(axis=1).sort_values(ascending=False).index[: args.top]
            df = df.loc[order]

    # Drop zero-variance rows
    df = df.loc[df.std(axis=1) > 0]
    n_features = len(df)
    if n_features < 2:
        fail(f"fewer than 2 variable features remain after filtering (got {n_features})")

    feat_names = df.index.tolist()
    x = df.to_numpy(dtype=float)  # features x samples
    rho_mat, p_mat = spearmanr(x, axis=1)
    if np.ndim(rho_mat) == 0:
        rho_mat = np.array([[1.0, float(rho_mat)], [float(rho_mat), 1.0]])
        p_mat = np.array([[0.0, float(p_mat)], [float(p_mat), 0.0]])

    pairs = []
    pvals = []
    for i in range(n_features):
        for j in range(i + 1, n_features):
            r = float(rho_mat[i, j])
            pv = float(p_mat[i, j])
            if np.isnan(r) or np.isnan(pv):
                continue
            pairs.append((i, j, r, pv))
            pvals.append(pv)

    qvals = bh_adjust(pvals)
    edges = []
    for (i, j, r, pv), qv in zip(pairs, qvals):
        if abs(r) >= args.min_abs_rho and qv <= args.fdr:
            edges.append(
                {
                    "source": feat_names[i],
                    "target": feat_names[j],
                    "rho": round(r, 4),
                    "abs_rho": round(abs(r), 4),
                    "direction": "positive" if r > 0 else "negative",
                    "p_value": float(pv),
                    "q_value": float(qv),
                }
            )
    edges.sort(key=lambda e: (-e["abs_rho"], e["q_value"], e["source"], e["target"]))

    node_stats = {
        f: {
            "feature": f,
            "degree": 0,
            "pos_edges": 0,
            "neg_edges": 0,
            "weighted_degree": 0.0,
            "mean_abundance": float(df.loc[f].mean()),
            "prevalence": float((df.loc[f] > 0).mean()),
        }
        for f in feat_names
    }
    for e in edges:
        for end in (e["source"], e["target"]):
            node_stats[end]["degree"] += 1
            node_stats[end]["weighted_degree"] += e["abs_rho"]
            if e["direction"] == "positive":
                node_stats[end]["pos_edges"] += 1
            else:
                node_stats[end]["neg_edges"] += 1

    connected_nodes = [v for v in node_stats.values() if v["degree"] > 0]
    connected_nodes.sort(key=lambda v: (-v["degree"], -v["weighted_degree"], v["feature"]))
    for v in connected_nodes:
        v["weighted_degree"] = round(v["weighted_degree"], 4)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    params = {
        "min_prevalence": args.min_prevalence,
        "top": args.top,
        "min_abs_rho": args.min_abs_rho,
        "fdr": args.fdr,
    }
    edges_df = pd.DataFrame(
        edges,
        columns=["source", "target", "rho", "abs_rho", "direction", "p_value", "q_value"],
    )
    edges_path = outdir / "network-edges.tsv"
    edges_df.to_csv(edges_path, sep="\t", index=False)

    nodes_df = pd.DataFrame(
        connected_nodes,
        columns=["feature", "degree", "pos_edges", "neg_edges", "weighted_degree", "mean_abundance", "prevalence"],
    )
    nodes_path = outdir / "network-nodes.tsv"
    nodes_df.to_csv(nodes_path, sep="\t", index=False)

    artifacts = [
        write_artifact(
            edges_path,
            "mg_network",
            "Spearman co-occurrence network edge table",
            skill="profiling-workflow",
            inputs=[args.input],
            parameters=params,
            shape=edges_df.shape,
        ),
        write_artifact(
            nodes_path,
            "mg_network",
            "Co-occurrence network node degree and hub table",
            skill="profiling-workflow",
            inputs=[args.input],
            parameters=params,
            shape=nodes_df.shape,
        ),
    ]
    outputs = {"edges": str(edges_path), "nodes": str(nodes_path)}

    # Render circular co-occurrence graph when edges exist
    if connected_nodes:
        plt = setup_matplotlib()
        fig, ax = plt.subplots(figsize=(7, 7))
        n_conn = len(connected_nodes)
        angles = np.linspace(0, 2 * np.pi, n_conn, endpoint=False)
        pos = {
            node["feature"]: (float(np.cos(ang)), float(np.sin(ang)))
            for node, ang in zip(connected_nodes, angles)
        }
        for e in edges:
            x0, y0 = pos[e["source"]]
            x1, y1 = pos[e["target"]]
            color = "#0EA5E9" if e["direction"] == "positive" else "#E11D48"
            lw = 0.8 + 2.5 * e["abs_rho"]
            ax.plot([x0, x1], [y0, y1], color=color, linewidth=lw, alpha=0.55, zorder=1)
        max_deg = max(n["degree"] for n in connected_nodes)
        for node in connected_nodes:
            fx, fy = pos[node["feature"]]
            size = 90 + 360 * (node["degree"] / max(1, max_deg))
            ax.scatter([fx], [fy], s=size, color="#1E293B", edgecolors="#38BDF8", linewidths=1.5, zorder=2)
            label = node["feature"]
            if len(label) > 26:
                label = label[:23] + "..."
            ax.text(fx * 1.14, fy * 1.14, label, fontsize=7, ha="center", va="center")
        ax.set_xlim(-1.4, 1.4)
        ax.set_ylim(-1.4, 1.4)
        ax.axis("off")
        ax.set_title(
            f"Spearman co-occurrence network (|rho|>={args.min_abs_rho}, q<={args.fdr}; "
            f"{n_conn} nodes, {len(edges)} edges)"
        )
        fig.tight_layout()
        fig_path = outdir / "cooccurrence-network.png"
        fig.savefig(fig_path, dpi=200)
        plt.close(fig)
        outputs["network_figure"] = str(fig_path)
        artifacts.append(
            write_artifact(
                fig_path,
                "mg_network",
                "Spearman co-occurrence network plot",
                skill="profiling-workflow",
                inputs=[args.input],
                parameters=params,
            )
        )

    n_pos = sum(1 for e in edges if e["direction"] == "positive")
    n_neg = sum(1 for e in edges if e["direction"] == "negative")
    emit(
        {
            "tool": "mg_network",
            "n_samples": int(n_samples),
            "n_features_tested": int(n_features),
            "n_pairs_tested": int(len(pairs)),
            "n_nodes_connected": int(len(connected_nodes)),
            "n_edges": int(len(edges)),
            "n_positive_edges": int(n_pos),
            "n_negative_edges": int(n_neg),
            "top_hubs": connected_nodes[:10],
            "outputs": outputs,
            "artifacts": artifacts,
            "caveat": "Spearman correlations on relative abundances can reflect compositional closure "
                      "effects; interpret negative edges cautiously and verify key hubs with "
                      "compositionally aware methods (e.g., SparCC/SPIEC-EASI) when needed.",
        }
    )


if __name__ == "__main__":
    main()
