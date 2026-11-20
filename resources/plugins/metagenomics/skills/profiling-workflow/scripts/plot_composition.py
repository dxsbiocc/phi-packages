#!/usr/bin/env python3
"""mg_plot: stacked composition barplots and taxa heatmaps from a matrix."""
import argparse
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from _common import emit, fail, read_matrix, read_metadata, setup_matplotlib, write_artifact


def cluster_order(vals):
    if vals.shape[0] < 2:
        return None
    return leaves_list(linkage(pdist(vals), method="average"))


def top_taxa(matrix, top):
    means = matrix.mean(axis=1).sort_values(ascending=False)
    keep = means.head(top).index
    other = matrix.loc[~matrix.index.isin(keep)].sum(axis=0)
    plotted = matrix.loc[keep]
    if (other > 0).any():
        plotted.loc["Other"] = other
    return plotted


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--kind", choices=["barplot", "heatmap"], default="barplot")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--metadata")
    ap.add_argument("--group_col", default="group")
    ap.add_argument("--trans", choices=["none", "log10"], default="none")
    ap.add_argument("--pseudocount", type=float, default=1e-6)
    ap.add_argument("--cluster", choices=["none", "rows", "cols", "both"], default="none")
    args = ap.parse_args()

    if args.cluster != "none" and args.kind != "heatmap":
        fail("--cluster only applies to --kind heatmap")

    matrix = read_matrix(args.input)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotted = top_taxa(matrix, args.top)
    plt = setup_matplotlib()
    artifacts = []
    applied = []
    inputs = [args.input] + ([args.metadata] if args.metadata else [])
    params = {"kind": args.kind, "top": args.top, "trans": args.trans, "cluster": args.cluster}

    if args.kind == "barplot":
        fig, ax = plt.subplots(figsize=(max(7, len(plotted.columns) * 0.6), 5))
        bottom = np.zeros(len(plotted.columns))
        cmap = plt.get_cmap("tab20")
        for i, (taxon, row) in enumerate(plotted.iterrows()):
            vals = row.to_numpy()
            ax.bar(plotted.columns, vals, bottom=bottom, label=taxon,
                   color=cmap(i % 20), width=0.8)
            bottom += vals
        ax.set_ylabel("abundance")
        ax.set_title(f"Composition (top {args.top} taxa)")
        ax.tick_params(axis="x", rotation=45)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.02, 1),
                  loc="upper left", fontsize=7)
        fig.tight_layout()
        p = outdir / "composition-barplot.png"
        fig.savefig(p, dpi=150, bbox_inches="tight")
        plt.close(fig)
        artifacts.append(write_artifact(p, "mg_plot", f"Composition barplot (top {args.top} taxa)",
                                        skill="profiling-workflow", inputs=inputs, parameters=params))

        if args.metadata:
            meta = read_metadata(args.metadata, matrix.columns.tolist())
            if args.group_col not in meta.columns:
                fail(f"metadata has no column {args.group_col}")
            groups = meta[args.group_col].astype(str)
            gmean = plotted.T.groupby(groups).mean().T
            fig, ax = plt.subplots(figsize=(max(5, len(gmean.columns) * 0.8), 5))
            bottom = np.zeros(len(gmean.columns))
            for i, (taxon, row) in enumerate(gmean.iterrows()):
                vals = row.to_numpy()
                ax.bar(gmean.columns, vals, bottom=bottom, label=taxon,
                       color=cmap(i % 20), width=0.6)
                bottom += vals
            ax.set_ylabel("mean abundance")
            ax.set_title(f"Composition by {args.group_col}")
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.02, 1),
                      loc="upper left", fontsize=7)
            fig.tight_layout()
            p = outdir / "composition-barplot-groups.png"
            fig.savefig(p, dpi=150, bbox_inches="tight")
            plt.close(fig)
            artifacts.append(write_artifact(p, "mg_plot", f"Composition by {args.group_col} (group means)",
                                            skill="profiling-workflow", inputs=inputs, parameters=params))
    else:
        vals = plotted.to_numpy()
        if args.trans == "log10":
            vals = np.log10(vals + args.pseudocount)
        if args.cluster in ("rows", "both"):
            order = cluster_order(vals)
            if order is not None:
                vals, plotted = vals[order], plotted.iloc[order]
                applied.append("rows")
        if args.cluster in ("cols", "both"):
            order = cluster_order(vals.T)
            if order is not None:
                vals, plotted = vals[:, order], plotted.iloc[:, order]
                applied.append("cols")
        fig, ax = plt.subplots(figsize=(max(7, len(plotted.columns) * 0.5),
                                        max(5, len(plotted) * 0.3)))
        im = ax.imshow(vals, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(plotted.columns)), plotted.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(plotted)), plotted.index, fontsize=7)
        fig.colorbar(im, ax=ax, label="log10(abundance)" if args.trans == "log10" else "abundance")
        suffix = " (clustered)" if applied else ""
        ax.set_title(f"Top {args.top} taxa{suffix}")
        fig.tight_layout()
        p = outdir / "taxa-heatmap.png"
        fig.savefig(p, dpi=150, bbox_inches="tight")
        plt.close(fig)
        artifacts.append(write_artifact(p, "mg_plot",
                                        f"Taxa heatmap (top {args.top}{', clustered' if applied else ''})",
                                        skill="profiling-workflow", inputs=inputs, parameters=params))

    emit({
        "kind": args.kind,
        "top": args.top,
        "samples": matrix.columns.tolist(),
        "taxa_plotted": int(len(plotted)),
        "cluster": {
            "mode": args.cluster,
            "method": "average",
            "metric": "euclidean",
            "applied": applied,
        },
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()