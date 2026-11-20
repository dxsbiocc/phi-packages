#!/usr/bin/env python3
"""mb_plot: clustered heatmap of a feature x sample matrix."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from _common import emit, fail, read_matrix, setup_matplotlib, write_artifact

SKILL = "lcms-workflow"


def order(vals):
    if vals.shape[0] < 2:
        return np.arange(vals.shape[0])
    return leaves_list(linkage(pdist(vals), method="average"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--trans", choices=["none", "log10"], default="none")
    ap.add_argument("--pseudocount", type=float, default=1e-6)
    ap.add_argument("--cluster", choices=["none", "rows", "cols", "both"],
                    default="none")
    args = ap.parse_args()

    meta, X = read_matrix(args.input)
    if X.shape[1] < 2:
        fail("heatmap needs at least two sample columns")
    X = X.dropna(how="all")
    if X.empty:
        fail("no features with observed values")
    var = X.var(axis=1).fillna(0)
    keep = var.sort_values(ascending=False).head(args.top).index
    M = X.loc[keep].fillna(0.0)
    if args.trans == "log10":
        if (M < 0).any().any():
            fail("log10 transform needs non-negative values")
        M = np.log10(M + args.pseudocount)
    vals = M.to_numpy(dtype=float)
    row_order = order(vals) if args.cluster in ("rows", "both") else np.arange(len(M))
    col_order = order(vals.T) if args.cluster in ("cols", "both") else np.arange(M.shape[1])
    vals = vals[row_order][:, col_order]
    rows = [str(M.index[i]) for i in row_order]
    cols = [str(M.columns[j]) for j in col_order]

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(cols) + 2),
                                    max(5, 0.22 * len(rows) + 1.5)))
    im = ax.imshow(vals, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90, fontsize=7)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=6)
    ax.set_title(f"Feature heatmap (top {len(rows)} by variance"
                 + (f", {args.trans}" if args.trans != "none" else "")
                 + (", clustered" if args.cluster != "none" else "") + ")")
    fig.colorbar(im, ax=ax, shrink=0.6)
    fig.tight_layout()
    fig_path = out / "feature-heatmap.png"
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)

    params = {"top": args.top, "trans": args.trans,
              "pseudocount": args.pseudocount, "cluster": args.cluster}
    artifacts = [
        write_artifact(fig_path, "mb_plot", "Feature heatmap",
                       skill=SKILL, inputs=[args.input], parameters=params),
    ]
    emit({
        "input": args.input,
        "kind": "heatmap",
        "n_features_shown": len(rows),
        "n_samples": len(cols),
        "trans": args.trans,
        "cluster": args.cluster,
        "output_dir": str(out),
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()