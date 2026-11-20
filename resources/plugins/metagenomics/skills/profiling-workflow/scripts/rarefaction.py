#!/usr/bin/env python3
"""mg_rarefaction: observed-richness rarefaction curves for a count matrix."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import emit, fail, read_matrix, read_metadata, setup_matplotlib


def subsample_richness(counts, depth, rng):
    """Richness after drawing depth individuals without replacement.

    Sequential hypergeometric draws over taxa: equivalent to subsampling the
    whole community, at O(n_taxa) per draw instead of O(depth).
    """
    remaining = depth
    total = int(counts.sum())
    richness = 0
    for c in counts:
        if remaining == 0:
            break
        c = int(c)
        if c == 0:
            continue
        k = int(rng.hypergeometric(c, total - c, remaining))
        if k > 0:
            richness += 1
        remaining -= k
        total -= c
    return richness


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--points", type=int, default=10)
    ap.add_argument("--iterations", type=int, default=10)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--metadata")
    ap.add_argument("--group_col", default="group")
    args = ap.parse_args()
    if args.points < 2:
        fail("--points must be >= 2")
    if args.iterations < 1:
        fail("--iterations must be >= 1")

    matrix = read_matrix(args.input)
    arr = matrix.to_numpy()
    if (arr < 0).any():
        fail("negative values in matrix; rarefaction needs non-negative counts")
    if not np.allclose(arr, np.round(arr), atol=1e-6):
        fail("rarefaction needs integer counts; fractional values (e.g. Bracken "
             "new_est_reads or relative abundances) cannot be subsampled — round "
             "estimates first or use mg_diversity observed richness instead")
    arr = np.round(arr).astype(np.int64)
    totals = arr.sum(axis=0)
    empty = matrix.columns[totals == 0].tolist()
    if empty:
        fail(f"samples with zero total counts: {empty}")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    max_depth = int(totals.max())
    min_depth = int(totals.min())
    grid = np.unique(np.maximum(1, np.round(np.linspace(max_depth / args.points, max_depth, args.points)).astype(int)))

    rng = np.random.default_rng(args.seed)
    records = []
    curves = {}
    for j, sample in enumerate(matrix.columns):
        counts = arr[:, j]
        total = int(totals[j])
        depths = np.unique(np.append(grid[grid <= total], total))
        means, sds = [], []
        for d in depths:
            draws = [subsample_richness(counts, int(d), rng) for _ in range(args.iterations)]
            means.append(float(np.mean(draws)))
            sds.append(float(np.std(draws)))
        curves[sample] = (depths, np.array(means), np.array(sds))
        for d, m, s in zip(depths, means, sds):
            records.append({"depth": int(d), "sample": sample,
                            "richness_mean": m, "richness_sd": s})

    table = pd.DataFrame(records)
    table_path = outdir / "rarefaction.tsv"
    table.to_csv(table_path, sep="\t", index=False)

    not_plateaued = []
    for sample, (depths, means, _) in curves.items():
        if len(means) >= 2 and means[-2] > 0 and (means[-1] - means[-2]) / means[-2] > 0.05:
            not_plateaued.append(sample)

    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(8, 5))
    groups = None
    if args.metadata:
        meta = read_metadata(args.metadata, matrix.columns.tolist())
        if args.group_col not in meta.columns:
            fail(f"metadata has no column {args.group_col}")
        groups = meta[args.group_col].astype(str)
    cmap = plt.get_cmap("tab10")
    if groups is not None:
        palette = {g: cmap(i % 10) for i, g in enumerate(sorted(groups.unique()))}
        for sample, (depths, means, sds) in curves.items():
            color = palette[groups[sample]]
            ax.plot(depths, means, color=color, lw=1.2)
            ax.fill_between(depths, means - sds, means + sds, color=color, alpha=0.15)
        handles = [plt.Line2D([], [], color=palette[g], label=g) for g in sorted(palette)]
        ax.legend(handles=handles, title=args.group_col, fontsize=8)
    else:
        for i, (sample, (depths, means, sds)) in enumerate(curves.items()):
            color = cmap(i % 10)
            ax.plot(depths, means, color=color, lw=1.2, label=sample)
            ax.fill_between(depths, means - sds, means + sds, color=color, alpha=0.15)
        ax.legend(fontsize=7)
    ax.axvline(min_depth, color="grey", ls="--", lw=1)
    ax.text(min_depth, ax.get_ylim()[1], f" min library size ({min_depth})",
            fontsize=7, color="grey", va="top")
    ax.set_xlabel("sequencing depth (reads)")
    ax.set_ylabel("observed richness (taxa)")
    ax.set_title(f"Rarefaction ({args.iterations} iterations, seed {args.seed})")
    fig.tight_layout()
    fig_path = outdir / "rarefaction-curves.png"
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)

    emit({
        "metric": "observed",
        "iterations": args.iterations,
        "seed": args.seed,
        "samples": matrix.columns.tolist(),
        "library_sizes": {s: int(t) for s, t in zip(matrix.columns, totals)},
        "min_library_size": min_depth,
        "max_depth": max_depth,
        "not_plateaued": not_plateaued,
        "output": str(table_path),
        "artifacts": [str(table_path), str(fig_path)],
    })


if __name__ == "__main__":
    main()