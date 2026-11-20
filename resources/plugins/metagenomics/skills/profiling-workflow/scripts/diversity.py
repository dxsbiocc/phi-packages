#!/usr/bin/env python3
"""mg_diversity: alpha/beta diversity, PCoA ordination and PERMANOVA."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from _common import emit, fail, read_matrix, read_metadata, setup_matplotlib, write_artifact


def alpha_indices(counts, values):
    x = np.asarray(counts, dtype=float)
    total = x.sum()
    out = {"observed": None, "chao1": None, "shannon": None, "simpson": None, "invsimpson": None}
    if total <= 0:
        return out
    p = x[x > 0] / total
    out["shannon"] = float(-(p * np.log(p)).sum())
    d = float((p * p).sum())
    out["simpson"] = 1.0 - d
    out["invsimpson"] = 1.0 / d if d > 0 else None
    if values == "counts":
        s_obs = int((x > 0).sum())
        f1 = int((x == 1).sum())
        f2 = int((x == 2).sum())
        out["observed"] = s_obs
        if f2 > 0:
            out["chao1"] = s_obs + f1 * f1 / (2.0 * f2)
        else:
            out["chao1"] = s_obs + f1 * (f1 - 1) / 2.0
    return out


def pcoa(dist):
    n = dist.shape[0]
    d2 = dist * dist
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos = eigvals > 1e-10
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    var = eigvals[pos] / eigvals[pos].sum() if pos.any() else np.array([])
    return coords, var


def permanova(dist, labels, permutations, seed):
    groups = np.unique(labels)
    a, n = len(groups), len(labels)
    if a < 2 or a >= n:
        return None
    d2 = dist * dist
    ss_total = d2.sum() / n
    ss_within = 0.0
    for g in groups:
        idx = np.where(labels == g)[0]
        ng = len(idx)
        if ng > 1:
            ss_within += d2[np.ix_(idx, idx)].sum() / ng
    ss_between = ss_total - ss_within
    if ss_within <= 0:
        return None
    f_obs = (ss_between / (a - 1)) / (ss_within / (n - a))
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(permutations):
        perm = rng.permutation(labels)
        ss_w = 0.0
        for g in groups:
            idx = np.where(perm == g)[0]
            ng = len(idx)
            if ng > 1:
                ss_w += d2[np.ix_(idx, idx)].sum() / ng
        ss_b = ss_total - ss_w
        f_perm = (ss_b / (a - 1)) / (ss_w / (n - a)) if ss_w > 0 else np.inf
        if f_perm >= f_obs:
            ge += 1
    return {
        "pseudo_f": float(f_obs),
        "p_value": float((1 + ge) / (1 + permutations)),
        "permutations": int(permutations),
        "groups": groups.tolist(),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--values", choices=["counts", "relative"], default="counts")
    ap.add_argument("--metadata")
    ap.add_argument("--group_col", default="group")
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--permutations", type=int, default=999)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    matrix = read_matrix(args.input)  # taxa x samples
    samples = matrix.columns.tolist()
    x = matrix.to_numpy().T  # samples x taxa
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = read_metadata(args.metadata, samples) if args.metadata else None
    labels = None
    if meta is not None:
        if args.group_col not in meta.columns:
            fail(f"metadata has no column {args.group_col}; columns: {meta.columns.tolist()}")
        labels = meta[args.group_col].astype(str).to_numpy()

    alpha = {s: alpha_indices(row, args.values) for s, row in zip(samples, x)}
    alpha_df = pd.DataFrame(alpha).T
    alpha_df.index.name = "sample"
    alpha_path = outdir / "alpha-diversity.tsv"
    alpha_df.to_csv(alpha_path, sep="\t")

    dist = squareform(pdist(x, metric="braycurtis"))
    dist_df = pd.DataFrame(dist, index=samples, columns=samples)
    dist_path = outdir / "bray-curtis.tsv"
    dist_df.to_csv(dist_path, sep="\t")

    coords, var = pcoa(dist)
    axes = [f"PCo{i + 1}" for i in range(coords.shape[1])]
    pcoa_df = pd.DataFrame(coords, index=samples, columns=axes)
    pcoa_path = outdir / "pcoa.tsv"
    pcoa_df.to_csv(pcoa_path, sep="\t")

    perm = permanova(dist, labels, args.permutations, args.seed) if labels is not None else None

    plt = setup_matplotlib()
    inputs = [args.input] + ([args.metadata] if args.metadata else [])
    params = {"values": args.values, "permutations": args.permutations, "seed": args.seed}
    artifacts = [
        write_artifact(alpha_path, "mg_diversity", "Alpha diversity indices",
                       skill="profiling-workflow", inputs=inputs, parameters=params,
                       shape=alpha_df.shape),
        write_artifact(dist_path, "mg_diversity", "Bray-Curtis distance matrix",
                       skill="profiling-workflow", inputs=inputs, parameters=params,
                       shape=dist_df.shape),
        write_artifact(pcoa_path, "mg_diversity", "PCoA coordinates",
                       skill="profiling-workflow", inputs=inputs, parameters=params,
                       shape=pcoa_df.shape),
    ]

    fig, ax = plt.subplots(figsize=(max(6, len(samples) * 0.5), 4))
    indices = [c for c in alpha_df.columns if alpha_df[c].notna().any()]
    data = [alpha_df[c].dropna().to_numpy() for c in indices]
    ax.boxplot(data, tick_labels=indices)
    ax.set_ylabel("index value")
    ax.set_title(f"Alpha diversity ({args.values})")
    fig.tight_layout()
    p = outdir / "alpha-boxplot.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    artifacts.append(write_artifact(p, "mg_diversity", f"Alpha diversity boxplot ({args.values})",
                                    skill="profiling-workflow", inputs=inputs, parameters=params))

    if coords.shape[1] >= 2:
        fig, ax = plt.subplots(figsize=(6, 5))
        if labels is not None:
            for g in np.unique(labels):
                m = labels == g
                ax.scatter(coords[m, 0], coords[m, 1], label=str(g), s=50, alpha=0.8)
            ax.legend()
        else:
            ax.scatter(coords[:, 0], coords[:, 1], s=50, alpha=0.8)
        for i, s in enumerate(samples):
            ax.annotate(s, (coords[i, 0], coords[i, 1]), fontsize=7)
        ax.set_xlabel(f"PCo1 ({var[0] * 100:.1f}%)")
        ax.set_ylabel(f"PCo2 ({var[1] * 100:.1f}%)")
        ax.set_title("PCoA (Bray-Curtis)")
        fig.tight_layout()
        p = outdir / "pcoa.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        artifacts.append(write_artifact(p, "mg_diversity", "PCoA ordination (Bray-Curtis)",
                                        skill="profiling-workflow", inputs=inputs, parameters=params))

    notes = []
    if args.values == "relative":
        notes.append("observed/chao1 skipped: relative abundances are not counts")
    emit({
        "samples": samples,
        "values": args.values,
        "alpha": alpha,
        "beta_method": "braycurtis",
        "pcoa_variance_explained": var.tolist(),
        "permanova": perm,
        "notes": notes,
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()