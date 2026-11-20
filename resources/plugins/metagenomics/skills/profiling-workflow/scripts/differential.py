#!/usr/bin/env python3
"""mg_differential: non-parametric differential abundance (Wilcoxon / Kruskal-Wallis)."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from _common import bh_fdr, emit, fail, read_matrix, read_metadata, setup_matplotlib


def wilcoxon_z(a, b):
    """Two-sided Mann-Whitney U with tie-corrected normal approximation."""
    n1, n2 = len(a), len(b)
    combined = np.concatenate([a, b])
    ranks = stats.rankdata(combined)
    r1 = ranks[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    _, counts = np.unique(combined, return_counts=True)
    n = n1 + n2
    tie_term = ((counts ** 3 - counts).sum()) / (n * (n - 1)) if n > 1 else 0.0
    sigma = np.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie_term))
    if sigma == 0:
        return 0.0, 1.0
    z = (u1 - mu) / sigma
    p = 2.0 * stats.norm.sf(abs(z))
    return float(z), float(min(p, 1.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--group_col", default="group")
    ap.add_argument("--case")
    ap.add_argument("--control")
    ap.add_argument("--values", choices=["counts", "relative"], default="counts")
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--fc_threshold", type=float, default=1.0)
    args = ap.parse_args()

    matrix = read_matrix(args.input)
    samples = matrix.columns.tolist()
    meta = read_metadata(args.metadata, samples)
    if args.group_col not in meta.columns:
        fail(f"metadata has no column {args.group_col}")
    labels = meta[args.group_col].astype(str)
    levels = sorted(labels.unique())
    if args.case or args.control:
        if not (args.case and args.control):
            fail("case and control must be given together")
        if args.case not in levels or args.control not in levels:
            fail(f"case/control not in group levels: {levels}")
        levels = [args.control, args.case]
        labels = labels[labels.isin(levels)]
        matrix = matrix[labels.index]
    if len(levels) < 2:
        fail(f"need at least 2 groups, got {levels}")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    pseudo = 0.5 if args.values == "counts" else 1e-6
    two_group = len(levels) == 2
    rows = []
    for taxon, row in matrix.iterrows():
        vals = row.to_numpy()
        groups = [vals[(labels == g).to_numpy()] for g in levels]
        means = [float(g.mean()) for g in groups]
        if all(m == 0 for m in means):
            continue
        if two_group:
            z, p = wilcoxon_z(groups[1], groups[0])
            log2fc = float(np.log2((means[1] + pseudo) / (means[0] + pseudo)))
            effect = abs(z) / np.sqrt(len(vals))
        else:
            stat, p = stats.kruskal(*groups)
            z, log2fc, effect = float(stat), None, None
        rows.append({
            "taxon": taxon,
            **{f"mean_{g}": m for g, m in zip(levels, means)},
            "log2fc": log2fc,
            "statistic": z,
            "effect_r": effect,
            "p_value": float(p),
        })
    if not rows:
        fail("no taxa with non-zero abundance")
    df = pd.DataFrame(rows).set_index("taxon")
    df["fdr"] = bh_fdr(df["p_value"].to_numpy())
    if two_group:
        df["significant"] = (df["fdr"] < args.fdr) & (df["log2fc"].abs() >= args.fc_threshold)
    else:
        df["significant"] = df["fdr"] < args.fdr
    df = df.sort_values("p_value")
    table_path = outdir / "differential-abundance.tsv"
    df.to_csv(table_path, sep="\t")

    plt = setup_matplotlib()
    artifacts = [str(table_path)]
    if two_group:
        fig, ax = plt.subplots(figsize=(6, 5))
        sig = df["significant"].to_numpy()
        x = df["log2fc"].to_numpy()
        y = -np.log10(df["fdr"].clip(lower=1e-300).to_numpy())
        ax.scatter(x[~sig], y[~sig], s=12, color="grey", alpha=0.6, label="ns")
        ax.scatter(x[sig], y[sig], s=18, color="crimson", alpha=0.8, label="significant")
        ax.axvline(args.fc_threshold, ls="--", lw=0.8, color="black")
        ax.axvline(-args.fc_threshold, ls="--", lw=0.8, color="black")
        ax.axhline(-np.log10(args.fdr), ls="--", lw=0.8, color="black")
        ax.set_xlabel(f"log2 fold change ({levels[1]} / {levels[0]})")
        ax.set_ylabel("-log10 FDR")
        ax.legend()
        fig.tight_layout()
        p = outdir / "volcano.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        artifacts.append(str(p))

    emit({
        "test": "mannwhitneyu" if two_group else "kruskal",
        "groups": levels,
        "values": args.values,
        "n_taxa": int(len(df)),
        "n_significant": int(df["significant"].sum()),
        "fdr_threshold": args.fdr,
        "fc_threshold": args.fc_threshold if two_group else None,
        "top_taxa": df.head(10).reset_index().to_dict(orient="records"),
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()