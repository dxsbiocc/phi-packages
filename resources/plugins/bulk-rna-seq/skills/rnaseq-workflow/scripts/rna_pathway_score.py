#!/usr/bin/env python3
"""rna_pathway_score: single-sample pathway activity scoring (mean gene z-score
per GMT pathway), optional group differential testing (Mann-Whitney U /
Kruskal-Wallis + BH FDR), and pathway-by-sample activity heatmap.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

from _common import bh_fdr, emit, fail, read_expression_matrix, read_gmt, read_metadata, setup_matplotlib, write_artifact


def to_log_scale(df):
    x = df.to_numpy(dtype=float)
    if (x < 0).any() or float(np.median(x.sum(axis=0))) <= 50.0:
        return x
    col_sums = np.where(x.sum(axis=0) > 0, x.sum(axis=0), 1.0)
    return np.log2((x / col_sums * 1e6) + 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Gene x sample expression matrix TSV/CSV.")
    ap.add_argument("--gmt", required=True, help="GMT gene-set file.")
    ap.add_argument("--metadata", help="Optional sample metadata for group comparison of pathway scores.")
    ap.add_argument("--sample_col", default=None)
    ap.add_argument("--group_col", default="condition")
    ap.add_argument("--min_set_size", type=int, default=3)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    df = read_expression_matrix(args.input)
    samples = df.columns.tolist()
    x_log = to_log_scale(df)
    # Per-gene z-score across samples
    mu = x_log.mean(axis=1, keepdims=True)
    sd = x_log.std(axis=1, keepdims=True)
    valid_var = (sd[:, 0] > 1e-12)
    z_mat = np.zeros_like(x_log)
    z_mat[valid_var] = (x_log[valid_var] - mu[valid_var]) / sd[valid_var]
    z_df = pd.DataFrame(z_mat, index=df.index, columns=samples)

    gene_sets = read_gmt(args.gmt)
    score_rows = {}
    gene_index = set(z_df.index)
    for gs_name, gs_genes in gene_sets.items():
        matched = [g for g in gs_genes if g in gene_index]
        if len(matched) < args.min_set_size:
            continue
        score_rows[gs_name] = z_df.loc[matched].mean(axis=0).to_numpy(dtype=float)

    if not score_rows:
        fail(f"no GMT pathways had at least {args.min_set_size} genes matching the expression matrix")

    score_df = pd.DataFrame.from_dict(score_rows, orient="index", columns=samples)
    score_df.index.name = "pathway"

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = [args.input, args.gmt] + ([args.metadata] if args.metadata else [])
    params = {"method": "mean_zscore", "min_set_size": args.min_set_size}

    score_path = outdir / "pathway-activity-scores.tsv"
    score_df.round(4).to_csv(score_path, sep="\t")
    artifacts = [
        write_artifact(score_path, "rna_pathway_score", "Single-sample pathway activity matrix (mean z-score)", inputs=inputs, parameters=params, shape=score_df.shape)
    ]
    outputs = {"scores_matrix": str(score_path)}

    diff_rows = []
    if args.metadata:
        meta = read_metadata(args.metadata, samples=samples, sample_col=args.sample_col)
        if args.group_col not in meta.columns:
            fail(f"metadata missing group_col {args.group_col!r}")
        groups = meta[args.group_col].astype(str).to_numpy()
        uniq = sorted(np.unique(groups).tolist())
        if len(uniq) >= 2:
            for pw in score_df.index:
                vals = score_df.loc[pw].to_numpy(dtype=float)
                grp_vals = [vals[groups == g] for g in uniq]
                if len(uniq) == 2:
                    stat, pval = mannwhitneyu(grp_vals[1], grp_vals[0], alternative="two-sided")
                    diff_mean = float(grp_vals[1].mean() - grp_vals[0].mean())
                    test_name = "mannwhitneyu"
                else:
                    stat, pval = kruskal(*grp_vals)
                    diff_mean = float(max(g.mean() for g in grp_vals) - min(g.mean() for g in grp_vals))
                    test_name = "kruskal"
                diff_rows.append(
                    {
                        "pathway": pw,
                        "test": test_name,
                        "mean_diff": round(diff_mean, 4),
                        "stat": round(float(stat), 4),
                        "p_value": float(pval),
                    }
                )
            diff_df = pd.DataFrame(diff_rows)
            diff_df["q_value"] = bh_fdr(diff_df["p_value"].to_numpy())
            diff_df = diff_df.sort_values(["q_value", "p_value", "pathway"]).reset_index(drop=True)
            diff_path = outdir / "pathway-group-tests.tsv"
            diff_df.to_csv(diff_path, sep="\t", index=False)
            outputs["group_tests"] = str(diff_path)
            artifacts.append(
                write_artifact(diff_path, "rna_pathway_score", "Group comparison of pathway activity scores", inputs=inputs, parameters=params, shape=diff_df.shape)
            )
            diff_rows = diff_df.to_dict(orient="records")

    # Heatmap of top 30 pathways by variance
    show_df = score_df
    if len(show_df) > 30:
        top_pw = show_df.var(axis=1).sort_values(ascending=False).index[:30]
        show_df = show_df.loc[top_pw]
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(samples) + 2), max(4, 0.35 * len(show_df) + 1.5)))
    vmax = max(1.5, float(np.percentile(np.abs(show_df.to_numpy()), 98)))
    im = ax.imshow(show_df.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(show_df)))
    ax.set_yticklabels([p[:34] for p in show_df.index], fontsize=7.5)
    ax.set_title("Single-sample pathway activity (mean gene z-score)")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Activity z-score")
    fig.tight_layout()
    hm_path = outdir / "pathway-activity-heatmap.png"
    fig.savefig(hm_path, dpi=180)
    plt.close(fig)
    outputs["heatmap_figure"] = str(hm_path)
    artifacts.append(write_artifact(hm_path, "rna_pathway_score", "Pathway activity heatmap", inputs=inputs, parameters=params))

    emit(
        {
            "tool": "rna_pathway_score",
            "n_samples": int(len(samples)),
            "n_pathways_scored": int(len(score_df)),
            "top_group_tests": diff_rows[:10],
            "outputs": outputs,
            "scientificContext": {
                str(score_path): {
                    "valueSemantics": "zscore",
                    "sampleIds": [str(s) for s in samples],
                }
            },
            "artifacts": artifacts,
        }
    )


if __name__ == "__main__":
    main()
