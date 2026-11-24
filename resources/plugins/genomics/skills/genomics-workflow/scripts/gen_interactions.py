#!/usr/bin/env python3
"""gen_interactions: somatic co-occurrence/exclusivity, 10 oncogenic pathways & cohort comparison."""

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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    NON_SYNONYMOUS_CLASSES,
    SANCHEZ_VEGA_PATHWAYS,
    bh_fdr,
    emit,
    fail,
    load_clinical,
    read_maf,
)


def binary_gene_sample_matrix(df: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    nonsyn = df[df["Variant_Classification"].isin(NON_SYNONYMOUS_CLASSES)]
    ct = pd.crosstab(nonsyn["Hugo_Symbol"], nonsyn["Tumor_Sample_Barcode"])
    ct = (ct > 0).astype(int)
    return ct.reindex(columns=samples, fill_value=0)


def compare_two_cohorts(
    mat1: pd.DataFrame, mat2: pd.DataFrame, name1: str, name2: str, min_mut: int
) -> pd.DataFrame:
    n1 = mat1.shape[1]
    n2 = mat2.shape[1]
    genes = sorted(set(mat1.index) | set(mat2.index))
    rows = []
    for g in genes:
        k1 = int(mat1.loc[g].sum()) if g in mat1.index else 0
        k2 = int(mat2.loc[g].sum()) if g in mat2.index else 0
        if max(k1, k2) < min_mut:
            continue
        table = np.array([[k1, n1 - k1], [k2, n2 - k2]], dtype=int)
        odds_r, pval = stats.fisher_exact(table, alternative="two-sided")
        rows.append(
            {
                "Hugo_Symbol": g,
                f"{name1}_mutated": k1,
                f"{name1}_total": n1,
                f"{name2}_mutated": k2,
                f"{name2}_total": n2,
                "odds_ratio": round(float(odds_r), 4) if np.isfinite(odds_r) else float("inf"),
                "pvalue": float(pval),
            }
        )
    if not rows:
        return pd.DataFrame(columns=["Hugo_Symbol", f"{name1}_mutated", f"{name2}_mutated", "odds_ratio", "pvalue", "padj"])
    res = pd.DataFrame(rows)
    res["padj"] = bh_fdr(res["pvalue"].to_numpy(dtype=float))
    return res.sort_values(["padj", "pvalue"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--compare_maf", default=None)
    parser.add_argument("--clinical", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default=None)
    parser.add_argument("--case", default=None)
    parser.add_argument("--control", default=None)
    parser.add_argument("--top_genes", type=int, default=20)
    parser.add_argument("--min_mut", type=int, default=3)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_maf(args.maf)
    all_samples = sorted(df["Tumor_Sample_Barcode"].unique())
    n_samples = len(all_samples)
    if n_samples < 2:
        fail("need >= 2 samples for interaction/pathway analysis")

    mat = binary_gene_sample_matrix(df, all_samples)
    gene_order = list(mat.sum(axis=1).sort_values(ascending=False).head(args.top_genes).index)

    # 1. Pairwise somatic interactions (co-occurrence / mutual exclusivity)
    pair_rows = []
    for i in range(len(gene_order)):
        g1 = gene_order[i]
        v1 = mat.loc[g1].to_numpy(dtype=int)
        for j in range(i + 1, len(gene_order)):
            g2 = gene_order[j]
            v2 = mat.loc[g2].to_numpy(dtype=int)
            n11 = int(np.sum((v1 == 1) & (v2 == 1)))
            n10 = int(np.sum((v1 == 1) & (v2 == 0)))
            n01 = int(np.sum((v1 == 0) & (v2 == 1)))
            n00 = int(np.sum((v1 == 0) & (v2 == 0)))
            odds_r, pval = stats.fisher_exact([[n11, n10], [n01, n00]], alternative="two-sided")
            event = "Co_Occurrence" if odds_r > 1.0 else "Mutually_Exclusive"
            pair_rows.append(
                {
                    "gene1": g1,
                    "gene2": g2,
                    "n11": n11,
                    "n10": n10,
                    "n01": n01,
                    "n00": n00,
                    "odds_ratio": round(float(odds_r), 4) if np.isfinite(odds_r) else float("inf"),
                    "event": event,
                    "pvalue": float(pval),
                }
            )
    inter_df = pd.DataFrame(pair_rows)
    if len(inter_df):
        inter_df["padj"] = bh_fdr(inter_df["pvalue"].to_numpy(dtype=float))
        inter_df = inter_df.sort_values(["padj", "pvalue"])
    else:
        inter_df = pd.DataFrame(columns=["gene1", "gene2", "odds_ratio", "event", "pvalue", "padj"])
    inter_path = os.path.join(args.output_dir, "somatic_interactions.tsv")
    inter_df.to_csv(inter_path, sep="\t", index=False)

    # 2. Sanchez-Vega (Cell 2018) 10 canonical oncogenic signaling pathways
    pw_rows = []
    for pw_name, members in SANCHEZ_VEGA_PATHWAYS.items():
        present = [g for g in members if g in mat.index]
        if present:
            altered_mask = (mat.loc[present].sum(axis=0) > 0).to_numpy()
            n_alt_samples = int(altered_mask.sum())
        else:
            n_alt_samples = 0
        pw_rows.append(
            {
                "pathway": pw_name,
                "n_pathway_genes": len(members),
                "n_affected_genes": len(present),
                "n_altered_samples": n_alt_samples,
                "altered_sample_fraction": round(n_alt_samples / n_samples, 4),
                "affected_genes": ";".join(present[:15]),
            }
        )
    pw_df = pd.DataFrame(pw_rows).sort_values("altered_sample_fraction", ascending=False)
    pw_path = os.path.join(args.output_dir, "oncogenic_pathways.tsv")
    pw_df.to_csv(pw_path, sep="\t", index=False)

    figures: list[str] = []

    # Plot oncogenic pathways barplot
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    pw_plot = pw_df.iloc[::-1]
    ax.barh(pw_plot["pathway"], pw_plot["altered_sample_fraction"] * 100.0, color="#F43F5E", alpha=0.85)
    ax.set_xlabel("Altered samples (%)")
    ax.set_title("Sanchez-Vega 10 Canonical Oncogenic Pathways")
    fig.tight_layout()
    pw_fig = os.path.join(args.output_dir, "oncogenic-pathways.png")
    fig.savefig(pw_fig, dpi=150)
    plt.close(fig)
    figures.append(pw_fig)

    # Plot pairwise interaction matrix (-log10 p signed by Co_Occurrence vs Mutually_Exclusive)
    if len(gene_order) >= 2:
        k = len(gene_order)
        grid = np.zeros((k, k), dtype=float)
        idx_map = {g: i for i, g in enumerate(gene_order)}
        for _, r in inter_df.iterrows():
            i, j = idx_map[r["gene1"]], idx_map[r["gene2"]]
            score = min(5.0, -np.log10(max(float(r["pvalue"]), 1e-10)))
            signed = score if r["event"] == "Co_Occurrence" else -score
            grid[i, j] = signed
            grid[j, i] = signed
        fig, ax = plt.subplots(figsize=(max(5.5, k * 0.42), max(5.0, k * 0.42)))
        im = ax.imshow(grid, cmap="BrBG", vmin=-3.5, vmax=3.5)
        ax.set_xticks(range(k))
        ax.set_yticks(range(k))
        ax.set_xticklabels(gene_order, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(gene_order, fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Signed -log10(p) (+Co-occ / -Excl)")
        ax.set_title("Somatic Interactions (Fisher's exact)")
        fig.tight_layout()
        inter_fig = os.path.join(args.output_dir, "somatic-interactions.png")
        fig.savefig(inter_fig, dpi=150)
        plt.close(fig)
        figures.append(inter_fig)

    # 3. Optional cohort or clinical-group differential mutation comparison
    cohort_cmp = {"skipped": "no compare_maf or clinical case/control contrast provided"}
    if args.compare_maf:
        df2 = read_maf(args.compare_maf)
        s2 = sorted(df2["Tumor_Sample_Barcode"].unique())
        mat2 = binary_gene_sample_matrix(df2, s2)
        cmp_df = compare_two_cohorts(mat, mat2, "cohort1", "cohort2", args.min_mut)
        cmp_path = os.path.join(args.output_dir, "cohort_comparison.tsv")
        cmp_df.to_csv(cmp_path, sep="\t", index=False)
        cohort_cmp = {
            "table": os.path.abspath(cmp_path),
            "mode": "two_maf",
            "n_genes_tested": int(len(cmp_df)),
            "n_significant": int((cmp_df["padj"] <= args.fdr).sum()) if len(cmp_df) else 0,
            "top": cmp_df.head(10).to_dict(orient="records"),
        }
    elif args.clinical and args.group_col and args.case and args.control:
        clin, scol = load_clinical(args.clinical, args.sample_col)
        if args.group_col not in clin.columns:
            fail(f"group_col '{args.group_col}' not found in clinical table")
        gmap = dict(zip(clin[scol], clin[args.group_col].astype(str)))
        s_case = [s for s in all_samples if gmap.get(s) == args.case]
        s_ctrl = [s for s in all_samples if gmap.get(s) == args.control]
        if len(s_case) < 2 or len(s_ctrl) < 2:
            fail(f"need >= 2 samples per clinical group; found {len(s_case)} for {args.case}, {len(s_ctrl)} for {args.control}")
        cmp_df = compare_two_cohorts(mat[s_case], mat[s_ctrl], args.case, args.control, args.min_mut)
        cmp_path = os.path.join(args.output_dir, "cohort_comparison.tsv")
        cmp_df.to_csv(cmp_path, sep="\t", index=False)
        cohort_cmp = {
            "table": os.path.abspath(cmp_path),
            "mode": "clinical_contrast",
            "contrast": f"{args.case}_vs_{args.control}",
            "n_case": len(s_case),
            "n_control": len(s_ctrl),
            "n_genes_tested": int(len(cmp_df)),
            "n_significant": int((cmp_df["padj"] <= args.fdr).sum()) if len(cmp_df) else 0,
            "top": cmp_df.head(10).to_dict(orient="records"),
        }

    emit(
        {
            "maf": os.path.abspath(args.maf),
            "interactions_table": os.path.abspath(inter_path),
            "pathways_table": os.path.abspath(pw_path),
            "n_samples": n_samples,
            "n_pairs_tested": int(len(inter_df)),
            "n_significant_pairs": int((inter_df["padj"] <= args.fdr).sum()) if len(inter_df) else 0,
            "top_interactions": inter_df.head(10).to_dict(orient="records"),
            "pathways": pw_df.to_dict(orient="records"),
            "cohort_comparison": cohort_cmp,
            "figures": [os.path.abspath(f) for f in figures],
        }
    )


if __name__ == "__main__":
    main()
