#!/usr/bin/env python3
"""gen_summary: cohort MAF summary, 6-class substitution spectrum & Ti/Tv, and waterfall oncoplot."""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    NON_SYNONYMOUS_CLASSES,
    SIX_CLASSES,
    TRANSITIONS,
    emit,
    fail,
    pyrimidine_sub,
    read_maf,
)

VC_COLORS = {
    "Missense_Mutation": "#22C55E",
    "Nonsense_Mutation": "#EF4444",
    "Frame_Shift_Del": "#3B82F6",
    "Frame_Shift_Ins": "#8B5CF6",
    "Splice_Site": "#F97316",
    "In_Frame_Del": "#06B6D4",
    "In_Frame_Ins": "#EC4899",
    "Translation_Start_Site": "#EAB308",
    "Nonstop_Mutation": "#14B8A6",
    "Multi_Hit": "#1E293B",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--top_genes", type=int, default=20)
    parser.add_argument("--include_silent", type=lambda x: str(x).lower() in ("true", "1", "yes"), default=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_maf(args.maf)
    all_samples = sorted(df["Tumor_Sample_Barcode"].unique())
    n_samples = len(all_samples)
    if n_samples == 0:
        fail("MAF contains 0 samples")

    nonsyn = df[df["Variant_Classification"].isin(NON_SYNONYMOUS_CLASSES)].copy()
    if len(nonsyn) == 0:
        fail("MAF contains 0 non-synonymous variants")

    # Per-gene summary
    gene_rows = []
    for gene, grp in nonsyn.groupby("Hugo_Symbol"):
        mut_samples = grp["Tumor_Sample_Barcode"].nunique()
        gene_rows.append(
            {
                "Hugo_Symbol": str(gene),
                "n_variants": int(len(grp)),
                "mutated_samples": int(mut_samples),
                "altered_fraction": round(mut_samples / n_samples, 4),
            }
        )
    gene_df = pd.DataFrame(gene_rows).sort_values(["mutated_samples", "n_variants"], ascending=False)
    gene_path = os.path.join(args.output_dir, "gene_summary.tsv")
    gene_df.to_csv(gene_path, sep="\t", index=False)

    # Per-sample summary
    samp_counts = nonsyn.groupby("Tumor_Sample_Barcode").size().reindex(all_samples, fill_value=0)
    samp_df = pd.DataFrame({"Tumor_Sample_Barcode": all_samples, "n_nonsyn_variants": samp_counts.to_numpy(dtype=int)})
    samp_path = os.path.join(args.output_dir, "sample_summary.tsv")
    samp_df.to_csv(samp_path, sep="\t", index=False)

    # 6-class substitution spectrum & Ti/Tv
    snv_src = df if args.include_silent else nonsyn
    snvs = snv_src[snv_src["Variant_Type"] == "SNP"].copy()
    subs = [
        pyrimidine_sub(r, a)
        for r, a in zip(snvs["Reference_Allele"], snvs["Tumor_Seq_Allele2"])
    ]
    snvs["sub6"] = subs
    snvs = snvs.dropna(subset=["sub6"])
    sub_counts = {c: int((snvs["sub6"] == c).sum()) for c in SIX_CLASSES}
    n_ti = sum(sub_counts[c] for c in TRANSITIONS)
    n_tv = sum(sub_counts[c] for c in SIX_CLASSES if c not in TRANSITIONS)
    titv_ratio = round(n_ti / n_tv, 4) if n_tv > 0 else None

    titv_df = pd.DataFrame([{"substitution": c, "count": sub_counts[c], "class": "Ti" if c in TRANSITIONS else "Tv"} for c in SIX_CLASSES])
    titv_path = os.path.join(args.output_dir, "titv_summary.tsv")
    titv_df.to_csv(titv_path, sep="\t", index=False)

    # Figure 1: MAF Summary Dashboard (Variant Classification + 6-class Ti/Tv + Top Genes)
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.2))
    vc_vc = nonsyn["Variant_Classification"].value_counts()
    axes[0].barh(vc_vc.index[::-1], vc_vc.values[::-1], color=[VC_COLORS.get(k, "#64748B") for k in vc_vc.index[::-1]])
    axes[0].set_title("Variant Classification")
    axes[0].set_xlabel("Count")

    sub_colors = ["#38BDF8" if c in TRANSITIONS else "#F59E0B" for c in SIX_CLASSES]
    axes[1].bar(SIX_CLASSES, [sub_counts[c] for c in SIX_CLASSES], color=sub_colors)
    axes[1].set_title(f"SNV Spectrum (Ti/Tv = {titv_ratio if titv_ratio is not None else 'NA'})")
    axes[1].set_ylabel("Count")

    top10 = gene_df.head(10).iloc[::-1]
    axes[2].barh(top10["Hugo_Symbol"], top10["altered_fraction"] * 100.0, color="#818CF8")
    axes[2].set_title("Top Mutated Genes")
    axes[2].set_xlabel("Mutated samples (%)")
    fig.tight_layout()
    summary_fig = os.path.join(args.output_dir, "maf-summary.png")
    fig.savefig(summary_fig, dpi=150)
    plt.close(fig)

    # Figure 2: Waterfall Oncoplot for top_genes
    top_g = list(gene_df.head(args.top_genes)["Hugo_Symbol"])
    cell_mat = pd.DataFrame("", index=top_g, columns=all_samples)
    sub_ns = nonsyn[nonsyn["Hugo_Symbol"].isin(top_g)]
    for (g, s), grp in sub_ns.groupby(["Hugo_Symbol", "Tumor_Sample_Barcode"]):
        classes = grp["Variant_Classification"].unique()
        cell_mat.loc[g, s] = "Multi_Hit" if len(grp) > 1 else str(classes[0])

    # Sort samples in classic waterfall order (binary presence across top_g)
    def sample_sort_key(s: str) -> tuple:
        return tuple(0 if cell_mat.loc[g, s] != "" else 1 for g in top_g) + (-int(samp_counts.loc[s]), s)

    ordered_samples = sorted(all_samples, key=sample_sort_key)
    cell_mat = cell_mat[ordered_samples]

    fig, ax = plt.subplots(figsize=(max(7.0, min(14.0, len(ordered_samples) * 0.28 + 3.0)), max(4.2, len(top_g) * 0.32)))
    ax.set_facecolor("#F1F5F9")
    present_types = set()
    for r_idx, g in enumerate(top_g):
        for c_idx, s in enumerate(ordered_samples):
            vc = cell_mat.loc[g, s]
            if vc:
                present_types.add(vc)
                rect = plt.Rectangle((c_idx, len(top_g) - 1 - r_idx), 0.88, 0.82, color=VC_COLORS.get(vc, "#475569"))
                ax.add_patch(rect)
    ax.set_xlim(0, len(ordered_samples))
    ax.set_ylim(0, len(top_g))
    ax.set_yticks(np.arange(len(top_g)) + 0.41)
    pct_map = dict(zip(gene_df["Hugo_Symbol"], gene_df["altered_fraction"] * 100.0))
    ax.set_yticklabels([f"{g} ({pct_map[g]:.0f}%)" for g in top_g[::-1]], fontsize=8)
    ax.set_xticks([])
    ax.set_xlabel(f"Samples (n={n_samples}, waterfall sorted)")
    ax.set_title(f"Somatic Oncoplot (Top {len(top_g)} genes)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=VC_COLORS.get(t, "#475569"), label=t) for t in sorted(present_types)]
    if handles:
        ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False, fontsize=7.5)
    fig.tight_layout()
    onco_fig = os.path.join(args.output_dir, "oncoplot.png")
    fig.savefig(onco_fig, dpi=150)
    plt.close(fig)

    emit(
        {
            "maf": os.path.abspath(args.maf),
            "gene_summary_table": os.path.abspath(gene_path),
            "sample_summary_table": os.path.abspath(samp_path),
            "titv_table": os.path.abspath(titv_path),
            "n_samples": n_samples,
            "n_nonsyn_variants": int(len(nonsyn)),
            "median_variants_per_sample": round(float(samp_counts.median()), 2),
            "titv_ratio": titv_ratio,
            "substitution_counts": sub_counts,
            "top_genes": gene_df.head(10).to_dict(orient="records"),
            "figures": [os.path.abspath(summary_fig), os.path.abspath(onco_fig)],
        }
    )


if __name__ == "__main__":
    main()
