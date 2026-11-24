#!/usr/bin/env python3
"""prot_ptm: PTM/phosphosite protein-normalized DE and Kinase-Substrate Enrichment Analysis (KSEA)."""

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
    bh_fdr,
    clean_sample_name,
    detect_format_and_columns,
    emit,
    fail,
    load_metadata,
    moderated_t_test,
    normalize_gene_id,
    read_table,
)


def infer_parent_protein(site_id: str) -> str:
    """Extract parent gene/protein symbol from site IDs like 'AKT1_S473' or 'AKT1:S473'."""
    text = str(site_id).strip()
    for sep in ("_", ":", "-", "|"):
        if sep in text:
            return text.split(sep)[0].strip()
    return text


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ptm_matrix", required=True)
    parser.add_argument("--protein_matrix", default=None)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--site_col", default=None)
    parser.add_argument("--protein_col", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default="condition")
    parser.add_argument("--case", required=True)
    parser.add_argument("--control", required=True)
    parser.add_argument("--kinase_table", default=None)
    parser.add_argument("--min_substrates", type=int, default=3)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--log2fc", type=float, default=1.0)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    ptm_df = read_table(args.ptm_matrix)
    lower_cols = {c.lower(): c for c in ptm_df.columns}
    inferred_site_col = args.site_col or lower_cols.get("site_id") or lower_cols.get("ptm_id") or ptm_df.columns[0]
    det = detect_format_and_columns(ptm_df, inferred_site_col)
    site_col = det["id_col"]
    sample_cols = det["sample_cols"]

    mat = ptm_df[sample_cols].apply(pd.to_numeric, errors="coerce")
    mat.columns = [clean_sample_name(c) for c in sample_cols]
    site_ids = ptm_df[site_col].astype(str).str.strip().to_numpy()

    if args.protein_col and args.protein_col in ptm_df.columns:
        parent_proteins = ptm_df[args.protein_col].astype(str).map(normalize_gene_id).to_numpy()
    elif "protein_id" in ptm_df.columns and site_col != "protein_id":
        parent_proteins = ptm_df["protein_id"].astype(str).map(normalize_gene_id).to_numpy()
    elif "gene" in ptm_df.columns and site_col != "gene":
        parent_proteins = ptm_df["gene"].astype(str).map(normalize_gene_id).to_numpy()
    else:
        parent_proteins = np.array([infer_parent_protein(s) for s in site_ids])

    warnings: list[str] = []
    protein_normalized = False
    n_sites_normalized = 0
    if args.protein_matrix:
        prot_df = read_table(args.protein_matrix)
        pdet = detect_format_and_columns(prot_df, None)
        pmat = prot_df[pdet["sample_cols"]].apply(pd.to_numeric, errors="coerce")
        pmat.columns = [clean_sample_name(c) for c in pdet["sample_cols"]]
        pmat.index = prot_df[pdet["id_col"]].astype(str).map(normalize_gene_id)
        pmat = pmat[~pmat.index.duplicated(keep="first")]
        shared_samples = [c for c in mat.columns if c in pmat.columns]
        if len(shared_samples) >= 2:
            norm_vals = mat[shared_samples].to_numpy(dtype=float).copy()
            for idx, prot in enumerate(parent_proteins):
                if prot in pmat.index:
                    norm_vals[idx, :] = norm_vals[idx, :] - pmat.loc[prot, shared_samples].to_numpy(dtype=float)
                    n_sites_normalized += 1
            mat = pd.DataFrame(norm_vals, columns=shared_samples)
            protein_normalized = True
        else:
            warnings.append("protein_matrix had < 2 overlapping samples with ptm_matrix; skipped protein normalization")
    else:
        warnings.append(
            "no protein_matrix provided: PTM fold-changes are NOT normalized for total protein abundance"
        )

    meta, scol = load_metadata(args.metadata, args.sample_col)
    if args.group_col not in meta.columns:
        fail(f"group_col '{args.group_col}' not found in metadata")
    meta_map = dict(zip(meta[scol], meta[args.group_col].astype(str)))

    case_cols = [c for c in mat.columns if meta_map.get(c) == args.case]
    ctrl_cols = [c for c in mat.columns if meta_map.get(c) == args.control]
    if len(case_cols) < 2 or len(ctrl_cols) < 2:
        fail(
            f"need >= 2 samples per group; found {len(case_cols)} for '{args.case}' and {len(ctrl_cols)} for '{args.control}'"
        )

    res = moderated_t_test(mat[case_cols].to_numpy(dtype=float), mat[ctrl_cols].to_numpy(dtype=float))
    direction = np.where(
        (res["padj"] <= args.fdr) & (res["log2fc"] >= args.log2fc),
        "up",
        np.where((res["padj"] <= args.fdr) & (res["log2fc"] <= -args.log2fc), "down", "ns"),
    )

    ptm_out = pd.DataFrame(
        {
            "site_id": site_ids,
            "parent_protein": parent_proteins,
            "protein_normalized": protein_normalized,
            "log2fc": np.round(res["log2fc"], 5),
            "stat": np.round(res["stat"], 5),
            "pvalue": res["pvalue"],
            "padj": res["padj"],
            "regulation": direction,
        }
    ).sort_values(["padj", "pvalue"], na_position="last")

    ptm_de_path = os.path.join(args.output_dir, "ptm_de_results.tsv")
    ptm_out.to_csv(ptm_de_path, sep="\t", index=False)

    figures: list[str] = []

    # Volcano plot
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    nlp = -np.log10(np.clip(ptm_out["pvalue"].to_numpy(dtype=float), 1e-300, 1.0))
    fc = ptm_out["log2fc"].to_numpy(dtype=float)
    reg = ptm_out["regulation"].to_numpy()
    ax.scatter(fc[reg == "ns"], nlp[reg == "ns"], s=14, color="#94A3B8", alpha=0.55)
    ax.scatter(fc[reg == "up"], nlp[reg == "up"], s=20, color="#F43F5E", alpha=0.85, label=f"up ({int((reg == 'up').sum())})")
    ax.scatter(fc[reg == "down"], nlp[reg == "down"], s=20, color="#2563EB", alpha=0.85, label=f"down ({int((reg == 'down').sum())})")
    ax.set_xlabel(f"PTM log2FC ({args.case} vs {args.control})")
    ax.set_ylabel("-log10(p-value)")
    title_suffix = "protein-normalized" if protein_normalized else "unnormalized"
    ax.set_title(f"PTM site differential analysis ({title_suffix})")
    ax.legend(frameon=False)
    fig.tight_layout()
    volc_path = os.path.join(args.output_dir, "ptm-volcano.png")
    fig.savefig(volc_path, dpi=150)
    plt.close(fig)
    figures.append(volc_path)

    # KSEA (Kinase-Substrate Enrichment Analysis)
    ksea_summary = {"skipped": "no kinase_table provided"}
    ksea_path = None
    if args.kinase_table:
        kt = read_table(args.kinase_table)
        kcols = {c.lower(): c for c in kt.columns}
        k_col = kcols.get("kinase") or kt.columns[0]
        s_col = kcols.get("site_id") or kcols.get("substrate") or kt.columns[1]
        fc_map = dict(zip(ptm_out["site_id"].astype(str), ptm_out["log2fc"].to_numpy(dtype=float)))
        all_fc = ptm_out["log2fc"].dropna().to_numpy(dtype=float)
        mean_all = float(np.mean(all_fc)) if len(all_fc) else 0.0
        sd_all = float(np.std(all_fc, ddof=1)) if len(all_fc) > 1 else 1.0
        sd_all = max(sd_all, 1e-6)

        ksea_rows = []
        for kinase, grp in kt.groupby(k_col):
            subs = list(dict.fromkeys(grp[s_col].astype(str).str.strip()))
            matched = [s for s in subs if s in fc_map and np.isfinite(fc_map[s])]
            m = len(matched)
            if m < args.min_substrates:
                continue
            sub_fcs = np.array([fc_map[s] for s in matched], dtype=float)
            mean_sub = float(np.mean(sub_fcs))
            zscore = (mean_sub - mean_all) * np.sqrt(m) / sd_all
            pval = float(2.0 * stats.norm.sf(abs(zscore)))
            ksea_rows.append(
                {
                    "kinase": str(kinase),
                    "n_substrates": m,
                    "mean_log2fc": round(mean_sub, 4),
                    "zscore": round(zscore, 4),
                    "pvalue": pval,
                    "substrates": ";".join(matched[:15]),
                }
            )
        if ksea_rows:
            ksea_df = pd.DataFrame(ksea_rows)
            ksea_df["padj"] = bh_fdr(ksea_df["pvalue"].to_numpy(dtype=float))
            ksea_df = ksea_df.sort_values(["padj", "pvalue"])
            ksea_path = os.path.join(args.output_dir, "ksea_results.tsv")
            ksea_df.to_csv(ksea_path, sep="\t", index=False)

            top_k = ksea_df.head(15).iloc[::-1]
            fig, ax = plt.subplots(figsize=(6.0, max(3.5, len(top_k) * 0.35)))
            colors = ["#EF4444" if z > 0 else "#3B82F6" for z in top_k["zscore"]]
            ax.barh(top_k["kinase"], top_k["zscore"], color=colors)
            ax.axvline(0.0, color="#334155", linewidth=0.9)
            ax.set_xlabel("KSEA z-score")
            ax.set_title("Kinase-Substrate Enrichment Analysis (KSEA)")
            fig.tight_layout()
            ksea_fig = os.path.join(args.output_dir, "ksea-barplot.png")
            fig.savefig(ksea_fig, dpi=150)
            plt.close(fig)
            figures.append(ksea_fig)

            ksea_summary = {
                "ksea_table": os.path.abspath(ksea_path),
                "n_kinases_tested": int(len(ksea_df)),
                "n_significant": int((ksea_df["padj"] <= args.fdr).sum()),
                "top_kinases": [
                    {
                        "kinase": str(r["kinase"]),
                        "n_substrates": int(r["n_substrates"]),
                        "zscore": round(float(r["zscore"]), 4),
                        "padj": float(r["padj"]),
                    }
                    for _, r in ksea_df.head(10).iterrows()
                ],
            }
        else:
            ksea_summary = {"skipped": f"no kinase had >= {args.min_substrates} quantified substrates"}

    emit(
        {
            "ptm_matrix": os.path.abspath(args.ptm_matrix),
            "ptm_de_table": os.path.abspath(ptm_de_path),
            "protein_normalized": protein_normalized,
            "n_sites_protein_normalized": n_sites_normalized,
            "contrast": f"{args.case}_vs_{args.control}",
            "n_sites_tested": int(ptm_out["pvalue"].notna().sum()),
            "n_up": int((ptm_out["regulation"] == "up").sum()),
            "n_down": int((ptm_out["regulation"] == "down").sum()),
            "ksea": ksea_summary,
            "figures": [os.path.abspath(f) for f in figures],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
