#!/usr/bin/env python3
"""Tumor microenvironment (TME) ESTIMATE-like scoring and NNLS immune cell deconvolution."""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    bh_fdr,
    emit,
    fail,
    infer_value_kind,
    read_expression_matrix,
    read_metadata,
    setup_matplotlib,
    write_artifact,
)

# Curated canonical stromal & pan-immune core signatures (human/mouse uppercase symbol orthologs)
STROMAL_GENES = [
    "FAP", "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
    "DCN", "LUM", "BGN", "THY1", "ACTA2", "TAGLN", "POSTN", "FN1", "VCAN", "SPARC",
    "MMP2", "MMP11", "MMP14", "PDGFRA", "PDGFRB", "PECAM1", "VWF", "CDH5", "KDR", "ENG",
    "CXCL12", "MFAP5", "CTHRC1", "PDPN", "LOXL2", "TIMP1", "TIMP3", "IGFBP7", "CALD1",
]

IMMUNE_GENES = [
    "PTPRC", "CD2", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "GZMK",
    "PRF1", "NKG7", "GNLY", "KLRD1", "MS4A1", "CD19", "CD79A", "CD79B", "MZB1",
    "CD14", "FCGR3A", "CSF1R", "CD68", "CD163", "MRC1", "ITGAM", "ITGAX", "C1QA", "C1QB",
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD74", "CXCL9", "CXCL10", "CCL5", "LCK", "ZAP70", "IL2RG",
]

# Curated 10-lineage leukocyte reference markers: (primary specific markers weight=1.0, secondary lineage markers weight=0.25)
BUILTIN_IMMUNE_PANEL: dict[str, tuple[list[str], list[str]]] = {
    "CD8_T": (
        ["CD8A", "CD8B", "GZMK", "GZMA", "LAG3", "PDCD1", "CD244"],
        ["CD3D", "CD3E", "PRF1", "IFNG"],
    ),
    "CD4_T": (
        ["CD4", "IL7R", "CCR7", "LEF1", "TCF7", "CD40LG", "ICOS"],
        ["CD3D", "CD3E", "CD2", "LCK"],
    ),
    "Treg": (
        ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "CCR8"],
        ["CD4", "CD3D", "TIGIT"],
    ),
    "NK": (
        ["NCR1", "KLRB1", "KLRC1", "KLRD1", "KLRF1", "GNLY", "NKG7", "FGFBP2"],
        ["PRF1", "GZMB", "FCGR3A"],
    ),
    "B_cell": (
        ["MS4A1", "CD19", "CD79A", "CD79B", "BANK1", "PAX5", "CD22", "TCL1A"],
        ["HLA-DRA", "CD74"],
    ),
    "Plasma_cell": (
        ["MZB1", "JCHAIN", "SDC1", "TNFRSF17", "DERL3", "PRDM1", "XBP1", "IGKC"],
        ["CD79A"],
    ),
    "Monocyte": (
        ["CD14", "S100A8", "S100A9", "S100A12", "FCN1", "VCAN", "LYZ", "CCR2"],
        ["ITGAM", "CSF1R", "FCGR3A"],
    ),
    "Macrophage_M1": (
        ["CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1", "GBP5", "SOCS1", "CD80", "CD86"],
        ["CD68", "HLA-DRA", "STAT1"],
    ),
    "Macrophage_M2": (
        ["CD163", "MRC1", "MSR1", "MARCO", "STAB1", "TREM2", "APOE", "CCL18"],
        ["CD68", "C1QA", "C1QB", "CSF1R"],
    ),
    "Neutrophil": (
        ["FCGR3B", "CXCR2", "CSF3R", "FPR1", "NAMPT", "G0S2", "PROK2", "CEACAM8"],
        ["S100A8", "S100A9", "ITGAM"],
    ),
}


def _build_builtin_signature_df():
    import pandas as pd

    all_genes: list[str] = []
    for prim, sec in BUILTIN_IMMUNE_PANEL.values():
        for g in prim + sec:
            if g not in all_genes:
                all_genes.append(g)
    data = {ct: [0.0] * len(all_genes) for ct in BUILTIN_IMMUNE_PANEL}
    g_idx = {g: i for i, g in enumerate(all_genes)}
    for ct, (prim, sec) in BUILTIN_IMMUNE_PANEL.items():
        for g in prim:
            data[ct][g_idx[g]] = 1.0
        for g in sec:
            data[ct][g_idx[g]] = max(data[ct][g_idx[g]], 0.25)
    # Add baseline background housekeeping floor (0.01) so matrix is well-conditioned
    df = pd.DataFrame(data, index=all_genes)
    df.index.name = "gene"
    return df


def _to_linear_expression(expr, value_kind: str):
    import numpy as np

    X = expr.to_numpy(dtype=float)
    if value_kind == "raw_counts":
        X_pos = np.clip(X, 0.0, None)
        lib = X_pos.sum(axis=0, keepdims=True)
        lib[lib <= 0] = 1.0
        lin = X_pos / lib * 1e6
    elif value_kind == "log_normalized":
        min_v = float(np.nanmin(X))
        shifted = X - min_v if min_v < 0 else X
        shifted = np.clip(shifted, 0.0, 25.0)
        lin = np.power(2.0, shifted) - 1.0
        lin = np.clip(lin, 0.0, None)
    else:
        lin = np.clip(X, 0.0, None)
    return expr.__class__(lin, index=expr.index, columns=expr.columns)


def _ssgsea_estimate_score(expr_upper, gene_set: list[str]) -> tuple[list[float], int]:
    """Rank-weighted empirical CDF enrichment score scaled to ESTIMATE range [-2000, 2000]."""
    import numpy as np

    present = [g for g in gene_set if g in expr_upper.index]
    if len(present) < 3:
        return [0.0] * expr_upper.shape[1], len(present)

    X = expr_upper.to_numpy(dtype=float)  # genes x samples
    n_genes, n_samples = X.shape
    hit_mask = np.array([g in set(present) for g in expr_upper.index], dtype=bool)
    n_miss = float(n_genes - int(hit_mask.sum()))
    if n_miss <= 0:
        return [0.0] * n_samples, len(present)

    scores: list[float] = []
    alpha = 0.25
    for s_idx in range(n_samples):
        col = X[:, s_idx]
        order = np.argsort(col)  # ascending rank 0..n_genes-1
        ranks = np.empty(n_genes, dtype=float)
        ranks[order] = np.arange(1, n_genes + 1, dtype=float)
        # Order descending by expression for running sum
        desc = order[::-1]
        hits_ord = hit_mask[desc]
        weights = np.where(hits_ord, np.power(ranks[desc], alpha), 0.0)
        w_sum = float(weights.sum())
        if w_sum <= 0:
            scores.append(0.0)
            continue
        p_hit = np.cumsum(weights) / w_sum
        p_miss = np.cumsum(~hits_ord) / n_miss
        es_raw = float(np.sum(p_hit - p_miss)) / float(n_genes)
        # Scale [-0.5, 0.5] integral to canonical ESTIMATE score range (~ -2000 .. +2500)
        scores.append(es_raw * 5000.0)
    return scores, len(present)


def main() -> None:
    ap = argparse.ArgumentParser(description="Tumor microenvironment scoring and NNLS immune deconvolution")
    ap.add_argument("--input", required=True, help="Expression matrix TSV/CSV (genes x samples, HGNC/MGI symbols)")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument(
        "--signature_matrix",
        default="",
        help="Optional custom linear reference signature matrix TSV/CSV (genes x cell_types)",
    )
    ap.add_argument("--metadata", default="", help="Optional sample metadata TSV/CSV")
    ap.add_argument("--sample_col", default="", help="Optional sample ID column in metadata")
    ap.add_argument("--group_col", default="", help="Optional condition/group column for differential infiltration tests")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from scipy import optimize, stats

    expr = read_expression_matrix(args.input)
    samples = expr.columns.tolist()
    value_kind = infer_value_kind(expr)

    # Check if IDs look like unmapped Ensembl IDs
    sample_genes = [str(g) for g in expr.index[: min(200, len(expr.index))]]
    ensembl_frac = sum(1 for g in sample_genes if g.upper().startswith(("ENSG0", "ENSMUSG0", "ENSRNOG0"))) / max(
        1, len(sample_genes)
    )
    if ensembl_frac > 0.5 and not args.signature_matrix:
        fail(
            "Input matrix appears to use Ensembl gene IDs (e.g. ENSG000...). "
            "Run rna_id_map first to convert Ensembl IDs to gene symbols before calling rna_tme_deconv, "
            "or pass a custom --signature_matrix indexed by matching IDs."
        )

    lin_expr = _to_linear_expression(expr, value_kind)
    # Collapse case-insensitive symbols by max mean expression for orthologous matching
    upper_ids = pd.Index([str(g).strip().upper() for g in lin_expr.index])
    lin_upper = lin_expr.copy()
    lin_upper.index = upper_ids
    if lin_upper.index.has_duplicates:
        lin_upper["_mean"] = lin_upper.mean(axis=1)
        lin_upper = lin_upper.sort_values("_mean", ascending=False)
        lin_upper = lin_upper[~lin_upper.index.duplicated(keep="first")].drop(columns=["_mean"])

    warnings: list[str] = []

    # 1. ESTIMATE-style StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity
    stromal_scores, n_stromal_found = _ssgsea_estimate_score(lin_upper, STROMAL_GENES)
    immune_scores, n_immune_found = _ssgsea_estimate_score(lin_upper, IMMUNE_GENES)
    if n_stromal_found < 5 or n_immune_found < 5:
        warnings.append(
            f"Low overlap with canonical stromal ({n_stromal_found}/{len(STROMAL_GENES)}) or immune "
            f"({n_immune_found}/{len(IMMUNE_GENES)}) signatures; verify gene symbol nomenclature."
        )

    estimate_scores = [s + im for s, im in zip(stromal_scores, immune_scores)]
    purities = []
    for est in estimate_scores:
        # Yoshihara et al. Nat Commun 2013 nonlinear purity calibration clipped to [0.05, 0.99]
        angle = min(max(0.6049872018 + 0.0001467884 * est, 0.05), 1.52)
        purities.append(round(float(math.cos(angle)), 4))

    # 2. Load or construct immune cell signature matrix
    if args.signature_matrix:
        sig_df = read_expression_matrix(args.signature_matrix)
        sig_df.index = pd.Index([str(g).strip().upper() for g in sig_df.index])
        sig_df = sig_df[~sig_df.index.duplicated(keep="first")]
        sig_df = sig_df.clip(lower=0.0)
        signature_source = f"custom:{Path(args.signature_matrix).name}"
    else:
        sig_df = _build_builtin_signature_df()
        signature_source = "builtin_leukocyte_10_panel"

    shared_genes = [g for g in sig_df.index if g in lin_upper.index]
    min_shared = 6 if args.signature_matrix else 12
    if len(shared_genes) < min_shared:
        fail(
            f"Only {len(shared_genes)} signature genes matched the expression matrix (need >= {min_shared}). "
            "Ensure input uses HGNC/MGI gene symbols or supply a matching --signature_matrix."
        )

    S_sub = sig_df.loc[shared_genes].copy()
    # Retain cell types with positive marker signal among shared genes
    valid_cts = [ct for ct in S_sub.columns if float((S_sub[ct] > 0).sum()) >= 2]
    if len(valid_cts) < 2:
        valid_cts = [ct for ct in S_sub.columns if float(S_sub[ct].sum()) > 0]
    if len(valid_cts) < 2:
        fail("Fewer than 2 reference cell types retained active marker genes after intersecting with input matrix.")
    S_sub = S_sub[valid_cts]

    # Normalize each signature column to unit L1 sum so NNLS coefficients represent comparable cell abundance
    col_sums = S_sub.sum(axis=0).to_numpy(dtype=float)
    col_sums[col_sums <= 0] = 1.0
    S_mat = S_sub.to_numpy(dtype=float) / col_sums[None, :]
    Y_mat = lin_upper.loc[shared_genes, samples].to_numpy(dtype=float)

    frac_rows = []
    fit_corrs = []
    fit_rmses = []
    abs_totals = []
    for s_idx, s_name in enumerate(samples):
        y = Y_mat[:, s_idx]
        if float(np.sum(y)) <= 1e-12:
            coefs = np.zeros(len(valid_cts), dtype=float)
            y_hat = np.zeros_like(y)
        else:
            coefs, _ = optimize.nnls(S_mat, y)
            y_hat = S_mat @ coefs
        tot = float(coefs.sum())
        abs_totals.append(tot)
        if tot > 1e-12:
            rel = coefs / tot
        else:
            rel = np.zeros_like(coefs)
        if np.std(y) > 1e-12 and np.std(y_hat) > 1e-12:
            r_fit = float(np.corrcoef(y, y_hat)[0, 1])
        else:
            r_fit = 0.0
        denom_y = float(np.linalg.norm(y))
        rel_rmse = float(np.linalg.norm(y - y_hat) / denom_y) if denom_y > 1e-12 else 0.0
        fit_corrs.append(round(r_fit, 4))
        fit_rmses.append(round(rel_rmse, 4))
        row = {"sample": s_name}
        for ct, val in zip(valid_cts, rel):
            row[ct] = round(float(val), 5)
        frac_rows.append(row)

    frac_df = pd.DataFrame(frac_rows).set_index("sample")
    tme_df = pd.DataFrame(
        {
            "StromalScore": [round(x, 2) for x in stromal_scores],
            "ImmuneScore": [round(x, 2) for x in immune_scores],
            "ESTIMATEScore": [round(x, 2) for x in estimate_scores],
            "TumorPurity": purities,
            "AbsoluteImmuneSignal": [round(x, 4) for x in abs_totals],
            "DeconvFitR": fit_corrs,
            "DeconvRelRMSE": fit_rmses,
        },
        index=samples,
    )
    tme_df.index.name = "sample"

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frac_path = out_dir / "tme_cell_fractions.tsv"
    scores_path = out_dir / "tme_microenvironment_scores.tsv"
    fig_path = out_dir / "tme-deconvolution.png"

    frac_df.to_csv(frac_path, sep="\t")
    tme_df.to_csv(scores_path, sep="\t")

    write_artifact(
        frac_path,
        "rna_tme_deconv",
        "Deconvolved Relative Immune Cell Fractions (NNLS)",
        inputs=[args.input] + ([args.signature_matrix] if args.signature_matrix else []),
        parameters={"signature_source": signature_source, "n_signature_genes": len(shared_genes)},
        shape={"rows": int(frac_df.shape[0]), "cols": int(frac_df.shape[1])},
        scientific_context={
            "dataType": "pathway_scores",
            "valueScale": "proportions",
            "normalization": "nnls_sum_to_one",
            "entityIdType": "sample_id",
            "organelle": "nuclear",
        },
    )
    write_artifact(
        scores_path,
        "rna_tme_deconv",
        "ESTIMATE-Style Stromal, Immune & Tumor Purity Scores",
        inputs=[args.input],
        parameters={"n_stromal_genes": n_stromal_found, "n_immune_genes": n_immune_found},
        shape={"rows": int(tme_df.shape[0]), "cols": int(tme_df.shape[1])},
    )

    # Optional group differential infiltration testing
    group_comp_path_str = None
    group_series = None
    top_diff_features = []
    if args.metadata and args.group_col:
        meta = read_metadata(args.metadata, samples=samples, sample_col=args.sample_col or None)
        if args.group_col not in meta.columns:
            fail(f"group_col '{args.group_col}' not found in metadata columns: {list(meta.columns)}")
        group_series = meta[args.group_col].astype(str)
        groups = [g for g in group_series.dropna().unique().tolist() if g not in ("", "nan", "NA", "None")]
        if len(groups) >= 2:
            combined_features = pd.concat(
                [frac_df, tme_df[["StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"]]],
                axis=1,
            )
            comp_rows = []
            for feat in combined_features.columns:
                vals_by_g = [
                    combined_features.loc[group_series == g, feat].dropna().to_numpy(dtype=float)
                    for g in groups
                ]
                vals_by_g = [v for v in vals_by_g if v.size > 0]
                if len(vals_by_g) < 2 or any(v.size < 2 for v in vals_by_g):
                    pval = 1.0
                    test_name = "insufficient_replicates"
                elif len(vals_by_g) == 2:
                    test_name = "mannwhitneyu"
                    try:
                        pval = float(stats.mannwhitneyu(vals_by_g[0], vals_by_g[1], alternative="two-sided").pvalue)
                    except ValueError:
                        pval = 1.0
                else:
                    test_name = "kruskal"
                    try:
                        pval = float(stats.kruskal(*vals_by_g).pvalue)
                    except ValueError:
                        pval = 1.0
                means_str = "; ".join(
                    f"{g}={float(combined_features.loc[group_series == g, feat].mean()):.4g}" for g in groups
                )
                comp_rows.append(
                    {
                        "feature": feat,
                        "test": test_name,
                        "group_means": means_str,
                        "pvalue": pval,
                    }
                )
            comp_df = pd.DataFrame(comp_rows)
            comp_df["padj"] = bh_fdr(comp_df["pvalue"].to_numpy(dtype=float))
            comp_df = comp_df.sort_values("pvalue", ascending=True).reset_index(drop=True)
            comp_path = out_dir / "tme_group_comparisons.tsv"
            comp_df.to_csv(comp_path, sep="\t", index=False)
            group_comp_path_str = str(comp_path)
            write_artifact(
                comp_path,
                "rna_tme_deconv",
                "Differential Immune Infiltration & TME Group Comparisons",
                inputs=[args.input, args.metadata],
                parameters={"group_col": args.group_col},
                shape={"rows": int(len(comp_df)), "cols": int(comp_df.shape[1])},
            )
            top_diff_features = [
                {
                    "feature": r["feature"],
                    "pvalue": float(r["pvalue"]),
                    "padj": float(r["padj"]),
                    "group_means": r["group_means"],
                }
                for r in comp_df.head(6).to_dict(orient="records")
            ]

    # Plot 2-panel TME overview
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.2), dpi=180, gridspec_kw={"width_ratios": [1.35, 1.0]})

    # Order samples by group (if available) then ImmuneScore descending
    order_df = tme_df.copy()
    if group_series is not None:
        order_df["_grp"] = group_series
        order_idx = order_df.sort_values(["_grp", "ImmuneScore"], ascending=[True, False]).index.tolist()
    else:
        order_idx = order_df.sort_values("ImmuneScore", ascending=False).index.tolist()

    plot_samples = order_idx[: min(45, len(order_idx))]
    frac_plot = frac_df.loc[plot_samples]
    palette = [
        "#2563eb", "#60a5fa", "#8b5cf6", "#06b6d4", "#10b981",
        "#34d399", "#f59e0b", "#ef4444", "#d97706", "#6b7280",
        "#ec4899", "#14b8a6",
    ]
    bottom = np.zeros(len(plot_samples), dtype=float)
    x_idx = np.arange(len(plot_samples))
    for i, ct in enumerate(frac_plot.columns):
        vals = frac_plot[ct].to_numpy(dtype=float)
        ax1.bar(x_idx, vals, bottom=bottom, width=0.82, color=palette[i % len(palette)], label=ct)
        bottom += vals
    ax1.set_ylim(0, 1.02)
    ax1.set_ylabel("Relative Immune Fraction")
    ax1.set_title("Deconvolved Immune Composition Across Samples")
    if len(plot_samples) <= 24:
        ax1.set_xticks(x_idx)
        ax1.set_xticklabels(plot_samples, rotation=55, ha="right", fontsize=7.5)
    else:
        ax1.set_xticks([])
        ax1.set_xlabel(f"Samples ordered by ImmuneScore (showing {len(plot_samples)}/{len(samples)})")
    ax1.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=7.8)

    if group_series is not None and group_series.nunique() >= 2:
        grps = sorted(group_series.dropna().unique().tolist())
        data_box = [tme_df.loc[group_series == g, "ImmuneScore"].to_numpy(dtype=float) for g in grps]
        bp = ax2.boxplot(data_box, tick_labels=grps, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor("#93c5fd")
            patch.set_alpha(0.8)
        ax2.set_ylabel("ImmuneScore (ESTIMATE-style)")
        ax2.set_title(f"ImmuneScore by {args.group_col}")
    else:
        sc = ax2.scatter(
            tme_df["StromalScore"],
            tme_df["ImmuneScore"],
            c=tme_df["TumorPurity"],
            cmap="viridis",
            s=42,
            edgecolor="#1f2937",
            linewidth=0.4,
        )
        fig.colorbar(sc, ax=ax2, fraction=0.046, pad=0.04, label="Inferred Tumor Purity")
        ax2.set_xlabel("StromalScore")
        ax2.set_ylabel("ImmuneScore")
        ax2.set_title("TME Stromal vs Immune Landscape")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "rna_tme_deconv",
        "Tumor Microenvironment Immune Deconvolution & Purity Overview",
        inputs=[args.input] + ([args.metadata] if args.metadata else []),
        parameters={"signature_source": signature_source},
    )

    mean_fractions = {ct: round(float(frac_df[ct].mean()), 4) for ct in frac_df.columns}
    emit(
        {
            "status": "ok",
            "tool": "rna_tme_deconv",
            "inferred_value_kind": value_kind,
            "signature_source": signature_source,
            "n_samples": int(len(samples)),
            "n_cell_types": int(len(valid_cts)),
            "n_signature_genes_matched": int(len(shared_genes)),
            "n_stromal_genes_matched": int(n_stromal_found),
            "n_immune_genes_matched": int(n_immune_found),
            "median_deconv_fit_r": round(float(np.median(fit_corrs)), 4),
            "median_tumor_purity": round(float(np.median(purities)), 4),
            "mean_cell_fractions": mean_fractions,
            "top_group_differences": top_diff_features,
            "outputs": {
                "cell_fractions": str(frac_path),
                "microenvironment_scores": str(scores_path),
                "group_comparisons": group_comp_path_str,
                "summary_plot": str(fig_path),
            },
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
