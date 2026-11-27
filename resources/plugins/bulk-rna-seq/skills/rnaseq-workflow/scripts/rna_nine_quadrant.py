#!/usr/bin/env python3
"""Cross-modal 9-quadrant integration between bulk RNA-seq DEGs and another omics layer
(ATAC-seq, ChIP-seq / CUT&Tag, proteomics, or DNA methylation).

All file paths, statistical metrics, and thresholds come strictly from the inputs.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import emit, fail, read_de_table, setup_matplotlib, write_artifact  # noqa: E402


QUADRANT_DEFS = [
    ("Q1_rna_up_other_up", 1, 1, "concordant_up", "RNA Up & {label} Up"),
    ("Q2_rna_ns_other_up", 0, 1, "other_only", "RNA NS & {label} Up"),
    ("Q3_rna_down_other_up", -1, 1, "discordant", "RNA Down & {label} Up"),
    ("Q4_rna_up_other_ns", 1, 0, "rna_only", "RNA Up & {label} NS"),
    ("Q5_rna_ns_other_ns", 0, 0, "ns", "Not Significant"),
    ("Q6_rna_down_other_ns", -1, 0, "rna_only", "RNA Down & {label} NS"),
    ("Q7_rna_up_other_down", 1, -1, "discordant", "RNA Up & {label} Down"),
    ("Q8_rna_ns_other_down", 0, -1, "other_only", "RNA NS & {label} Down"),
    ("Q9_rna_down_other_down", -1, -1, "concordant_down", "RNA Down & {label} Down"),
]

QUADRANT_LOOKUP = {(r, o): (q, c) for q, r, o, c, _ in QUADRANT_DEFS}


def _resolve_col(df, explicit, candidates, kind, required=True):
    lower_map = {str(c).strip().lower(): str(c) for c in df.columns}
    if explicit:
        key = explicit.strip()
        if key in df.columns:
            return key
        if key.lower() in lower_map:
            return lower_map[key.lower()]
        fail(f"specified --other_{kind}_col {explicit!r} not found in columns: {list(df.columns)}")
    for cand in candidates:
        if cand in lower_map:
            return lower_map[cand]
    if required:
        fail(f"could not auto-detect {kind} column in other_de; columns: {list(df.columns)}")
    return None


def _read_other_table(args):
    path = Path(args.other_de)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    raw = pd.read_csv(path, sep=sep)
    if raw.empty:
        fail(f"other_de table is empty: {path}")
    raw.columns = [str(c).strip() for c in raw.columns]

    gene_col = _resolve_col(
        raw,
        args.other_gene_col,
        [
            "gene",
            "gene_name",
            "symbol",
            "gene_symbol",
            "external_gene_name",
            "annotation_gene",
            "nearest_gene",
            "target",
            "id",
        ],
        "gene",
        required=False,
    )
    if gene_col is None:
        gene_col = raw.columns[0]

    lfc_col = _resolve_col(
        raw,
        args.other_lfc_col,
        [
            "log2fc",
            "log2foldchange",
            "logfc",
            "fold",
            "log2ratio",
            "diff",
            "meth.diff",
            "delta",
            "effect_size",
            "stat",
        ],
        "lfc",
        required=True,
    )
    padj_col = _resolve_col(
        raw,
        args.other_padj_col,
        [
            "padj",
            "fdr",
            "adj.p.val",
            "qvalue",
            "p.adjust",
            "pvalue",
            "pval",
            "p.value",
            "p_value",
        ],
        "padj",
        required=bool(args.other_padj_col),
    )
    dist_col = _resolve_col(
        raw,
        args.other_dist_col,
        ["distancetotss", "distance_to_tss", "dist_to_tss", "tss_distance", "distance"],
        "dist",
        required=bool(args.other_dist_col),
    )

    df = pd.DataFrame(
        {
            "gene": raw[gene_col].astype(str).str.strip(),
            "other_log2fc": pd.to_numeric(raw[lfc_col], errors="coerce"),
            "other_padj": pd.to_numeric(raw[padj_col], errors="coerce") if padj_col else 0.0,
        }
    )
    if dist_col:
        df["other_dist"] = pd.to_numeric(raw[dist_col], errors="coerce").abs()
    else:
        df["other_dist"] = np.nan

    df = df.dropna(subset=["other_log2fc"])
    df = df[(df["gene"] != "") & (df["gene"].str.lower() != "nan")]
    df["other_padj"] = df["other_padj"].fillna(1.0)
    if df.empty:
        fail(f"no valid numeric rows parsed from other_de: {path}")

    used_strategy = args.peak_agg
    if used_strategy == "nearest_tss" and df["other_dist"].notna().sum() == 0:
        used_strategy = "max_abs"

    if used_strategy == "mean":
        agg = (
            df.groupby("gene", as_index=False)
            .agg(
                other_log2fc=("other_log2fc", "mean"),
                other_padj=("other_padj", "min"),
                other_n_features=("other_log2fc", "size"),
            )
        )
    elif used_strategy == "nearest_tss":
        df = df.assign(
            _dist_rank=df["other_dist"].fillna(np.inf),
            _abs_lfc=-df["other_log2fc"].abs(),
        )
        counts = df.groupby("gene")["other_log2fc"].size().rename("other_n_features")
        picked = (
            df.sort_values(["_dist_rank", "_abs_lfc", "other_padj"])
            .drop_duplicates(subset=["gene"], keep="first")
            .loc[:, ["gene", "other_log2fc", "other_padj"]]
        )
        agg = picked.merge(counts.reset_index(), on="gene", how="left")
    else:
        df = df.assign(_abs_lfc=-df["other_log2fc"].abs())
        counts = df.groupby("gene")["other_log2fc"].size().rename("other_n_features")
        picked = (
            df.sort_values(["_abs_lfc", "other_padj"])
            .drop_duplicates(subset=["gene"], keep="first")
            .loc[:, ["gene", "other_log2fc", "other_padj"]]
        )
        agg = picked.merge(counts.reset_index(), on="gene", how="left")

    return agg, used_strategy, {
        "gene_col": gene_col,
        "lfc_col": lfc_col,
        "padj_col": padj_col,
        "dist_col": dist_col,
    }


def _classify_dir(lfc_series, padj_series, lfc_thr, padj_thr):
    sig = padj_series <= padj_thr
    out = np.zeros(len(lfc_series), dtype=int)
    out[(sig) & (lfc_series >= lfc_thr)] = 1
    out[(sig) & (lfc_series <= -lfc_thr)] = -1
    return out


def _top_records(df, limit=10):
    cols = ["gene", "rna_log2fc", "rna_padj", "other_log2fc", "other_padj", "quadrant"]
    records = []
    for _, row in df.head(limit).iterrows():
        records.append(
            {
                "gene": str(row["gene"]),
                "rna_log2fc": round(float(row["rna_log2fc"]), 4),
                "rna_padj": float(row["rna_padj"]),
                "other_log2fc": round(float(row["other_log2fc"]), 4),
                "other_padj": float(row["other_padj"]),
                "quadrant": str(row["quadrant"]),
            }
        )
    return records


def _plot_nine_quadrant(df, args, r_val, rho_val, fisher_p, counts, out_png):
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(8.0, 7.0), dpi=180)

    palette = {
        "concordant_up": ("#c0392b", "Concordant Up", 26, 0.82),
        "concordant_down": ("#2980b9", "Concordant Down", 26, 0.82),
        "discordant": ("#8e44ad", "Discordant", 24, 0.78),
        "rna_only": ("#e67e22", "RNA Only", 14, 0.45),
        "other_only": ("#16a085", f"{args.other_label} Only", 14, 0.45),
        "ns": ("#bdc3c7", "Not Significant", 9, 0.30),
    }

    draw_order = ["ns", "rna_only", "other_only", "discordant", "concordant_down", "concordant_up"]
    for cls_name in draw_order:
        sub = df[df["concordance_class"] == cls_name]
        color, label, size, alpha = palette[cls_name]
        ax.scatter(
            sub["other_log2fc"],
            sub["rna_log2fc"],
            c=color,
            s=size,
            alpha=alpha,
            edgecolors="none",
            label=f"{label} (n={len(sub)})",
        )

    ax.axvline(args.other_lfc, color="#555555", linestyle="--", linewidth=0.9)
    ax.axvline(-args.other_lfc, color="#555555", linestyle="--", linewidth=0.9)
    ax.axhline(args.rna_lfc, color="#555555", linestyle="--", linewidth=0.9)
    ax.axhline(-args.rna_lfc, color="#555555", linestyle="--", linewidth=0.9)

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_hi = max(xlim[1], args.other_lfc * 1.6)
    x_lo = min(xlim[0], -args.other_lfc * 1.6)
    y_hi = max(ylim[1], args.rna_lfc * 1.6)
    y_lo = min(ylim[0], -args.rna_lfc * 1.6)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)

    corner_ann = [
        (0.97, 0.96, "right", "top", f"Q1 Up/Up: {counts['Q1_rna_up_other_up']}", "#c0392b"),
        (0.03, 0.96, "left", "top", f"Q7 RNA+/Other-: {counts['Q7_rna_up_other_down']}", "#8e44ad"),
        (0.97, 0.04, "right", "bottom", f"Q3 RNA-/Other+: {counts['Q3_rna_down_other_up']}", "#8e44ad"),
        (0.03, 0.04, "left", "bottom", f"Q9 Down/Down: {counts['Q9_rna_down_other_down']}", "#2980b9"),
    ]
    for tx, ty, ha, va, text, col in corner_ann:
        ax.text(
            tx,
            ty,
            text,
            transform=ax.transAxes,
            ha=ha,
            va=va,
            fontsize=8.5,
            fontweight="bold",
            color=col,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.8, "edgecolor": col, "linewidth": 0.6},
        )

    if args.top_labels > 0:
        priority = df[df["concordance_class"].isin(["concordant_up", "concordant_down", "discordant"])]
        if priority.empty:
            priority = df[df["concordance_class"] != "ns"]
        to_label = priority.sort_values("combined_score", ascending=False).head(args.top_labels)
        for idx, (_, row) in enumerate(to_label.iterrows()):
            dx = 6 if idx % 2 == 0 else -6
            dy = 5 if idx % 3 != 0 else -8
            ax.annotate(
                str(row["gene"]),
                xy=(float(row["other_log2fc"]), float(row["rna_log2fc"])),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=7.5,
                alpha=0.9,
            )

    ax.set_xlabel(f"{args.other_label} log2FC / Effect Size", fontsize=11)
    ax.set_ylabel("RNA-seq log2FC", fontsize=11)
    ax.set_title(
        f"9-Quadrant Integration: RNA-seq vs {args.other_label}\n"
        f"n={len(df)} | Pearson r={r_val:.3f} | Spearman rho={rho_val:.3f} | Fisher concordance p={fisher_p:.2e}",
        fontsize=11,
        pad=10,
    )
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Cross-modal 9-quadrant DE integration.")
    parser.add_argument("--rna_de", required=True, help="RNA-seq DE table TSV/CSV.")
    parser.add_argument("--other_de", required=True, help="Second omics DE/DA table TSV/CSV.")
    parser.add_argument("--other_label", default="Other", help="Label for second modality (default Other).")
    parser.add_argument("--other_gene_col", default=None, help="Gene symbol column in other_de.")
    parser.add_argument("--other_lfc_col", default=None, help="Effect size / log2FC column in other_de.")
    parser.add_argument("--other_padj_col", default=None, help="Adjusted p-value column in other_de.")
    parser.add_argument("--other_dist_col", default=None, help="Optional TSS distance column for nearest_tss.")
    parser.add_argument("--rna_lfc", type=float, default=1.0, help="Abs log2FC threshold for RNA DEGs.")
    parser.add_argument("--rna_padj", type=float, default=0.05, help="Adjusted p-value threshold for RNA DEGs.")
    parser.add_argument("--other_lfc", type=float, default=0.585, help="Abs log2FC threshold for second modality.")
    parser.add_argument("--other_padj", type=float, default=0.05, help="Adjusted p-value threshold for second modality.")
    parser.add_argument(
        "--peak_agg",
        choices=["max_abs", "mean", "nearest_tss"],
        default="max_abs",
        help="Aggregation strategy when multiple features in other_de map to one gene.",
    )
    parser.add_argument("--top_labels", type=int, default=12, help="Number of extreme genes to label on the plot.")
    parser.add_argument("--output_dir", required=True, help="Directory for tables and figure.")
    args = parser.parse_args()

    rna_df = read_de_table(args.rna_de).rename(
        columns={"log2fc": "rna_log2fc", "pvalue": "rna_pvalue", "padj": "rna_padj"}
    )
    rna_df["rna_padj"] = pd.to_numeric(rna_df["rna_padj"], errors="coerce").fillna(1.0)
    rna_df = rna_df.dropna(subset=["rna_log2fc"])

    other_df, used_agg, resolved_cols = _read_other_table(args)

    merged = rna_df[["gene", "rna_log2fc", "rna_padj"]].merge(other_df, on="gene", how="inner")
    if len(merged) < 3:
        fail(f"fewer than 3 shared genes between rna_de and other_de ({len(merged)} matched)")

    rna_dir = _classify_dir(merged["rna_log2fc"], merged["rna_padj"], args.rna_lfc, args.rna_padj)
    other_dir = _classify_dir(merged["other_log2fc"], merged["other_padj"], args.other_lfc, args.other_padj)

    quadrants = []
    classes = []
    for rd, od in zip(rna_dir, other_dir):
        q, c = QUADRANT_LOOKUP[(int(rd), int(od))]
        quadrants.append(q)
        classes.append(c)

    merged["rna_dir"] = rna_dir
    merged["other_dir"] = other_dir
    merged["quadrant"] = quadrants
    merged["concordance_class"] = classes
    merged["combined_score"] = merged["rna_log2fc"].abs() + merged["other_log2fc"].abs()

    class_priority = {
        "concordant_up": 0,
        "concordant_down": 1,
        "discordant": 2,
        "rna_only": 3,
        "other_only": 4,
        "ns": 5,
    }
    merged["_prio"] = merged["concordance_class"].map(class_priority)
    merged = merged.sort_values(["_prio", "combined_score"], ascending=[True, False]).drop(columns=["_prio"])

    pearson_r, pearson_p = stats.pearsonr(merged["rna_log2fc"], merged["other_log2fc"])
    spearman_rho, spearman_p = stats.spearmanr(merged["rna_log2fc"], merged["other_log2fc"])
    if not np.isfinite(pearson_r):
        pearson_r, pearson_p = 0.0, 1.0
    if not np.isfinite(spearman_rho):
        spearman_rho, spearman_p = 0.0, 1.0

    q_counts = {q: int((merged["quadrant"] == q).sum()) for q, _, _, _, _ in QUADRANT_DEFS}
    n_uu = q_counts["Q1_rna_up_other_up"]
    n_ud = q_counts["Q7_rna_up_other_down"]
    n_du = q_counts["Q3_rna_down_other_up"]
    n_dd = q_counts["Q9_rna_down_other_down"]
    table_2x2 = [[n_uu, n_ud], [n_du, n_dd]]
    if (n_uu + n_ud > 0) and (n_du + n_dd > 0) and (n_uu + n_du > 0) and (n_ud + n_dd > 0):
        fisher_or, fisher_p = stats.fisher_exact(table_2x2, alternative="two-sided")
        if not np.isfinite(fisher_or):
            fisher_or = float("inf")
    else:
        fisher_or, fisher_p = 1.0, 1.0

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genes_tsv = out_dir / "nine-quadrant-genes.tsv"
    summary_tsv = out_dir / "nine-quadrant-summary.tsv"
    plot_png = out_dir / "nine-quadrant-plot.png"

    export_cols = [
        "gene",
        "rna_log2fc",
        "rna_padj",
        "other_log2fc",
        "other_padj",
        "other_n_features",
        "quadrant",
        "concordance_class",
        "combined_score",
    ]
    merged[export_cols].to_csv(genes_tsv, sep="\t", index=False)

    summary_rows = []
    n_total = len(merged)
    for q, rd, od, cls_name, interp in QUADRANT_DEFS:
        sub = merged[merged["quadrant"] == q]
        summary_rows.append(
            {
                "quadrant": q,
                "rna_state": {1: "up", -1: "down", 0: "ns"}[rd],
                "other_state": {1: "up", -1: "down", 0: "ns"}[od],
                "concordance_class": cls_name,
                "interpretation": interp.format(label=args.other_label),
                "n_genes": int(len(sub)),
                "fraction": round(float(len(sub) / n_total), 6),
                "mean_rna_log2fc": round(float(sub["rna_log2fc"].mean()), 4) if len(sub) else 0.0,
                "mean_other_log2fc": round(float(sub["other_log2fc"].mean()), 4) if len(sub) else 0.0,
                "top_genes": ",".join(sub.sort_values("combined_score", ascending=False)["gene"].head(5).tolist()),
            }
        )
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(summary_tsv, sep="\t", index=False)

    _plot_nine_quadrant(merged, args, float(pearson_r), float(spearman_rho), float(fisher_p), q_counts, plot_png)

    sci_ctx = {
        "modality": "bulk_rna_seq_multiomics_integration",
        "otherModalityLabel": args.other_label,
        "peakAggregation": used_agg,
        "valueSemantics": "cross_modal_log2fc",
        "nSharedFeatures": int(n_total),
        "pearsonR": round(float(pearson_r), 4),
        "spearmanRho": round(float(spearman_rho), 4),
        "fisherConcordanceP": float(fisher_p),
        "thresholds": {
            "rna_lfc": args.rna_lfc,
            "rna_padj": args.rna_padj,
            "other_lfc": args.other_lfc,
            "other_padj": args.other_padj,
        },
        "resolvedOtherColumns": resolved_cols,
    }

    params = {
        "other_label": args.other_label,
        "peak_agg": used_agg,
        "rna_lfc": args.rna_lfc,
        "rna_padj": args.rna_padj,
        "other_lfc": args.other_lfc,
        "other_padj": args.other_padj,
        "top_labels": args.top_labels,
    }
    inputs = [args.rna_de, args.other_de]
    artifacts = [
        write_artifact(
            genes_tsv,
            "rna_nine_quadrant",
            f"9-Quadrant Gene Table (RNA vs {args.other_label})",
            inputs=inputs,
            parameters=params,
            shape=(len(merged), len(export_cols)),
            scientific_context=sci_ctx,
        ),
        write_artifact(
            summary_tsv,
            "rna_nine_quadrant",
            f"9-Quadrant Summary (RNA vs {args.other_label})",
            inputs=inputs,
            parameters=params,
            shape=(len(summary_df), len(summary_df.columns)),
            scientific_context=sci_ctx,
        ),
        write_artifact(
            plot_png,
            "rna_nine_quadrant",
            f"9-Quadrant Scatter Plot (RNA vs {args.other_label})",
            inputs=inputs,
            parameters=params,
            scientific_context=sci_ctx,
        ),
    ]

    emit(
        {
            "tool": "rna_nine_quadrant",
            "other_label": args.other_label,
            "peak_agg": used_agg,
            "thresholds": sci_ctx["thresholds"],
            "n_shared_features": int(n_total),
            "pearson_r": round(float(pearson_r), 6),
            "pearson_p": float(pearson_p),
            "spearman_rho": round(float(spearman_rho), 6),
            "spearman_p": float(spearman_p),
            "fisher_concordance_or": float(fisher_or) if np.isfinite(fisher_or) else 999999.0,
            "fisher_concordance_p": float(fisher_p),
            "quadrant_counts": q_counts,
            "top_concordant_up": _top_records(merged[merged["concordance_class"] == "concordant_up"]),
            "top_concordant_down": _top_records(merged[merged["concordance_class"] == "concordant_down"]),
            "top_discordant": _top_records(merged[merged["concordance_class"] == "discordant"]),
            "outputs": {
                "genes_tsv": str(genes_tsv),
                "summary_tsv": str(summary_tsv),
                "plot_png": str(plot_png),
            },
            "scientificContext": sci_ctx,
            "artifacts": artifacts,
        }
    )


if __name__ == "__main__":
    main()
