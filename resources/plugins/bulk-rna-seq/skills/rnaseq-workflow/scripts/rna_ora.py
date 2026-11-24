#!/usr/bin/env python3
"""rna_ora: offline hypergeometric Over-Representation Analysis (ORA) from a
DE table (or explicit gene list + background universe) and a GMT gene-set file.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from _common import bh_fdr, emit, fail, read_de_table, read_gmt, setup_matplotlib, write_artifact


def _read_gene_list(path):
    p = Path(path)
    if not p.is_file():
        fail(f"gene list file not found: {p}")
    genes = []
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            genes.append(s.split("\t")[0].split(",")[0].strip())
    return sorted(set(g for g in genes if g))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de_table", help="Full DESeq2/edgeR/limma table (supplies both foreground DEGs and background universe).")
    ap.add_argument("--genes", help="Explicit foreground gene list (one gene per line) if --de_table is omitted.")
    ap.add_argument("--background", help="Explicit background universe gene list (required when --genes is used without --de_table).")
    ap.add_argument("--gmt", required=True, help="GMT gene-set file.")
    ap.add_argument("--direction", choices=["both", "up", "down"], default="both")
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--log2fc", type=float, default=1.0)
    ap.add_argument("--min_set_size", type=int, default=5)
    ap.add_argument("--max_set_size", type=int, default=500)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    if args.de_table:
        de = read_de_table(args.de_table)
        universe = set(de["gene"].tolist())
        if args.background:
            universe = set(_read_gene_list(args.background))
        if args.direction == "up":
            fg = set(de.loc[(de["padj"] < args.padj) & (de["log2fc"] >= args.log2fc), "gene"])
        elif args.direction == "down":
            fg = set(de.loc[(de["padj"] < args.padj) & (de["log2fc"] <= -args.log2fc), "gene"])
        else:
            fg = set(de.loc[(de["padj"] < args.padj) & (de["log2fc"].abs() >= args.log2fc), "gene"])
        inputs = [args.de_table, args.gmt] + ([args.background] if args.background else [])
    elif args.genes and args.background:
        fg = set(_read_gene_list(args.genes))
        universe = set(_read_gene_list(args.background))
        inputs = [args.genes, args.background, args.gmt]
    else:
        fail("provide --de_table, or both --genes and --background (explicit universe is mandatory for ORA)")

    fg = fg & universe
    n_univ = len(universe)
    n_fg = len(fg)
    if n_univ < 10:
        fail(f"background universe too small ({n_univ} genes)")
    if n_fg == 0:
        fail("zero foreground genes matched the criteria and background universe")

    gene_sets = read_gmt(args.gmt)
    rows = []
    for gs_name, gs_genes in gene_sets.items():
        gs_in_u = set(gs_genes) & universe
        m_size = len(gs_in_u)
        if m_size < args.min_set_size or m_size > args.max_set_size:
            continue
        overlap = sorted(fg & gs_in_u)
        k = len(overlap)
        pval = float(hypergeom.sf(k - 1, n_univ, m_size, n_fg)) if k > 0 else 1.0
        expected = (n_fg * m_size) / float(n_univ)
        fold_enrich = (k / expected) if expected > 0 else 0.0
        rows.append(
            {
                "gene_set": gs_name,
                "overlap_count": int(k),
                "set_size_in_universe": int(m_size),
                "gene_ratio": round(k / float(n_fg), 4),
                "fold_enrichment": round(fold_enrich, 4),
                "p_value": pval,
                "overlap_genes": ";".join(overlap),
            }
        )

    if not rows:
        fail(
            f"no GMT gene sets had between {args.min_set_size} and {args.max_set_size} genes overlapping the universe"
        )

    res_df = pd.DataFrame(rows)
    res_df["q_value"] = bh_fdr(res_df["p_value"].to_numpy())
    res_df["significant"] = res_df["q_value"] <= args.fdr
    res_df = res_df.sort_values(["q_value", "p_value", "gene_set"]).reset_index(drop=True)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    params = {
        "direction": args.direction,
        "padj": args.padj,
        "log2fc": args.log2fc,
        "min_set_size": args.min_set_size,
        "max_set_size": args.max_set_size,
        "fdr": args.fdr,
    }
    tsv_path = outdir / "ora-enrichment.tsv"
    res_df.to_csv(tsv_path, sep="\t", index=False)

    # Dotplot of top 15 gene sets
    top_df = res_df.head(15).iloc[::-1]
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(7.5, max(4, 0.35 * len(top_df) + 1.5)))
    y_pos = np.arange(len(top_df))
    neg_log_q = -np.log10(np.maximum(top_df["q_value"].to_numpy(dtype=float), 1e-300))
    sizes = 40 + 18 * top_df["overlap_count"].to_numpy(dtype=float)
    sc = ax.scatter(top_df["fold_enrichment"], y_pos, s=sizes, c=neg_log_q, cmap="viridis", alpha=0.9)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([s[:40] for s in top_df["gene_set"]], fontsize=8)
    ax.set_xlabel("Fold enrichment")
    ax.set_title(f"ORA pathway enrichment ({args.direction}, n_fg={n_fg}, universe={n_univ})")
    cbar = fig.colorbar(sc, ax=ax, shrink=0.8)
    cbar.set_label("-log10(q-value)")
    fig.tight_layout()
    fig_path = outdir / "ora-dotplot.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    artifacts = [
        write_artifact(tsv_path, "rna_ora", "Hypergeometric ORA enrichment table", inputs=inputs, parameters=params, shape=res_df.shape),
        write_artifact(fig_path, "rna_ora", "ORA enrichment dotplot", inputs=inputs, parameters=params),
    ]

    emit(
        {
            "tool": "rna_ora",
            "direction": args.direction,
            "n_foreground": int(n_fg),
            "n_universe": int(n_univ),
            "n_sets_tested": int(len(res_df)),
            "n_significant": int(res_df["significant"].sum()),
            "top_pathways": res_df.head(10).drop(columns=["overlap_genes"]).to_dict(orient="records"),
            "outputs": {"table": str(tsv_path), "dotplot_figure": str(fig_path)},
            "artifacts": artifacts,
            "caveat": "ORA depends on foreground cutoff choices and background universe; complement with threshold-free rna_gsea when pathway shifts are subtle.",
        }
    )


if __name__ == "__main__":
    main()
