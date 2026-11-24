#!/usr/bin/env python3
"""rna_gsea: offline deterministic pre-ranked GSEA (Subramanian et al. 2005
weighted Kolmogorov-Smirnov enrichment score, gene-set permutation NES/p/BH FDR,
leading-edge genes, and enrichment running-sum curves) in pure numpy/scipy.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import bh_fdr, emit, fail, read_de_table, read_gmt, setup_matplotlib, write_artifact


def compute_es_curve(weights, hit_mask):
    n = len(weights)
    n_hit = int(hit_mask.sum())
    n_miss = n - n_hit
    if n_hit == 0 or n_miss == 0:
        return 0.0, np.zeros(n, dtype=float), 0
    hit_w = weights * hit_mask
    nr = float(hit_w.sum())
    if nr <= 0.0:
        hit_inc = hit_mask.astype(float) / n_hit
    else:
        hit_inc = hit_w / nr
    miss_dec = (~hit_mask).astype(float) / n_miss
    running = np.cumsum(hit_inc - miss_dec)
    max_pos = float(running.max())
    min_neg = float(running.min())
    if abs(max_pos) >= abs(min_neg):
        peak_idx = int(np.argmax(running))
        es = max_pos
    else:
        peak_idx = int(np.argmin(running))
        es = min_neg
    return es, running, peak_idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de_table", required=True, help="Full unfiltered DE table (DESeq2/edgeR/limma).")
    ap.add_argument("--gmt", required=True, help="GMT gene-set file.")
    ap.add_argument("--rank_by", choices=["stat", "log2fc"], default="stat")
    ap.add_argument("--weight", type=float, default=1.0, help="GSEA weighting exponent p (default 1.0).")
    ap.add_argument("--min_set_size", type=int, default=5)
    ap.add_argument("--max_set_size", type=int, default=500)
    ap.add_argument("--permutations", type=int, default=500)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    de = read_de_table(args.de_table)
    de = de.dropna(subset=[args.rank_by]).sort_values(args.rank_by, ascending=False).reset_index(drop=True)
    if len(de) < 20:
        fail(f"GSEA requires at least 20 ranked genes (got {len(de)})")

    genes = de["gene"].to_numpy()
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    scores = de[args.rank_by].to_numpy(dtype=float)
    weights = np.abs(scores) ** args.weight
    n_genes = len(genes)

    gene_sets = read_gmt(args.gmt)
    rng = np.random.default_rng(args.seed)

    rows = []
    curves = {}
    for gs_name, gs_genes in gene_sets.items():
        hit_idx = [gene_to_idx[g] for g in gs_genes if g in gene_to_idx]
        m = len(hit_idx)
        if m < args.min_set_size or m > args.max_set_size:
            continue
        hit_mask = np.zeros(n_genes, dtype=bool)
        hit_mask[hit_idx] = True
        es_obs, running, peak_idx = compute_es_curve(weights, hit_mask)

        # Leading-edge genes
        if es_obs >= 0:
            le_idx = [i for i in hit_idx if i <= peak_idx]
        else:
            le_idx = [i for i in hit_idx if i >= peak_idx]
        le_genes = [str(genes[i]) for i in sorted(le_idx)]

        # Seeded gene-set permutations (exact same set size m)
        perm_es = np.empty(args.permutations, dtype=float)
        for p_i in range(args.permutations):
            chosen = rng.choice(n_genes, size=m, replace=False)
            pmask = np.zeros(n_genes, dtype=bool)
            pmask[chosen] = True
            perm_es[p_i], _, _ = compute_es_curve(weights, pmask)

        if es_obs >= 0:
            pos_null = perm_es[perm_es >= 0]
            mean_pos = float(pos_null.mean()) if len(pos_null) > 0 else 1.0
            nes = es_obs / mean_pos if mean_pos > 1e-12 else 0.0
            pval = float((1 + int((pos_null >= es_obs).sum())) / (1 + len(pos_null)))
        else:
            neg_null = perm_es[perm_es <= 0]
            mean_neg = float(np.abs(neg_null).mean()) if len(neg_null) > 0 else 1.0
            nes = es_obs / mean_neg if mean_neg > 1e-12 else 0.0
            pval = float((1 + int((neg_null <= es_obs).sum())) / (1 + len(neg_null)))

        rows.append(
            {
                "gene_set": gs_name,
                "set_size": int(m),
                "es": round(float(es_obs), 4),
                "nes": round(float(nes), 4),
                "p_value": pval,
                "n_leading_edge": int(len(le_genes)),
                "leading_edge_genes": ";".join(le_genes),
            }
        )
        curves[gs_name] = (running, hit_mask, nes, pval)

    if not rows:
        fail(f"no GMT gene sets matched between {args.min_set_size} and {args.max_set_size} genes")

    res_df = pd.DataFrame(rows)
    res_df["q_value"] = bh_fdr(res_df["p_value"].to_numpy())
    res_df["significant"] = res_df["q_value"] <= args.fdr
    res_df = res_df.sort_values(["q_value", "p_value", "gene_set"]).reset_index(drop=True)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = [args.de_table, args.gmt]
    params = {
        "rank_by": args.rank_by,
        "weight": args.weight,
        "permutations": args.permutations,
        "fdr": args.fdr,
        "seed": args.seed,
    }
    tsv_path = outdir / "gsea-results.tsv"
    res_df.to_csv(tsv_path, sep="\t", index=False)

    # Plot top up to 4 enriched running-sum curves
    top_sets = res_df.head(4)["gene_set"].tolist()
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    palette = ["#E11D48", "#0EA5E9", "#10B981", "#F59E0B"]
    for idx_s, gs in enumerate(top_sets):
        running, _, nes, _ = curves[gs]
        q_val = float(res_df.loc[res_df["gene_set"] == gs, "q_value"].iloc[0])
        ax.plot(running, color=palette[idx_s % len(palette)], lw=2.0, label=f"{gs[:30]} (NES={nes:.2f}, q={q_val:.3g})")
    ax.axhline(0.0, color="#475569", lw=1.0, ls="--")
    ax.set_xlabel(f"Gene rank (sorted by {args.rank_by})")
    ax.set_ylabel("Running Enrichment Score (ES)")
    ax.set_title("Pre-ranked GSEA enrichment curves (top pathways)")
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig_path = outdir / "gsea-curves.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    artifacts = [
        write_artifact(tsv_path, "rna_gsea", "Pre-ranked GSEA enrichment table", inputs=inputs, parameters=params, shape=res_df.shape),
        write_artifact(fig_path, "rna_gsea", "Pre-ranked GSEA running-sum curves", inputs=inputs, parameters=params),
    ]

    emit(
        {
            "tool": "rna_gsea",
            "n_genes_ranked": int(n_genes),
            "rank_by": args.rank_by,
            "n_sets_tested": int(len(res_df)),
            "n_significant": int(res_df["significant"].sum()),
            "top_pathways": res_df.head(10).drop(columns=["leading_edge_genes"]).to_dict(orient="records"),
            "outputs": {"table": str(tsv_path), "curves_figure": str(fig_path)},
            "artifacts": artifacts,
            "caveat": "Gene-set permutation GSEA assumes gene independence and can yield slightly optimistic p-values compared to phenotype permutation on large cohorts.",
        }
    )


if __name__ == "__main__":
    main()
