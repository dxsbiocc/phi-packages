#!/usr/bin/env python3
"""prot_enrichment: offline detected-proteome ORA and deterministic pre-ranked GSEA."""

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
from _common import bh_fdr, emit, fail, normalize_gene_id, read_gmt, read_table, write_artifact


def run_ora(
    foreground: set[str],
    universe: set[str],
    gmt: dict[str, list[str]],
    min_size: int,
    max_size: int,
) -> pd.DataFrame:
    M = len(universe)
    N = len(foreground)
    rows = []
    for pathway, genes in gmt.items():
        gset = set(genes) & universe
        K = len(gset)
        if K < min_size or K > max_size:
            continue
        hits = sorted(foreground & gset)
        k = len(hits)
        pval = float(stats.hypergeom.sf(k - 1, M, K, N)) if k > 0 else 1.0
        expected = (K * N) / M if M > 0 else 0.0
        enrichment_ratio = (k / expected) if expected > 0 else 0.0
        rows.append(
            {
                "pathway": pathway,
                "overlap": k,
                "set_size": K,
                "enrichment_ratio": round(enrichment_ratio, 4),
                "pvalue": pval,
                "genes": ";".join(hits),
            }
        )
    if not rows:
        return pd.DataFrame(columns=["pathway", "overlap", "set_size", "enrichment_ratio", "pvalue", "padj", "genes"])
    df = pd.DataFrame(rows)
    df["padj"] = bh_fdr(df["pvalue"].to_numpy(dtype=float))
    return df.sort_values(["padj", "pvalue"])


def gsea_es_and_le(
    ranked_genes: list[str],
    ranked_scores: np.ndarray,
    gene_set: set[str],
    weight: float = 1.0,
) -> tuple[float, list[str]]:
    n = len(ranked_genes)
    in_set = np.array([g in gene_set for g in ranked_genes], dtype=bool)
    n_hit = int(in_set.sum())
    n_miss = n - n_hit
    if n_hit == 0 or n_miss == 0:
        return 0.0, []
    weights = np.abs(ranked_scores[in_set]) ** weight
    sum_w = float(weights.sum())
    if sum_w <= 0:
        weights = np.ones(n_hit, dtype=float)
        sum_w = float(n_hit)
    increments = np.full(n, -1.0 / n_miss, dtype=float)
    increments[in_set] = weights / sum_w
    running = np.cumsum(increments)
    max_idx = int(np.argmax(running))
    min_idx = int(np.argmin(running))
    if abs(running[max_idx]) >= abs(running[min_idx]):
        es = float(running[max_idx])
        le = [ranked_genes[i] for i in range(max_idx + 1) if in_set[i]]
    else:
        es = float(running[min_idx])
        le = [ranked_genes[i] for i in range(min_idx, n) if in_set[i]]
    return es, le


def run_gsea(
    ranked_genes: list[str],
    ranked_scores: np.ndarray,
    gmt: dict[str, list[str]],
    min_size: int,
    max_size: int,
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n = len(ranked_genes)
    universe = set(ranked_genes)
    abs_scores = np.abs(ranked_scores)

    rows = []
    for pathway, genes in gmt.items():
        gset = set(genes) & universe
        k = len(gset)
        if k < min_size or k > max_size:
            continue
        es, le = gsea_es_and_le(ranked_genes, ranked_scores, gset)
        null_es = np.empty(permutations, dtype=float)
        miss_step = -1.0 / (n - k)
        for p_idx in range(permutations):
            hit_idx = rng.choice(n, size=k, replace=False)
            w = abs_scores[hit_idx]
            sw = float(w.sum())
            if sw <= 0:
                w = np.ones(k, dtype=float)
                sw = float(k)
            inc = np.full(n, miss_step, dtype=float)
            inc[hit_idx] = w / sw
            run = np.cumsum(inc)
            hi, lo = float(run.max()), float(run.min())
            null_es[p_idx] = hi if abs(hi) >= abs(lo) else lo

        if es >= 0:
            pos_null = null_es[null_es >= 0]
            mean_pos = float(pos_null.mean()) if len(pos_null) else 1.0
            nes = es / max(mean_pos, 1e-6)
            pval = float((np.sum(pos_null >= es) + 1) / (len(pos_null) + 1)) if len(pos_null) else 1.0
        else:
            neg_null = null_es[null_es < 0]
            mean_neg = float(np.abs(neg_null).mean()) if len(neg_null) else 1.0
            nes = es / max(mean_neg, 1e-6)
            pval = float((np.sum(neg_null <= es) + 1) / (len(neg_null) + 1)) if len(neg_null) else 1.0

        rows.append(
            {
                "pathway": pathway,
                "set_size": k,
                "es": round(es, 4),
                "nes": round(nes, 4),
                "pvalue": pval,
                "leading_edge_size": len(le),
                "leading_edge": ";".join(le[:20]),
            }
        )
    if not rows:
        return pd.DataFrame(columns=["pathway", "set_size", "es", "nes", "pvalue", "padj", "leading_edge_size", "leading_edge"])
    df = pd.DataFrame(rows)
    df["padj"] = bh_fdr(df["pvalue"].to_numpy(dtype=float))
    return df.sort_values(["padj", "pvalue"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--de_table", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--mode", choices=["both", "ora", "gsea"], default="both")
    parser.add_argument("--direction", choices=["both", "up", "down"], default="both")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--log2fc", type=float, default=1.0)
    parser.add_argument("--min_set_size", type=int, default=5)
    parser.add_argument("--max_set_size", type=int, default=500)
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    de = read_table(args.de_table)
    id_col = "protein_id" if "protein_id" in de.columns else de.columns[0]
    de = de.copy()
    de["_gene"] = de[id_col].astype(str).map(normalize_gene_id)
    de = de.dropna(subset=["log2fc"]).drop_duplicates("_gene")
    if len(de) < 5:
        fail("DE table has fewer than 5 valid proteins")

    gmt = read_gmt(args.gmt)
    universe = set(de["_gene"])

    figures: list[str] = []
    artifacts: list[str] = []
    in_files = (args.de_table, args.gmt)
    params_meta = {"mode": args.mode, "direction": args.direction, "fdr": args.fdr}
    ora_summary = {"skipped": "mode=gsea"}
    gsea_summary = {"skipped": "mode=ora"}

    if args.mode in ("both", "ora"):
        padj_col = "padj" if "padj" in de.columns else "pvalue"
        sig = de[de[padj_col] <= args.fdr]
        if args.direction == "up":
            fg = set(sig[sig["log2fc"] >= args.log2fc]["_gene"])
        elif args.direction == "down":
            fg = set(sig[sig["log2fc"] <= -args.log2fc]["_gene"])
        else:
            fg = set(sig[sig["log2fc"].abs() >= args.log2fc]["_gene"])

        ora_df = run_ora(fg, universe, gmt, args.min_set_size, args.max_set_size)
        ora_path = os.path.join(args.output_dir, "ora_results.tsv")
        ora_df.to_csv(ora_path, sep="\t", index=False)
        artifacts.append(
            write_artifact(
                ora_path,
                "prot_enrichment",
                "Proteomics over-representation analysis (ORA) results",
                inputs=in_files,
                parameters=params_meta,
                shape=ora_df.shape,
            )
        )
        ora_summary = {
            "ora_table": os.path.abspath(ora_path),
            "universe_size": len(universe),
            "foreground_size": len(fg),
            "n_pathways_tested": int(len(ora_df)),
            "n_significant": int((ora_df["padj"] <= args.fdr).sum()) if len(ora_df) else 0,
            "top": [
                {
                    "pathway": str(r["pathway"]),
                    "overlap": int(r["overlap"]),
                    "set_size": int(r["set_size"]),
                    "padj": float(r["padj"]),
                }
                for _, r in ora_df.head(8).iterrows()
            ],
        }

    if args.mode in ("both", "gsea"):
        if "stat" in de.columns and de["stat"].notna().sum() >= 5:
            rank_series = de["stat"].astype(float)
        else:
            p = np.clip(de["pvalue"].astype(float), 1e-300, 1.0)
            rank_series = -np.log10(p) * np.sign(de["log2fc"].astype(float))
        de_rank = de.assign(_rank=rank_series).dropna(subset=["_rank"]).sort_values("_rank", ascending=False)
        gsea_df = run_gsea(
            list(de_rank["_gene"]),
            de_rank["_rank"].to_numpy(dtype=float),
            gmt,
            args.min_set_size,
            args.max_set_size,
            args.permutations,
            args.seed,
        )
        gsea_path = os.path.join(args.output_dir, "gsea_results.tsv")
        gsea_df.to_csv(gsea_path, sep="\t", index=False)
        artifacts.append(
            write_artifact(
                gsea_path,
                "prot_enrichment",
                "Proteomics pre-ranked GSEA results",
                inputs=in_files,
                parameters=params_meta,
                shape=gsea_df.shape,
            )
        )
        gsea_summary = {
            "gsea_table": os.path.abspath(gsea_path),
            "n_ranked_proteins": int(len(de_rank)),
            "n_pathways_tested": int(len(gsea_df)),
            "n_significant": int((gsea_df["padj"] <= args.fdr).sum()) if len(gsea_df) else 0,
            "top": [
                {
                    "pathway": str(r["pathway"]),
                    "nes": float(r["nes"]),
                    "padj": float(r["padj"]),
                }
                for _, r in gsea_df.head(8).iterrows()
            ],
        }
        if len(gsea_df) > 0:
            top_p = gsea_df.head(12).iloc[::-1]
            fig, ax = plt.subplots(figsize=(6.2, max(3.5, len(top_p) * 0.35)))
            colors = ["#EF4444" if v >= 0 else "#3B82F6" for v in top_p["nes"]]
            ax.barh(top_p["pathway"], top_p["nes"], color=colors)
            ax.axvline(0.0, color="#334155", linewidth=0.9)
            ax.set_xlabel("Normalized Enrichment Score (NES)")
            ax.set_title("Proteomics pre-ranked GSEA")
            fig.tight_layout()
            fig_path = os.path.join(args.output_dir, "enrichment-dotplot.png")
            fig.savefig(fig_path, dpi=150)
            plt.close(fig)
            figures.append(fig_path)
            artifacts.insert(
                0,
                write_artifact(
                    fig_path,
                    "prot_enrichment",
                    "Proteomics pre-ranked GSEA NES barplot",
                    inputs=in_files,
                    parameters=params_meta,
                ),
            )

    emit(
        {
            "de_table": os.path.abspath(args.de_table),
            "gmt": os.path.abspath(args.gmt),
            "mode": args.mode,
            "ora": ora_summary,
            "gsea": gsea_summary,
            "figures": [os.path.abspath(f) for f in figures],
            "artifacts": artifacts,
        }
    )


if __name__ == "__main__":
    main()
