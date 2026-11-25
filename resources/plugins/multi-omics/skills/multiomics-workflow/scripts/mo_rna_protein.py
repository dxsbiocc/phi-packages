#!/usr/bin/env python3
"""mo_rna_protein: Transcriptome-Proteome 9-quadrant differential concordance
(buffering vs post-transcriptional regulation) and/or across-sample gene-wise
mRNA-protein Spearman correlation with BH FDR.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import bh_fdr, emit, fail, read_matrix, read_table, write_artifact

FC_CANDS = ("log2foldchange", "log2fc", "logfc", "lfc", "diff")
FDR_CANDS = ("padj", "fdr", "adj.p.val", "qvalue", "q_val", "p_val_adj")
P_CANDS = ("pvalue", "pval", "p.value", "p_val")
GENE_CANDS = ("gene_symbol", "genes", "gene", "symbol", "gene_name", "hugo_symbol", "gene_id", "protein", "id")
VERSION_RE = re.compile(r"^(ENS[A-Z]*G\d+)\.\d+$", re.IGNORECASE)


def clean_gene(g: str, gene_map: dict[str, str] | None) -> str:
    s = str(g).strip()
    if ";" in s:
        s = s.split(";", 1)[0].strip()
    m = VERSION_RE.match(s)
    if m:
        s = m.group(1)
    if gene_map and s in gene_map:
        return gene_map[s]
    return s


def parse_de_standardized(path: str, gene_map: dict[str, str] | None, label: str) -> pd.DataFrame:
    df = read_table(path)
    lower = {str(c).strip().lower(): c for c in df.columns}
    gcol = next((lower[c] for c in GENE_CANDS if c in lower), df.columns[0])
    fccol = next((lower[c] for c in FC_CANDS if c in lower), None)
    fdrcol = next((lower[c] for c in FDR_CANDS if c in lower), None)
    if fdrcol is None:
        fdrcol = next((lower[c] for c in P_CANDS if c in lower), None)
    if fccol is None or fdrcol is None:
        fail(f"{label} DE table '{path}' must contain log2FC and padj/FDR/pvalue columns; found {list(df.columns)}")

    out = pd.DataFrame(
        {
            "gene": [clean_gene(v, gene_map) for v in df[gcol]],
            "log2fc": pd.to_numeric(df[fccol], errors="coerce"),
            "fdr": pd.to_numeric(df[fdrcol], errors="coerce"),
        }
    )
    out = out.dropna(subset=["gene", "log2fc", "fdr"])
    out = out[out["gene"] != ""]
    # Keep smallest FDR on duplicate symbols
    out = out.sort_values("fdr", ascending=True).drop_duplicates(subset=["gene"], keep="first").set_index("gene")
    return out


def classify_quadrant(rna_fc: float, rna_q: float, prot_fc: float, prot_q: float, r_cut: float, p_cut: float, q_cut: float) -> str:
    r_sig = 1 if (rna_q < q_cut and rna_fc >= r_cut) else (-1 if (rna_q < q_cut and rna_fc <= -r_cut) else 0)
    p_sig = 1 if (prot_q < q_cut and prot_fc >= p_cut) else (-1 if (prot_q < q_cut and prot_fc <= -p_cut) else 0)
    if r_sig == 1 and p_sig == 1:
        return "concordant_up"
    if r_sig == -1 and p_sig == -1:
        return "concordant_down"
    if r_sig == 1 and p_sig == 0:
        return "rna_only_up_buffered"
    if r_sig == -1 and p_sig == 0:
        return "rna_only_down_buffered"
    if r_sig == 0 and p_sig == 1:
        return "protein_only_up_post_tx"
    if r_sig == 0 and p_sig == -1:
        return "protein_only_down_post_tx"
    if r_sig == 1 and p_sig == -1:
        return "discordant_rna_up_prot_down"
    if r_sig == -1 and p_sig == 1:
        return "discordant_rna_down_prot_up"
    return "unchanged"


QUAD_COLORS = {
    "concordant_up": "#E11D48",
    "concordant_down": "#2563EB",
    "rna_only_up_buffered": "#F59E0B",
    "rna_only_down_buffered": "#38BDF8",
    "protein_only_up_post_tx": "#10B981",
    "protein_only_down_post_tx": "#14B8A6",
    "discordant_rna_up_prot_down": "#8B5CF6",
    "discordant_rna_down_prot_up": "#A855F7",
    "unchanged": "#CBD5E1",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna_de", default=None)
    parser.add_argument("--prot_de", default=None)
    parser.add_argument("--rna_matrix", default=None)
    parser.add_argument("--prot_matrix", default=None)
    parser.add_argument("--gene_map", default=None)
    parser.add_argument("--rna_fc_cutoff", type=float, default=1.0)
    parser.add_argument("--prot_fc_cutoff", type=float, default=0.585)
    parser.add_argument("--fdr_cutoff", type=float, default=0.05)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    has_de = bool(args.rna_de and args.prot_de)
    has_mat = bool(args.rna_matrix and args.prot_matrix)
    if not has_de and not has_mat:
        fail("provide (--rna_de + --prot_de) and/or (--rna_matrix + --prot_matrix)")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    gmap: dict[str, str] | None = None
    if args.gene_map:
        mdf = read_table(args.gene_map)
        if mdf.shape[1] < 2:
            fail("gene_map must have at least 2 columns (source_id, gene_symbol)")
        gmap = {
            str(k).strip(): str(v).strip()
            for k, v in zip(mdf.iloc[:, 0], mdf.iloc[:, 1])
            if str(k).strip() and str(v).strip()
        }

    in_files = [args.rna_de, args.prot_de, args.rna_matrix, args.prot_matrix, args.gene_map]
    params_meta = {
        "rna_fc_cutoff": args.rna_fc_cutoff,
        "prot_fc_cutoff": args.prot_fc_cutoff,
        "fdr_cutoff": args.fdr_cutoff,
    }

    de_summary = None
    sample_summary = None
    outputs: dict[str, str | None] = {"nine_quadrant_table": None, "genewise_correlation_table": None}
    figures: list[str] = []
    artifacts: list[str] = []
    sci_ctx: dict[str, dict] = {}
    warnings: list[str] = []

    if has_de:
        r_df = parse_de_standardized(args.rna_de, None, "RNA")
        p_df = parse_de_standardized(args.prot_de, gmap, "Protein")
        shared_genes = sorted(set(r_df.index) & set(p_df.index))
        if len(shared_genes) < 5:
            fail(
                f"only {len(shared_genes)} shared gene IDs between RNA DE and Protein DE; "
                "pass --gene_map or run rna_id_map / prot_inspect to align gene symbols"
            )

        merged = pd.DataFrame(
            {
                "gene": shared_genes,
                "rna_log2fc": r_df.loc[shared_genes, "log2fc"].to_numpy(dtype=float),
                "rna_fdr": r_df.loc[shared_genes, "fdr"].to_numpy(dtype=float),
                "prot_log2fc": p_df.loc[shared_genes, "log2fc"].to_numpy(dtype=float),
                "prot_fdr": p_df.loc[shared_genes, "fdr"].to_numpy(dtype=float),
            }
        )
        merged["quadrant"] = [
            classify_quadrant(
                rf, rq, pf, pq, args.rna_fc_cutoff, args.prot_fc_cutoff, args.fdr_cutoff
            )
            for rf, rq, pf, pq in zip(
                merged["rna_log2fc"], merged["rna_fdr"], merged["prot_log2fc"], merged["prot_fdr"]
            )
        ]

        sp_r, sp_p = stats.spearmanr(merged["rna_log2fc"], merged["prot_log2fc"])
        pe_r, pe_p = stats.pearsonr(merged["rna_log2fc"], merged["prot_log2fc"])

        q_table_path = outdir / "rna_protein_9quadrant.tsv"
        merged.to_csv(q_table_path, sep="\t", index=False)
        outputs["nine_quadrant_table"] = str(q_table_path)
        sci_ctx[str(q_table_path)] = {"valueSemantics": "de_stats"}

        # Plot 9-quadrant scatter
        fig, ax = plt.subplots(figsize=(6.8, 6.0))
        unch = merged[merged["quadrant"] == "unchanged"]
        if len(unch) > 0:
            ax.scatter(
                unch["rna_log2fc"], unch["prot_log2fc"], s=14, alpha=0.35, color=QUAD_COLORS["unchanged"], label=f"unchanged ({len(unch)})"
            )
        for qname, color in QUAD_COLORS.items():
            if qname == "unchanged":
                continue
            sub_q = merged[merged["quadrant"] == qname]
            if len(sub_q) > 0:
                ax.scatter(
                    sub_q["rna_log2fc"], sub_q["prot_log2fc"], s=26, alpha=0.85, color=color, label=f"{qname} ({len(sub_q)})"
                )
        ax.axvline(args.rna_fc_cutoff, color="#64748B", linestyle="--", linewidth=0.9)
        ax.axvline(-args.rna_fc_cutoff, color="#64748B", linestyle="--", linewidth=0.9)
        ax.axhline(args.prot_fc_cutoff, color="#64748B", linestyle="--", linewidth=0.9)
        ax.axhline(-args.prot_fc_cutoff, color="#64748B", linestyle="--", linewidth=0.9)
        ax.set_xlabel("RNA log2 Fold Change")
        ax.set_ylabel("Protein log2 Fold Change")
        ax.set_title(
            f"RNA–Protein 9-Quadrant Concordance (n={len(shared_genes)}, ρ={float(sp_r):.3f}, p={float(sp_p):.2e})",
            fontsize=10,
        )
        ax.legend(frameon=False, fontsize=7.5, loc="best")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        q_fig = outdir / "nine-quadrant.png"
        fig.savefig(q_fig, dpi=160)
        plt.close(fig)

        figures.append(str(q_fig))
        artifacts.append(write_artifact(q_fig, "mo_rna_protein", "RNA-Protein 9-quadrant scatter plot", inputs=in_files, parameters=params_meta))
        artifacts.append(write_artifact(q_table_path, "mo_rna_protein", "RNA-Protein 9-quadrant regulatory table", inputs=in_files, parameters=params_meta, shape=merged.shape))

        q_counts = {str(k): int(v) for k, v in merged["quadrant"].value_counts().items()}
        de_summary = {
            "n_shared_genes": len(shared_genes),
            "spearman_r": round(float(sp_r), 4),
            "spearman_pvalue": float(sp_p),
            "pearson_r": round(float(pe_r), 4),
            "pearson_pvalue": float(pe_p),
            "quadrant_counts": q_counts,
        }

    if has_mat:
        r_mat = read_matrix(args.rna_matrix)
        r_mat.index = [clean_gene(g, None) for g in r_mat.index]
        r_mat = r_mat.groupby(level=0).mean()

        p_mat = read_matrix(args.prot_matrix)
        p_mat.index = [clean_gene(g, gmap) for g in p_mat.index]
        p_mat = p_mat.groupby(level=0).mean()

        shared_samples = [s for s in r_mat.columns if s in set(p_mat.columns)]
        shared_m_genes = sorted(set(r_mat.index) & set(p_mat.index))
        if len(shared_samples) < 5:
            fail(f"across-sample correlation requires >= 5 matched samples, found {len(shared_samples)}")
        if len(shared_m_genes) < 3:
            fail(f"only {len(shared_m_genes)} shared genes between RNA and Protein matrices")

        r_arr = r_mat.loc[shared_m_genes, shared_samples].to_numpy(dtype=float)
        p_arr = p_mat.loc[shared_m_genes, shared_samples].to_numpy(dtype=float)

        gw_rows = []
        for idx_g, gname in enumerate(shared_m_genes):
            rv = r_arr[idx_g]
            pv = p_arr[idx_g]
            ok = np.isfinite(rv) & np.isfinite(pv)
            n_ok = int(ok.sum())
            if n_ok < 5 or np.std(rv[ok]) < 1e-12 or np.std(pv[ok]) < 1e-12:
                continue
            rho, pval = stats.spearmanr(rv[ok], pv[ok])
            gw_rows.append(
                {
                    "gene": gname,
                    "n_samples": n_ok,
                    "spearman_r": round(float(rho), 4),
                    "pvalue": float(pval),
                }
            )
        if not gw_rows:
            warnings.append("0 genes had >=5 valid non-constant paired observations for across-sample correlation")
        else:
            gw_df = pd.DataFrame(gw_rows)
            gw_df["fdr"] = bh_fdr(gw_df["pvalue"].to_numpy(dtype=float))
            gw_df = gw_df.sort_values("pvalue", ascending=True)
            gw_path = outdir / "genewise_mrna_protein_correlation.tsv"
            gw_df.to_csv(gw_path, sep="\t", index=False)
            outputs["genewise_correlation_table"] = str(gw_path)

            fig2, ax2 = plt.subplots(figsize=(6.0, 4.2))
            ax2.hist(gw_df["spearman_r"], bins=25, color="#0EA5E9", alpha=0.8, edgecolor="white")
            med_r = float(gw_df["spearman_r"].median())
            ax2.axvline(med_r, color="#E11D48", linestyle="--", linewidth=1.5, label=f"Median ρ = {med_r:.3f}")
            ax2.axvline(0.0, color="#64748B", linestyle=":", linewidth=1.0)
            ax2.set_xlabel("Across-sample mRNA–Protein Spearman ρ")
            ax2.set_ylabel("Number of genes")
            ax2.set_title(f"Gene-wise mRNA–Protein Correlation ({len(shared_samples)} samples)", fontsize=10)
            ax2.legend(frameon=False)
            ax2.spines["top"].set_visible(False)
            ax2.spines["right"].set_visible(False)
            fig2.tight_layout()
            gw_fig = outdir / "genewise-cor-hist.png"
            fig2.savefig(gw_fig, dpi=160)
            plt.close(fig2)

            figures.append(str(gw_fig))
            artifacts.insert(1 if has_de else 0, write_artifact(gw_fig, "mo_rna_protein", "Across-sample gene-wise mRNA-protein correlation distribution", inputs=in_files, parameters=params_meta))
            artifacts.append(write_artifact(gw_path, "mo_rna_protein", "Across-sample gene-wise mRNA-protein correlation table", inputs=in_files, parameters=params_meta, shape=gw_df.shape))

            sample_summary = {
                "n_shared_samples": len(shared_samples),
                "n_genes_tested": int(len(gw_df)),
                "median_spearman_r": round(med_r, 4),
                "n_positive_fdr05": int(((gw_df["fdr"] < args.fdr_cutoff) & (gw_df["spearman_r"] > 0)).sum()),
                "n_negative_fdr05": int(((gw_df["fdr"] < args.fdr_cutoff) & (gw_df["spearman_r"] < 0)).sum()),
            }

    emit(
        {
            "tool": "mo_rna_protein",
            "de_concordance": de_summary,
            "sample_level_correlation": sample_summary,
            "outputs": outputs,
            "figures": figures,
            **({"scientificContext": sci_ctx} if sci_ctx else {}),
            "artifacts": artifacts[:8],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
