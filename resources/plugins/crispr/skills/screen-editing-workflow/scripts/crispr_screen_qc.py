#!/usr/bin/env python3
"""Pooled CRISPR screen count matrix QC (Gini, dropout, DESeq size factors, PCA, and NNMD/ROC-AUC)."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import deseq_size_factors, emit, fail, pick_col, read_table, setup_matplotlib, write_artifact


def _gini(x) -> float:
    import numpy as np

    arr = np.sort(np.maximum(np.asarray(x, dtype=float), 0.0))
    n = arr.size
    total = float(arr.sum())
    if n == 0 or total <= 0:
        return 0.0
    idx = np.arange(1, n + 1, dtype=float)
    return float(np.sum((2.0 * idx - n - 1.0) * arr) / (n * total))


def _load_gene_set(path: str | Path) -> set[str]:
    p = Path(path)
    if not p.is_file():
        fail(f"Gene list file not found: {path}")
    genes: set[str] = set()
    with open(p, "rt", encoding="utf-8", errors="replace") as fh:
        for i, raw in enumerate(fh):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            tok = line.split("\t")[0].split(",")[0].strip()
            if i == 0 and tok.lower() in ("gene", "gene_id", "symbol", "genes"):
                continue
            if tok:
                genes.add(tok)
    return genes


def _pick_reference_sample(sample_cols: list[str], meta_path: str) -> str:
    if meta_path:
        mdf = read_table(meta_path).copy()
        sc = pick_col(list(mdf.columns), "", ["sample", "sample_id", "name"], required=False)
        if sc:
            for col_cand in ("timepoint", "condition", "group", "arm"):
                tc = pick_col(list(mdf.columns), "", [col_cand], required=False)
                if tc:
                    for _, row in mdf.iterrows():
                        s_name = str(row[sc]).strip()
                        val = str(row[tc]).strip().lower()
                        if s_name in sample_cols and any(
                            k in val for k in ("t0", "day0", "d0", "plasmid", "baseline", "input", "control", "ctrl")
                        ):
                            return s_name
    for s in sample_cols:
        sl = s.lower()
        if any(k in sl for k in ("t0", "day0", "d0", "plasmid", "baseline", "input")):
            return s
    return sample_cols[0]


def main() -> None:
    ap = argparse.ArgumentParser(description="Pooled CRISPR screen count matrix QC")
    ap.add_argument("--counts", required=True, help="sgRNA count matrix TSV/CSV")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--sgrna_col", default="")
    ap.add_argument("--gene_col", default="")
    ap.add_argument("--metadata", default="")
    ap.add_argument("--essential_genes", default="")
    ap.add_argument("--nonessential_genes", default="")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA
    from sklearn.metrics import roc_auc_score

    df = read_table(args.counts).copy()
    cols = list(df.columns)
    sg_c = pick_col(cols, args.sgrna_col, ["sgRNA", "sgrna", "guide_id", "id"])
    gn_c = pick_col(cols, args.gene_col, ["Gene", "gene", "gene_id", "target", "symbol"])

    sample_cols = [c for c in cols if c not in (sg_c, gn_c)]
    if not sample_cols:
        fail("No sample count columns detected in counts matrix.")

    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float, copy=True)
    if (mat < 0).any():
        fail("Negative read counts encountered in count matrix.")

    sf = deseq_size_factors(mat)
    norm_mat = mat / sf[None, :]
    log_norm = np.log2(norm_mat + 1.0)

    n_sgrnas = int(len(df))
    n_genes = int(df[gn_c].astype(str).nunique())
    warnings: list[str] = []

    # PCA & correlation
    if len(sample_cols) >= 2:
        corr_mat = np.corrcoef(log_norm.T)
        if np.ndim(corr_mat) == 0:
            corr_mat = np.array([[1.0]])
    else:
        corr_mat = np.array([[1.0]])

    if len(sample_cols) >= 2 and n_sgrnas >= 2:
        n_comp = min(2, len(sample_cols), n_sgrnas)
        coords = PCA(n_components=n_comp).fit_transform(log_norm.T)
    else:
        coords = np.zeros((len(sample_cols), 2), dtype=float)

    sample_rows: list[dict] = []
    for i, s in enumerate(sample_cols):
        col_vec = mat[:, i]
        total_reads = int(round(float(col_vec.sum())))
        med_reads = float(np.median(col_vec))
        mean_reads = float(np.mean(col_vec))
        zero_cnt = int((col_vec == 0).sum())
        dropout_pct = round(100.0 * zero_cnt / max(1, n_sgrnas), 3)
        gini = round(_gini(col_vec), 4)
        pc1 = round(float(coords[i, 0]), 4)
        pc2 = round(float(coords[i, 1]), 4) if coords.shape[1] > 1 else 0.0

        flags: list[str] = []
        if med_reads < 100.0:
            flags.append("low_coverage_lt100")
        elif med_reads < 300.0:
            flags.append("caution_coverage_lt300")
        if gini > 0.30:
            flags.append("severe_bottlenecking_gini_gt0.30")
        elif gini > 0.25:
            flags.append("bottlenecking_gini_gt0.25")
        elif gini > 0.15:
            flags.append("elevated_gini_gt0.15")

        if flags:
            warnings.append(f"Sample '{s}': {', '.join(flags)} (median={med_reads:.1f}, gini={gini:.3f}).")

        sample_rows.append(
            {
                "sample": s,
                "total_reads": total_reads,
                "median_reads_per_sgrna": round(med_reads, 2),
                "mean_reads_per_sgrna": round(mean_reads, 2),
                "zero_count_sgrnas": zero_cnt,
                "zero_dropout_pct": dropout_pct,
                "gini_index": gini,
                "deseq_size_factor": round(float(sf[i]), 4),
                "pc1": pc1,
                "pc2": pc2,
                "qc_flags": ";".join(flags) if flags else "pass",
            }
        )

    ref_sample = _pick_reference_sample(sample_cols, args.metadata)
    ref_idx = sample_cols.index(ref_sample)
    nnmd_rows: list[dict] = []

    if args.essential_genes and args.nonessential_genes:
        ess_set = _load_gene_set(args.essential_genes)
        noness_set = _load_gene_set(args.nonessential_genes)
        genes_series = df[gn_c].astype(str).str.strip()
        for i, s in enumerate(sample_cols):
            if i == ref_idx and len(sample_cols) > 1:
                continue
            lfc = np.log2((norm_mat[:, i] + 1.0) / (norm_mat[:, ref_idx] + 1.0))
            gdf = pd.DataFrame({"gene": genes_series, "lfc": lfc}).groupby("gene", as_index=False)["lfc"].mean()
            ess_vals = gdf[gdf["gene"].isin(ess_set)]["lfc"].to_numpy(dtype=float, copy=True)
            non_vals = gdf[gdf["gene"].isin(noness_set)]["lfc"].to_numpy(dtype=float, copy=True)
            if ess_vals.size >= 3 and non_vals.size >= 3:
                sd_non = float(np.std(non_vals, ddof=1))
                nnmd = float((np.mean(ess_vals) - np.mean(non_vals)) / max(sd_non, 1e-6))
                y_true = np.concatenate([np.ones(ess_vals.size), np.zeros(non_vals.size)])
                scores = np.concatenate([-ess_vals, -non_vals])
                auc = float(roc_auc_score(y_true, scores))
                nnmd_rows.append(
                    {
                        "sample": s,
                        "reference_sample": ref_sample,
                        "n_essential_genes": int(ess_vals.size),
                        "n_nonessential_genes": int(non_vals.size),
                        "mean_ess_lfc": round(float(np.mean(ess_vals)), 4),
                        "mean_noness_lfc": round(float(np.mean(non_vals)), 4),
                        "nnmd": round(nnmd, 4),
                        "essentiality_roc_auc": round(auc, 4),
                        "strong_separation": bool(nnmd < -1.25),
                    }
                )
                if nnmd >= -1.25:
                    warnings.append(
                        f"Weak essentiality separation in '{s}' vs '{ref_sample}' (NNMD={nnmd:.2f} >= -1.25, AUC={auc:.3f})."
                    )

    qc_df = pd.DataFrame(sample_rows)
    if nnmd_rows:
        nnmd_df = pd.DataFrame(nnmd_rows)[["sample", "nnmd", "essentiality_roc_auc"]]
        qc_df = qc_df.merge(nnmd_df, on="sample", how="left")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / "screen_sample_qc.tsv"
    fig_path = out_dir / "screen_qc_dashboard.png"
    qc_df.to_csv(tsv_path, sep="\t", index=False)

    write_artifact(
        tsv_path,
        "crispr_screen_qc",
        "Pooled CRISPR Screen Sample QC & Representation Metrics",
        inputs=[args.counts, args.metadata, args.essential_genes, args.nonessential_genes],
        parameters={"reference_sample": ref_sample},
        shape=(int(len(qc_df)), int(qc_df.shape[1])),
    )

    # 3-panel QC dashboard
    plt = setup_matplotlib()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14.8, 4.5), dpi=180)
    palette = ["#0d9488", "#2563eb", "#d97706", "#7c3aed", "#dc2626", "#059669", "#db2777", "#4f46e5"]

    # Panel A: Lorenz curves
    for i, s in enumerate(sample_cols):
        vals = np.sort(np.maximum(mat[:, i], 0.0))
        tot = float(vals.sum())
        if tot > 0:
            cum_y = np.concatenate([[0.0], np.cumsum(vals) / tot])
        else:
            cum_y = np.linspace(0.0, 1.0, len(vals) + 1)
        cum_x = np.linspace(0.0, 1.0, len(cum_y))
        ax1.plot(cum_x, cum_y, lw=1.7, color=palette[i % len(palette)], label=f"{s} (G={sample_rows[i]['gini_index']:.3f})")
    ax1.plot([0, 1], [0, 1], "--", color="#64748b", lw=1.0, label="Uniform")
    ax1.set_xlabel("Cumulative Fraction of sgRNAs")
    ax1.set_ylabel("Cumulative Fraction of Reads")
    ax1.set_title("Lorenz Curve & Gini Skew")
    ax1.legend(frameon=False, fontsize=7, loc="upper left")

    # Panel B: Zero-dropout % & median coverage
    x_pos = np.arange(len(sample_cols))
    dropouts = [r["zero_dropout_pct"] for r in sample_rows]
    medians = [r["median_reads_per_sgrna"] for r in sample_rows]
    ax2.bar(x_pos, dropouts, color="#f59e0b", width=0.55)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(sample_cols, rotation=25, ha="right", fontsize=8)
    ax2.set_ylabel("Zero-Count Dropout (%)")
    ax2.set_title("Dropout & Median sgRNA Depth")
    ax2b = ax2.twinx()
    ax2b.plot(x_pos, medians, "o-", color="#0f766e", lw=1.5, ms=5)
    ax2b.axhline(100.0, color="#dc2626", ls=":", lw=0.9)
    ax2b.set_ylabel("Median Reads / sgRNA", color="#0f766e")

    # Panel C: Sample Correlation Heatmap or PCA
    if len(sample_cols) >= 2:
        im = ax3.imshow(corr_mat, vmin=max(0.0, float(np.min(corr_mat)) - 0.05), vmax=1.0, cmap="YlGnBu")
        ax3.set_xticks(range(len(sample_cols)))
        ax3.set_yticks(range(len(sample_cols)))
        ax3.set_xticklabels(sample_cols, rotation=30, ha="right", fontsize=7.5)
        ax3.set_yticklabels(sample_cols, fontsize=7.5)
        for r_i in range(len(sample_cols)):
            for c_i in range(len(sample_cols)):
                ax3.text(c_i, r_i, f"{corr_mat[r_i, c_i]:.2f}", ha="center", va="center", fontsize=7, color="#0f172a")
        fig.colorbar(im, ax=ax3, fraction=0.046, pad=0.04, label="Pearson r")
        ax3.set_title("log2(Norm+1) Sample Correlation")
    else:
        ax3.text(0.5, 0.5, "Single sample\n(PCA/Corr skipped)", ha="center", va="center")
        ax3.axis("off")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "crispr_screen_qc",
        "Pooled CRISPR Screen QC Dashboard",
        inputs=[args.counts],
        parameters={"reference_sample": ref_sample},
    )

    emit(
        {
            "status": "ok",
            "tool": "crispr_screen_qc",
            "n_sgrnas": n_sgrnas,
            "n_genes": n_genes,
            "n_samples": len(sample_cols),
            "reference_sample": ref_sample,
            "samples": sample_rows,
            "nnmd_summaries": nnmd_rows,
            "outputs": {
                "screen_sample_qc": str(tsv_path),
                "qc_dashboard": str(fig_path),
            },
            "artifacts": [str(tsv_path), str(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
