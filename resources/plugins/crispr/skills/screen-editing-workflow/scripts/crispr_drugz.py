#!/usr/bin/env python3
"""crispr_drugz: Chemogenetic Drug-vs-Control pooled CRISPR screen analysis
(Collimore & Hart 2019 DrugZ algorithm): RPM normalization, abundance-ordered
empirical-Bayes sliding-window variance estimation, gene-level sumZ / normZ,
and BH FDR for synthetic lethality (sensitizers) & drug resistance (suppressors).
"""
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
    parse_col_list,
    pick_col,
    read_table,
    setup_matplotlib,
    write_artifact,
)


def _eb_sliding_sd(ctrl_abund, delta_fc, window_size: int) -> "np.ndarray":
    import numpy as np
    import pandas as pd

    n = len(delta_fc)
    win = max(25, min(window_size, max(25, n // 3)))
    order = np.argsort(-ctrl_abund, kind="mergesort")
    fc_sorted = delta_fc[order]

    # Robust rolling MAD-based variance + empirical-Bayes shrinkage
    s_ser = pd.Series(fc_sorted)
    roll_med = s_ser.rolling(window=win, center=True, min_periods=max(5, win // 4)).median().to_numpy()
    abs_dev = np.abs(fc_sorted - roll_med)
    roll_mad = (
        pd.Series(abs_dev)
        .rolling(window=win, center=True, min_periods=max(5, win // 4))
        .median()
        .bfill()
        .ffill()
        .to_numpy()
    )
    global_sd = max(float(np.median(abs_dev)) * 1.4826, 0.08)
    eb_sd_sorted = np.maximum(0.05, 0.75 * (roll_mad * 1.4826) + 0.25 * global_sd)
    # Enforce monotonic non-decreasing variance from high to low control read counts (Hart DrugZ prior)
    eb_sd_sorted = np.maximum.accumulate(eb_sd_sorted)

    eb_sd = np.empty(n, dtype=float)
    eb_sd[order] = eb_sd_sorted
    return eb_sd


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--drug_cols", action="append", required=True)
    ap.add_argument("--control_cols", action="append", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--sgrna_col", default="")
    ap.add_argument("--gene_col", default="")
    ap.add_argument("--unpaired", action="store_true", default=False)
    ap.add_argument("--pseudocount", type=float, default=5.0)
    ap.add_argument("--window_size", type=int, default=500)
    ap.add_argument("--remove_genes", action="append", default=None)
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from scipy import stats

    drug_cols = parse_col_list(args.drug_cols)
    ctrl_cols = parse_col_list(args.control_cols)
    rem_genes = set(parse_col_list(args.remove_genes))
    if not drug_cols or not ctrl_cols:
        fail("Both --drug_cols and --control_cols must be non-empty.")
    if args.pseudocount <= 0:
        fail("--pseudocount must be > 0.")

    df = read_table(args.counts).copy()
    cols = list(df.columns)
    sg_c = pick_col(cols, args.sgrna_col, ["sgRNA", "sgrna", "guide_id", "id"])
    gn_c = pick_col(cols, args.gene_col, ["Gene", "gene", "gene_id", "target", "symbol"])
    for c in drug_cols + ctrl_cols:
        if c not in cols:
            fail(f"Sample column '{c}' not found in count matrix columns: {cols}")

    warnings: list[str] = []
    if rem_genes:
        before_n = len(df)
        df = df[~df[gn_c].astype(str).str.strip().isin(rem_genes)].reset_index(drop=True)
        warnings.append(f"Excluded {before_n - len(df)} sgRNAs targeting {len(rem_genes)} specified core-essential gene(s).")

    if len(df) < 30:
        fail("Count matrix requires >= 30 sgRNAs for empirical-Bayes DrugZ variance estimation.")

    use_paired = (not args.unpaired) and (len(drug_cols) == len(ctrl_cols))
    if (not args.unpaired) and (len(drug_cols) != len(ctrl_cols)):
        warnings.append(
            f"Unequal drug ({len(drug_cols)}) and control ({len(ctrl_cols)}) columns; switching automatically to unpaired mode."
        )

    # Normalize each sample to 1e7 total reads (DrugZ standard 10M RPM normalization)
    all_cols = ctrl_cols + drug_cols
    raw_mat = df[all_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float, copy=True)
    col_sums = np.maximum(raw_mat.sum(axis=0), 1.0)
    rpm_mat = raw_mat * (1e7 / col_sums[None, :])

    n_ctrl = len(ctrl_cols)
    ctrl_rpm = rpm_mat[:, :n_ctrl]
    drug_rpm = rpm_mat[:, n_ctrl:]
    p0 = float(args.pseudocount)

    if use_paired:
        mode = "paired"
        n_rep = len(drug_cols)
        z_rep_mat = np.zeros((len(df), n_rep), dtype=float)
        fc_rep_mat = np.zeros((len(df), n_rep), dtype=float)
        eb_rep_mat = np.zeros((len(df), n_rep), dtype=float)
        for r in range(n_rep):
            fc_r = np.log2((drug_rpm[:, r] + p0) / (ctrl_rpm[:, r] + p0))
            fc_r = fc_r - float(np.median(fc_r))  # center per-replicate median LFC at 0
            eb_sd = _eb_sliding_sd(ctrl_rpm[:, r], fc_r, args.window_size)
            fc_rep_mat[:, r] = fc_r
            eb_rep_mat[:, r] = eb_sd
            z_rep_mat[:, r] = fc_r / eb_sd
        mean_fc = np.mean(fc_rep_mat, axis=1)
        mean_eb_sd = np.mean(eb_rep_mat, axis=1)
        mean_ctrl = np.mean(ctrl_rpm, axis=1)
    else:
        mode = "unpaired"
        n_rep = 1
        mean_ctrl = np.mean(ctrl_rpm, axis=1)
        mean_drug = np.mean(drug_rpm, axis=1)
        fc = np.log2((mean_drug + p0) / (mean_ctrl + p0))
        fc = fc - float(np.median(fc))
        mean_eb_sd = _eb_sliding_sd(mean_ctrl, fc, args.window_size)
        mean_fc = fc
        z_rep_mat = (fc / mean_eb_sd)[:, None]

    sg_df = pd.DataFrame(
        {
            "sgRNA": df[sg_c].astype(str).str.strip(),
            "gene": df[gn_c].astype(str).str.strip(),
            "ctrl_mean_rpm": np.round(mean_ctrl, 2),
            "centered_log2fc": np.round(mean_fc, 4),
            "eb_std": np.round(mean_eb_sd, 4),
            "mean_z": np.round(np.mean(z_rep_mat, axis=1), 4),
        }
    )

    gene_rows: list[dict] = []
    for gene, idx in sg_df.groupby("gene", sort=False).groups.items():
        idx_arr = np.asarray(list(idx), dtype=int)
        z_sub = z_rep_mat[idx_arr]
        n_obs = int(z_sub.size)
        sum_z = float(np.sum(z_sub))
        raw_norm_z = sum_z / math.sqrt(max(1, n_obs))
        gene_rows.append(
            {
                "gene": gene,
                "n_sgrnas": int(idx_arr.size),
                "n_observations": n_obs,
                "mean_log2fc": round(float(np.mean(mean_fc[idx_arr])), 4),
                "sumZ": round(sum_z, 4),
                "_raw_normZ": raw_norm_z,
            }
        )

    gene_df = pd.DataFrame(gene_rows)
    raw_nz = gene_df["_raw_normZ"].to_numpy(dtype=float, copy=True)
    # DrugZ empirical null standardization via median and MAD
    med_nz = float(np.median(raw_nz))
    mad_nz = max(float(np.median(np.abs(raw_nz - med_nz))) * 1.4826, 0.5)
    norm_z = (raw_nz - med_nz) / mad_nz

    p_synth = stats.norm.cdf(norm_z)  # negative normZ = drug sensitization / synthetic lethality
    p_supp = stats.norm.sf(norm_z)   # positive normZ = drug resistance / suppressor

    gene_df["normZ"] = np.round(norm_z, 4)
    gene_df["p_synth"] = p_synth
    gene_df["fdr_synth"] = bh_fdr(p_synth)
    gene_df["rank_synth"] = gene_df["normZ"].rank(method="min", ascending=True).astype(int)
    gene_df["p_supp"] = p_supp
    gene_df["fdr_supp"] = bh_fdr(p_supp)
    gene_df["rank_supp"] = gene_df["normZ"].rank(method="min", ascending=False).astype(int)
    gene_df = gene_df.drop(columns=["_raw_normZ"]).sort_values("normZ", ascending=True).reset_index(drop=True)

    n_syn = int((gene_df["fdr_synth"] < 0.05).sum())
    n_sup = int((gene_df["fdr_supp"] < 0.05).sum())
    if n_syn == 0 and n_sup == 0:
        warnings.append("No synthetic-lethal or suppressor genes reached BH FDR < 0.05.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sg_path = out_dir / "drugz_sgrna_zscores.tsv"
    gn_path = out_dir / "drugz_gene_summary.tsv"
    fig_path = out_dir / "drugz_rank_waterfall.png"

    sg_df.to_csv(sg_path, sep="\t", index=False)
    gene_df.to_csv(gn_path, sep="\t", index=False)

    write_artifact(
        sg_path,
        "crispr_drugz",
        "DrugZ sgRNA-Level Empirical-Bayes Z-Scores",
        inputs=[args.counts],
        parameters={"mode": mode, "pseudocount": p0, "window_size": args.window_size},
        shape=(int(len(sg_df)), int(sg_df.shape[1])),
    )
    write_artifact(
        gn_path,
        "crispr_drugz",
        "DrugZ Gene-Level normZ Synthetic Lethality & Resistance Summary",
        inputs=[args.counts],
        parameters={"mode": mode, "pseudocount": p0},
        shape=(int(len(gene_df)), int(gene_df.shape[1])),
    )

    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.8), dpi=180)
    ranks = np.arange(1, len(gene_df) + 1)
    nz_vals = gene_df["normZ"].to_numpy()
    cols_pt = np.where(gene_df["fdr_synth"] < 0.05, "#2563eb", np.where(gene_df["fdr_supp"] < 0.05, "#dc2626", "#cbd5e1"))
    ax1.scatter(ranks, nz_vals, c=cols_pt, s=16, alpha=0.88)
    ax1.axhline(0.0, color="#64748b", ls="--", lw=0.9)
    for _, row in gene_df[gene_df["fdr_synth"] < 0.05].head(5).iterrows():
        ax1.annotate(row["gene"], (row["rank_synth"], row["normZ"]), fontsize=7.5, xytext=(6, -2), textcoords="offset points", color="#1e3a8a")
    for _, row in gene_df[gene_df["fdr_supp"] < 0.05].tail(5).iterrows():
        r_idx = int(len(gene_df) - row["rank_supp"] + 1)
        ax1.annotate(row["gene"], (r_idx, row["normZ"]), fontsize=7.5, xytext=(-36, 2), textcoords="offset points", color="#991b1b")
    ax1.set_xlabel("Gene Rank (Sensitizers → Suppressors)")
    ax1.set_ylabel("DrugZ Normalized Z-Score (normZ)")
    ax1.set_title(f"DrugZ Chemogenetic Waterfall (SL={n_syn}, Res={n_sup})")

    ax2.scatter(np.log10(mean_ctrl + 1.0), np.abs(mean_fc), s=8, c="#94a3b8", alpha=0.4, label="|log2FC|")
    ord_c = np.argsort(mean_ctrl)
    ax2.plot(np.log10(mean_ctrl[ord_c] + 1.0), mean_eb_sd[ord_c], color="#0d9488", lw=2.0, label="EB smoothed SD")
    ax2.set_xlabel("log10(Control RPM + 1)")
    ax2.set_ylabel("Dispersion / |log2FC|")
    ax2.set_title("Empirical-Bayes Sliding-Window Variance")
    ax2.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "crispr_drugz",
        "DrugZ Chemogenetic Rank Waterfall & EB Variance Profile",
        inputs=[args.counts],
        parameters={"mode": mode},
    )

    top_sl = gene_df.sort_values("normZ", ascending=True).head(15).to_dict(orient="records")
    top_sup = gene_df.sort_values("normZ", ascending=False).head(15).to_dict(orient="records")

    emit(
        {
            "status": "ok",
            "tool": "crispr_drugz",
            "mode": mode,
            "pseudocount": p0,
            "window_size": args.window_size,
            "n_sgrnas": int(len(sg_df)),
            "n_genes": int(len(gene_df)),
            "n_replicates": n_rep,
            "n_synthetic_lethal_fdr05": n_syn,
            "n_suppressor_fdr05": n_sup,
            "top_synthetic_lethal": top_sl,
            "top_suppressors": top_sup,
            "outputs": {
                "drugz_sgrna_zscores": str(sg_path),
                "drugz_gene_summary": str(gn_path),
                "drugz_waterfall_plot": str(fig_path),
            },
            "artifacts": [str(sg_path), str(gn_path), str(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
