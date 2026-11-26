#!/usr/bin/env python3
"""Gene-level positive & negative selection hit calling via alpha-RRA with permutation calibration & BH FDR."""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    bh_fdr,
    deseq_size_factors,
    emit,
    fail,
    parse_col_list,
    pick_col,
    read_table,
    setup_matplotlib,
    write_artifact,
)


def _rra_rho(u_sorted, alpha: float) -> float:
    from scipy.special import betainc

    k = len(u_sorted)
    if k == 0:
        return 1.0
    best = 1.0
    for j, u in enumerate(u_sorted, start=1):
        if u > alpha:
            break
        u_clip = min(max(float(u), 1e-12), 1.0 - 1e-12)
        p = float(betainc(j, k - j + 1, u_clip))
        if p < best:
            best = p
    return best


def _calibrate_rra_pvalue(rho_obs: float, null_sorted: "np.ndarray", k: int, n_perm: int) -> float:
    import numpy as np

    n_le = int(np.searchsorted(null_sorted, rho_obs, side="right"))
    if n_le >= 5:
        return (1.0 + n_le) / (1.0 + n_perm)
    # In the extreme tail beyond permutation resolution, calibrate P(rho <= t) ~ c_k * t
    # using the lower 1%-5% empirical null quantile ratio, bounded by [1, k].
    ratios = []
    for q in (0.01, 0.02, 0.05):
        idx = max(0, min(len(null_sorted) - 1, int(q * len(null_sorted))))
        val = float(null_sorted[idx])
        if 0.0 < val < 1.0:
            ratios.append(q / val)
    c_k = float(np.clip(np.median(ratios) if ratios else float(k), 1.0, float(max(1, k))))
    p_tail = min(5.0 / (1.0 + n_perm), max(1e-300, c_k * float(rho_obs)))
    return float(p_tail)


def _null_rra_sorted(k: int, n_sgrnas: int, alpha: float, n_perm: int, rng) -> "np.ndarray":
    import numpy as np
    from scipy.special import betainc

    # Sample k uniform rank quantiles per permutation
    ranks = rng.uniform(1.0 / n_sgrnas, 1.0, size=(n_perm, k))
    ranks.sort(axis=1)
    j_arr = np.arange(1, k + 1, dtype=float)[None, :]
    a_arr = j_arr
    b_arr = (k - j_arr + 1.0)
    p_mat = betainc(a_arr, b_arr, np.clip(ranks, 1e-12, 1.0 - 1e-12))
    p_masked = np.where(ranks <= alpha, p_mat, 1.0)
    rhos = np.min(p_masked, axis=1)
    rhos.sort()
    return rhos


def main() -> None:
    ap = argparse.ArgumentParser(description="Call gene hits via alpha-RRA from CRISPR screen counts")
    ap.add_argument("--counts", required=True, help="sgRNA count matrix TSV/CSV")
    ap.add_argument("--treatment_cols", action="append", required=True)
    ap.add_argument("--control_cols", action="append", required=True)
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--sgrna_col", default="")
    ap.add_argument("--gene_col", default="")
    ap.add_argument("--paired", action="store_true", default=False)
    ap.add_argument("--norm_method", choices=["median", "control_sgrna", "total"], default="median")
    ap.add_argument("--control_gene_prefix", default="NTC,NonTargeting")
    ap.add_argument("--alpha_cutoff", type=float, default=0.10)
    ap.add_argument("--n_permutations", type=int, default=2000)
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from scipy import stats as sp_stats

    trt_cols = parse_col_list(args.treatment_cols)
    ctrl_cols = parse_col_list(args.control_cols)
    if not trt_cols or not ctrl_cols:
        fail("Both --treatment_cols and --control_cols must be non-empty.")
    if args.paired and len(trt_cols) != len(ctrl_cols):
        fail(f"Paired mode requires equal treatment and control replicate counts ({len(trt_cols)} vs {len(ctrl_cols)}).")
    if not (0.0 < args.alpha_cutoff <= 1.0):
        fail("--alpha_cutoff must be in (0, 1].")

    df = read_table(args.counts).copy()
    cols = list(df.columns)
    sg_c = pick_col(cols, args.sgrna_col, ["sgRNA", "sgrna", "guide_id", "id"])
    gn_c = pick_col(cols, args.gene_col, ["Gene", "gene", "gene_id", "target", "symbol"])

    for c in trt_cols + ctrl_cols:
        if c not in cols:
            fail(f"Sample column '{c}' not found in counts matrix columns: {cols}")

    all_cols = ctrl_cols + trt_cols
    raw_mat = df[all_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float, copy=True)
    n_ctrl = len(ctrl_cols)
    n_trt = len(trt_cols)
    warnings: list[str] = []

    norm_used = args.norm_method
    if norm_used == "control_sgrna":
        prefixes = [p.strip().lower() for p in args.control_gene_prefix.split(",") if p.strip()]
        genes_low = df[gn_c].astype(str).str.strip().str.lower()
        ctrl_mask = genes_low.map(lambda g: any(g.startswith(p) or p in g for p in prefixes)).to_numpy()
        if int(ctrl_mask.sum()) >= 5:
            sf = deseq_size_factors(raw_mat[ctrl_mask])
        else:
            warnings.append(
                f"Fewer than 5 control sgRNAs matched prefixes '{args.control_gene_prefix}'; falling back to median-of-ratios normalization."
            )
            sf = deseq_size_factors(raw_mat)
            norm_used = "median"
    elif norm_used == "total":
        totals = raw_mat.sum(axis=0)
        med_tot = float(np.median(totals[totals > 0])) if np.any(totals > 0) else 1.0
        sf = np.where(totals > 0, totals / med_tot, 1.0)
    else:
        sf = deseq_size_factors(raw_mat)

    norm_mat = raw_mat / sf[None, :]
    ctrl_mat = norm_mat[:, :n_ctrl]
    trt_mat = norm_mat[:, n_ctrl:]

    mu_c = np.mean(ctrl_mat, axis=1)
    mu_t = np.mean(trt_mat, axis=1)

    if args.paired:
        pair_lfc = np.log2((trt_mat + 0.5) / (ctrl_mat + 0.5))
        lfc = np.mean(pair_lfc, axis=1)
        obs_var = np.var(pair_lfc, axis=1, ddof=1) if n_trt >= 2 else np.full(len(df), 0.15)
        # Moderated log2FC variance
        mean_abund = 0.5 * (mu_c + mu_t)
        trend_var = np.maximum(0.04, 1.2 / (mean_abund + 4.0) + 0.03)
        mod_var = 0.35 * obs_var + 0.65 * trend_var
        se = np.sqrt(mod_var / max(1, n_trt))
        z_sgrna = lfc / np.maximum(se, 1e-4)
    else:
        lfc = np.log2((mu_t + 0.5) / (mu_c + 0.5))
        if n_ctrl >= 2 and n_trt >= 2:
            var_c = np.var(ctrl_mat, axis=1, ddof=1)
            var_t = np.var(trt_mat, axis=1, ddof=1)
            pooled_var = 0.5 * (var_c + var_t)
        elif n_ctrl >= 2:
            pooled_var = np.var(ctrl_mat, axis=1, ddof=1)
        elif n_trt >= 2:
            pooled_var = np.var(trt_mat, axis=1, ddof=1)
        else:
            pooled_var = np.maximum(mu_c, 1.0)

        # Fit empirical mean-variance trend in log space
        mu_all = np.maximum(0.5 * (mu_c + mu_t), 0.5)
        fit_mask = (mu_all > 1.0) & (pooled_var > 0)
        if int(fit_mask.sum()) >= 20:
            slope, intercept = np.polyfit(np.log(mu_all[fit_mask]), np.log(pooled_var[fit_mask]), 1)
            slope = float(np.clip(slope, 0.8, 2.2))
            trend_var = np.exp(intercept + slope * np.log(mu_all))
        else:
            trend_var = mu_all + 0.05 * (mu_all ** 2)

        # Shrink observed variance toward mean-variance trend with Poisson floor
        mod_var = np.maximum(mu_all + 0.5, 0.3 * pooled_var + 0.7 * trend_var)
        se_diff = np.sqrt(mod_var * (1.0 / max(1, n_ctrl) + 1.0 / max(1, n_trt)))
        z_sgrna = (mu_t - mu_c) / np.maximum(se_diff, 1e-4)

    p_neg_sg = sp_stats.norm.cdf(z_sgrna)
    p_pos_sg = sp_stats.norm.sf(z_sgrna)

    M = len(df)
    # Ranks 1..M (most negative z -> smallest depletion rank; most positive z -> smallest enrichment rank)
    neg_order = np.argsort(z_sgrna, kind="mergesort")
    pos_order = np.argsort(-z_sgrna, kind="mergesort")
    u_neg = np.empty(M, dtype=float)
    u_pos = np.empty(M, dtype=float)
    u_neg[neg_order] = np.arange(1, M + 1, dtype=float) / M
    u_pos[pos_order] = np.arange(1, M + 1, dtype=float) / M

    sgrna_df = pd.DataFrame(
        {
            "sgRNA": df[sg_c].astype(str).str.strip(),
            "gene": df[gn_c].astype(str).str.strip(),
            "control_mean": np.round(mu_c, 3),
            "treatment_mean": np.round(mu_t, 3),
            "log2fc": np.round(lfc, 4),
            "z_score": np.round(z_sgrna, 4),
            "p_neg": p_neg_sg,
            "p_pos": p_pos_sg,
            "neg_rank_u": np.round(u_neg, 6),
            "pos_rank_u": np.round(u_pos, 6),
        }
    )

    # Group by gene and run alpha-RRA + permutation null stratified by k
    rng = np.random.default_rng(42)
    null_cache: dict[int, np.ndarray] = {}
    gene_rows: list[dict] = []

    for gene, idx in sgrna_df.groupby("gene", sort=False).groups.items():
        idx_arr = np.asarray(list(idx), dtype=int)
        k = int(idx_arr.size)
        if k not in null_cache:
            null_cache[k] = _null_rra_sorted(k, M, args.alpha_cutoff, args.n_permutations, rng)
        null_sorted = null_cache[k]

        un = np.sort(u_neg[idx_arr])
        up = np.sort(u_pos[idx_arr])
        rho_neg = _rra_rho(un, args.alpha_cutoff)
        rho_pos = _rra_rho(up, args.alpha_cutoff)

        p_neg_g = _calibrate_rra_pvalue(rho_neg, null_sorted, k, args.n_permutations)
        p_pos_g = _calibrate_rra_pvalue(rho_pos, null_sorted, k, args.n_permutations)

        g_lfc = lfc[idx_arr]
        g_z = z_sgrna[idx_arr]
        stouffer_z = float(np.sum(g_z) / math.sqrt(k))
        neg_good = int(np.sum((u_neg[idx_arr] <= args.alpha_cutoff) & (g_lfc < 0)))
        pos_good = int(np.sum((u_pos[idx_arr] <= args.alpha_cutoff) & (g_lfc > 0)))

        gene_rows.append(
            {
                "gene": gene,
                "n_sgrnas": k,
                "mean_lfc": round(float(np.mean(g_lfc)), 4),
                "median_lfc": round(float(np.median(g_lfc)), 4),
                "stouffer_z": round(stouffer_z, 4),
                "neg_rra_score": rho_neg,
                "neg_pvalue": p_neg_g,
                "neg_good_sgrnas": neg_good,
                "pos_rra_score": rho_pos,
                "pos_pvalue": p_pos_g,
                "pos_good_sgrnas": pos_good,
            }
        )

    gene_df = pd.DataFrame(gene_rows)
    gene_df["neg_fdr"] = bh_fdr(gene_df["neg_pvalue"].to_numpy())
    gene_df["pos_fdr"] = bh_fdr(gene_df["pos_pvalue"].to_numpy())
    gene_df["neg_rank"] = gene_df["neg_rra_score"].rank(method="min").astype(int)
    gene_df["pos_rank"] = gene_df["pos_rra_score"].rank(method="min").astype(int)
    gene_df = gene_df.sort_values(["neg_fdr", "neg_rra_score", "median_lfc"], ascending=[True, True, True])

    n_dep_05 = int((gene_df["neg_fdr"] < 0.05).sum())
    n_enr_05 = int((gene_df["pos_fdr"] < 0.05).sum())
    if n_dep_05 == 0 and n_enr_05 == 0:
        warnings.append("No genes reached BH FDR < 0.05 for positive or negative selection.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sg_path = out_dir / "sgrna_summary.tsv"
    gn_path = out_dir / "gene_hits_summary.tsv"
    fig_path = out_dir / "screen_hits_volcano_rank.png"

    sgrna_df.to_csv(sg_path, sep="\t", index=False)
    gene_df.to_csv(gn_path, sep="\t", index=False)

    write_artifact(
        sg_path,
        "crispr_screen_hits",
        "sgRNA-Level Normalized Counts, log2FC & Moderated Z-Scores",
        inputs=[args.counts],
        parameters={"norm_method": norm_used, "paired": bool(args.paired)},
        shape=(int(len(sgrna_df)), int(sgrna_df.shape[1])),
    )
    write_artifact(
        gn_path,
        "crispr_screen_hits",
        "Gene-Level alpha-RRA Positive & Negative Selection Results",
        inputs=[args.counts],
        parameters={"alpha_cutoff": args.alpha_cutoff, "n_permutations": args.n_permutations},
        shape=(int(len(gene_df)), int(gene_df.shape[1])),
    )

    # 2-panel figure: Volcano + Gene Rank Waterfall
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.8), dpi=180)

    min_fdr = np.minimum(gene_df["neg_fdr"].to_numpy(), gene_df["pos_fdr"].to_numpy())
    neg_log10_fdr = -np.log10(np.clip(min_fdr, 1e-10, 1.0))
    lfc_g = gene_df["median_lfc"].to_numpy()
    is_dep = gene_df["neg_fdr"].to_numpy() < 0.05
    is_enr = gene_df["pos_fdr"].to_numpy() < 0.05
    is_bg = ~(is_dep | is_enr)

    ax1.scatter(lfc_g[is_bg], neg_log10_fdr[is_bg], s=14, c="#cbd5e1", alpha=0.7, label="NS")
    if is_dep.any():
        ax1.scatter(lfc_g[is_dep], neg_log10_fdr[is_dep], s=26, c="#2563eb", alpha=0.9, label=f"Depleted (FDR<0.05: {n_dep_05})")
    if is_enr.any():
        ax1.scatter(lfc_g[is_enr], neg_log10_fdr[is_enr], s=26, c="#dc2626", alpha=0.9, label=f"Enriched (FDR<0.05: {n_enr_05})")
    ax1.axhline(-math.log10(0.05), color="#64748b", ls="--", lw=0.9)
    ax1.axvline(0.0, color="#94a3b8", ls=":", lw=0.8)

    for _, row in gene_df[gene_df["neg_fdr"] < 0.05].head(6).iterrows():
        y_v = -math.log10(max(float(row["neg_fdr"]), 1e-10))
        ax1.annotate(row["gene"], (row["median_lfc"], y_v), fontsize=7, xytext=(-4, 3), textcoords="offset points", color="#1e3a8a")
    for _, row in gene_df.sort_values(["pos_fdr", "pos_rra_score"]).head(6).iterrows():
        if row["pos_fdr"] < 0.05:
            y_v = -math.log10(max(float(row["pos_fdr"]), 1e-10))
            ax1.annotate(row["gene"], (row["median_lfc"], y_v), fontsize=7, xytext=(4, 3), textcoords="offset points", color="#991b1b")

    ax1.set_xlabel("Gene Median log2(Treatment / Control)")
    ax1.set_ylabel("-log10(min BH FDR)")
    ax1.set_title("CRISPR Screen Volcano (alpha-RRA)")
    ax1.legend(frameon=False, fontsize=7.5, loc="upper center")

    wf = gene_df.sort_values("stouffer_z", ascending=True).reset_index(drop=True)
    ranks = np.arange(1, len(wf) + 1)
    wf_z = wf["stouffer_z"].to_numpy()
    cols_wf = np.where(wf["neg_fdr"] < 0.05, "#2563eb", np.where(wf["pos_fdr"] < 0.05, "#dc2626", "#94a3b8"))
    ax2.scatter(ranks, wf_z, s=14, c=cols_wf, alpha=0.85)
    ax2.axhline(0.0, color="#64748b", ls="--", lw=0.9)
    for _, row in wf.head(4).iterrows():
        ax2.annotate(row["gene"], (row.name + 1, row["stouffer_z"]), fontsize=7, xytext=(4, 0), textcoords="offset points", color="#1e3a8a")
    for _, row in wf.tail(4).iterrows():
        ax2.annotate(row["gene"], (row.name + 1, row["stouffer_z"]), fontsize=7, xytext=(-28, -4), textcoords="offset points", color="#991b1b")
    ax2.set_xlabel("Gene Rank (Depleted → Enriched)")
    ax2.set_ylabel("Stouffer's Combined Z-Score")
    ax2.set_title("Gene Selection Waterfall")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "crispr_screen_hits",
        "CRISPR Screen Volcano & Gene Rank Waterfall",
        inputs=[args.counts],
        parameters={"alpha_cutoff": args.alpha_cutoff},
    )

    top_dep = gene_df.sort_values(["neg_fdr", "neg_rra_score", "median_lfc"]).head(15).to_dict(orient="records")
    top_enr = gene_df.sort_values(["pos_fdr", "pos_rra_score", "median_lfc"], ascending=[True, True, False]).head(15).to_dict(orient="records")

    emit(
        {
            "status": "ok",
            "tool": "crispr_screen_hits",
            "norm_method": norm_used,
            "paired": bool(args.paired),
            "alpha_cutoff": args.alpha_cutoff,
            "n_permutations": args.n_permutations,
            "n_sgrnas": M,
            "n_genes": int(len(gene_df)),
            "n_depleted_fdr05": n_dep_05,
            "n_enriched_fdr05": n_enr_05,
            "top_depleted_genes": top_dep,
            "top_enriched_genes": top_enr,
            "outputs": {
                "sgrna_summary": str(sg_path),
                "gene_hits_summary": str(gn_path),
                "hits_plot": str(fig_path),
            },
            "artifacts": [str(sg_path), str(gn_path), str(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
