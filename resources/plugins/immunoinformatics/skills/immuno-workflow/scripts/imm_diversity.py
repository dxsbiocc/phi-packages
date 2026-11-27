#!/usr/bin/env python3
"""imm_diversity: Compute repertoire diversity (Shannon bits/nats, Pielou
clonality, Gini-Simpson, Inverse Simpson, D50), clonal homeostasis space tiers,
V-J gene pairing usage matrix, optional Mann-Whitney / Kruskal-Wallis + BH FDR
group tests, and a 2-panel diversity/homeostasis figure.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    bh_fdr,
    emit,
    fail,
    harmonize_repertoire_file,
    normalize_gene_base,
    read_table,
    resolve_column,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)


def _compute_sample_metrics(counts: np.ndarray, sample_id: str) -> dict[str, Any]:
    counts = np.sort(counts[counts > 0])[::-1]
    richness = int(counts.size)
    total = float(np.sum(counts))
    if richness == 0 or total <= 0:
        return {
            "sample_id": sample_id,
            "total_reads_umis": 0.0,
            "richness_clonotypes": 0,
            "shannon_entropy_nats": 0.0,
            "shannon_entropy_bits": 0.0,
            "pielou_clonality": np.nan,
            "simpson_dominance": np.nan,
            "gini_simpson": np.nan,
            "inverse_simpson": np.nan,
            "d50_count": 0,
            "d50_pct": np.nan,
            "homeostasis_hyperexpanded_pct": 0.0,
            "homeostasis_large_pct": 0.0,
            "homeostasis_medium_pct": 0.0,
            "homeostasis_small_pct": 0.0,
        }

    p = counts / total
    h_nats = float(-np.sum(p * np.log(p)))
    h_bits = float(h_nats / math.log(2.0))
    clonality = float(1.0 - h_nats / math.log(richness)) if richness > 1 else 0.0
    clonality = float(np.clip(clonality, 0.0, 1.0))

    simpson_dom = float(np.sum(p ** 2))
    gini_simp = float(1.0 - simpson_dom)
    inv_simp = float(1.0 / simpson_dom) if simpson_dom > 0 else np.nan

    cum_p = np.cumsum(p)
    d50_idx = int(np.searchsorted(cum_p, 0.50, side="left")) + 1
    d50_idx = min(d50_idx, richness)
    d50_pct = float(100.0 * d50_idx / richness)

    hyper = float(100.0 * np.sum(p[p >= 0.01]))
    large = float(100.0 * np.sum(p[(p >= 0.001) & (p < 0.01)]))
    medium = float(100.0 * np.sum(p[(p >= 0.0001) & (p < 0.001)]))
    small = float(100.0 * np.sum(p[p < 0.0001]))

    return {
        "sample_id": sample_id,
        "total_reads_umis": round(total, 2),
        "richness_clonotypes": richness,
        "shannon_entropy_nats": round(h_nats, 5),
        "shannon_entropy_bits": round(h_bits, 5),
        "pielou_clonality": round(clonality, 5),
        "simpson_dominance": round(simpson_dom, 6),
        "gini_simpson": round(gini_simp, 6),
        "inverse_simpson": round(inv_simp, 3),
        "d50_count": d50_idx,
        "d50_pct": round(d50_pct, 3),
        "homeostasis_hyperexpanded_pct": round(hyper, 3),
        "homeostasis_large_pct": round(large, 3),
        "homeostasis_medium_pct": round(medium, 3),
        "homeostasis_small_pct": round(small, 3),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repertoire", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--group_col", default=None)
    ap.add_argument("--locus", default="TRB")
    args = ap.parse_args()

    df, detected_fmt = harmonize_repertoire_file(args.repertoire, fmt="auto")
    warnings: list[str] = []

    target_locus = args.locus.strip().upper()
    if target_locus != "ALL":
        df = df[df["locus"] == target_locus].copy()
    if df.empty:
        fail(f"No records found for locus '{args.locus}'.")

    prod_df = df[(df["productive"]) & (df["cdr3_aa"].astype(str).str.len() > 0) & (df["duplicate_count"] > 0)].copy()
    if prod_df.empty:
        fail(f"No productive clonotypes with positive counts remain for locus '{args.locus}'.")

    prod_df["v_gene"] = prod_df["v_call"].map(normalize_gene_base).replace("", "UNKNOWN_V")
    prod_df["j_gene"] = prod_df["j_call"].map(normalize_gene_base).replace("", "UNKNOWN_J")

    # Aggregate clonotypes by (sample_id, cdr3_aa, v_gene) for diversity
    clono_df = (
        prod_df.groupby(["sample_id", "cdr3_aa", "v_gene"], as_index=False, sort=False)["duplicate_count"]
        .sum()
    )

    sample_order = list(dict.fromkeys(prod_df["sample_id"].astype(str).tolist()))
    div_rows: list[dict[str, Any]] = []
    for sid in sample_order:
        sub = clono_df[clono_df["sample_id"] == sid]
        cnts = sub["duplicate_count"].to_numpy(dtype=float, copy=True)
        div_rows.append(_compute_sample_metrics(cnts, sid))
    div_df = pd.DataFrame(div_rows)

    # V-J pairing matrix
    vj_counts = (
        prod_df.groupby(["sample_id", "v_gene", "j_gene"], as_index=False)["duplicate_count"]
        .sum()
    )
    sample_totals = vj_counts.groupby("sample_id")["duplicate_count"].transform("sum")
    vj_counts["freq"] = np.where(sample_totals > 0, vj_counts["duplicate_count"] / sample_totals, 0.0)
    vj_piv = vj_counts.pivot_table(
        index=["v_gene", "j_gene"],
        columns="sample_id",
        values="freq",
        fill_value=0.0,
    )
    vj_tot = vj_counts.groupby(["v_gene", "j_gene"])["duplicate_count"].sum()
    vj_matrix = vj_piv.copy()
    vj_matrix.columns = [f"{c}_freq" for c in vj_matrix.columns]
    vj_matrix.insert(0, "mean_frequency", vj_piv.mean(axis=1).round(6))
    vj_matrix.insert(0, "total_count", vj_tot.reindex(vj_matrix.index).fillna(0.0).round(2))
    vj_matrix = vj_matrix.reset_index().sort_values("total_count", ascending=False).reset_index(drop=True)

    # Optional metadata + group comparison
    group_col_used: str | None = args.group_col
    group_comp_rows: list[dict[str, Any]] = []
    if args.metadata:
        meta_df = read_table(args.metadata)
        ms_col = resolve_column(meta_df, ["sample_id", "sample", "repertoire_id", "subject_id"], "sample_id")
        if group_col_used is None:
            group_col_used = resolve_column(
                meta_df,
                ["group", "condition", "response", "arm", "treatment", "cohort"],
                "group",
                required=False,
            )
            if group_col_used is None:
                fail("Metadata provided without --group_col and no standard group column found.")
        else:
            group_col_used = resolve_column(meta_df, [group_col_used], "group_col")

        meta_sub = meta_df[[ms_col, group_col_used]].copy()
        meta_sub.columns = ["sample_id", "group"]
        meta_sub["sample_id"] = meta_sub["sample_id"].astype(str).str.strip()
        meta_sub["group"] = meta_sub["group"].astype(str).str.strip()
        meta_sub = meta_sub.drop_duplicates(subset=["sample_id"])
        div_df = div_df.merge(meta_sub, on="sample_id", how="left")

        valid_grp = div_df.dropna(subset=["group"]).copy()
        valid_grp = valid_grp[valid_grp["group"].astype(str).str.lower() != "nan"]
        groups = [g for g, sub in valid_grp.groupby("group", sort=False) if len(sub) >= 1]
        if len(groups) < 2:
            warnings.append("Fewer than 2 groups available in metadata for statistical testing.")
        else:
            rep_counts = [int((valid_grp["group"] == g).sum()) for g in groups]
            if min(rep_counts) < 2:
                warnings.append(
                    f"Some groups have <2 biological samples ({dict(zip(groups, rep_counts))}); skipping group hypothesis tests."
                )
            else:
                metrics_to_test = [
                    "pielou_clonality",
                    "shannon_entropy_bits",
                    "gini_simpson",
                    "inverse_simpson",
                    "d50_pct",
                    "homeostasis_hyperexpanded_pct",
                ]
                raw_pvals: list[float] = []
                for metric in metrics_to_test:
                    arrays = [
                        valid_grp.loc[valid_grp["group"] == g, metric].dropna().to_numpy(dtype=float, copy=True)
                        for g in groups
                    ]
                    if any(a.size < 2 for a in arrays):
                        continue
                    grp_means = {g: round(float(np.mean(a)), 5) for g, a in zip(groups, arrays)}
                    if len(groups) == 2:
                        stat, pval = stats.mannwhitneyu(arrays[0], arrays[1], alternative="two-sided")
                        test_name = "mann_whitney_u"
                    else:
                        stat, pval = stats.kruskal(*arrays)
                        test_name = "kruskal_wallis"
                    raw_pvals.append(float(pval))
                    group_comp_rows.append(
                        {
                            "metric": metric,
                            "test": test_name,
                            "n_groups": len(groups),
                            "groups": ",".join(groups),
                            "group_means": json_inline(grp_means),
                            "statistic": round(float(stat), 5),
                            "p_value": float(pval),
                        }
                    )
                if group_comp_rows:
                    fdrs = bh_fdr([r["p_value"] for r in group_comp_rows])
                    for r, q in zip(group_comp_rows, fdrs):
                        r["fdr"] = round(q, 6)
                        r["p_value"] = round(r["p_value"], 6)
                        r["significant_fdr_05"] = bool(q < 0.05)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    params = {
        "locus": target_locus,
        "group_col": group_col_used,
        "detected_format": detected_fmt,
    }
    sci_ctx = {"locus": target_locus, "n_samples": int(len(div_df))}

    div_path = write_tsv(
        div_df,
        out_dir / "diversity_metrics.tsv",
        tool="imm_diversity",
        title="Per-Sample Repertoire Diversity & Clonal Homeostasis Metrics",
        inputs=[args.repertoire, args.metadata],
        parameters=params,
        scientific_context=sci_ctx,
    )
    vj_path = write_tsv(
        vj_matrix,
        out_dir / "vj_usage_matrix.tsv",
        tool="imm_diversity",
        title="V-J Gene Pairing Usage Matrix",
        inputs=[args.repertoire],
        parameters=params,
        scientific_context=sci_ctx,
    )
    artifacts = [div_path, vj_path]

    comp_path: str | None = None
    if group_comp_rows:
        comp_df = pd.DataFrame(group_comp_rows)
        comp_path = write_tsv(
            comp_df,
            out_dir / "diversity_group_comparisons.tsv",
            tool="imm_diversity",
            title="Group Statistical Comparisons of Repertoire Diversity",
            inputs=[args.repertoire, args.metadata],
            parameters=params,
            scientific_context=sci_ctx,
        )
        artifacts.append(comp_path)

    # 2-panel figure: Homeostasis + Clonality vs Shannon
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.8))

    sids = div_df["sample_id"].tolist()
    x = np.arange(len(sids))
    h_hyp = div_df["homeostasis_hyperexpanded_pct"].to_numpy(dtype=float, copy=True)
    h_lrg = div_df["homeostasis_large_pct"].to_numpy(dtype=float, copy=True)
    h_med = div_df["homeostasis_medium_pct"].to_numpy(dtype=float, copy=True)
    h_sml = div_df["homeostasis_small_pct"].to_numpy(dtype=float, copy=True)

    ax1.bar(x, h_hyp, color="#b91c1c", label="Hyperexpanded (>=1%)")
    ax1.bar(x, h_lrg, bottom=h_hyp, color="#f59e0b", label="Large (0.1–1%)")
    ax1.bar(x, h_med, bottom=h_hyp + h_lrg, color="#0ea5e9", label="Medium (0.01–0.1%)")
    ax1.bar(x, h_sml, bottom=h_hyp + h_lrg + h_med, color="#0f766e", label="Small (<0.01%)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(sids, rotation=30, ha="right")
    ax1.set_ylim(0, 105)
    ax1.set_ylabel("Clonal Space (%)")
    ax1.set_title(f"Clonal Homeostasis Space ({target_locus})")
    ax1.legend(frameon=False, loc="upper right", fontsize=8)

    shannon_vals = div_df["shannon_entropy_bits"].to_numpy(dtype=float, copy=True)
    clonality_vals = div_df["pielou_clonality"].fillna(0.0).to_numpy(dtype=float, copy=True)
    if "group" in div_df.columns and div_df["group"].notna().any():
        uniq_grps = [str(g) for g in div_df["group"].fillna("NA").unique()]
        pal = plt.cm.Set1(np.linspace(0, 1, max(3, len(uniq_grps))))
        for gi, gname in enumerate(uniq_grps):
            mask = div_df["group"].fillna("NA").astype(str) == gname
            ax2.scatter(
                shannon_vals[mask],
                clonality_vals[mask],
                s=75,
                color=pal[gi % len(pal)],
                alpha=0.88,
                edgecolors="#1e293b",
                linewidth=0.7,
                label=gname,
            )
        ax2.legend(frameon=False, title=group_col_used or "Group")
    else:
        ax2.scatter(shannon_vals, clonality_vals, s=75, color="#0f766e", alpha=0.88, edgecolors="#1e293b", linewidth=0.7)

    for idx_pt, sid in enumerate(sids):
        ax2.annotate(sid, (shannon_vals[idx_pt], clonality_vals[idx_pt]), textcoords="offset points", xytext=(5, 4), fontsize=8)
    ax2.set_xlabel("Shannon Entropy (bits)")
    ax2.set_ylabel("Pielou Clonality (1 - H / ln S)")
    ax2.set_title(f"Diversity vs. Clonality ({target_locus})")
    ax2.grid(alpha=0.25, linestyle="--")

    fig.tight_layout()
    fig_file = out_dir / "repertoire_diversity_homeostasis.png"
    fig.savefig(fig_file, bbox_inches="tight")
    plt.close(fig)
    fig_path = write_artifact(
        fig_file,
        tool="imm_diversity",
        title="Repertoire Clonal Homeostasis & Diversity Landscape",
        inputs=[args.repertoire, args.metadata],
        parameters=params,
        scientific_context=sci_ctx,
    )
    artifacts.append(fig_path)

    emit(
        {
            "status": "ok",
            "tool": "imm_diversity",
            "locus": target_locus,
            "group_col": group_col_used,
            "n_samples": int(len(div_df)),
            "samples": div_df.to_dict(orient="records"),
            "group_comparisons": group_comp_rows,
            "outputs": {
                "diversity_metrics": div_path,
                "vj_usage_matrix": vj_path,
                "diversity_group_comparisons": comp_path,
                "repertoire_diversity_homeostasis": fig_path,
            },
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


def json_inline(obj: dict[str, Any]) -> str:
    import json

    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"))


if __name__ == "__main__":
    main()
