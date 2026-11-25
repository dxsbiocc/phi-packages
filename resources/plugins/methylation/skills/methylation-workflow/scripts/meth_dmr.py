#!/usr/bin/env python3
"""meth_dmr: Empirical-Bayes moderated linear model on logit M-values (with
optional confounder covariates) + biological delta-beta DMC calling + sliding
coordinate-based Differentially Methylated Region (DMR) aggregation (Stouffer Z).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    beta_to_m,
    bh_fdr,
    emit,
    fail,
    m_to_beta,
    moderated_ols_contrast,
    parse_cpg_coord,
    read_methylation_matrix,
    read_table,
    write_artifact,
)


def load_probe_coords(path: str) -> dict[str, tuple[str, int]]:
    df = read_table(path)
    if df.shape[1] < 3:
        fail("probe_coords file must contain at least 3 columns: probe_id, chrom, pos (or BED chrom, start, end, probe_id)")
    lower = {str(c).strip().lower(): i for i, c in enumerate(df.columns)}
    if "probe_id" in lower or "id" in lower or "cpg_id" in lower:
        id_i = lower.get("probe_id", lower.get("cpg_id", lower.get("id", 0)))
        ch_i = lower.get("chrom", lower.get("chr", lower.get("chromosome", 1)))
        po_i = lower.get("pos", lower.get("start", lower.get("mapinfo", 2)))
    elif df.shape[1] >= 4 and str(df.iloc[0, 0]).startswith("chr"):
        ch_i, po_i, id_i = 0, 1, 3
    else:
        id_i, ch_i, po_i = 0, 1, 2
    out: dict[str, tuple[str, int]] = {}
    for _, row in df.iterrows():
        try:
            pid = str(row.iloc[id_i]).strip()
            ch = str(row.iloc[ch_i]).strip()
            pos = int(float(row.iloc[po_i]))
            if pid and ch:
                out[pid] = (ch, pos)
        except (ValueError, TypeError):
            continue
    return out


def call_dmrs(
    dmc_df: pd.DataFrame, max_dist: int, min_cpgs: int, delta_beta_cut: float, fdr_cut: float
) -> pd.DataFrame:
    """Chain spatially adjacent CpGs with seed pvalue < 0.05 and consistent sign(delta_beta)."""
    sub = dmc_df.dropna(subset=["chrom", "pos", "pvalue", "delta_beta"]).copy()
    if sub.empty:
        return pd.DataFrame()
    sub["pos"] = sub["pos"].astype(int)
    sub = sub.sort_values(["chrom", "pos"]).reset_index(drop=True)

    # Seed CpGs: nominal p < 0.05 and |delta_beta| >= 0.5 * delta_beta_cut
    seeds = sub[(sub["pvalue"] < 0.05) & (sub["delta_beta"].abs() >= 0.5 * delta_beta_cut)].copy()
    if len(seeds) < min_cpgs:
        return pd.DataFrame()

    clusters: list[pd.DataFrame] = []
    cur_idx: list[int] = []
    prev_ch = None
    prev_pos = -10**9
    prev_sign = 0

    for idx, row in seeds.iterrows():
        ch = str(row["chrom"])
        pos = int(row["pos"])
        sgn = 1 if float(row["delta_beta"]) > 0 else -1
        if ch == prev_ch and sgn == prev_sign and (pos - prev_pos) <= max_dist:
            cur_idx.append(idx)
        else:
            if len(cur_idx) >= min_cpgs:
                clusters.append(seeds.loc[cur_idx])
            cur_idx = [idx]
        prev_ch, prev_pos, prev_sign = ch, pos, sgn
    if len(cur_idx) >= min_cpgs:
        clusters.append(seeds.loc[cur_idx])

    if not clusters:
        return pd.DataFrame()

    dmr_rows = []
    for i, cdf in enumerate(clusters, start=1):
        ch = str(cdf["chrom"].iloc[0])
        start = int(cdf["pos"].min())
        end = int(cdf["pos"].max())
        pvals = np.clip(cdf["pvalue"].to_numpy(dtype=float), 1e-300, 1.0 - 1e-15)
        signs = np.sign(cdf["delta_beta"].to_numpy(dtype=float))
        # Signed Stouffer's Z across regional CpGs
        z_i = stats.norm.isf(pvals / 2.0) * signs
        z_comb = float(np.sum(z_i) / np.sqrt(len(z_i)))
        p_comb = float(2.0 * stats.norm.sf(abs(z_comb)))
        mean_db = float(cdf["delta_beta"].mean())
        dmr_rows.append(
            {
                "dmr_id": f"DMR_{i:04d}",
                "chrom": ch,
                "start": start,
                "end": max(end, start + 1),
                "width_bp": max(1, end - start + 1),
                "n_cpgs": int(len(cdf)),
                "mean_beta_case": round(float(cdf["mean_beta_case"].mean()), 4),
                "mean_beta_control": round(float(cdf["mean_beta_control"].mean()), 4),
                "mean_delta_beta": round(mean_db, 4),
                "direction": "hypermethylated" if mean_db > 0 else "hypomethylated",
                "stouffer_z": round(z_comb, 4),
                "region_pvalue": p_comb,
            }
        )

    out = pd.DataFrame(dmr_rows)
    out["region_fdr"] = bh_fdr(out["region_pvalue"].to_numpy(dtype=float))
    out["significant"] = (out["region_fdr"] < fdr_cut) & (out["mean_delta_beta"].abs() >= delta_beta_cut)
    return out.sort_values(["region_pvalue", "mean_delta_beta"], ascending=[True, False]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--control", required=True)
    parser.add_argument("--covariates", action="extend", nargs="+", default=[])
    parser.add_argument("--probe_coords", default=None)
    parser.add_argument("--delta_beta", type=float, default=0.10)
    parser.add_argument("--fdr_cutoff", type=float, default=0.05)
    parser.add_argument("--max_dist", type=int, default=500)
    parser.add_argument("--min_cpgs", type=int, default=3)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    covs = [c.strip() for item in args.covariates for c in str(item).split(",") if c.strip()]
    mat_df, scale = read_methylation_matrix(args.matrix)
    meta = read_table(args.metadata)
    scol = args.sample_col or str(meta.columns[0])
    if scol not in meta.columns or args.group_col not in meta.columns:
        fail(f"sample_col '{scol}' or group_col '{args.group_col}' missing from metadata")

    meta[scol] = meta[scol].astype(str).str.strip()
    meta[args.group_col] = meta[args.group_col].astype(str).str.strip()
    sub_meta = meta[
        meta[scol].isin(set(mat_df.columns)) & meta[args.group_col].isin({args.case, args.control})
    ].drop_duplicates(subset=[scol])

    case_samples = sub_meta[sub_meta[args.group_col] == args.case][scol].tolist()
    ctrl_samples = sub_meta[sub_meta[args.group_col] == args.control][scol].tolist()
    if len(case_samples) < 2 or len(ctrl_samples) < 2:
        fail(
            f"meth_dmr requires >= 2 biological replicates per group; found {args.case}={len(case_samples)}, {args.control}={len(ctrl_samples)}"
        )

    ordered_samples = ctrl_samples + case_samples
    meta_aligned = sub_meta.set_index(scol).loc[ordered_samples]

    # Build OLS design X: [Intercept, is_case, ...covariates]
    design_cols = [
        np.ones(len(ordered_samples), dtype=float),
        (meta_aligned[args.group_col] == args.case).to_numpy(dtype=float),
    ]
    for cov in covs:
        if cov not in meta_aligned.columns:
            fail(f"covariate '{cov}' not found in metadata")
        s = meta_aligned[cov]
        num = pd.to_numeric(s, errors="coerce")
        if num.notna().all() and num.nunique() > 2:
            v = num.to_numpy(dtype=float)
            design_cols.append((v - v.mean()) / max(v.std(ddof=1), 1e-12))
        else:
            cats = s.fillna("MISSING").astype(str)
            for lvl in sorted(cats.unique())[1:]:
                design_cols.append((cats == lvl).to_numpy(dtype=float))
    X_design = np.column_stack(design_cols)

    raw_sub = mat_df[ordered_samples].to_numpy(dtype=float)
    beta_sub = raw_sub if scale == "beta" else m_to_beta(raw_sub)
    m_sub = beta_to_m(beta_sub) if scale == "beta" else raw_sub

    fit = moderated_ols_contrast(m_sub, X_design, col_idx=1)

    n_ctrl = len(ctrl_samples)
    mean_ctrl = np.nanmean(beta_sub[:, :n_ctrl], axis=1)
    mean_case = np.nanmean(beta_sub[:, n_ctrl:], axis=1)
    delta_b = mean_case - mean_ctrl

    pcoords = load_probe_coords(args.probe_coords) if args.probe_coords else {}
    chroms: list[str | None] = []
    positions: list[int | None] = []
    for cid in mat_df.index:
        parsed = parse_cpg_coord(cid)
        if parsed is not None:
            chroms.append(parsed[0])
            positions.append(parsed[1])
        elif cid in pcoords:
            chroms.append(pcoords[cid][0])
            positions.append(pcoords[cid][1])
        else:
            chroms.append(None)
            positions.append(None)

    dmc_df = pd.DataFrame(
        {
            "cpg_id": mat_df.index,
            "chrom": chroms,
            "pos": positions,
            "mean_beta_case": np.round(mean_case, 5),
            "mean_beta_control": np.round(mean_ctrl, 5),
            "delta_beta": np.round(delta_b, 5),
            "log2fc": np.round(fit["logFC_M"], 5),
            "t": np.round(fit["t"], 4),
            "pvalue": fit["pvalue"],
            "padj": fit["fdr"],
        }
    )
    dmc_df = dmc_df.dropna(subset=["pvalue"]).reset_index(drop=True)
    dmc_df["status"] = np.where(
        (dmc_df["padj"] < args.fdr_cutoff) & (dmc_df["delta_beta"] >= args.delta_beta),
        "hyper",
        np.where(
            (dmc_df["padj"] < args.fdr_cutoff) & (dmc_df["delta_beta"] <= -args.delta_beta),
            "hypo",
            "ns",
        ),
    )

    dmr_df = call_dmrs(dmc_df, args.max_dist, args.min_cpgs, args.delta_beta, args.fdr_cutoff)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    dmc_path = outdir / "dmc_results.tsv"
    dmc_df.sort_values("pvalue", ascending=True).to_csv(dmc_path, sep="\t", index=False)

    dmr_path = outdir / "dmr_regions.tsv"
    if dmr_df.empty:
        pd.DataFrame(
            columns=[
                "dmr_id", "chrom", "start", "end", "width_bp", "n_cpgs",
                "mean_beta_case", "mean_beta_control", "mean_delta_beta",
                "direction", "stouffer_z", "region_pvalue", "region_fdr", "significant",
            ]
        ).to_csv(dmr_path, sep="\t", index=False)
    else:
        dmr_df.to_csv(dmr_path, sep="\t", index=False)

    # Volcano plot (x = delta_beta, y = -log10(pvalue))
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    nlp = -np.log10(np.clip(dmc_df["pvalue"].to_numpy(dtype=float), 1e-300, 1.0))
    db_v = dmc_df["delta_beta"].to_numpy(dtype=float)
    st_v = dmc_df["status"].to_numpy()
    for st_lbl, col, lab in (
        ("ns", "#CBD5E1", f"NS ({int((st_v == 'ns').sum())})"),
        ("hypo", "#2563EB", f"Hypomethylated ({int((st_v == 'hypo').sum())})"),
        ("hyper", "#E11D48", f"Hypermethylated ({int((st_v == 'hyper').sum())})"),
    ):
        m = st_v == st_lbl
        if np.any(m):
            ax.scatter(db_v[m], nlp[m], s=14 if st_lbl == "ns" else 24, alpha=0.45 if st_lbl == "ns" else 0.85, color=col, label=lab)
    ax.axvline(args.delta_beta, color="#64748B", linestyle="--", linewidth=0.9)
    ax.axvline(-args.delta_beta, color="#64748B", linestyle="--", linewidth=0.9)
    ax.set_xlabel(f"Δβ ({args.case} − {args.control})")
    ax.set_ylabel("−log10(moderated M-value P)")
    ax.set_title(f"Differential Methylation Volcano ({args.case} vs {args.control})", fontsize=10)
    ax.legend(frameon=False, fontsize=8.5, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig_path = outdir / "meth-dmr.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    in_files = (args.matrix, args.metadata, args.probe_coords)
    params_meta = {
        "contrast": f"{args.case}_vs_{args.control}",
        "covariates": covs,
        "delta_beta": args.delta_beta,
        "fdr_cutoff": args.fdr_cutoff,
        "max_dist": args.max_dist,
        "min_cpgs": args.min_cpgs,
    }
    artifacts = [
        write_artifact(fig_path, "meth_dmr", f"Differential methylation Δβ volcano ({args.case} vs {args.control})", inputs=in_files, parameters=params_meta),
        write_artifact(dmr_path, "meth_dmr", "Called Differentially Methylated Regions (DMRs)", inputs=in_files, parameters=params_meta, shape=(len(dmr_df), 14)),
        write_artifact(dmc_path, "meth_dmr", "Full CpG moderated M-value & Δβ DMC results", inputs=in_files, parameters=params_meta, shape=dmc_df.shape),
    ]

    warnings: list[str] = []
    if dmc_df["chrom"].isna().all():
        warnings.append(
            "no genomic coordinates available on CpG IDs: DMR aggregation skipped (pass --probe_coords for cg* array probes)"
        )

    n_sig_dmr = int(dmr_df["significant"].sum()) if not dmr_df.empty else 0
    n_hyper_dmr = int(((dmr_df["significant"]) & (dmr_df["mean_delta_beta"] > 0)).sum()) if not dmr_df.empty else 0
    n_hypo_dmr = int(((dmr_df["significant"]) & (dmr_df["mean_delta_beta"] < 0)).sum()) if not dmr_df.empty else 0

    emit(
        {
            "tool": "meth_dmr",
            "contrast": f"{args.case}_vs_{args.control}",
            "n_samples": len(ordered_samples),
            "covariates_adjusted": covs,
            "n_cpgs_tested": int(len(dmc_df)),
            "prior_df_d0": fit["d0"],
            "dmc_summary": {
                "fdr_cutoff": args.fdr_cutoff,
                "delta_beta_cutoff": args.delta_beta,
                "n_hyper_dmcs": int((dmc_df["status"] == "hyper").sum()),
                "n_hypo_dmcs": int((dmc_df["status"] == "hypo").sum()),
            },
            "dmr_summary": {
                "n_candidate_dmrs": int(len(dmr_df)),
                "n_significant_dmrs": n_sig_dmr,
                "n_hyper_dmrs": n_hyper_dmr,
                "n_hypo_dmrs": n_hypo_dmr,
                "top_dmrs": dmr_df.head(10).to_dict(orient="records") if not dmr_df.empty else [],
            },
            "outputs": {
                "dmc_table": str(dmc_path),
                "dmr_table": str(dmr_path),
            },
            "figures": [str(fig_path)],
            "scientificContext": {
                str(dmc_path): {"valueSemantics": "de_stats"},
                str(dmr_path): {"valueSemantics": "de_stats"},
            },
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
