#!/usr/bin/env python3
"""MIQE-compliant RT-qPCR relative quantification, reference stability, and group statistics."""
from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
from scipy import stats

from _common import (
    bh_fdr,
    emit,
    fail,
    parse_str_list,
    read_table,
    resolve_column,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="MIQE-compliant RT-qPCR relative quantification (wa_qpcr_quant)")
    ap.add_argument("--cq_table", required=True, help="TSV/CSV with sample, group, gene, cq")
    ap.add_argument("--reference_genes", action="append", required=True, help="Reference gene name(s)")
    ap.add_argument("--control_group", required=True, help="Calibrator / control condition name")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--standard_curve", default=None, help="Optional standard curve TSV/CSV (gene, dilution_or_copy, cq)")
    ap.add_argument("--efficiencies", action="append", default=None, help="Optional GENE=E efficiency overrides")
    ap.add_argument("--ntc_group", default="NTC", help="Group or sample label for NTC wells (default NTC)")
    return ap.parse_args()


def _estimate_efficiencies(
    genes: list[str],
    sc_path: str | None,
    eff_overrides: list[str],
    warnings: list[str],
) -> tuple[dict[str, float], dict[str, dict[str, Any]]]:
    eff_map: dict[str, float] = {g: 2.0 for g in genes}
    sc_details: dict[str, dict[str, Any]] = {}

    if sc_path:
        sc_df = read_table(sc_path)
        g_col = resolve_column(sc_df, ["gene", "target", "assay"], "standard_curve gene")
        q_col = resolve_column(
            sc_df,
            ["dilution_or_copy", "dilution", "copy", "copies", "quantity", "conc", "concentration"],
            "standard_curve dilution_or_copy",
        )
        c_col = resolve_column(sc_df, ["cq", "ct", "mean_cq"], "standard_curve cq")
        sc_df = sc_df.copy()
        sc_df["_gene"] = sc_df[g_col].astype(str).str.strip()
        sc_df["_qty"] = pd.to_numeric(sc_df[q_col], errors="coerce")
        sc_df["_cq"] = pd.to_numeric(sc_df[c_col], errors="coerce")

        for gene, sub in sc_df.groupby("_gene", sort=False):
            valid = sub[(sub["_qty"] > 0) & np.isfinite(sub["_qty"]) & np.isfinite(sub["_cq"])]
            if len(valid) < 3 or valid["_qty"].nunique() < 2:
                warnings.append(f"Standard curve for {gene} has <3 points or <2 dilutions; defaulting E=2.0.")
                continue
            log_q = np.log10(valid["_qty"].to_numpy(copy=True, dtype=float))
            cq_v = valid["_cq"].to_numpy(copy=True, dtype=float)
            reg = stats.linregress(log_q, cq_v)
            slope = float(reg.slope)
            r2 = float(reg.rvalue ** 2)
            if slope >= -1e-6:
                warnings.append(f"Non-negative standard curve slope ({slope:.3f}) for {gene}; defaulting E=2.0.")
                continue
            eff = float(10.0 ** (-1.0 / slope))
            eff_pct = (eff - 1.0) * 100.0
            eff_map[gene] = eff
            sc_details[gene] = {
                "efficiency": round(eff, 4),
                "efficiency_pct": round(eff_pct, 2),
                "slope": round(slope, 4),
                "r2": round(r2, 4),
                "n_points": int(len(valid)),
                "miqe_efficiency_ok": bool(1.90 <= eff <= 2.10 and r2 >= 0.98),
            }
            if eff < 1.90 or eff > 2.10:
                warnings.append(
                    f"MIQE efficiency warning for {gene}: E={eff:.3f} ({eff_pct:.1f}%) outside 90%-110% [1.90, 2.10]."
                )
            if r2 < 0.98:
                warnings.append(f"MIQE standard curve linearity warning for {gene}: R^2={r2:.4f} < 0.98.")

    for spec in eff_overrides:
        if "=" not in spec:
            fail(f"Invalid --efficiencies entry '{spec}'; expected GENE=VALUE (e.g. GAPDH=1.98)")
        g_raw, v_raw = spec.split("=", 1)
        g_name = g_raw.strip()
        try:
            eff_val = float(v_raw.strip())
        except ValueError:
            fail(f"Non-numeric efficiency value in '{spec}'")
        if eff_val <= 1.0 or eff_val > 3.0:
            fail(f"Invalid amplification efficiency {eff_val} for {g_name}; expected in (1.0, 3.0]")
        eff_map[g_name] = eff_val
        sc_details[g_name] = {
            "efficiency": round(eff_val, 4),
            "efficiency_pct": round((eff_val - 1.0) * 100.0, 2),
            "source": "manual_override",
            "miqe_efficiency_ok": bool(1.90 <= eff_val <= 2.10),
        }
        if eff_val < 1.90 or eff_val > 2.10:
            warnings.append(
                f"MIQE efficiency warning for {g_name}: override E={eff_val:.3f} outside 90%-110% [1.90, 2.10]."
            )

    return eff_map, sc_details


def _check_ntc_wells(
    ntc_df: pd.DataFrame,
    sample_df: pd.DataFrame,
    warnings: list[str],
) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    if ntc_df.empty:
        return checks
    all_genes = sorted(set(sample_df["_gene"].unique()) | set(ntc_df["_gene"].unique()))
    for gene in all_genes:
        n_sub = ntc_df[ntc_df["_gene"] == gene]
        if n_sub.empty:
            continue
        s_sub = sample_df[(sample_df["_gene"] == gene) & np.isfinite(sample_df["_cq"])]
        min_sample_cq = float(s_sub["_cq"].min()) if not s_sub.empty else float("nan")
        finite_ntc = n_sub["_cq"].to_numpy(copy=True, dtype=float)
        finite_ntc = finite_ntc[np.isfinite(finite_ntc)]
        if finite_ntc.size == 0:
            checks.append(
                {
                    "gene": gene,
                    "n_ntc_wells": int(len(n_sub)),
                    "min_ntc_cq": None,
                    "min_sample_cq": min_sample_cq if math.isfinite(min_sample_cq) else None,
                    "delta_cq_ntc_minus_sample": None,
                    "contaminated": False,
                    "status": "clean_undetermined",
                }
            )
            continue
        min_ntc_cq = float(np.min(finite_ntc))
        delta = (min_ntc_cq - min_sample_cq) if math.isfinite(min_sample_cq) else float("nan")
        contaminated = bool(min_ntc_cq < 35.0 or (math.isfinite(delta) and delta < 5.0))
        status = "contaminated" if contaminated else "trace_late_amplification"
        if contaminated:
            warnings.append(
                f"NTC contamination flagged for {gene}: min NTC Cq={min_ntc_cq:.2f}"
                + (f" (delta vs sample={delta:.2f} < 5.0)" if math.isfinite(delta) else " (< 35.0)")
            )
        checks.append(
            {
                "gene": gene,
                "n_ntc_wells": int(len(n_sub)),
                "min_ntc_cq": round(min_ntc_cq, 3),
                "min_sample_cq": round(min_sample_cq, 3) if math.isfinite(min_sample_cq) else None,
                "delta_cq_ntc_minus_sample": round(delta, 3) if math.isfinite(delta) else None,
                "contaminated": contaminated,
                "status": status,
            }
        )
    return checks


def _compute_reference_stability(
    sample_meta: pd.DataFrame,
    ref_pivot: pd.DataFrame,
    ref_genes: list[str],
    eff_map: dict[str, float],
    control_group: str,
    warnings: list[str],
) -> pd.DataFrame:
    """Compute Vandesompele geNorm M_j and clean-room Andersen (2004) rho_j stability."""
    samples = ref_pivot.index.tolist()
    groups = sample_meta.set_index("_sample").loc[samples, "_group"].to_numpy(copy=True)
    n_samples = len(samples)
    n_refs = len(ref_genes)

    # Convert Cq to log2 relative quantity centered on control_group mean
    ctrl_mask = groups == control_group
    y_mat = np.zeros((n_refs, n_samples), dtype=float)
    for j, gene in enumerate(ref_genes):
        cqs = ref_pivot[gene].to_numpy(copy=True, dtype=float)
        ctrl_mean = float(np.mean(cqs[ctrl_mask]))
        eff = float(eff_map.get(gene, 2.0))
        y_mat[j, :] = (ctrl_mean - cqs) * math.log2(eff)

    # Vandesompele et al. (2002) geNorm M_j
    genorm_m = np.full(n_refs, np.nan, dtype=float)
    if n_refs >= 2:
        for j in range(n_refs):
            pairwise_sd = []
            for k in range(n_refs):
                if j == k:
                    continue
                diff = y_mat[j, :] - y_mat[k, :]
                pairwise_sd.append(float(np.std(diff, ddof=1)) if n_samples > 1 else 0.0)
            genorm_m[j] = float(np.mean(pairwise_sd))

    # Andersen et al. (2004) NormFinder intra- and inter-group variance stability model
    unique_groups = list(dict.fromkeys(groups.tolist()))
    n_groups = len(unique_groups)
    andersen_rho = np.full(n_refs, np.nan, dtype=float)
    intra_sd_arr = np.full(n_refs, np.nan, dtype=float)
    inter_dev_arr = np.full(n_refs, np.nan, dtype=float)

    if n_refs >= 2 and n_groups >= 2:
        # Mean across genes per sample and grand mean
        sample_panel_mean = np.mean(y_mat, axis=0)  # (S,)
        gene_grand_mean = np.mean(y_mat, axis=1)  # (R,)
        grand_mean = float(np.mean(y_mat))

        d_jg = np.zeros((n_refs, n_groups), dtype=float)
        s2_jg = np.zeros((n_refs, n_groups), dtype=float)
        n_g = np.zeros(n_groups, dtype=float)

        for g_idx, grp in enumerate(unique_groups):
            idx = np.where(groups == grp)[0]
            ng = len(idx)
            n_g[g_idx] = float(ng)
            sub_y = y_mat[:, idx]  # (R, ng)
            grp_gene_mean = np.mean(sub_y, axis=1)
            grp_panel_mean = float(np.mean(sample_panel_mean[idx]))
            d_jg[:, g_idx] = (grp_gene_mean - grp_panel_mean) - (gene_grand_mean - grand_mean)

            if ng > 1:
                resid = (sub_y - grp_gene_mean[:, None]) - (sample_panel_mean[idx][None, :] - grp_panel_mean)
                # Correction factor (1 - 1/R) for panel centering
                denom = max(1e-6, (ng - 1.0) * (1.0 - 1.0 / n_refs))
                s2_jg[:, g_idx] = np.maximum(np.sum(resid ** 2, axis=1) / denom, 1e-8)
            else:
                s2_jg[:, g_idx] = 1e-4

        ms_inter = float(np.sum(d_jg ** 2) / max(1, (n_refs - 1) * (n_groups - 1)))
        mean_samp_err = float(np.mean(s2_jg / np.maximum(n_g[None, :], 1.0)))
        tau2 = max(0.0, ms_inter - mean_samp_err)

        for j in range(n_refs):
            se2_jg = s2_jg[j, :] / np.maximum(n_g, 1.0)
            if tau2 > 0.0:
                shrink = tau2 / (tau2 + se2_jg)
                bias_jg = shrink * d_jg[j, :]
                post_var_jg = shrink * se2_jg
            else:
                bias_jg = np.zeros(n_groups, dtype=float)
                post_var_jg = np.zeros(n_groups, dtype=float)
            rho_j = float(np.mean(np.abs(bias_jg) + np.sqrt(post_var_jg + se2_jg)))
            andersen_rho[j] = rho_j
            intra_sd_arr[j] = float(math.sqrt(np.mean(s2_jg[j, :])))
            inter_dev_arr[j] = float(math.sqrt(np.mean(d_jg[j, :] ** 2)))
    elif n_refs >= 2:
        # Single group fallback: intra-group residual SD
        sample_panel_mean = np.mean(y_mat, axis=0)
        for j in range(n_refs):
            resid = (y_mat[j, :] - np.mean(y_mat[j, :])) - (sample_panel_mean - np.mean(sample_panel_mean))
            sd_j = float(np.std(resid, ddof=1)) if n_samples > 1 else 0.0
            andersen_rho[j] = sd_j
            intra_sd_arr[j] = sd_j
            inter_dev_arr[j] = 0.0

    rows = []
    for j, gene in enumerate(ref_genes):
        m_val = float(genorm_m[j]) if math.isfinite(genorm_m[j]) else None
        r_val = float(andersen_rho[j]) if math.isfinite(andersen_rho[j]) else None
        if m_val is not None and m_val > 1.0:
            warnings.append(f"Unstable reference gene {gene}: geNorm M={m_val:.3f} > 1.0.")
        rows.append(
            {
                "gene": gene,
                "efficiency": round(float(eff_map.get(gene, 2.0)), 4),
                "genorm_m": round(m_val, 4) if m_val is not None else None,
                "andersen_rho": round(r_val, 4) if r_val is not None else None,
                "intra_group_sd": round(float(intra_sd_arr[j]), 4) if math.isfinite(intra_sd_arr[j]) else None,
                "inter_group_rms_dev": round(float(inter_dev_arr[j]), 4) if math.isfinite(inter_dev_arr[j]) else None,
                "miqe_stability_tier": (
                    "single_ref"
                    if m_val is None
                    else ("optimal_M<0.5" if m_val < 0.5 else ("acceptable_M<1.0" if m_val < 1.0 else "unstable_M>=1.0"))
                ),
            }
        )

    stab_df = pd.DataFrame(rows)
    if n_refs >= 2:
        stab_df = stab_df.sort_values(["genorm_m", "andersen_rho"], ascending=[True, True], na_position="last").reset_index(drop=True)
    stab_df["stability_rank"] = np.arange(1, len(stab_df) + 1)
    return stab_df


def _run_group_statistics(
    quant_df: pd.DataFrame,
    target_genes: list[str],
    control_group: str,
) -> pd.DataFrame:
    """Test biological-replicate log2 fold changes vs control_group (never pseudoreplicate wells)."""
    rows: list[dict[str, Any]] = []
    all_groups = list(dict.fromkeys(quant_df["group"].astype(str).tolist()))
    test_groups = [g for g in all_groups if g != control_group]

    for gene in target_genes:
        sub = quant_df[quant_df["gene"] == gene]
        ctrl_vals = sub[sub["group"] == control_group]["log2_fold_change"].to_numpy(copy=True, dtype=float)
        ctrl_vals = ctrl_vals[np.isfinite(ctrl_vals)]
        ctrl_fc = sub[sub["group"] == control_group]["pfaffl_fold_change"].to_numpy(copy=True, dtype=float)
        ctrl_fc = ctrl_fc[np.isfinite(ctrl_fc)]

        group_arrays: list[np.ndarray] = []
        group_names: list[str] = []
        for grp in all_groups:
            vals = sub[sub["group"] == grp]["log2_fold_change"].to_numpy(copy=True, dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size > 0:
                group_arrays.append(vals)
                group_names.append(grp)

        anova_p = float("nan")
        tukey_res = None
        if len(group_arrays) >= 3 and all(arr.size >= 2 for arr in group_arrays):
            try:
                _, anova_p = stats.f_oneway(*group_arrays)
                anova_p = float(anova_p)
            except Exception:
                anova_p = float("nan")
            try:
                tukey_res = stats.tukey_hsd(*group_arrays)
            except Exception:
                tukey_res = None

        ctrl_idx = group_names.index(control_group) if control_group in group_names else -1

        for grp in test_groups:
            grp_vals = sub[sub["group"] == grp]["log2_fold_change"].to_numpy(copy=True, dtype=float)
            grp_vals = grp_vals[np.isfinite(grp_vals)]
            grp_fc = sub[sub["group"] == grp]["pfaffl_fold_change"].to_numpy(copy=True, dtype=float)
            grp_fc = grp_fc[np.isfinite(grp_fc)]

            if grp_vals.size < 2 or ctrl_vals.size < 2:
                welch_t, welch_p = float("nan"), float("nan")
            else:
                res_t = stats.ttest_ind(grp_vals, ctrl_vals, equal_var=False)
                welch_t = float(res_t.statistic)
                welch_p = float(res_t.pvalue)

            posthoc_p = welch_p
            test_method = "welch_t"
            if len(group_arrays) >= 3:
                test_method = "anova_tukey_hsd"
                if tukey_res is not None and ctrl_idx >= 0 and grp in group_names:
                    g_idx = group_names.index(grp)
                    posthoc_p = float(tukey_res.pvalue[g_idx, ctrl_idx])

            mean_log2_fc = float(np.mean(grp_vals) - np.mean(ctrl_vals)) if (grp_vals.size and ctrl_vals.size) else float("nan")
            mean_fc = float(np.mean(grp_fc)) if grp_fc.size else float("nan")
            sd_fc = float(np.std(grp_fc, ddof=1)) if grp_fc.size > 1 else 0.0
            ctrl_mean_fc = float(np.mean(ctrl_fc)) if ctrl_fc.size else float("nan")

            rows.append(
                {
                    "gene": gene,
                    "group": grp,
                    "control_group": control_group,
                    "n_biol_group": int(grp_vals.size),
                    "n_biol_control": int(ctrl_vals.size),
                    "control_mean_fc": round(ctrl_mean_fc, 4) if math.isfinite(ctrl_mean_fc) else None,
                    "group_mean_fc": round(mean_fc, 4) if math.isfinite(mean_fc) else None,
                    "group_sd_fc": round(sd_fc, 4),
                    "mean_log2_fc": round(mean_log2_fc, 4) if math.isfinite(mean_log2_fc) else None,
                    "test_method": test_method,
                    "anova_p": round(anova_p, 6) if math.isfinite(anova_p) else None,
                    "welch_t": round(welch_t, 4) if math.isfinite(welch_t) else None,
                    "welch_p": welch_p if math.isfinite(welch_p) else None,
                    "p_value": posthoc_p if math.isfinite(posthoc_p) else None,
                }
            )

    comp_df = pd.DataFrame(rows)
    if not comp_df.empty:
        pvals = comp_df["p_value"].fillna(1.0).to_numpy(copy=True, dtype=float)
        comp_df["bh_fdr"] = bh_fdr(pvals)
        comp_df["significant_fdr_05"] = comp_df["bh_fdr"] < 0.05
    return comp_df


def _plot_summary(
    quant_df: pd.DataFrame,
    stab_df: pd.DataFrame,
    tech_df: pd.DataFrame,
    target_genes: list[str],
    control_group: str,
    out_png: Path,
) -> None:
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.8), gridspec_kw={"width_ratios": [1.35, 1.0]})

    groups = [control_group] + [g for g in dict.fromkeys(quant_df["group"].tolist()) if g != control_group]
    palette = ["#64748b", "#0f766e", "#0284c7", "#d97706", "#dc2626", "#7c3aed"]
    x_base = np.arange(len(target_genes), dtype=float)
    width = min(0.75 / max(1, len(groups)), 0.28)
    rng = np.random.default_rng(42)

    for g_idx, grp in enumerate(groups):
        offset = (g_idx - (len(groups) - 1) / 2.0) * width
        color = palette[g_idx % len(palette)]
        means = []
        sds = []
        for t_idx, gene in enumerate(target_genes):
            vals = quant_df[(quant_df["gene"] == gene) & (quant_df["group"] == grp)]["pfaffl_fold_change"].to_numpy(
                copy=True, dtype=float
            )
            vals = vals[np.isfinite(vals)]
            means.append(float(np.mean(vals)) if vals.size else 0.0)
            sds.append(float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0)
            if vals.size:
                jitter = rng.uniform(-width * 0.18, width * 0.18, size=vals.size)
                ax1.scatter(
                    np.full(vals.size, x_base[t_idx] + offset) + jitter,
                    vals,
                    color="#0f172a",
                    s=18,
                    zorder=4,
                    alpha=0.8,
                )
        ax1.bar(
            x_base + offset,
            means,
            width=width * 0.9,
            yerr=sds,
            capsize=3,
            color=color,
            alpha=0.78,
            label=grp,
            zorder=2,
        )

    ax1.axhline(1.0, color="#475569", linestyle="--", linewidth=1.0, alpha=0.7)
    ax1.set_xticks(x_base)
    ax1.set_xticklabels(target_genes, rotation=20 if len(target_genes) > 4 else 0)
    ax1.set_ylabel("Relative Expression (Pfaffl Fold Change)")
    ax1.set_title("A. Target Gene Expression (Biological Replicates)")
    ax1.legend(frameon=False, title="Group")

    # Panel B: Reference stability (geNorm M & Andersen rho) + Technical SD QC summary
    ref_names = stab_df["gene"].astype(str).tolist()
    y_pos = np.arange(len(ref_names), dtype=float)
    m_vals = stab_df["genorm_m"].fillna(0.0).to_numpy(copy=True, dtype=float)
    r_vals = stab_df["andersen_rho"].fillna(0.0).to_numpy(copy=True, dtype=float)

    ax2.barh(y_pos - 0.17, m_vals, height=0.32, color="#0f766e", alpha=0.85, label="geNorm M")
    ax2.barh(y_pos + 0.17, r_vals, height=0.32, color="#0284c7", alpha=0.85, label="Andersen rho")
    ax2.axvline(0.5, color="#16a34a", linestyle=":", linewidth=1.1, label="M=0.5 optimal")
    ax2.axvline(1.0, color="#dc2626", linestyle="--", linewidth=1.1, label="M=1.0 unstable")
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(ref_names)
    ax2.invert_yaxis()
    ax2.set_xlabel("Stability Metric (Lower = More Stable)")
    high_sd = int((tech_df["sd_cq"] > 0.5).sum())
    ax2.set_title(f"B. Reference Stability & Tech QC (SD>0.5: {high_sd}/{len(tech_df)})")
    ax2.legend(frameon=False, loc="lower right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    ref_genes = parse_str_list(args.reference_genes)
    if not ref_genes:
        fail("At least one reference gene must be supplied via --reference_genes")
    eff_overrides = parse_str_list(args.efficiencies)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = read_table(args.cq_table)
    s_col = resolve_column(df, ["sample", "sample_id", "bio_rep", "replicate"], "sample")
    g_col = resolve_column(df, ["group", "condition", "treatment"], "group")
    t_col = resolve_column(df, ["gene", "target", "assay"], "gene")
    c_col = resolve_column(df, ["cq", "ct", "cp"], "cq")

    work = df[[s_col, g_col, t_col, c_col]].copy()
    work.columns = ["_sample", "_group", "_gene", "_cq_raw"]
    work["_sample"] = work["_sample"].astype(str).str.strip()
    work["_group"] = work["_group"].astype(str).str.strip()
    work["_gene"] = work["_gene"].astype(str).str.strip()
    work["_cq"] = pd.to_numeric(work["_cq_raw"], errors="coerce")

    ntc_tok = args.ntc_group.strip().lower()
    is_ntc = (work["_group"].str.lower() == ntc_tok) | (work["_sample"].str.lower() == ntc_tok)
    ntc_df = work[is_ntc].copy()
    sample_df = work[~is_ntc].copy()

    if sample_df.empty:
        fail("No non-NTC biological sample rows found in cq_table.")

    warnings: list[str] = []
    ntc_checks = _check_ntc_wells(ntc_df, sample_df, warnings)

    n_undetermined = int((~np.isfinite(sample_df["_cq"])).sum())
    if n_undetermined > 0:
        warnings.append(f"Dropped {n_undetermined} non-NTC wells with Undetermined/NaN Cq values.")
    sample_df = sample_df[np.isfinite(sample_df["_cq"])].copy()
    if sample_df.empty:
        fail("No finite Cq values remain in biological sample rows.")

    # Verify sample -> group uniqueness
    grp_per_sample = sample_df.groupby("_sample")["_group"].nunique()
    ambig = grp_per_sample[grp_per_sample > 1]
    if not ambig.empty:
        fail(f"Sample(s) assigned to multiple groups: {ambig.index.tolist()[:5]}")

    present_genes = list(dict.fromkeys(sample_df["_gene"].tolist()))
    missing_refs = [r for r in ref_genes if r not in present_genes]
    if missing_refs:
        fail(f"Reference gene(s) not found in cq_table: {missing_refs}")
    target_genes = [g for g in present_genes if g not in ref_genes]
    if not target_genes:
        fail("No target genes remain after excluding reference_genes.")

    if args.control_group not in set(sample_df["_group"].unique()):
        fail(f"control_group '{args.control_group}' not found in groups: {sorted(sample_df['_group'].unique())}")

    eff_map, sc_details = _estimate_efficiencies(present_genes, args.standard_curve, eff_overrides, warnings)

    # Aggregate technical replicates per (sample, group, gene)
    tech_df = (
        sample_df.groupby(["_sample", "_group", "_gene"], sort=False)["_cq"]
        .agg(n_tech="count", mean_cq="mean", sd_cq=lambda x: float(np.std(x.to_numpy(copy=True, dtype=float), ddof=1)) if len(x) > 1 else 0.0)
        .reset_index()
    )
    tech_df.rename(columns={"_sample": "sample", "_group": "group", "_gene": "gene"}, inplace=True)

    high_sd = tech_df[tech_df["sd_cq"] > 0.5]
    high_sd_records = []
    for _, row in high_sd.iterrows():
        rec = {
            "sample": row["sample"],
            "group": row["group"],
            "gene": row["gene"],
            "n_tech": int(row["n_tech"]),
            "mean_cq": round(float(row["mean_cq"]), 3),
            "sd_cq": round(float(row["sd_cq"]), 3),
        }
        high_sd_records.append(rec)
        warnings.append(
            f"MIQE technical replicate precision warning: {row['sample']} / {row['gene']} has sd_cq={row['sd_cq']:.3f} > 0.5."
        )

    pivot = tech_df.pivot(index="sample", columns="gene", values="mean_cq")
    sample_meta = sample_df[["_sample", "_group"]].drop_duplicates().copy()
    valid_samples = pivot[ref_genes].dropna().index.tolist()
    if not valid_samples:
        fail("No samples have finite Cq values across all specified reference_genes.")

    ref_pivot = pivot.loc[valid_samples, ref_genes]
    sample_meta_valid = sample_meta[sample_meta["_sample"].isin(valid_samples)].copy()
    stab_df = _compute_reference_stability(
        sample_meta_valid, ref_pivot, ref_genes, eff_map, args.control_group, warnings
    )

    # Calibrator mean Cq in control_group for every gene
    ctrl_samples = set(sample_meta_valid[sample_meta_valid["_group"] == args.control_group]["_sample"].tolist())
    ctrl_mean_cq: dict[str, float] = {}
    for g in present_genes:
        vals = pivot.loc[pivot.index.isin(ctrl_samples), g].dropna().to_numpy(copy=True, dtype=float)
        if vals.size == 0:
            fail(f"Control group '{args.control_group}' has no valid measurements for gene '{g}'.")
        ctrl_mean_cq[g] = float(np.mean(vals))

    # Normalization factor per biological replicate (geometric mean of E_r^(ctrl_mean_r - Cq_r,s))
    grp_lookup = dict(zip(sample_meta_valid["_sample"], sample_meta_valid["_group"]))
    nf_pfaffl: dict[str, float] = {}
    mean_ref_dcq: dict[str, float] = {}
    for s in valid_samples:
        log2_nf = 0.0
        dcq_sum = 0.0
        for r in ref_genes:
            dcq_r = ctrl_mean_cq[r] - float(ref_pivot.loc[s, r])
            log2_nf += dcq_r * math.log2(float(eff_map[r]))
            dcq_sum += dcq_r
        nf_pfaffl[s] = 2.0 ** (log2_nf / len(ref_genes))
        mean_ref_dcq[s] = dcq_sum / len(ref_genes)

    # Compute sample-level Pfaffl & 2^-ddCq quantities
    quant_rows = []
    tech_lookup = tech_df.set_index(["sample", "gene"])
    for s in valid_samples:
        grp = grp_lookup[s]
        for t in target_genes:
            if (s, t) not in tech_lookup.index:
                continue
            row_t = tech_lookup.loc[(s, t)]
            mean_cq_t = float(row_t["mean_cq"])
            sd_cq_t = float(row_t["sd_cq"])
            n_tech_t = int(row_t["n_tech"])
            eff_t = float(eff_map[t])
            dcq_t = ctrl_mean_cq[t] - mean_cq_t
            pfaffl_fc = (eff_t ** dcq_t) / nf_pfaffl[s]
            ddcq = -(dcq_t - mean_ref_dcq[s])  # standard ddCq = (Cq_t - Cq_ref)_sample - (Cq_t - Cq_ref)_ctrl
            ddcq_fc = 2.0 ** (-ddcq)
            quant_rows.append(
                {
                    "sample": s,
                    "group": grp,
                    "gene": t,
                    "n_tech": n_tech_t,
                    "mean_cq": round(mean_cq_t, 4),
                    "sd_cq": round(sd_cq_t, 4),
                    "high_sd_flag": bool(sd_cq_t > 0.5),
                    "efficiency": round(eff_t, 4),
                    "norm_factor_geom": round(nf_pfaffl[s], 4),
                    "ddcq": round(ddcq, 4),
                    "ddcq_fold_change": round(ddcq_fc, 6),
                    "pfaffl_fold_change": round(pfaffl_fc, 6),
                    "log2_fold_change": round(math.log2(max(pfaffl_fc, 1e-12)), 6),
                }
            )

    quant_df = pd.DataFrame(quant_rows)
    if quant_df.empty:
        fail("No valid target-gene measurements computed across samples.")

    comp_df = _run_group_statistics(quant_df, target_genes, args.control_group)

    inputs = [args.cq_table] + ([args.standard_curve] if args.standard_curve else [])
    params = {
        "reference_genes": ref_genes,
        "control_group": args.control_group,
        "ntc_group": args.ntc_group,
    }
    p_quant = write_tsv(
        quant_df,
        out_dir / "qpcr_sample_quantities.tsv",
        tool="wa_qpcr_quant",
        title="RT-qPCR Biological Replicate Normalized Quantities (Pfaffl & 2^-ddCq)",
        inputs=inputs,
        parameters=params,
    )
    p_comp = write_tsv(
        comp_df,
        out_dir / "qpcr_group_comparisons.tsv",
        tool="wa_qpcr_quant",
        title="RT-qPCR Group Statistical Comparisons vs Control (Biological Replicates)",
        inputs=inputs,
        parameters=params,
    )
    p_stab = write_tsv(
        stab_df,
        out_dir / "qpcr_reference_stability.tsv",
        tool="wa_qpcr_quant",
        title="RT-qPCR Reference Gene Stability Ranking (geNorm M & Andersen rho)",
        inputs=inputs,
        parameters=params,
    )
    fig_path = out_dir / "qpcr_expression_summary.png"
    _plot_summary(quant_df, stab_df, tech_df, target_genes, args.control_group, fig_path)
    p_fig = write_artifact(
        fig_path,
        tool="wa_qpcr_quant",
        title="RT-qPCR Relative Expression & MIQE Reference Stability Summary",
        inputs=inputs,
        parameters=params,
    )

    emit(
        {
            "status": "ok",
            "tool": "wa_qpcr_quant",
            "control_group": args.control_group,
            "n_samples": len(valid_samples),
            "n_target_genes": len(target_genes),
            "n_reference_genes": len(ref_genes),
            "efficiencies": {g: round(float(eff_map[g]), 4) for g in present_genes},
            "standard_curve_qc": sc_details,
            "reference_stability": stab_df.to_dict(orient="records"),
            "comparisons": comp_df.to_dict(orient="records"),
            "qc": {
                "n_tech_groups": int(len(tech_df)),
                "n_high_sd_replicates": len(high_sd_records),
                "high_sd_replicates": high_sd_records,
                "ntc_checks": ntc_checks,
                "n_contaminated_ntc_genes": int(sum(1 for c in ntc_checks if c["contaminated"])),
            },
            "outputs": {
                "qpcr_sample_quantities": p_quant,
                "qpcr_group_comparisons": p_comp,
                "qpcr_reference_stability": p_stab,
                "qpcr_expression_summary": p_fig,
            },
            "artifacts": [p_quant, p_comp, p_stab, p_fig],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
