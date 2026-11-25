#!/usr/bin/env python3
"""mr_estimate: Two-Sample Mendelian Randomization suite (IVW random/fixed, Wald
ratio, MR-Egger + intercept pleiotropy, Weighted Median, Simple & Weighted Mode,
Cochran's Q, outlier & Leave-One-Out sensitivity, Steiger directionality, and
4-panel diagnostic plots).
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
from _common import emit, fail, write_artifact


def ivw_fit(bx: np.ndarray, by: np.ndarray, sey: np.ndarray, random_effects: bool = True) -> dict:
    w = (bx * bx) / (sey * sey)
    ratio = by / bx
    b_ivw = float(np.sum(w * ratio) / np.sum(w))
    se_fe = float(1.0 / np.sqrt(np.sum(w)))
    k = len(bx)
    if k == 1:
        pval = float(2.0 * stats.norm.sf(abs(b_ivw / se_fe)))
        return {"beta": b_ivw, "se": se_fe, "pvalue": pval, "q": 0.0, "q_df": 0, "q_pval": None, "i2": 0.0}
    q = float(np.sum(w * (ratio - b_ivw) ** 2))
    df = k - 1
    q_pval = float(stats.chi2.sf(q, df=df))
    i2 = max(0.0, 100.0 * (q - df) / max(q, 1e-12))
    sigma = max(1.0, np.sqrt(q / df)) if random_effects else 1.0
    se = se_fe * float(sigma)
    pval = float(2.0 * stats.norm.sf(abs(b_ivw / se)))
    return {"beta": b_ivw, "se": se, "pvalue": pval, "q": q, "q_df": df, "q_pval": q_pval, "i2": i2}


def egger_fit(bx: np.ndarray, by: np.ndarray, sey: np.ndarray) -> dict:
    # Orient bx > 0
    sgn = np.where(bx >= 0, 1.0, -1.0)
    bx_o = bx * sgn
    by_o = by * sgn
    w = 1.0 / (sey * sey)
    X = np.column_stack([np.ones_like(bx_o), bx_o])
    XtW = X.T * w
    cov_unscaled = np.linalg.pinv(XtW @ X)
    theta = cov_unscaled @ (XtW @ by_o)
    resid = by_o - X @ theta
    df = len(bx) - 2
    q_prime = float(np.sum(w * resid * resid))
    q_pval = float(stats.chi2.sf(q_prime, df=df))
    sigma2 = max(1.0, q_prime / df)
    se_vec = np.sqrt(np.diag(cov_unscaled) * sigma2)
    b_int, b_slope = float(theta[0]), float(theta[1])
    se_int, se_slope = float(se_vec[0]), float(se_vec[1])
    p_int = float(2.0 * stats.t.sf(abs(b_int / max(se_int, 1e-15)), df=df))
    p_slope = float(2.0 * stats.t.sf(abs(b_slope / max(se_slope, 1e-15)), df=df))
    return {
        "beta": b_slope,
        "se": se_slope,
        "pvalue": p_slope,
        "intercept": b_int,
        "intercept_se": se_int,
        "intercept_pvalue": p_int,
        "q_prime": q_prime,
        "q_df": df,
        "q_pvalue": q_pval,
    }


def weighted_median_point(ratios: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(ratios)
    r_s = ratios[order]
    w_s = weights[order] / np.sum(weights)
    cum = np.cumsum(w_s) - 0.5 * w_s
    return float(np.interp(0.5, cum, r_s))


def mode_point(ratios: np.ndarray, weights: np.ndarray) -> float:
    sd = float(np.std(ratios, ddof=1)) if len(ratios) > 1 else 1.0
    iqr = float(np.subtract(*np.percentile(ratios, [75, 25])))
    mad = float(np.median(np.abs(ratios - np.median(ratios)))) * 1.4826
    s_scale = min(sd, iqr / 1.349 if iqr > 1e-9 else sd)
    if s_scale < 1e-9:
        s_scale = max(mad, 1e-3)
    bw = 0.9 * s_scale / (len(ratios) ** 0.2)
    grid = np.linspace(float(ratios.min()) - 2 * bw, float(ratios.max()) + 2 * bw, 512)
    # Vectorized Gaussian kernel density with normalized weights
    w_norm = weights / np.sum(weights)
    diffs = (grid[:, None] - ratios[None, :]) / max(bw, 1e-9)
    dens = np.exp(-0.5 * diffs * diffs) @ w_norm
    return float(grid[int(np.argmax(dens))])


def fmt_method_row(name: str, n_snp: int, beta: float, se: float, pval: float, binary: bool) -> dict:
    ci_lo = beta - 1.959964 * se
    ci_hi = beta + 1.959964 * se
    row = {
        "method": name,
        "nsnp": n_snp,
        "beta": round(beta, 5),
        "se": round(se, 5),
        "ci_low": round(ci_lo, 5),
        "ci_high": round(ci_hi, 5),
        "pvalue": float(pval),
    }
    if binary:
        row["or"] = round(float(np.exp(beta)), 4)
        row["or_ci_low"] = round(float(np.exp(ci_lo)), 4)
        row["or_ci_high"] = round(float(np.exp(ci_hi)), 4)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized", required=True)
    parser.add_argument("--exposure_name", default="Exposure")
    parser.add_argument("--outcome_name", default="Outcome")
    parser.add_argument(
        "--outcome_binary",
        nargs="?",
        const=True,
        type=lambda x: str(x).lower() not in ("false", "0", "no"),
        default=False,
    )
    parser.add_argument("--n_boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.harmonized, sep="," if args.harmonized.lower().endswith(".csv") else "\t")
    req = {"SNP", "beta_exp", "se_exp", "beta_out", "se_out"}
    if not req.issubset(set(df.columns)):
        fail(f"harmonized table missing required columns {req}; found {list(df.columns)}")
    if df.empty:
        fail("harmonized table has 0 SNPs")

    snps = df["SNP"].astype(str).tolist()
    bx = df["beta_exp"].to_numpy(dtype=float)
    sex = df["se_exp"].to_numpy(dtype=float)
    by = df["beta_out"].to_numpy(dtype=float)
    sey = df["se_out"].to_numpy(dtype=float)
    k = len(snps)

    methods: list[dict] = []
    if k == 1:
        w_fit = ivw_fit(bx, by, sey, random_effects=False)
        methods.append(fmt_method_row("Wald ratio", 1, w_fit["beta"], w_fit["se"], w_fit["pvalue"], args.outcome_binary))
        het_info = {"cochran_q": None, "q_df": 0, "q_pvalue": None, "i2_pct": None}
        pleio_info = {"egger_intercept": None, "intercept_se": None, "intercept_pvalue": None}
        b_ivw_primary = w_fit["beta"]
        b_egger, int_egger = None, 0.0
        b_wm = None
    else:
        ivw_re = ivw_fit(bx, by, sey, random_effects=True)
        ivw_fe = ivw_fit(bx, by, sey, random_effects=False)
        methods.append(fmt_method_row("Inverse variance weighted (MRE)", k, ivw_re["beta"], ivw_re["se"], ivw_re["pvalue"], args.outcome_binary))
        methods.append(fmt_method_row("Inverse variance weighted (FE)", k, ivw_fe["beta"], ivw_fe["se"], ivw_fe["pvalue"], args.outcome_binary))
        b_ivw_primary = ivw_re["beta"]
        het_info = {
            "cochran_q": round(ivw_re["q"], 4),
            "q_df": ivw_re["q_df"],
            "q_pvalue": ivw_re["q_pval"],
            "i2_pct": round(ivw_re["i2"], 2),
        }

        if k >= 3:
            eg = egger_fit(bx, by, sey)
            methods.append(fmt_method_row("MR Egger", k, eg["beta"], eg["se"], eg["pvalue"], args.outcome_binary))
            pleio_info = {
                "egger_intercept": round(eg["intercept"], 5),
                "intercept_se": round(eg["intercept_se"], 5),
                "intercept_pvalue": eg["intercept_pvalue"],
            }
            b_egger, int_egger = eg["beta"], eg["intercept"]

            # Parametric bootstrap for Weighted Median, Simple Mode, Weighted Mode
            ratios = by / bx
            se_ratio = np.sqrt((sey / bx) ** 2 + ((by * sex) / (bx * bx)) ** 2)
            w_iv = 1.0 / (se_ratio * se_ratio)

            b_wm = weighted_median_point(ratios, w_iv)
            b_sm = mode_point(ratios, np.ones_like(w_iv))
            b_wmode = mode_point(ratios, w_iv)

            rng = np.random.default_rng(args.seed)
            bx_boot = rng.normal(bx[None, :], sex[None, :], size=(args.n_boot, k))
            bx_boot = np.where(np.abs(bx_boot) < 1e-6, np.sign(bx)[None, :] * 1e-6, bx_boot)
            by_boot = rng.normal(by[None, :], sey[None, :], size=(args.n_boot, k))
            r_boot = by_boot / bx_boot
            w_boot = (bx_boot * bx_boot) / (sey[None, :] ** 2)

            wm_b = np.array([weighted_median_point(r_boot[i], w_boot[i]) for i in range(args.n_boot)])
            sm_b = np.array([mode_point(r_boot[i], np.ones(k)) for i in range(args.n_boot)])
            wmo_b = np.array([mode_point(r_boot[i], w_boot[i]) for i in range(args.n_boot)])

            for lbl, pt, boots in (
                ("Weighted median", b_wm, wm_b),
                ("Weighted mode", b_wmode, wmo_b),
                ("Simple mode", b_sm, sm_b),
            ):
                se_b = max(float(np.std(boots, ddof=1)), 1e-9)
                pv_b = float(2.0 * stats.norm.sf(abs(pt / se_b)))
                methods.append(fmt_method_row(lbl, k, pt, se_b, pv_b, args.outcome_binary))
        else:
            pleio_info = {"egger_intercept": None, "intercept_se": None, "intercept_pvalue": None}
            b_egger, int_egger, b_wm = None, 0.0, None

    # Per-SNP Wald ratios, residual outlier check, and Leave-One-Out
    ratios_obs = by / bx
    se_ratios_obs = sey / np.abs(bx)
    outlier_snps: list[str] = []
    loo_rows = []
    for idx in range(k):
        p_single = float(2.0 * stats.norm.sf(abs(ratios_obs[idx] / se_ratios_obs[idx])))
        q_i = ((ratios_obs[idx] - b_ivw_primary) / se_ratios_obs[idx]) ** 2
        p_out = float(stats.chi2.sf(q_i, df=1))
        is_out = bool(k >= 3 and p_out < (0.05 / k))
        if is_out:
            outlier_snps.append(snps[idx])

        if k >= 2:
            mask = np.arange(k) != idx
            loo_fit = ivw_fit(bx[mask], by[mask], sey[mask], random_effects=(k - 1 >= 2))
            loo_b, loo_se, loo_p = loo_fit["beta"], loo_fit["se"], loo_fit["pvalue"]
        else:
            loo_b, loo_se, loo_p = np.nan, np.nan, np.nan

        loo_rows.append(
            {
                "SNP": snps[idx],
                "wald_beta": round(float(ratios_obs[idx]), 5),
                "wald_se": round(float(se_ratios_obs[idx]), 5),
                "wald_pvalue": p_single,
                "loo_ivw_beta": round(float(loo_b), 5) if np.isfinite(loo_b) else None,
                "loo_ivw_se": round(float(loo_se), 5) if np.isfinite(loo_se) else None,
                "loo_ivw_pvalue": float(loo_p) if np.isfinite(loo_p) else None,
                "bonferroni_outlier": is_out,
            }
        )

    # Global Steiger directionality test if r2_exp & r2_out present
    steiger_summary = None
    if {"r2_exp", "r2_out", "n_exp", "n_out"}.issubset(set(df.columns)) and df["r2_exp"].notna().all() and df["r2_out"].notna().all():
        r2_x = float(df["r2_exp"].sum())
        r2_y = float(df["r2_out"].sum())
        nx = float(df["n_exp"].median())
        ny = float(df["n_out"].median())
        rx = np.sqrt(min(0.9999, max(0.0, r2_x)))
        ry = np.sqrt(min(0.9999, max(0.0, r2_y)))
        zx, zy = np.arctanh(rx), np.arctanh(ry)
        se_diff = np.sqrt(1.0 / max(1.0, nx - 3.0) + 1.0 / max(1.0, ny - 3.0))
        z_st = (zx - zy) / se_diff
        steiger_summary = {
            "total_r2_exposure": round(r2_x, 5),
            "total_r2_outcome": round(r2_y, 5),
            "correct_causal_direction": bool(r2_x > r2_y),
            "steiger_pvalue": float(2.0 * stats.norm.sf(abs(z_st))),
        }

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    res_df = pd.DataFrame(methods)
    res_path = outdir / "mr_results.tsv"
    res_df.to_csv(res_path, sep="\t", index=False)

    loo_df = pd.DataFrame(loo_rows)
    loo_path = outdir / "mr_singlesnp_loo.tsv"
    loo_df.to_csv(loo_path, sep="\t", index=False)

    # 4-Panel Diagnostic Figure (Scatter, Method Forest, Leave-One-Out, Funnel)
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 9.6))
    ax_sc, ax_fo = axes[0, 0], axes[0, 1]
    ax_lo, ax_fu = axes[1, 0], axes[1, 1]

    # (A) Scatter plot
    ax_sc.errorbar(bx, by, xerr=sex, yerr=sey, fmt="o", color="#0F172A", ecolor="#94A3B8", elinewidth=0.9, markersize=4.5, alpha=0.85)
    x_line = np.linspace(0.0, float(np.max(bx)) * 1.12, 100)
    ax_sc.plot(x_line, b_ivw_primary * x_line, color="#2563EB", linewidth=1.8, label=f"IVW (β={b_ivw_primary:.3f})")
    if b_wm is not None:
        ax_sc.plot(x_line, b_wm * x_line, color="#10B981", linewidth=1.6, linestyle="--", label=f"Weighted Median (β={b_wm:.3f})")
    if b_egger is not None:
        ax_sc.plot(x_line, int_egger + b_egger * x_line, color="#E11D48", linewidth=1.6, linestyle=":", label=f"MR-Egger (β={b_egger:.3f})")
    ax_sc.axhline(0.0, color="#CBD5E1", linewidth=0.8)
    ax_sc.set_xlabel(f"SNP effect on {args.exposure_name}")
    ax_sc.set_ylabel(f"SNP effect on {args.outcome_name}")
    ax_sc.set_title("MR Effect Scatter Plot", fontsize=10)
    ax_sc.legend(frameon=False, fontsize=8, loc="best")
    ax_sc.spines["top"].set_visible(False)
    ax_sc.spines["right"].set_visible(False)

    # (B) Method comparison forest plot
    m_rev = methods[::-1]
    ys = np.arange(len(m_rev))
    ax_fo.axvline(0.0, color="#64748B", linestyle="--", linewidth=0.9)
    for i, m in enumerate(m_rev):
        col = "#2563EB" if m["pvalue"] < 0.05 else "#64748B"
        ax_fo.plot([m["ci_low"], m["ci_high"]], [i, i], color=col, linewidth=2.0)
        ax_fo.scatter([m["beta"]], [i], color=col, s=42, zorder=4)
    ax_fo.set_yticks(ys)
    ax_fo.set_yticklabels([m["method"] for m in m_rev], fontsize=8.5)
    ax_fo.set_xlabel("Causal Effect β (95% CI)")
    ax_fo.set_title("MR Estimator Triangulation", fontsize=10)
    ax_fo.spines["top"].set_visible(False)
    ax_fo.spines["right"].set_visible(False)

    # (C) Leave-one-out plot (top 25 if >25 SNPs)
    if k >= 2:
        show_loo = loo_rows[:25][::-1]
        yl = np.arange(len(show_loo))
        ax_lo.axvline(b_ivw_primary, color="#E11D48", linestyle=":", linewidth=1.1, label="All-SNP IVW")
        ax_lo.axvline(0.0, color="#94A3B8", linestyle="--", linewidth=0.8)
        for i, lr in enumerate(show_loo):
            lb, lse = lr["loo_ivw_beta"], lr["loo_ivw_se"]
            if lb is not None and lse is not None:
                ax_lo.plot([lb - 1.96 * lse, lb + 1.96 * lse], [i, i], color="#0F172A", linewidth=1.4)
                ax_lo.scatter([lb], [i], color="#0F172A", s=24)
        ax_lo.set_yticks(yl)
        ax_lo.set_yticklabels([f"Drop {r['SNP']}" for r in show_loo], fontsize=7.2)
        ax_lo.set_xlabel("Leave-One-Out IVW β (95% CI)")
        ax_lo.set_title("Leave-One-Out Sensitivity", fontsize=10)
    else:
        ax_lo.text(0.5, 0.5, "Single-SNP Wald Ratio\n(LOO requires ≥ 2 SNPs)", ha="center", va="center")
    ax_lo.spines["top"].set_visible(False)
    ax_lo.spines["right"].set_visible(False)

    # (D) Funnel plot (precision 1/SE_IV vs Wald ratio)
    prec = 1.0 / np.maximum(se_ratios_obs, 1e-12)
    ax_fu.scatter(ratios_obs, prec, s=28, color="#0F172A", alpha=0.8)
    ax_fu.axvline(b_ivw_primary, color="#2563EB", linewidth=1.5, label="IVW")
    if b_egger is not None:
        ax_fu.axvline(b_egger, color="#E11D48", linestyle=":", linewidth=1.4, label="MR-Egger")
    ax_fu.set_xlabel("Per-SNP Wald Ratio (β_Y / β_X)")
    ax_fu.set_ylabel("Instrument Precision (1 / SE_IV)")
    ax_fu.set_title("Funnel Asymmetry Plot", fontsize=10)
    ax_fu.legend(frameon=False, fontsize=8, loc="best")
    ax_fu.spines["top"].set_visible(False)
    ax_fu.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = outdir / "mr-diagnostics.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {"outcome_binary": args.outcome_binary, "n_boot": args.n_boot, "seed": args.seed}
    artifacts = [
        write_artifact(fig_path, "mr_estimate", f"MR 4-panel diagnostics ({args.exposure_name} -> {args.outcome_name})", inputs=(args.harmonized,), parameters=params_meta),
        write_artifact(res_path, "mr_estimate", "Multi-method Mendelian randomization summary table", inputs=(args.harmonized,), parameters=params_meta, shape=res_df.shape),
        write_artifact(loo_path, "mr_estimate", "Single-SNP Wald ratios and Leave-One-Out table", inputs=(args.harmonized,), parameters=params_meta, shape=loo_df.shape),
    ]

    warnings: list[str] = []
    if pleio_info["intercept_pvalue"] is not None and pleio_info["intercept_pvalue"] < 0.05:
        warnings.append(
            f"MR-Egger intercept p = {pleio_info['intercept_pvalue']:.3g} (< 0.05): significant directional horizontal pleiotropy detected"
        )
    if het_info["q_pvalue"] is not None and het_info["q_pvalue"] < 0.05:
        warnings.append(
            f"Cochran's Q heterogeneity p = {het_info['q_pvalue']:.3g} (I² = {het_info['i2_pct']}%): rely on multiplicative random-effects IVW & Weighted Median and inspect outlier_snps"
        )
    if outlier_snps:
        warnings.append(f"Bonferroni residual outlier SNP(s) detected: {outlier_snps}")
    if steiger_summary and not steiger_summary["correct_causal_direction"]:
        warnings.append("Global Steiger test indicates instruments explain more variance in outcome than exposure (potential reverse causation)")

    emit(
        {
            "tool": "mr_estimate",
            "exposure": args.exposure_name,
            "outcome": args.outcome_name,
            "n_snps": k,
            "methods": methods,
            "pleiotropy": pleio_info,
            "heterogeneity": het_info,
            "steiger_directionality": steiger_summary,
            "outlier_snps": outlier_snps,
            "outputs": {
                "results_table": str(res_path),
                "singlesnp_loo_table": str(loo_path),
            },
            "figures": [str(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
