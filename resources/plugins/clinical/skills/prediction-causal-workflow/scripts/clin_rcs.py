#!/usr/bin/env python3
"""clin_rcs: Harrell's Restricted Cubic Spline (RCS) regression (3, 4, or 5 knots)
for linear, logistic, and Cox PH models with Overall & Non-Linearity Wald chi2 tests
and delta-method 95% CI dose-response curve.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pred_common import (
    emit,
    encode_covariates,
    fail,
    fit_cox_weighted,
    fit_logistic_weighted,
    parse_binary_series,
    read_table,
    write_artifact,
)

HARRELL_PERCENTILES = {
    3: [10.0, 50.0, 90.0],
    4: [5.0, 35.0, 65.0, 95.0],
    5: [5.0, 27.5, 50.0, 72.5, 95.0],
}


def rcs_basis(x: np.ndarray, knots: np.ndarray) -> np.ndarray:
    """Harrell's restricted cubic spline basis (K - 1 columns)."""
    x_arr = np.asarray(x, dtype=float)
    k = len(knots)
    t1 = float(knots[0])
    tk = float(knots[-1])
    tk1 = float(knots[-2])
    denom = (tk - t1) ** 2
    if denom <= 1e-18 or abs(tk - tk1) <= 1e-12:
        fail("degenerate spline knots (exposure variance too small or heavily tied)")

    cols = [x_arr]
    for j in range(k - 2):
        tj = float(knots[j])
        term = (
            np.maximum(x_arr - tj, 0.0) ** 3
            - ((tk - tj) / (tk - tk1)) * (np.maximum(x_arr - tk1, 0.0) ** 3)
            + ((tk1 - tj) / (tk - tk1)) * (np.maximum(x_arr - tk, 0.0) ** 3)
        ) / denom
        cols.append(term)
    return np.column_stack(cols)


def fit_ols(X_no_int: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    n = len(y)
    X = np.column_stack([np.ones(n, dtype=float), X_no_int])
    p = X.shape[1]
    XtX = X.T @ X
    inv_xtx = np.linalg.pinv(XtX + 1e-12 * np.eye(p))
    beta = inv_xtx @ (X.T @ y)
    resid = y - X @ beta
    df_resid = max(1, n - p)
    sigma2 = float(np.sum(resid ** 2) / df_resid)
    cov = sigma2 * inv_xtx
    return {"beta": beta, "cov": cov}


def wald_test(beta_sub: np.ndarray, cov_sub: np.ndarray) -> tuple[float, int, float]:
    df = len(beta_sub)
    inv_c = np.linalg.pinv(cov_sub + 1e-12 * np.eye(df))
    chi2 = float(max(0.0, beta_sub.T @ inv_c @ beta_sub))
    pval = float(stats.chi2.sf(chi2, df=df))
    return chi2, df, pval


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--exposure_col", required=True)
    parser.add_argument("--outcome_col", required=True)
    parser.add_argument("--model_type", choices=["auto", "logistic", "cox", "linear"], default="auto")
    parser.add_argument("--time_col", default=None)
    parser.add_argument("--covariates", action="extend", nargs="+", default=[])
    parser.add_argument("--n_knots", type=int, choices=[3, 4, 5], default=4)
    parser.add_argument("--ref_value", type=float, default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.clinical)
    for req in (args.exposure_col, args.outcome_col):
        if req not in df.columns:
            fail(f"column '{req}' not found in clinical table: {list(df.columns)}")

    covs = [c.strip() for item in args.covariates for c in str(item).split(",") if c.strip() and c.strip() != args.exposure_col]
    for c in covs:
        if c not in df.columns:
            fail(f"covariate '{c}' not found in clinical table: {list(df.columns)}")

    exp_s = pd.to_numeric(df[args.exposure_col], errors="coerce")
    valid = exp_s.notna()

    mtype = args.model_type
    if mtype == "auto":
        if args.time_col:
            mtype = "cox"
        else:
            s_out = df.loc[valid, args.outcome_col].dropna()
            num_out = pd.to_numeric(s_out, errors="coerce")
            if s_out.astype(str).str.strip().nunique() == 2:
                mtype = "logistic"
            elif num_out.notna().sum() >= 0.85 * len(s_out):
                mtype = "linear"
            else:
                fail(f"cannot infer model_type for outcome '{args.outcome_col}'; specify --model_type")

    if mtype == "cox":
        if not args.time_col or args.time_col not in df.columns:
            fail("time_col is required and must exist in clinical table when model_type is cox")
        tm_s = pd.to_numeric(df[args.time_col], errors="coerce")
        ev_s, _ = parse_binary_series(df[args.outcome_col], args.outcome_col)
        valid = valid & tm_s.notna() & (tm_s > 0) & ev_s.notna()
    elif mtype == "logistic":
        y_s, _ = parse_binary_series(df[args.outcome_col], args.outcome_col)
        valid = valid & y_s.notna()
    else:
        y_s = pd.to_numeric(df[args.outcome_col], errors="coerce")
        valid = valid & y_s.notna()

    cov_df = None
    if covs:
        cov_df, _, _ = encode_covariates(df, covs)
        valid = valid & cov_df.notna().all(axis=1)

    warnings: list[str] = []
    n_drop = int((~valid).sum())
    if n_drop > 0:
        warnings.append(f"excluded {n_drop} record(s) with missing exposure, outcome, or covariates")

    x = exp_s.loc[valid].to_numpy(dtype=float, copy=True)
    n = len(x)
    if n < 20:
        fail(f"RCS requires >= 20 complete observations (got {n})")

    pcts = HARRELL_PERCENTILES[args.n_knots]
    knots = np.percentile(x, pcts)
    # Handle tied quantile knots gracefully
    uniq_k = np.unique(np.round(knots, 8))
    if len(uniq_k) < 3:
        fail("exposure has too few distinct values for restricted cubic spline knots")
    if len(uniq_k) < len(knots):
        warnings.append(f"exposure ties reduced effective knots from {len(knots)} to {len(uniq_k)}")
        knots = uniq_k

    k_eff = len(knots)
    B = rcs_basis(x, knots)
    n_spline_cols = B.shape[1]

    if cov_df is not None:
        X_cov = cov_df.loc[valid].to_numpy(dtype=float, copy=True)
        X_design = np.column_stack([B, X_cov])
    else:
        X_design = B

    n_events = None
    if mtype == "linear":
        y_arr = y_s.loc[valid].to_numpy(dtype=float, copy=True)
        fit = fit_ols(X_design, y_arr)
        beta_sp = fit["beta"][1 : 1 + n_spline_cols]
        cov_sp = fit["cov"][1 : 1 + n_spline_cols, 1 : 1 + n_spline_cols]
        effect_measure = "delta_y"
        ref_line = 0.0
    elif mtype == "logistic":
        y_arr = y_s.loc[valid].to_numpy(dtype=float, copy=True)
        n_events = int(np.sum(y_arr == 1))
        if n_events < 6 or (n - n_events) < 6:
            fail(f"logistic RCS needs >= 6 positives and >= 6 negatives (got {n_events}, {n - n_events})")
        fit = fit_logistic_weighted(X_design, y_arr, l2_penalty=0.0)
        beta_sp = fit["beta"][1 : 1 + n_spline_cols]
        cov_sp = fit["cov_model"][1 : 1 + n_spline_cols, 1 : 1 + n_spline_cols]
        effect_measure = "OR"
        ref_line = 1.0
    else:
        tm_arr = tm_s.loc[valid].to_numpy(dtype=float, copy=True)
        ev_arr = ev_s.loc[valid].to_numpy(dtype=int, copy=True)
        n_events = int(np.sum(ev_arr == 1))
        if n_events < 6:
            fail(f"Cox RCS needs >= 6 events (got {n_events})")
        fit = fit_cox_weighted(X_design, tm_arr, ev_arr)
        beta_sp = fit["beta"][:n_spline_cols]
        cov_sp = fit["cov_model"][:n_spline_cols, :n_spline_cols]
        effect_measure = "HR"
        ref_line = 1.0

    chi2_all, df_all, p_all = wald_test(beta_sp, cov_sp)
    chi2_nl, df_nl, p_nl = wald_test(beta_sp[1:], cov_sp[1:, 1:])

    ref_val = float(args.ref_value) if args.ref_value is not None else float(np.median(x))
    x_lo = float(np.percentile(x, 1.0))
    x_hi = float(np.percentile(x, 99.0))
    if ref_val < x_lo:
        x_lo = ref_val
    if ref_val > x_hi:
        x_hi = ref_val

    x_grid = np.linspace(x_lo, x_hi, 200)
    B_grid = rcs_basis(x_grid, knots)
    B_ref = rcs_basis(np.array([ref_val], dtype=float), knots)[0]
    dB = B_grid - B_ref[None, :]

    eta_diff = dB @ beta_sp
    se_diff = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", dB, cov_sp, dB), 0.0))

    if mtype == "linear":
        est = eta_diff
        ci_lo = eta_diff - 1.959964 * se_diff
        ci_hi = eta_diff + 1.959964 * se_diff
    else:
        est = np.exp(np.clip(eta_diff, -15.0, 15.0))
        ci_lo = np.exp(np.clip(eta_diff - 1.959964 * se_diff, -15.0, 15.0))
        ci_hi = np.exp(np.clip(eta_diff + 1.959964 * se_diff, -15.0, 15.0))

    # Inflection / extremum check on interior P5..P95 when p_nonlinear < 0.05
    inflection = None
    if p_nl < 0.05:
        interior = (x_grid >= knots[0]) & (x_grid <= knots[-1])
        xg_in = x_grid[interior]
        est_in = est[interior]
        if len(est_in) >= 5:
            min_i = int(np.argmin(est_in))
            max_i = int(np.argmax(est_in))
            if 0 < min_i < len(est_in) - 1:
                inflection = {"type": "nadir", "exposure": round(float(xg_in[min_i]), 4), "effect": round(float(est_in[min_i]), 4)}
            elif 0 < max_i < len(est_in) - 1:
                inflection = {"type": "peak", "exposure": round(float(xg_in[max_i]), 4), "effect": round(float(est_in[max_i]), 4)}

    curve_df = pd.DataFrame(
        {
            "exposure": np.round(x_grid, 5),
            "effect_measure": effect_measure,
            "estimate": np.round(est, 5),
            "ci_low": np.round(ci_lo, 5),
            "ci_high": np.round(ci_hi, 5),
            "se_linear_pred": np.round(se_diff, 5),
        }
    )
    curve_path = os.path.join(args.output_dir, "rcs_curve.tsv")
    curve_df.to_csv(curve_path, sep="\t", index=False)

    # Dose-response figure with background histogram
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax_hist = ax.twinx()
    x_in = x[(x >= x_lo) & (x <= x_hi)]
    ax_hist.hist(x_in, bins=28, color="#CBD5E1", alpha=0.45, edgecolor="none")
    ax_hist.set_ylim(0, ax_hist.get_ylim()[1] * 3.4)
    ax_hist.set_yticks([])
    ax_hist.spines["top"].set_visible(False)
    ax_hist.spines["right"].set_visible(False)

    ax.set_zorder(ax_hist.get_zorder() + 1)
    ax.patch.set_visible(False)
    ax.axhline(ref_line, color="#475569", linestyle="--", linewidth=1.1)
    ax.axvline(ref_val, color="#94A3B8", linestyle=":", linewidth=1.1, label=f"Ref = {ref_val:.3g}")
    ax.fill_between(x_grid, ci_lo, ci_hi, color="#2563EB", alpha=0.18, label="95% CI")
    ax.plot(x_grid, est, color="#2563EB", linewidth=2.3, label=f"RCS ({k_eff} knots)")
    if inflection:
        ax.scatter([inflection["exposure"]], [inflection["effect"]], color="#E11D48", s=55, zorder=5, label=f"{inflection['type']} @ {inflection['exposure']:.3g}")

    ax.set_xlabel(args.exposure_col)
    y_lbl = "Predicted Difference (Delta Y)" if mtype == "linear" else f"{effect_measure} (95% CI)"
    ax.set_ylabel(y_lbl)
    ax.set_title(
        f"Restricted Cubic Spline ({mtype.capitalize()}: {args.exposure_col} -> {args.outcome_col})\n"
        f"P_overall = {p_all:.3e}  |  P_non-linear = {p_nl:.3e}",
        fontsize=10,
    )
    ax.legend(frameon=False, loc="best", fontsize=8.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "rcs_dose_response.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {
        "exposure_col": args.exposure_col,
        "outcome_col": args.outcome_col,
        "model_type": mtype,
        "n_knots": k_eff,
        "ref_value": ref_val,
        "covariates": covs,
    }
    artifacts = [
        write_artifact(fig_path, "clin_rcs", f"RCS dose-response ({args.exposure_col})", inputs=(args.clinical,), parameters=params_meta),
        write_artifact(curve_path, "clin_rcs", "RCS dose-response curve estimates", inputs=(args.clinical,), parameters=params_meta, shape=curve_df.shape),
    ]

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "exposure_col": args.exposure_col,
            "outcome_col": args.outcome_col,
            "model_type": mtype,
            "effect_measure": effect_measure,
            "n_samples": n,
            "n_events": n_events,
            "n_knots": k_eff,
            "knots": [round(float(kv), 5) for kv in knots],
            "ref_value": round(ref_val, 5),
            "chisq_overall": round(chi2_all, 4),
            "df_overall": df_all,
            "p_overall": p_all,
            "chisq_nonlinear": round(chi2_nl, 4),
            "df_nonlinear": df_nl,
            "p_nonlinear": p_nl,
            "nonlinear_inflection": inflection,
            "curve_table": os.path.abspath(curve_path),
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
