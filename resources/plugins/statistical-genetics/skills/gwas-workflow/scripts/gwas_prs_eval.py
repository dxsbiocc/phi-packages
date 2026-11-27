#!/usr/bin/env python3
"""gwas_prs_eval: Standardize a Polygenic Risk Score (prs_z), compare Clinical
Covariates Only vs Covariates + PRS across binary (logistic), quantitative (OLS),
or time-to-event (Cox PH) outcomes, compute incremental delta-R2 / Nagelkerke /
Lee (2012) liability R2 / LRT P / delta-AUC / delta-C-index, and evaluate OR/HR
risk stratification across PRS quantiles relative to the Middle 40-60% reference.
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
from scipy import optimize, special, stats
from sklearn.metrics import roc_auc_score, roc_curve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    parse_str_list,
    read_table,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)


def _prepare_covariates(df: pd.DataFrame, cov_cols: list[str]) -> tuple[np.ndarray, list[str]]:
    """One-hot encode categorical covariates and standardize numeric columns."""
    if not cov_cols:
        return np.empty((len(df), 0), dtype=float), []
    sub = df[cov_cols].copy()
    num_parts: list[pd.DataFrame] = []
    for col in cov_cols:
        num_s = pd.to_numeric(sub[col], errors="coerce")
        if num_s.notna().sum() == len(sub):
            num_parts.append(num_s.to_frame(col))
        else:
            dummies = pd.get_dummies(sub[col].astype(str), prefix=col, drop_first=True, dtype=float)
            if not dummies.empty:
                num_parts.append(dummies)
    if not num_parts:
        return np.empty((len(df), 0), dtype=float), []
    cov_df = pd.concat(num_parts, axis=1)
    arr = cov_df.to_numpy(dtype=float, copy=True)
    sds = np.std(arr, axis=0, ddof=1)
    keep = sds > 1e-12
    arr = arr[:, keep]
    kept_names = [c for c, k in zip(cov_df.columns, keep) if k]
    if arr.shape[1] > 0:
        arr = (arr - np.mean(arr, axis=0)) / np.std(arr, axis=0, ddof=1)
    return arr, kept_names


def _fit_ols(X: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    """Fit OLS y = X @ beta with intercept as first column of X."""
    n, p = X.shape
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    rss = float(np.sum(resid * resid))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - (rss / tss) if tss > 0 else 0.0
    df_res = max(1, n - p)
    sigma2 = rss / df_res
    xtx_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(np.diag(xtx_inv) * sigma2, 1e-18))
    ll = -0.5 * n * (math.log(2.0 * math.pi * max(rss / n, 1e-18)) + 1.0)
    return {"beta": beta, "se": se, "r2": r2, "ll": ll, "pred": X @ beta}


def _fit_logistic(X: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    """Fit unpenalized logistic regression via damped Newton-Raphson (IRLS)."""
    n, p = X.shape
    beta = np.zeros(p, dtype=float)
    y_mean = float(np.clip(np.mean(y), 1e-6, 1.0 - 1e-6))
    beta[0] = math.log(y_mean / (1.0 - y_mean))

    def _ll(b_vec: np.ndarray) -> float:
        eta = np.clip(X @ b_vec, -30.0, 30.0)
        return float(np.sum(y * eta - np.logaddexp(0.0, eta)))

    cur_ll = _ll(beta)
    for _ in range(60):
        eta = np.clip(X @ beta, -30.0, 30.0)
        prob = special.expit(eta)
        w = np.maximum(prob * (1.0 - prob), 1e-9)
        grad = X.T @ (y - prob)
        hess = X.T @ (w[:, None] * X) + 1e-8 * np.eye(p)
        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(hess, grad, rcond=None)[0]
        if float(np.max(np.abs(step))) < 1e-8:
            break
        alpha = 1.0
        while alpha > 1e-4:
            cand = beta + alpha * step
            cand_ll = _ll(cand)
            if cand_ll >= cur_ll - 1e-12:
                beta = cand
                cur_ll = cand_ll
                break
            alpha *= 0.5

    eta = np.clip(X @ beta, -30.0, 30.0)
    prob = special.expit(eta)
    w = np.maximum(prob * (1.0 - prob), 1e-9)
    info = X.T @ (w[:, None] * X) + 1e-8 * np.eye(p)
    cov = np.linalg.pinv(info)
    se = np.sqrt(np.maximum(np.diag(cov), 1e-18))
    return {"beta": beta, "se": se, "ll": cur_ll, "prob": prob}


def _cox_breslow_eval(X_sorted: np.ndarray, e_sorted: np.ndarray, group_ends: np.ndarray, beta: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Compute Breslow log-partial-likelihood, gradient, and Fisher info on time-descending sorted data."""
    n, p = X_sorted.shape
    if p == 0:
        # Null model beta = 0
        eta = np.zeros(n, dtype=float)
    else:
        eta = np.clip(X_sorted @ beta, -30.0, 30.0)
    exp_eta = np.exp(eta)
    s0_all = np.cumsum(exp_eta)

    if p == 0:
        ll = 0.0
        prev_end = 0
        for g_end in group_ends:
            d_k = int(np.sum(e_sorted[prev_end:g_end]))
            if d_k > 0:
                ll += float(np.sum(eta[prev_end:g_end][e_sorted[prev_end:g_end] == 1])) - d_k * math.log(max(float(s0_all[g_end - 1]), 1e-300))
            prev_end = g_end
        return ll, np.zeros(0), np.zeros((0, 0))

    wx = exp_eta[:, None] * X_sorted
    s1_all = np.cumsum(wx, axis=0)
    # Accumulate outer products for Fisher information
    wxx = exp_eta[:, None, None] * (X_sorted[:, :, None] * X_sorted[:, None, :])
    s2_all = np.cumsum(wxx, axis=0)

    ll = 0.0
    grad = np.zeros(p, dtype=float)
    info = np.zeros((p, p), dtype=float)

    prev_end = 0
    for g_end in group_ends:
        ev_slice = e_sorted[prev_end:g_end] == 1
        d_k = int(np.sum(ev_slice))
        if d_k > 0:
            s0 = max(float(s0_all[g_end - 1]), 1e-300)
            s1 = s1_all[g_end - 1]
            s2 = s2_all[g_end - 1]
            x_bar = s1 / s0
            ll += float(np.sum(eta[prev_end:g_end][ev_slice])) - d_k * math.log(s0)
            grad += np.sum(X_sorted[prev_end:g_end][ev_slice], axis=0) - d_k * x_bar
            info += d_k * ((s2 / s0) - np.outer(x_bar, x_bar))
        prev_end = g_end

    return ll, grad, info


def _fit_cox(X: np.ndarray, time: np.ndarray, event: np.ndarray) -> dict[str, Any]:
    """Fit Cox Proportional Hazards model via Breslow partial likelihood Newton-Raphson."""
    order = np.argsort(-time, kind="mergesort")
    t_s = time[order]
    e_s = event[order].astype(int)
    X_s = X[order]

    # Identify end index (exclusive) of each tied-time block in descending time order
    change = np.where(np.diff(t_s) != 0)[0] + 1
    group_ends = np.append(change, len(t_s))

    ll_null, _, _ = _cox_breslow_eval(X_s[:, :0], e_s, group_ends, np.zeros(0))
    n, p = X_s.shape
    if p == 0:
        return {"beta": np.zeros(0), "se": np.zeros(0), "ll": ll_null, "ll_null": ll_null, "lp": np.zeros(n)}

    beta = np.zeros(p, dtype=float)
    cur_ll, grad, info = _cox_breslow_eval(X_s, e_s, group_ends, beta)
    for _ in range(50):
        hess = info + 1e-8 * np.eye(p)
        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(hess, grad, rcond=None)[0]
        if float(np.max(np.abs(step))) < 1e-8:
            break
        alpha = 1.0
        while alpha > 1e-4:
            cand = beta + alpha * step
            cand_ll, cand_g, cand_i = _cox_breslow_eval(X_s, e_s, group_ends, cand)
            if cand_ll >= cur_ll - 1e-12:
                beta = cand
                cur_ll, grad, info = cand_ll, cand_g, cand_i
                break
            alpha *= 0.5

    cov = np.linalg.pinv(info + 1e-8 * np.eye(p))
    se = np.sqrt(np.maximum(np.diag(cov), 1e-18))
    lp = X @ beta
    return {"beta": beta, "se": se, "ll": cur_ll, "ll_null": ll_null, "lp": lp}


def _harrell_c_index(time: np.ndarray, event: np.ndarray, risk_score: np.ndarray) -> float:
    """Compute Harrell's concordance index (C-index) efficiently in chunks."""
    ev_idx = np.where(event == 1)[0]
    if ev_idx.size == 0:
        return 0.5
    concordant = 0.0
    permissible = 0.0
    chunk = 500
    for start in range(0, len(ev_idx), chunk):
        i_idx = ev_idx[start : start + chunk]
        ti = time[i_idx][:, None]
        ri = risk_score[i_idx][:, None]
        # Compares t_i < t_j or tied time where j is censored
        comp = (ti < time[None, :]) | ((ti == time[None, :]) & (event[None, :] == 0))
        comp[np.arange(len(i_idx)), i_idx] = False
        n_comp = float(np.sum(comp))
        if n_comp == 0:
            continue
        permissible += n_comp
        concordant += float(np.sum(comp & (risk_score[None, :] < ri))) + 0.5 * float(np.sum(comp & (risk_score[None, :] == ri)))
    return float(concordant / permissible) if permissible > 0 else 0.5


def _nagelkerke_r2(ll_model: float, ll_null: float, n: int) -> float:
    cs = 1.0 - math.exp(min(0.0, (2.0 / n) * (ll_null - ll_model)))
    denom = 1.0 - math.exp(min(0.0, (2.0 / n) * ll_null))
    return float(cs / denom) if denom > 1e-12 else 0.0


def _lee2012_liability_r2(r2_obs: float, sample_prev: float, pop_prev: float) -> float:
    """Transform observed-scale binary R2 to liability scale via Lee et al. (2012) Eq. 15."""
    if r2_obs <= 0.0:
        return 0.0
    k = pop_prev
    p = sample_prev
    t = float(stats.norm.ppf(1.0 - k))
    z = float(stats.norm.pdf(t))
    m = z / k
    c_fac = (k * (1.0 - k) / (z * z)) * (k * (1.0 - k) / (p * (1.0 - p)))
    theta = m * ((p - k) / (1.0 - k)) * (m * ((p - k) / (1.0 - k)) - t)
    return float((r2_obs * c_fac) / (1.0 + r2_obs * theta * c_fac))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--prs_col", required=True)
    ap.add_argument("--outcome_col", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--covariates", nargs="+", action="extend", default=None)
    ap.add_argument("--time_col", default=None)
    ap.add_argument("--pop_prev", type=float, default=None)
    args = ap.parse_args()

    df = read_table(args.cohort)
    for req_col in [args.prs_col, args.outcome_col]:
        if req_col not in df.columns:
            fail(f"Required column '{req_col}' not found in cohort table; available: {list(df.columns)}")

    cov_req = parse_str_list(args.covariates)
    for c in cov_req:
        if c not in df.columns:
            fail(f"Requested covariate '{c}' not found in cohort table; available: {list(df.columns)}")
        if c in {args.prs_col, args.outcome_col, args.time_col}:
            fail(f"Covariate '{c}' overlaps with prs_col / outcome_col / time_col")

    keep_cols = [args.prs_col, args.outcome_col] + cov_req
    if args.time_col:
        if args.time_col not in df.columns:
            fail(f"Survival time column '{args.time_col}' not found in cohort table")
        keep_cols.append(args.time_col)

    work = df[keep_cols].copy()
    work[args.prs_col] = pd.to_numeric(work[args.prs_col], errors="coerce")
    work[args.outcome_col] = pd.to_numeric(work[args.outcome_col], errors="coerce")
    if args.time_col:
        work[args.time_col] = pd.to_numeric(work[args.time_col], errors="coerce")
    work = work.dropna().reset_index(drop=True)

    if args.time_col:
        work = work[work[args.time_col] > 0].reset_index(drop=True)
    n_samples = int(len(work))
    if n_samples < 30:
        fail(f"Too few complete samples ({n_samples}) after dropping missing values (need >= 30)")

    prs_raw = work[args.prs_col].to_numpy(dtype=float, copy=True)
    prs_sd = float(np.std(prs_raw, ddof=1))
    if prs_sd <= 1e-12:
        fail(f"PRS column '{args.prs_col}' has zero variance")
    prs_z = (prs_raw - float(np.mean(prs_raw))) / prs_sd

    y_raw = work[args.outcome_col].to_numpy(dtype=float, copy=True)
    uniq_y = set(np.unique(y_raw).tolist())

    if args.time_col:
        if not uniq_y.issubset({0.0, 1.0}) or len(uniq_y) < 2:
            fail(f"Survival mode requires binary 0/1 event indicator in '{args.outcome_col}'")
        mode = "survival"
        time_arr = work[args.time_col].to_numpy(dtype=float, copy=True)
        n_events = int(np.sum(y_raw == 1.0))
    elif uniq_y == {0.0, 1.0}:
        mode = "binary"
        time_arr = None
        n_events = int(np.sum(y_raw == 1.0))
    else:
        mode = "quantitative"
        time_arr = None
        n_events = None

    if args.pop_prev is not None and not (0.0 < args.pop_prev < 1.0):
        fail("--pop_prev must be strictly between 0 and 1")

    X_cov, cov_names_used = _prepare_covariates(work, cov_req)
    warnings: list[str] = []

    # Fit Covariates Only vs Covariates + PRS
    if mode == "survival":
        assert time_arr is not None
        base_fit = _fit_cox(X_cov, time_arr, y_raw)
        X_full = np.column_stack([X_cov, prs_z]) if X_cov.shape[1] > 0 else prs_z[:, None]
        full_fit = _fit_cox(X_full, time_arr, y_raw)

        beta_sd = float(full_fit["beta"][-1])
        se_sd = float(full_fit["se"][-1])
        z_sd = beta_sd / se_sd
        p_sd = float(max(2.0 * stats.norm.sf(abs(z_sd)), 1e-300))
        eff_measure = "HR"
        eff_est = math.exp(beta_sd)
        eff_lo = math.exp(beta_sd - 1.959963984540054 * se_sd)
        eff_hi = math.exp(beta_sd + 1.959963984540054 * se_sd)

        lrt_chi2 = max(0.0, 2.0 * (full_fit["ll"] - base_fit["ll"]))
        lrt_p = float(max(stats.chi2.sf(lrt_chi2, df=1), 1e-300))
        c_base = _harrell_c_index(time_arr, y_raw, base_fit["lp"]) if X_cov.shape[1] > 0 else 0.5
        c_full = _harrell_c_index(time_arr, y_raw, full_fit["lp"])
        r2_base = _nagelkerke_r2(base_fit["ll"], base_fit["ll_null"], n_samples)
        r2_full = _nagelkerke_r2(full_fit["ll"], full_fit["ll_null"], n_samples)

        model_comp: dict[str, Any] = {
            "metric_r2": "Nagelkerke_Cox_R2",
            "r2_base": round(r2_base, 6),
            "r2_full": round(r2_full, 6),
            "delta_r2": round(r2_full - r2_base, 6),
            "r2_liability_delta": None,
            "discrimination_metric": "Harrell_C_index",
            "discrimination_base": round(c_base, 5),
            "discrimination_full": round(c_full, 5),
            "delta_discrimination": round(c_full - c_base, 5),
            "lrt_chi2": round(lrt_chi2, 4),
            "lrt_p_value": lrt_p,
        }
        comp_rows = [
            {"model": "Covariates_Only", "n_predictors": X_cov.shape[1], "r2": r2_base, "discrimination": c_base, "log_likelihood": base_fit["ll"], "lrt_p_vs_base": np.nan},
            {"model": "Covariates_Plus_PRS", "n_predictors": X_full.shape[1], "r2": r2_full, "discrimination": c_full, "log_likelihood": full_fit["ll"], "lrt_p_vs_base": lrt_p},
        ]

    elif mode == "binary":
        ones = np.ones((n_samples, 1), dtype=float)
        X_null = ones
        X_base = np.column_stack([ones, X_cov]) if X_cov.shape[1] > 0 else ones
        X_full = np.column_stack([ones, X_cov, prs_z]) if X_cov.shape[1] > 0 else np.column_stack([ones, prs_z])

        null_fit = _fit_logistic(X_null, y_raw)
        base_fit = null_fit if X_cov.shape[1] == 0 else _fit_logistic(X_base, y_raw)
        full_fit = _fit_logistic(X_full, y_raw)

        beta_sd = float(full_fit["beta"][-1])
        se_sd = float(full_fit["se"][-1])
        z_sd = beta_sd / se_sd
        p_sd = float(max(2.0 * stats.norm.sf(abs(z_sd)), 1e-300))
        eff_measure = "OR"
        eff_est = math.exp(beta_sd)
        eff_lo = math.exp(beta_sd - 1.959963984540054 * se_sd)
        eff_hi = math.exp(beta_sd + 1.959963984540054 * se_sd)

        lrt_chi2 = max(0.0, 2.0 * (full_fit["ll"] - base_fit["ll"]))
        lrt_p = float(max(stats.chi2.sf(lrt_chi2, df=1), 1e-300))
        r2_base = _nagelkerke_r2(base_fit["ll"], null_fit["ll"], n_samples)
        r2_full = _nagelkerke_r2(full_fit["ll"], null_fit["ll"], n_samples)
        auc_base = float(roc_auc_score(y_raw, base_fit["prob"])) if X_cov.shape[1] > 0 else 0.5
        auc_full = float(roc_auc_score(y_raw, full_fit["prob"]))

        r2_liab_delta: float | None = None
        if args.pop_prev is not None:
            ols_base = _fit_ols(X_base, y_raw)
            ols_full = _fit_ols(X_full, y_raw)
            delta_r2_obs = max(0.0, float(ols_full["r2"] - ols_base["r2"]))
            samp_prev = float(np.mean(y_raw))
            r2_liab_delta = _lee2012_liability_r2(delta_r2_obs, samp_prev, float(args.pop_prev))

        model_comp = {
            "metric_r2": "Nagelkerke_R2",
            "r2_base": round(r2_base, 6),
            "r2_full": round(r2_full, 6),
            "delta_r2": round(r2_full - r2_base, 6),
            "r2_liability_delta": round(r2_liab_delta, 6) if r2_liab_delta is not None else None,
            "discrimination_metric": "AUC",
            "discrimination_base": round(auc_base, 5),
            "discrimination_full": round(auc_full, 5),
            "delta_discrimination": round(auc_full - auc_base, 5),
            "lrt_chi2": round(lrt_chi2, 4),
            "lrt_p_value": lrt_p,
        }
        comp_rows = [
            {"model": "Covariates_Only", "n_predictors": X_cov.shape[1], "r2": r2_base, "discrimination": auc_base, "log_likelihood": base_fit["ll"], "lrt_p_vs_base": np.nan},
            {"model": "Covariates_Plus_PRS", "n_predictors": X_cov.shape[1] + 1, "r2": r2_full, "discrimination": auc_full, "log_likelihood": full_fit["ll"], "lrt_p_vs_base": lrt_p},
        ]

    else:  # quantitative
        ones = np.ones((n_samples, 1), dtype=float)
        X_base = np.column_stack([ones, X_cov]) if X_cov.shape[1] > 0 else ones
        X_full = np.column_stack([ones, X_cov, prs_z]) if X_cov.shape[1] > 0 else np.column_stack([ones, prs_z])
        base_fit = _fit_ols(X_base, y_raw)
        full_fit = _fit_ols(X_full, y_raw)

        beta_sd = float(full_fit["beta"][-1])
        se_sd = float(full_fit["se"][-1])
        t_stat = beta_sd / se_sd
        p_sd = float(max(2.0 * stats.t.sf(abs(t_stat), df=max(1, n_samples - X_full.shape[1])), 1e-300))
        eff_measure = "BETA"
        eff_est = beta_sd
        eff_lo = beta_sd - 1.959963984540054 * se_sd
        eff_hi = beta_sd + 1.959963984540054 * se_sd

        lrt_chi2 = max(0.0, 2.0 * (full_fit["ll"] - base_fit["ll"]))
        lrt_p = float(max(stats.chi2.sf(lrt_chi2, df=1), 1e-300))
        r2_base = float(base_fit["r2"])
        r2_full = float(full_fit["r2"])

        model_comp = {
            "metric_r2": "OLS_R2",
            "r2_base": round(r2_base, 6),
            "r2_full": round(r2_full, 6),
            "delta_r2": round(r2_full - r2_base, 6),
            "r2_liability_delta": None,
            "discrimination_metric": "none",
            "discrimination_base": None,
            "discrimination_full": None,
            "delta_discrimination": None,
            "lrt_chi2": round(lrt_chi2, 4),
            "lrt_p_value": lrt_p,
        }
        comp_rows = [
            {"model": "Covariates_Only", "n_predictors": X_cov.shape[1], "r2": r2_base, "discrimination": np.nan, "log_likelihood": base_fit["ll"], "lrt_p_vs_base": np.nan},
            {"model": "Covariates_Plus_PRS", "n_predictors": X_cov.shape[1] + 1, "r2": r2_full, "discrimination": np.nan, "log_likelihood": full_fit["ll"], "lrt_p_vs_base": lrt_p},
        ]

    # Quantile risk stratification vs Middle 40-60% reference
    q10, q40, q60, q80, q90, q95 = np.quantile(prs_z, [0.10, 0.40, 0.60, 0.80, 0.90, 0.95])
    strata_defs = [
        ("Bottom 10%", prs_z < q10, False),
        ("Middle 40-60%", (prs_z >= q40) & (prs_z <= q60), True),
        ("Top 20%", prs_z >= q80, False),
        ("Top 10%", prs_z >= q90, False),
        ("Top 5%", prs_z >= q95, False),
    ]
    ref_mask = strata_defs[1][1]

    quant_rows: list[dict[str, Any]] = []
    for label, mask, is_ref in strata_defs:
        n_g = int(np.sum(mask))
        ev_g = int(np.sum(y_raw[mask] == 1.0)) if mode in {"binary", "survival"} else None
        rate_g = float(np.mean(y_raw[mask])) if n_g > 0 else np.nan
        if is_ref:
            est_val = 1.0 if mode in {"binary", "survival"} else 0.0
            ci_lo = est_val
            ci_hi = est_val
            pval = 1.0
        else:
            pair_mask = mask | ref_mask
            is_tgt = mask[pair_mask].astype(float)
            y_sub = y_raw[pair_mask]
            cov_sub = X_cov[pair_mask] if X_cov.shape[1] > 0 else np.empty((int(np.sum(pair_mask)), 0), dtype=float)
            if mode == "survival":
                assert time_arr is not None
                t_sub = time_arr[pair_mask]
                X_sub = np.column_stack([cov_sub, is_tgt]) if cov_sub.shape[1] > 0 else is_tgt[:, None]
                fit_s = _fit_cox(X_sub, t_sub, y_sub)
                b_s, se_s = float(fit_s["beta"][-1]), float(fit_s["se"][-1])
                est_val = math.exp(b_s)
                ci_lo = math.exp(b_s - 1.959963984540054 * se_s)
                ci_hi = math.exp(b_s + 1.959963984540054 * se_s)
                pval = float(max(2.0 * stats.norm.sf(abs(b_s / se_s)), 1e-300))
            elif mode == "binary":
                ones_s = np.ones((len(y_sub), 1), dtype=float)
                X_sub = np.column_stack([ones_s, cov_sub, is_tgt]) if cov_sub.shape[1] > 0 else np.column_stack([ones_s, is_tgt])
                fit_s = _fit_logistic(X_sub, y_sub)
                b_s, se_s = float(fit_s["beta"][-1]), float(fit_s["se"][-1])
                est_val = math.exp(b_s)
                ci_lo = math.exp(b_s - 1.959963984540054 * se_s)
                ci_hi = math.exp(b_s + 1.959963984540054 * se_s)
                pval = float(max(2.0 * stats.norm.sf(abs(b_s / se_s)), 1e-300))
            else:
                ones_s = np.ones((len(y_sub), 1), dtype=float)
                X_sub = np.column_stack([ones_s, cov_sub, is_tgt]) if cov_sub.shape[1] > 0 else np.column_stack([ones_s, is_tgt])
                fit_s = _fit_ols(X_sub, y_sub)
                b_s, se_s = float(fit_s["beta"][-1]), float(fit_s["se"][-1])
                est_val, ci_lo, ci_hi = b_s, b_s - 1.959963984540054 * se_s, b_s + 1.959963984540054 * se_s
                pval = float(max(2.0 * stats.norm.sf(abs(b_s / se_s)), 1e-300))

        quant_rows.append(
            {
                "stratum": label,
                "is_reference": is_ref,
                "n": n_g,
                "n_events": ev_g if ev_g is not None else np.nan,
                "outcome_mean_or_rate": round(rate_g, 5) if np.isfinite(rate_g) else np.nan,
                "effect_measure": eff_measure if mode in {"binary", "survival"} else "BETA_DIFF",
                "estimate": round(est_val, 4),
                "ci_lower": round(ci_lo, 4),
                "ci_upper": round(ci_hi, 4),
                "p_value": pval,
            }
        )

    quant_df = pd.DataFrame(quant_rows)
    comp_df = pd.DataFrame(comp_rows)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    quant_path = outdir / "prs_quantile_risk.tsv"
    comp_path = outdir / "prs_model_comparison.tsv"
    fig_path = outdir / "prs_stratification_plot.png"

    # 2-panel PRS stratification + discrimination/decile figure
    plt = setup_matplotlib()
    fig, (ax_q, ax_r) = plt.subplots(1, 2, figsize=(11.8, 4.6))

    y_pos = np.arange(len(quant_rows))[::-1]
    ests = [r["estimate"] for r in quant_rows]
    los = [max(0.0, r["estimate"] - r["ci_lower"]) for r in quant_rows]
    his = [max(0.0, r["ci_upper"] - r["estimate"]) for r in quant_rows]
    cols_q = ["#0284C7", "#64748B", "#F59E0B", "#EA580C", "#E11D48"]
    for yp, est_v, lo_v, hi_v, col_v in zip(y_pos, ests, los, his, cols_q):
        ax_q.errorbar(est_v, yp, xerr=[[lo_v], [hi_v]], fmt="o", color=col_v, ecolor=col_v, capsize=3.5, markersize=7)
    ref_line = 1.0 if mode in {"binary", "survival"} else 0.0
    ax_q.axvline(ref_line, color="#94A3B8", linestyle="--", linewidth=1.0)
    ax_q.set_yticks(y_pos)
    ax_q.set_yticklabels([r["stratum"] for r in quant_rows])
    ax_q.set_xlabel(f"{eff_measure} vs Middle 40–60% (95% CI)" if mode in {"binary", "survival"} else "Mean Difference vs Middle 40–60% (95% CI)")
    ax_q.set_title(f"PRS Risk Stratification ({eff_measure}/SD = {eff_est:.2f})", fontsize=10)

    if mode == "binary":
        if X_cov.shape[1] > 0:
            fpr_b, tpr_b, _ = roc_curve(y_raw, base_fit["prob"])
            ax_r.plot(fpr_b, tpr_b, color="#64748B", linewidth=1.6, label=f"Covariates Only (AUC={model_comp['discrimination_base']:.3f})")
        fpr_f, tpr_f, _ = roc_curve(y_raw, full_fit["prob"])
        ax_r.plot(fpr_f, tpr_f, color="#1D4ED8", linewidth=2.0, label=f"Covariates + PRS (AUC={model_comp['discrimination_full']:.3f})")
        ax_r.plot([0, 1], [0, 1], color="#CBD5E1", linestyle="--", linewidth=1.0)
        ax_r.set_xlabel("False Positive Rate")
        ax_r.set_ylabel("True Positive Rate")
        ax_r.set_title(f"ROC Comparison (ΔAUC = {model_comp['delta_discrimination']:+.3f})", fontsize=10)
        ax_r.legend(frameon=False, loc="lower right")
    elif mode == "survival":
        assert time_arr is not None
        # Plot Nelson-Aalen / empirical cumulative incidence curves for Bottom 10%, Middle 40-60%, Top 10%
        for lbl, m_s, col_s in [
            ("Bottom 10%", strata_defs[0][1], "#0284C7"),
            ("Middle 40-60%", strata_defs[1][1], "#64748B"),
            ("Top 10%", strata_defs[3][1], "#E11D48"),
        ]:
            t_sub = time_arr[m_s]
            e_sub = y_raw[m_s]
            ord_s = np.argsort(t_sub)
            t_o = t_sub[ord_s]
            e_o = e_sub[ord_s]
            n_at_risk = np.arange(len(t_o), 0, -1, dtype=float)
            surv = np.cumprod(1.0 - e_o / np.maximum(n_at_risk, 1.0))
            ax_r.step(np.concatenate([[0.0], t_o]), np.concatenate([[0.0], 1.0 - surv]), where="post", color=col_s, linewidth=1.8, label=lbl)
        ax_r.set_xlabel("Follow-Up Time")
        ax_r.set_ylabel("Cumulative Event Incidence")
        ax_r.set_title(f"Survival Stratification (ΔC = {model_comp['delta_discrimination']:+.3f})", fontsize=10)
        ax_r.legend(frameon=False, loc="upper left")
    else:
        dec_edges = np.quantile(prs_z, np.linspace(0.0, 1.0, 11))
        dec_means, dec_se = [], []
        for di in range(10):
            if di < 9:
                dm = (prs_z >= dec_edges[di]) & (prs_z < dec_edges[di + 1])
            else:
                dm = (prs_z >= dec_edges[di]) & (prs_z <= dec_edges[di + 1])
            dec_means.append(float(np.mean(y_raw[dm])))
            dec_se.append(float(np.std(y_raw[dm], ddof=1) / math.sqrt(max(1, np.sum(dm)))))
        ax_r.errorbar(np.arange(1, 11), dec_means, yerr=dec_se, fmt="o-", color="#1D4ED8", capsize=3)
        ax_r.set_xlabel("PRS Decile (1 = lowest, 10 = highest)")
        ax_r.set_ylabel(f"Mean {args.outcome_col}")
        ax_r.set_title(f"PRS Decile Trend (ΔR² = {model_comp['delta_r2']:.4f})", fontsize=10)

    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    params_meta = {
        "prs_col": args.prs_col,
        "outcome_col": args.outcome_col,
        "covariates": cov_names_used,
        "time_col": args.time_col,
        "pop_prev": args.pop_prev,
        "outcome_mode": mode,
    }
    artifacts = [
        write_tsv(quant_df, quant_path, "gwas_prs_eval", "PRS quantile risk stratification vs Middle 40-60%", inputs=(args.cohort,), parameters=params_meta),
        write_tsv(comp_df, comp_path, "gwas_prs_eval", "PRS nested model comparison (Covariates Only vs +PRS)", inputs=(args.cohort,), parameters=params_meta),
        write_artifact(fig_path, "gwas_prs_eval", "PRS quantile risk stratification and discrimination plot", inputs=(args.cohort,), parameters=params_meta),
    ]

    emit(
        {
            "tool": "gwas_prs_eval",
            "outcome_mode": mode,
            "n_samples": n_samples,
            "n_events": n_events,
            "prs_col": args.prs_col,
            "outcome_col": args.outcome_col,
            "time_col": args.time_col,
            "pop_prev": args.pop_prev,
            "covariates_used": cov_names_used,
            "prs_effect_measure": eff_measure,
            "prs_effect_per_sd": round(eff_est, 4),
            "prs_ci_lower": round(eff_lo, 4),
            "prs_ci_upper": round(eff_hi, 4),
            "prs_p_value": p_sd,
            "model_comparison": model_comp,
            "quantile_risk": quant_rows,
            "outputs": {
                "prs_quantile_risk": str(quant_path),
                "prs_model_comparison": str(comp_path),
            },
            "figures": [str(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
