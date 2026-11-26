#!/usr/bin/env python3
"""clin_psm_iptw: Logistic propensity score positivity audit, greedy 1:k caliper
matching (PSM), stabilized truncated IPTW with Effective Sample Size (ESS), Love-plot
SMD covariate balance audit, and causal treatment effect estimation.
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
from sklearn.linear_model import LogisticRegression

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


def weighted_smd(x: np.ndarray, t: np.ndarray, w: np.ndarray | None = None) -> float:
    if w is None:
        w = np.ones_like(t, dtype=float)
    m1 = t == 1
    m0 = t == 0
    w1 = w[m1]
    w0 = w[m0]
    x1 = x[m1]
    x0 = x[m0]
    s1 = float(np.sum(w1))
    s0 = float(np.sum(w0))
    if s1 <= 1e-12 or s0 <= 1e-12:
        return 0.0
    mu1 = float(np.sum(w1 * x1) / s1)
    mu0 = float(np.sum(w0 * x0) / s0)
    # Weighted variance (population/reliability formula)
    v1 = float(np.sum(w1 * ((x1 - mu1) ** 2)) / s1)
    v0 = float(np.sum(w0 * ((x0 - mu0) ** 2)) / s0)
    denom = math.sqrt(0.5 * (v1 + v0))
    if denom < 1e-12:
        return 0.0
    return float((mu1 - mu0) / denom)


def greedy_psm_1_to_k(
    logit_ps: np.ndarray, treat: np.ndarray, caliper_sd: float, ratio: int, seed: int = 42
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Greedy 1:k nearest-neighbor matching without replacement on logit(PS)."""
    n = len(treat)
    sd_lps = float(np.std(logit_ps, ddof=1))
    caliper = caliper_sd * sd_lps if sd_lps > 1e-12 else float("inf")

    rng = np.random.default_rng(seed)
    treated_idx = np.flatnonzero(treat == 1).copy()
    rng.shuffle(treated_idx)
    available_ctrl = set(np.flatnonzero(treat == 0).tolist())

    match_stratum = np.full(n, -1, dtype=int)
    match_weight = np.zeros(n, dtype=float)
    matched_mask = np.zeros(n, dtype=bool)
    stratum_id = 0

    for ti in treated_idx:
        if not available_ctrl:
            break
        ctrl_arr = np.fromiter(available_ctrl, dtype=int)
        dists = np.abs(logit_ps[ctrl_arr] - logit_ps[ti])
        within = np.flatnonzero(dists <= caliper)
        if len(within) == 0:
            continue
        order = within[np.argsort(dists[within])[:ratio]]
        chosen = ctrl_arr[order]
        stratum_id += 1
        match_stratum[ti] = stratum_id
        match_weight[ti] = 1.0
        matched_mask[ti] = True
        w_ctrl = 1.0 / float(len(chosen))
        for ci in chosen:
            available_ctrl.remove(int(ci))
            match_stratum[ci] = stratum_id
            match_weight[ci] = w_ctrl
            matched_mask[ci] = True

    return matched_mask, match_stratum, match_weight


def estimate_continuous_effect(y: np.ndarray, t: np.ndarray, w: np.ndarray) -> dict[str, Any]:
    """Weighted difference in means (WLS of y ~ 1 + t) with Huber-White sandwich SE."""
    n = len(y)
    X = np.column_stack([np.ones(n, dtype=float), t.astype(float)])
    XtWX = (X * w[:, None]).T @ X
    XtWy = (X * w[:, None]).T @ y
    beta = np.linalg.solve(XtWX + 1e-10 * np.eye(2), XtWy)
    resid = y - X @ beta
    score = X * (w * resid)[:, None]
    inv_h = np.linalg.pinv(XtWX + 1e-10 * np.eye(2))
    cov = inv_h @ (score.T @ score) @ inv_h
    est = float(beta[1])
    se = float(math.sqrt(max(cov[1, 1], 1e-18)))
    z = est / se
    p = float(2.0 * stats.norm.sf(abs(z)))
    return {
        "estimand": "mean_difference",
        "estimate": round(est, 5),
        "se": round(se, 5),
        "ci_low": round(est - 1.959964 * se, 5),
        "ci_high": round(est + 1.959964 * se, 5),
        "pvalue": p,
    }


def estimate_binary_effect(y: np.ndarray, t: np.ndarray, w: np.ndarray) -> dict[str, Any]:
    """Weighted logistic regression OR + weighted RD & RR with sandwich SE."""
    fit = fit_logistic_weighted(t[:, None].astype(float), y, weights=w, l2_penalty=0.0)
    beta = float(fit["beta"][1])
    se = float(fit["se_robust"][1])
    z = beta / se
    p = float(2.0 * stats.norm.sf(abs(z)))

    w1 = w[t == 1]
    w0 = w[t == 0]
    p1 = float(np.sum(w1 * y[t == 1]) / max(1e-12, float(np.sum(w1))))
    p0 = float(np.sum(w0 * y[t == 0]) / max(1e-12, float(np.sum(w0))))
    rd = p1 - p0
    rr = (p1 / p0) if p0 > 1e-9 else None
    return {
        "estimand": "odds_ratio",
        "or": round(float(math.exp(np.clip(beta, -20.0, 20.0))), 4),
        "ci_low": round(float(math.exp(np.clip(beta - 1.959964 * se, -20.0, 20.0))), 4),
        "ci_high": round(float(math.exp(np.clip(beta + 1.959964 * se, -20.0, 20.0))), 4),
        "risk_treated": round(p1, 4),
        "risk_control": round(p0, 4),
        "risk_difference": round(rd, 4),
        "risk_ratio": round(rr, 4) if rr is not None else None,
        "pvalue": p,
    }


def estimate_survival_effect(times: np.ndarray, events: np.ndarray, t: np.ndarray, w: np.ndarray) -> dict[str, Any]:
    fit = fit_cox_weighted(t[:, None].astype(float), times, events, weights=w)
    beta = float(fit["beta"][0])
    se = float(fit["se_robust"][0])
    z = beta / se
    p = float(2.0 * stats.norm.sf(abs(z)))
    return {
        "estimand": "hazard_ratio",
        "hr": round(float(math.exp(np.clip(beta, -20.0, 20.0))), 4),
        "ci_low": round(float(math.exp(np.clip(beta - 1.959964 * se, -20.0, 20.0))), 4),
        "ci_high": round(float(math.exp(np.clip(beta + 1.959964 * se, -20.0, 20.0))), 4),
        "pvalue": p,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--treatment_col", required=True)
    parser.add_argument("--covariates", action="extend", nargs="+", default=[], required=True)
    parser.add_argument("--treated_label", default=None)
    parser.add_argument("--outcome_col", default=None)
    parser.add_argument("--outcome_type", choices=["auto", "binary", "continuous", "survival"], default="auto")
    parser.add_argument("--time_col", default=None)
    parser.add_argument("--method", choices=["both", "psm", "iptw"], default="both")
    parser.add_argument("--caliper_sd", type=float, default=0.2)
    parser.add_argument("--ratio", type=int, default=1)
    parser.add_argument("--truncate_pct", type=float, default=0.01)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.clinical)
    if args.treatment_col not in df.columns:
        fail(f"treatment_col '{args.treatment_col}' not found in columns: {list(df.columns)}")

    covs = [c.strip() for item in args.covariates for c in str(item).split(",") if c.strip()]
    if not covs:
        fail("at least one covariate is required")
    for c in covs:
        if c not in df.columns:
            fail(f"covariate '{c}' not found in columns: {list(df.columns)}")

    t_ser, treat_map = parse_binary_series(df[args.treatment_col], args.treatment_col, args.treated_label)
    enc_df, _, _ = encode_covariates(df, covs)
    valid = t_ser.notna() & enc_df.notna().all(axis=1)

    # Validate optional outcome columns up front
    o_type = args.outcome_type
    if args.outcome_col:
        if args.outcome_col not in df.columns:
            fail(f"outcome_col '{args.outcome_col}' not found in columns")
        if o_type == "auto":
            if args.time_col:
                o_type = "survival"
            else:
                s_o = df.loc[valid, args.outcome_col].dropna()
                num_o = pd.to_numeric(s_o, errors="coerce")
                if s_o.astype(str).str.strip().nunique() == 2:
                    o_type = "binary"
                elif num_o.notna().sum() >= 0.85 * len(s_o):
                    o_type = "continuous"
                else:
                    fail(f"cannot infer outcome_type for '{args.outcome_col}'; specify --outcome_type")
        if o_type == "survival":
            if not args.time_col or args.time_col not in df.columns:
                fail("time_col is required and must exist in clinical table when outcome_type is survival")
            tm_s = pd.to_numeric(df[args.time_col], errors="coerce")
            ev_s, _ = parse_binary_series(df[args.outcome_col], args.outcome_col)
            valid = valid & tm_s.notna() & (tm_s > 0) & ev_s.notna()
        elif o_type == "binary":
            y_bin_s, _ = parse_binary_series(df[args.outcome_col], args.outcome_col)
            valid = valid & y_bin_s.notna()
        else:
            y_cont_s = pd.to_numeric(df[args.outcome_col], errors="coerce")
            valid = valid & y_cont_s.notna()

    warnings: list[str] = []
    n_drop = int((~valid).sum())
    if n_drop > 0:
        warnings.append(f"excluded {n_drop} record(s) with missing treatment, covariates, or outcome")

    sub = df.loc[valid].copy()
    t_arr = t_ser.loc[valid].to_numpy(dtype=int, copy=True)
    X_cov = enc_df.loc[valid].to_numpy(dtype=float, copy=True)
    cov_names = list(enc_df.columns)

    n = len(t_arr)
    n_treated = int(np.sum(t_arr == 1))
    n_control = int(np.sum(t_arr == 0))
    if n_treated < 5 or n_control < 5:
        fail(f"need >= 5 treated and >= 5 control subjects (got treated={n_treated}, control={n_control})")

    # Fit logistic propensity score model
    mu = X_cov.mean(axis=0)
    sd = np.where(X_cov.std(axis=0, ddof=0) > 1e-12, X_cov.std(axis=0, ddof=0), 1.0)
    X_std = (X_cov - mu) / sd
    ps_clf = LogisticRegression(C=100.0, max_iter=1000)
    ps_clf.fit(X_std, t_arr)
    ps = np.clip(ps_clf.predict_proba(X_std)[:, 1], 1e-5, 1.0 - 1e-5)
    logit_ps = np.log(ps / (1.0 - ps))

    # Positivity / common support audit
    ps1 = ps[t_arr == 1]
    ps0 = ps[t_arr == 0]
    ov_lo = float(max(np.min(ps1), np.min(ps0)))
    ov_hi = float(min(np.max(ps1), np.max(ps0)))
    n_ext_low = int(np.sum(ps < 0.05))
    n_ext_high = int(np.sum(ps > 0.95))
    n_outside = int(np.sum((ps < ov_lo) | (ps > ov_hi)))
    positivity_audit = {
        "ps_treated_min": round(float(np.min(ps1)), 4),
        "ps_treated_max": round(float(np.max(ps1)), 4),
        "ps_control_min": round(float(np.min(ps0)), 4),
        "ps_control_max": round(float(np.max(ps0)), 4),
        "common_support_range": [round(ov_lo, 4), round(ov_hi, 4)],
        "n_outside_common_support": n_outside,
        "n_extreme_ps_lt05": n_ext_low,
        "n_extreme_ps_gt95": n_ext_high,
    }
    if (n_ext_low + n_ext_high) / float(n) > 0.10:
        warnings.append(
            f"{n_ext_low + n_ext_high}/{n} subjects have extreme propensity scores (<0.05 or >0.95); "
            "positivity is strained — check common support and stabilized weight truncation"
        )

    # PSM
    psm_summary = None
    matched_mask = None
    match_weight = None
    matched_path = None
    if args.method in ("both", "psm"):
        matched_mask, match_stratum, match_weight = greedy_psm_1_to_k(
            logit_ps, t_arr, args.caliper_sd, args.ratio
        )
        n_m_tr = int(np.sum(matched_mask & (t_arr == 1)))
        n_m_ct = int(np.sum(matched_mask & (t_arr == 0)))
        if n_m_tr < 3 or n_m_ct < 3:
            warnings.append(f"PSM retained too few matched pairs (treated={n_m_tr}, control={n_m_ct}); try wider caliper_sd")
        matched_df = sub.loc[matched_mask].copy()
        matched_df["propensity_score"] = np.round(ps[matched_mask], 6)
        matched_df["logit_ps"] = np.round(logit_ps[matched_mask], 6)
        matched_df["match_stratum"] = match_stratum[matched_mask]
        matched_df["match_weight"] = np.round(match_weight[matched_mask], 5)
        matched_path = os.path.join(args.output_dir, "psm_matched_cohort.tsv")
        matched_df.to_csv(matched_path, sep="\t", index=False)
        psm_summary = {
            "caliper_sd": args.caliper_sd,
            "ratio": args.ratio,
            "n_matched_treated": n_m_tr,
            "n_matched_control": n_m_ct,
            "treated_match_rate": round(n_m_tr / float(n_treated), 4),
        }

    # Stabilized IPTW
    iptw_summary = None
    w_iptw = None
    if args.method in ("both", "iptw"):
        p_t = n_treated / float(n)
        w_raw = np.where(t_arr == 1, p_t / ps, (1.0 - p_t) / (1.0 - ps))
        if 0.0 < args.truncate_pct < 0.25:
            q_lo, q_hi = np.quantile(w_raw, [args.truncate_pct, 1.0 - args.truncate_pct])
            w_iptw = np.clip(w_raw, float(q_lo), float(q_hi))
        else:
            w_iptw = w_raw
        w1 = w_iptw[t_arr == 1]
        w0 = w_iptw[t_arr == 0]
        ess_tr = float((np.sum(w1) ** 2) / max(1e-12, float(np.sum(w1 ** 2))))
        ess_ct = float((np.sum(w0) ** 2) / max(1e-12, float(np.sum(w0 ** 2))))
        iptw_summary = {
            "stabilized": True,
            "truncate_pct": args.truncate_pct,
            "weight_min": round(float(np.min(w_iptw)), 4),
            "weight_max": round(float(np.max(w_iptw)), 4),
            "ess_treated": round(ess_tr, 2),
            "ess_control": round(ess_ct, 2),
        }

    # Covariate balance audit
    bal_rows: list[dict[str, Any]] = []
    imbal_unadj = []
    imbal_psm = []
    imbal_iptw = []
    for j, cname in enumerate(cov_names):
        xj = X_cov[:, j]
        smd_u = weighted_smd(xj, t_arr, None)
        if abs(smd_u) >= 0.10:
            imbal_unadj.append(cname)
        row: dict[str, Any] = {
            "covariate": cname,
            "smd_unadjusted": round(smd_u, 4),
            "abs_smd_unadjusted": round(abs(smd_u), 4),
        }
        if matched_mask is not None and match_weight is not None and int(matched_mask.sum()) >= 4:
            smd_p = weighted_smd(xj[matched_mask], t_arr[matched_mask], match_weight[matched_mask])
            row["smd_psm"] = round(smd_p, 4)
            row["abs_smd_psm"] = round(abs(smd_p), 4)
            if abs(smd_p) >= 0.10:
                imbal_psm.append(cname)
        else:
            row["smd_psm"] = None
            row["abs_smd_psm"] = None

        if w_iptw is not None:
            smd_w = weighted_smd(xj, t_arr, w_iptw)
            row["smd_iptw"] = round(smd_w, 4)
            row["abs_smd_iptw"] = round(abs(smd_w), 4)
            if abs(smd_w) >= 0.10:
                imbal_iptw.append(cname)
        else:
            row["smd_iptw"] = None
            row["abs_smd_iptw"] = None
        bal_rows.append(row)

    if imbal_psm:
        warnings.append(f"PSM residual imbalance |SMD| >= 0.10 for: {imbal_psm}")
    if imbal_iptw:
        warnings.append(f"IPTW residual imbalance |SMD| >= 0.10 for: {imbal_iptw}")

    bal_df = pd.DataFrame(bal_rows)
    bal_path = os.path.join(args.output_dir, "covariate_balance_smd.tsv")
    bal_df.to_csv(bal_path, sep="\t", index=False)

    # Causal treatment effects
    treatment_effects = None
    if args.outcome_col:
        treatment_effects = {"outcome_col": args.outcome_col, "outcome_type": o_type}
        if o_type == "continuous":
            y_c = pd.to_numeric(sub[args.outcome_col], errors="coerce").to_numpy(dtype=float, copy=True)
            treatment_effects["unadjusted"] = estimate_continuous_effect(y_c, t_arr, np.ones(n))
            if matched_mask is not None and match_weight is not None and int(matched_mask.sum()) >= 6:
                treatment_effects["psm"] = estimate_continuous_effect(y_c[matched_mask], t_arr[matched_mask], match_weight[matched_mask])
            if w_iptw is not None:
                treatment_effects["iptw"] = estimate_continuous_effect(y_c, t_arr, w_iptw)
        elif o_type == "binary":
            y_b_s, _ = parse_binary_series(sub[args.outcome_col], args.outcome_col)
            y_b = y_b_s.to_numpy(dtype=float, copy=True)
            treatment_effects["unadjusted"] = estimate_binary_effect(y_b, t_arr, np.ones(n))
            if matched_mask is not None and match_weight is not None and int(matched_mask.sum()) >= 6:
                treatment_effects["psm"] = estimate_binary_effect(y_b[matched_mask], t_arr[matched_mask], match_weight[matched_mask])
            if w_iptw is not None:
                treatment_effects["iptw"] = estimate_binary_effect(y_b, t_arr, w_iptw)
        else:
            tm_arr = pd.to_numeric(sub[args.time_col], errors="coerce").to_numpy(dtype=float, copy=True)
            ev_s, _ = parse_binary_series(sub[args.outcome_col], args.outcome_col)
            ev_arr = ev_s.to_numpy(dtype=int, copy=True)
            treatment_effects["unadjusted"] = estimate_survival_effect(tm_arr, ev_arr, t_arr, np.ones(n))
            if matched_mask is not None and match_weight is not None and int(matched_mask.sum()) >= 6 and int(ev_arr[matched_mask].sum()) >= 2:
                treatment_effects["psm"] = estimate_survival_effect(tm_arr[matched_mask], ev_arr[matched_mask], t_arr[matched_mask], match_weight[matched_mask])
            if w_iptw is not None:
                treatment_effects["iptw"] = estimate_survival_effect(tm_arr, ev_arr, t_arr, w_iptw)

    # 2-panel Diagnostic Figure: Love plot + Propensity Score Mirror Histogram
    fig, (ax_love, ax_ps) = plt.subplots(1, 2, figsize=(12.0, max(4.4, 0.35 * len(cov_names) + 2.2)))

    rev_rows = bal_rows[::-1]
    y_pos = np.arange(len(rev_rows))
    ax_love.axvline(0.10, color="#DC2626", linestyle="--", linewidth=1.1, label="Balance threshold (|SMD| = 0.10)")
    ax_love.scatter([r["abs_smd_unadjusted"] for r in rev_rows], y_pos, color="#64748B", marker="o", s=46, label="Unadjusted")
    if any(r["abs_smd_psm"] is not None for r in rev_rows):
        ax_love.scatter(
            [r["abs_smd_psm"] if r["abs_smd_psm"] is not None else np.nan for r in rev_rows],
            y_pos,
            color="#2563EB",
            marker="^",
            s=52,
            label="PSM Matched",
        )
    if any(r["abs_smd_iptw"] is not None for r in rev_rows):
        ax_love.scatter(
            [r["abs_smd_iptw"] if r["abs_smd_iptw"] is not None else np.nan for r in rev_rows],
            y_pos,
            color="#10B981",
            marker="s",
            s=46,
            label="Stabilized IPTW",
        )
    ax_love.set_yticks(y_pos)
    ax_love.set_yticklabels([r["covariate"] for r in rev_rows], fontsize=8.2)
    ax_love.set_xlabel("Absolute Standardized Mean Difference (|SMD|)")
    ax_love.set_title("Covariate Balance Love Plot", fontsize=10)
    ax_love.legend(frameon=False, loc="lower right", fontsize=8)
    ax_love.spines["top"].set_visible(False)
    ax_love.spines["right"].set_visible(False)

    # Mirror histogram of propensity scores
    bins_ps = np.linspace(0.0, 1.0, 21)
    h1, _ = np.histogram(ps1, bins=bins_ps, density=True)
    h0, _ = np.histogram(ps0, bins=bins_ps, density=True)
    centers = 0.5 * (bins_ps[:-1] + bins_ps[1:])
    ax_ps.bar(centers, h1, width=0.044, color="#E11D48", alpha=0.72, label=f"Treated (n={n_treated})")
    ax_ps.bar(centers, -h0, width=0.044, color="#2563EB", alpha=0.72, label=f"Control (n={n_control})")
    ax_ps.axhline(0.0, color="#0F172A", linewidth=0.9)
    ax_ps.set_xlim(0.0, 1.0)
    ax_ps.set_xlabel("Estimated Propensity Score e(X)")
    ax_ps.set_ylabel("Density (+ Treated / - Control)")
    ax_ps.set_title("Propensity Score Overlap & Positivity", fontsize=10)
    ax_ps.legend(frameon=False, loc="upper right", fontsize=8)
    ax_ps.spines["top"].set_visible(False)
    ax_ps.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "psm_iptw_diagnostics.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {
        "treatment_col": args.treatment_col,
        "covariates": covs,
        "method": args.method,
        "caliper_sd": args.caliper_sd,
        "ratio": args.ratio,
    }
    artifacts = [
        write_artifact(fig_path, "clin_psm_iptw", "PSM/IPTW Love plot & PS overlap", inputs=(args.clinical,), parameters=params_meta),
        write_artifact(bal_path, "clin_psm_iptw", "Covariate balance SMD table", inputs=(args.clinical,), parameters=params_meta, shape=bal_df.shape),
    ]
    if matched_path:
        artifacts.append(
            write_artifact(matched_path, "clin_psm_iptw", "PSM matched cohort table", inputs=(args.clinical,), parameters=params_meta, shape=matched_df.shape)
        )

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "treatment_col": args.treatment_col,
            "treatment_mapping": treat_map,
            "method": args.method,
            "n_samples": n,
            "n_treated": n_treated,
            "n_control": n_control,
            "positivity_audit": positivity_audit,
            "psm_summary": psm_summary,
            "iptw_summary": iptw_summary,
            "covariate_balance": bal_rows,
            "treatment_effects": treatment_effects,
            "balance_table": os.path.abspath(bal_path),
            "matched_cohort_table": os.path.abspath(matched_path) if matched_path else None,
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
