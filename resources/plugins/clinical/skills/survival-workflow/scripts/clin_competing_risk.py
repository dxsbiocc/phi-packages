#!/usr/bin/env python3
"""clin_competing_risk: Aalen-Johansen Cumulative Incidence Functions (CIF),
naive 1-KM overestimation audit, Gray's subdistribution test, and Fine-Gray IPCW
subdistribution proportional-hazards regression (sHR).
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
from _common import (
    bh_fdr,
    emit,
    fail,
    infer_time_unit,
    kaplan_meier,
    read_table,
    write_artifact,
)

PALETTE = ["#2563EB", "#E11D48", "#10B981", "#F59E0B", "#8B5CF6", "#06B6D4"]

CENSOR_TOKENS = {
    "0", "0.0", "false", "f", "no", "n", "alive", "living", "censored",
    "censor", "diseasefree", "disease-free", "0:living", "0:diseasefree", "0:censored",
}


def _norm_token(val: Any) -> str:
    s = str(val).strip()
    try:
        f = float(s)
        if math.isfinite(f) and f == int(f):
            return str(int(f))
    except ValueError:
        pass
    return s.lower()


def parse_competing_status(
    series: pd.Series,
    primary_code: str,
    competing_codes: list[str],
    col_name: str,
) -> tuple[pd.Series, dict[str, int], dict[str, str]]:
    """Map status column to {0: censored, 1: primary event, 2: competing event, NaN: missing}."""
    p_norm = _norm_token(primary_code)
    c_norms = {_norm_token(c) for c in competing_codes}
    explicit_comp = len(c_norms) > 0

    mapped: list[float] = []
    counts: dict[str, int] = {}
    role_map: dict[str, str] = {}
    unknown: set[str] = set()

    for val in series:
        if pd.isna(val):
            mapped.append(np.nan)
            continue
        raw = str(val).strip()
        if not raw or raw.lower() in ("na", "nan", "null", "none", "."):
            mapped.append(np.nan)
            continue
        counts[raw] = counts.get(raw, 0) + 1
        tok = _norm_token(raw)
        if tok == p_norm:
            mapped.append(1.0)
            role_map[raw] = "primary"
        elif tok in c_norms:
            mapped.append(2.0)
            role_map[raw] = "competing"
        elif tok in CENSOR_TOKENS:
            mapped.append(0.0)
            role_map[raw] = "censored"
        elif not explicit_comp:
            mapped.append(2.0)
            role_map[raw] = "competing"
        else:
            unknown.add(raw)
            mapped.append(np.nan)

    if unknown:
        fail(
            f"unrecognized status values in '{col_name}' with explicit competing_codes: {sorted(unknown)[:8]}"
        )
    return pd.Series(mapped, index=series.index, dtype=float), counts, role_map


def _loglog_cif_ci(cif_vals: np.ndarray, var_vals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lo = np.zeros_like(cif_vals, dtype=float)
    hi = np.zeros_like(cif_vals, dtype=float)
    for idx, (f_val, v_val) in enumerate(zip(cif_vals, var_vals)):
        if f_val <= 0.0:
            lo[idx], hi[idx] = 0.0, 0.0
        elif f_val >= 1.0:
            lo[idx], hi[idx] = 1.0, 1.0
        elif v_val <= 0.0 or not math.isfinite(v_val):
            lo[idx], hi[idx] = float(f_val), float(f_val)
        else:
            se_f = math.sqrt(v_val)
            log_f = math.log(max(f_val, 1e-15))
            se_theta = se_f / (f_val * abs(log_f))
            fac = math.exp(1.959964 * se_theta)
            lo[idx] = float(np.clip(f_val ** fac, 0.0, 1.0))
            hi[idx] = float(np.clip(f_val ** (1.0 / fac), 0.0, 1.0))
    return lo, hi


def aalen_johansen_cif(times: np.ndarray, status: np.ndarray) -> dict[str, Any]:
    """Aalen-Johansen CIF for primary (1) and competing (2) events with delta-method 95% CIs."""
    t = np.asarray(times, dtype=float)
    st = np.asarray(status, dtype=int)
    order = np.argsort(t)
    t = t[order]
    st = st[order]
    n_total = len(t)
    uniq_t = np.unique(t)

    curve_t = [0.0]
    curve_s = [1.0]
    curve_f1 = [0.0]
    curve_f2 = [0.0]
    curve_nr = [n_total]
    curve_d1 = [0]
    curve_d2 = [0]
    curve_c = [0]

    # Store event-time increments for Marubini-Valsecchi / Aalen delta-method variance
    ev_n: list[int] = []
    ev_d: list[int] = []
    ev_d1: list[int] = []
    ev_d2: list[int] = []
    ev_s_prev: list[float] = []
    ev_f1_at: list[float] = []
    ev_f2_at: list[float] = []

    surv = 1.0
    f1 = 0.0
    f2 = 0.0
    var1_list = [0.0]
    var2_list = [0.0]

    for ut in uniq_t:
        at_risk = int(np.sum(t >= ut))
        d1 = int(np.sum((t == ut) & (st == 1)))
        d2 = int(np.sum((t == ut) & (st == 2)))
        c_i = int(np.sum((t == ut) & (st == 0)))
        d_all = d1 + d2

        if d_all > 0 and at_risk > 0:
            s_prev = surv
            f1 += s_prev * (d1 / float(at_risk))
            f2 += s_prev * (d2 / float(at_risk))
            surv *= 1.0 - (d_all / float(at_risk))

            ev_n.append(at_risk)
            ev_d.append(d_all)
            ev_d1.append(d1)
            ev_d2.append(d2)
            ev_s_prev.append(s_prev)
            ev_f1_at.append(f1)
            ev_f2_at.append(f2)

        # Compute variance at current time point over all event times <= ut
        v1 = 0.0
        v2 = 0.0
        for nj, dj, d1j, d2j, spj, f1j, f2j in zip(
            ev_n, ev_d, ev_d1, ev_d2, ev_s_prev, ev_f1_at, ev_f2_at
        ):
            gw = (dj / float(nj * (nj - dj))) if nj > dj else 0.0
            diff1 = f1 - f1j
            v1 += (
                (diff1 ** 2) * gw
                + (spj ** 2) * (d1j * (nj - d1j)) / float(nj ** 3)
                - 2.0 * diff1 * spj * (d1j / float(nj ** 2))
            )
            diff2 = f2 - f2j
            v2 += (
                (diff2 ** 2) * gw
                + (spj ** 2) * (d2j * (nj - d2j)) / float(nj ** 3)
                - 2.0 * diff2 * spj * (d2j / float(nj ** 2))
            )
        var1_list.append(max(0.0, v1))
        var2_list.append(max(0.0, v2))

        curve_t.append(float(ut))
        curve_s.append(float(surv))
        curve_f1.append(float(np.clip(f1, 0.0, 1.0)))
        curve_f2.append(float(np.clip(f2, 0.0, 1.0)))
        curve_nr.append(at_risk)
        curve_d1.append(d1)
        curve_d2.append(d2)
        curve_c.append(c_i)

    arr_f1 = np.array(curve_f1, dtype=float)
    arr_f2 = np.array(curve_f2, dtype=float)
    lo1, hi1 = _loglog_cif_ci(arr_f1, np.array(var1_list, dtype=float))
    lo2, hi2 = _loglog_cif_ci(arr_f2, np.array(var2_list, dtype=float))

    return {
        "n": n_total,
        "n_primary": int(np.sum(st == 1)),
        "n_competing": int(np.sum(st == 2)),
        "n_censored": int(np.sum(st == 0)),
        "time": np.array(curve_t, dtype=float),
        "surv": np.array(curve_s, dtype=float),
        "cif1": arr_f1,
        "cif1_ci_low": lo1,
        "cif1_ci_high": hi1,
        "cif2": arr_f2,
        "cif2_ci_low": lo2,
        "cif2_ci_high": hi2,
        "n_risk": np.array(curve_nr, dtype=int),
        "n_ev1": np.array(curve_d1, dtype=int),
        "n_ev2": np.array(curve_d2, dtype=int),
        "n_censor": np.array(curve_c, dtype=int),
    }


def cif_at(aj: dict[str, Any], horizon: float) -> tuple[float, float, float, float, float, float, int]:
    t = aj["time"]
    idx = int(np.searchsorted(t, horizon, side="right") - 1)
    idx = max(0, min(idx, len(t) - 1))
    r_idx = int(np.searchsorted(t[1:], horizon, side="left")) + 1
    n_risk = int(aj["n_risk"][r_idx]) if r_idx < len(aj["n_risk"]) else 0
    return (
        float(aj["cif1"][idx]),
        float(aj["cif1_ci_low"][idx]),
        float(aj["cif1_ci_high"][idx]),
        float(aj["cif2"][idx]),
        float(aj["cif2_ci_low"][idx]),
        float(aj["cif2_ci_high"][idx]),
        n_risk,
    )


def _censoring_km_lookup(times: np.ndarray, status: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    km_c = kaplan_meier(times, (status == 0).astype(int))
    return km_c["time"], np.maximum(km_c["surv"], 0.02)


def _g_left(g_times: np.ndarray, g_surv: np.ndarray, u: float) -> float:
    idx = int(np.searchsorted(g_times, u, side="left") - 1)
    idx = max(0, min(idx, len(g_surv) - 1))
    return float(g_surv[idx])


def _build_ipcw_risk_sets(
    times: np.ndarray, status: np.ndarray
) -> tuple[np.ndarray, list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    g_times, g_surv = _censoring_km_lookup(times, status)
    g_at_i = np.array([_g_left(g_times, g_surv, float(ti)) for ti in times], dtype=float)
    uniq_ev1 = np.unique(times[status == 1])

    risk_sets: list[np.ndarray] = []
    event_sets: list[np.ndarray] = []
    weights_list: list[np.ndarray] = []
    for ut in uniq_ev1:
        g_ut = _g_left(g_times, g_surv, float(ut))
        mask_at_risk = times >= ut
        mask_comp_prior = (times < ut) & (status == 2)
        r_idx = np.flatnonzero(mask_at_risk | mask_comp_prior)
        w = np.where(times[r_idx] >= ut, 1.0, np.clip(g_ut / g_at_i[r_idx], 0.0, 1.0))
        d_idx = np.flatnonzero((times == ut) & (status == 1))
        risk_sets.append(r_idx)
        event_sets.append(d_idx)
        weights_list.append(w)
    return uniq_ev1, risk_sets, event_sets, weights_list


def grays_test_ipcw(times: np.ndarray, status: np.ndarray, groups: np.ndarray) -> dict[str, Any]:
    """Subdistribution score test (Fine-Gray score test at beta=0) across k groups."""
    t = np.asarray(times, dtype=float)
    st = np.asarray(status, dtype=int)
    g = np.asarray(groups, dtype=str)
    levels = sorted(set(g.tolist()))
    k = len(levels)
    if k < 2 or int(np.sum(st == 1)) < 2:
        return {"chisq": None, "df": None, "pvalue": None, "levels": levels}

    g_idx = np.array([levels.index(val) for val in g], dtype=int)
    Z = np.zeros((len(t), k - 1), dtype=float)
    for i, gi in enumerate(g_idx):
        if gi > 0:
            Z[i, gi - 1] = 1.0

    _, risk_sets, event_sets, weights_list = _build_ipcw_risk_sets(t, st)
    u_subj = np.zeros((len(t), k - 1), dtype=float)

    for r_idx, d_idx, w in zip(risk_sets, event_sets, weights_list):
        d_m = len(d_idx)
        w_sum = float(np.sum(w))
        if w_sum <= 1e-12:
            continue
        z_bar = (w[:, None] * Z[r_idx]).sum(axis=0) / w_sum
        u_subj[d_idx] += Z[d_idx] - z_bar[None, :]
        u_subj[r_idx] -= (d_m * (w / w_sum))[:, None] * (Z[r_idx] - z_bar[None, :])

    u_tot = u_subj.sum(axis=0)
    v_rob = u_subj.T @ u_subj
    chisq = float(max(0.0, u_tot.T @ np.linalg.pinv(v_rob) @ u_tot))
    df = k - 1
    pval = float(stats.chi2.sf(chisq, df=df))
    return {"chisq": round(chisq, 4), "df": df, "pvalue": pval, "levels": levels}


def fit_fine_gray(
    X_raw: np.ndarray, times: np.ndarray, status: np.ndarray, max_iter: int = 60, tol: float = 1e-7
) -> dict[str, Any]:
    """Fine-Gray IPCW subdistribution PH regression with robust sandwich variance."""
    X = np.asarray(X_raw, dtype=float)
    t = np.asarray(times, dtype=float)
    st = np.asarray(status, dtype=int)
    n, p = X.shape
    n_ev1 = int(np.sum(st == 1))
    n_ev2 = int(np.sum(st == 2))
    if n_ev1 < 2:
        fail("Fine-Gray regression requires >= 2 primary events")

    means = X.mean(axis=0)
    sds = X.std(axis=0, ddof=0)
    sds_safe = np.where(sds > 1e-12, sds, 1.0)
    Z = (X - means) / sds_safe

    _, risk_sets, event_sets, weights_list = _build_ipcw_risk_sets(t, st)

    def nll_grad_hess(beta_z: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        eta = np.clip(Z @ beta_z, -50.0, 50.0)
        nll = 0.0
        grad = np.zeros(p, dtype=float)
        hess = np.zeros((p, p), dtype=float)
        for r_idx, d_idx, w in zip(risk_sets, event_sets, weights_list):
            d_m = len(d_idx)
            eta_r = eta[r_idx]
            shift = float(np.max(eta_r))
            ew = w * np.exp(eta_r - shift)
            s0 = float(np.sum(ew))
            if s0 <= 1e-15:
                continue
            w_norm = ew / s0
            Z_r = Z[r_idx]
            x_bar = w_norm @ Z_r
            nll -= float(np.sum(eta[d_idx])) - d_m * (shift + math.log(s0))
            grad -= np.sum(Z[d_idx], axis=0) - d_m * x_bar
            Z_c = Z_r - x_bar[None, :]
            hess += d_m * ((Z_c * w_norm[:, None]).T @ Z_c)
        return nll, grad, hess

    beta_z = np.zeros(p, dtype=float)
    converged = False
    nll_cur, grad_cur, hess_cur = nll_grad_hess(beta_z)

    for _ in range(max_iter):
        reg_hess = hess_cur + 1e-9 * np.eye(p)
        try:
            step = np.linalg.solve(reg_hess, grad_cur)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(reg_hess) @ grad_cur
        alpha = 1.0
        accepted = False
        for _ls in range(20):
            cand = beta_z - alpha * step
            if np.max(np.abs(cand)) <= 25.0:
                nll_new, grad_new, hess_new = nll_grad_hess(cand)
                if math.isfinite(nll_new) and nll_new <= nll_cur + 1e-10:
                    beta_z = cand
                    nll_cur, grad_cur, hess_cur = nll_new, grad_new, hess_new
                    accepted = True
                    break
            alpha *= 0.5
        if not accepted:
            break
        if float(np.max(np.abs(alpha * step))) < tol:
            converged = True
            break

    # Robust Lin-Wei / Fine-Gray score-residual sandwich covariance
    eta_hat = np.clip(Z @ beta_z, -50.0, 50.0)
    u_subj = np.zeros((n, p), dtype=float)
    for r_idx, d_idx, w in zip(risk_sets, event_sets, weights_list):
        d_m = len(d_idx)
        ew = w * np.exp(eta_hat[r_idx] - float(np.max(eta_hat[r_idx])))
        s0 = float(np.sum(ew))
        if s0 <= 1e-15:
            continue
        w_norm = ew / s0
        Z_r = Z[r_idx]
        x_bar = w_norm @ Z_r
        u_subj[d_idx] += Z[d_idx] - x_bar[None, :]
        u_subj[r_idx] -= (d_m * w_norm)[:, None] * (Z_r - x_bar[None, :])

    inv_info = np.linalg.pinv(hess_cur + 1e-9 * np.eye(p))
    meat_z = u_subj.T @ u_subj
    cov_z_rob = inv_info @ meat_z @ inv_info

    inv_sd = 1.0 / sds_safe
    beta = beta_z * inv_sd
    cov_rob = cov_z_rob * np.outer(inv_sd, inv_sd)
    se = np.sqrt(np.maximum(np.diag(cov_rob), 1e-18))
    z_stat = beta / se
    pvals = 2.0 * stats.norm.sf(np.abs(z_stat))

    return {
        "beta": beta,
        "se": se,
        "z": z_stat,
        "pvalue": pvals,
        "shr": np.exp(np.clip(beta, -20.0, 20.0)),
        "ci_low": np.exp(np.clip(beta - 1.959964 * se, -20.0, 20.0)),
        "ci_high": np.exp(np.clip(beta + 1.959964 * se, -20.0, 20.0)),
        "converged": converged,
        "n_samples": n,
        "n_primary_events": n_ev1,
        "n_competing_events": n_ev2,
        "n_parameters": p,
    }


def encode_covariates(df: pd.DataFrame, covariates: list[str]) -> tuple[pd.DataFrame, list[str]]:
    parts = []
    notes = []
    for cov in covariates:
        s = df[cov]
        num = pd.to_numeric(s, errors="coerce")
        is_cont = num.notna().sum() >= 0.85 * max(1, s.notna().sum()) and num.dropna().nunique() > 2
        if is_cont:
            parts.append(num.astype(float).rename(cov))
        else:
            cats = s.where(s.notna() & (s.astype(str).str.strip() != ""), np.nan).astype(object)
            levels = sorted(set(str(x).strip() for x in cats.dropna()))
            if len(levels) < 2:
                fail(f"categorical covariate '{cov}' has < 2 non-missing levels")
            ref = levels[0]
            notes.append(f"{cov}: reference='{ref}'")
            for lvl in levels[1:]:
                ind = pd.Series(
                    np.where(cats.isna(), np.nan, (cats.astype(str).str.strip() == lvl).astype(float)),
                    index=df.index,
                    name=f"{cov} [{lvl} vs {ref}]",
                )
                parts.append(ind)
    return pd.concat(parts, axis=1), notes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--time_col", required=True)
    parser.add_argument("--status_col", required=True)
    parser.add_argument("--primary_code", default="1")
    parser.add_argument("--competing_codes", action="extend", nargs="+", default=[])
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default=None)
    parser.add_argument("--covariates", action="extend", nargs="+", default=[])
    parser.add_argument("--landmarks", action="extend", nargs="+", type=float, default=[])
    parser.add_argument("--time_label", default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    clin = read_table(args.clinical)
    for req_col in (args.time_col, args.status_col):
        if req_col not in clin.columns:
            fail(f"column '{req_col}' not found in clinical table: {list(clin.columns)}")

    comp_codes = [c.strip() for item in args.competing_codes for c in str(item).split(",") if c.strip()]
    covs = [c.strip() for item in args.covariates for c in str(item).split(",") if c.strip()]

    t_s = pd.to_numeric(clin[args.time_col], errors="coerce")
    st_s, raw_counts, role_map = parse_competing_status(
        clin[args.status_col], args.primary_code, comp_codes, args.status_col
    )
    valid = t_s.notna() & (t_s > 0) & st_s.notna()
    warnings: list[str] = []
    n_nonpos = int((t_s.notna() & (t_s <= 0)).sum())
    if n_nonpos > 0:
        warnings.append(f"excluded {n_nonpos} record(s) with non-positive follow-up time")

    sub = clin.loc[valid].copy()
    sub["_time"] = t_s.loc[valid].to_numpy(dtype=float, copy=True)
    sub["_status"] = st_s.loc[valid].to_numpy(dtype=int, copy=True)

    if len(sub) < 5:
        fail("fewer than 5 valid records for competing risks analysis")
    n_ev1 = int((sub["_status"] == 1).sum())
    n_ev2 = int((sub["_status"] == 2).sum())
    n_cens = int((sub["_status"] == 0).sum())
    if n_ev1 < 1:
        fail(f"0 primary events found for primary_code='{args.primary_code}'")
    if n_ev2 == 0:
        warnings.append("0 competing events observed; Aalen-Johansen CIF coincides with 1 - KM")

    t_arr = sub["_time"].to_numpy(dtype=float, copy=True)
    st_arr = sub["_status"].to_numpy(dtype=int, copy=True)
    aj_overall = aalen_johansen_cif(t_arr, st_arr)

    # Naive 1 - KM treating competing events as censored
    km_naive = kaplan_meier(t_arr, (st_arr == 1).astype(int))
    naive_1mkm_on_aj = []
    for tm in aj_overall["time"]:
        idx = int(np.searchsorted(km_naive["time"], tm, side="right") - 1)
        idx = max(0, min(idx, len(km_naive["time"]) - 1))
        naive_1mkm_on_aj.append(1.0 - float(km_naive["surv"][idx]))
    naive_arr = np.array(naive_1mkm_on_aj, dtype=float)
    overest = np.maximum(0.0, naive_arr - aj_overall["cif1"])
    max_idx = int(np.argmax(overest))
    bias_summary = {
        "max_overestimation": round(float(overest[max_idx]), 5),
        "max_overestimation_time": round(float(aj_overall["time"][max_idx]), 4),
        "final_naive_1_minus_km": round(float(naive_arr[-1]), 5),
        "final_aalen_johansen_cif1": round(float(aj_overall["cif1"][-1]), 5),
        "final_overestimation": round(float(overest[-1]), 5),
    }

    cif_rows: list[dict[str, Any]] = []
    group_summaries: list[dict[str, Any]] = []
    group_aj: dict[str, dict[str, Any]] = {}
    grays = None

    if args.group_col:
        if args.group_col not in sub.columns:
            fail(f"group_col '{args.group_col}' not found in clinical table")
        sub = sub[sub[args.group_col].notna() & (sub[args.group_col].astype(str).str.strip() != "")].copy()
        sub["_grp"] = sub[args.group_col].astype(str).str.strip()
        levels = sorted(sub["_grp"].unique().tolist())
        if len(levels) < 2:
            fail(f"group_col '{args.group_col}' has < 2 groups")
        grays = grays_test_ipcw(
            sub["_time"].to_numpy(dtype=float, copy=True),
            sub["_status"].to_numpy(dtype=int, copy=True),
            sub["_grp"].to_numpy(dtype=str, copy=True),
        )
        for g in levels:
            sg = sub[sub["_grp"] == g]
            aj_g = aalen_johansen_cif(
                sg["_time"].to_numpy(dtype=float, copy=True),
                sg["_status"].to_numpy(dtype=int, copy=True),
            )
            group_aj[g] = aj_g
            group_summaries.append(
                {
                    "group": g,
                    "n": aj_g["n"],
                    "n_primary": aj_g["n_primary"],
                    "n_competing": aj_g["n_competing"],
                    "n_censored": aj_g["n_censored"],
                    "final_cif_primary": round(float(aj_g["cif1"][-1]), 4),
                    "final_cif_competing": round(float(aj_g["cif2"][-1]), 4),
                }
            )
            for idx_t in range(len(aj_g["time"])):
                cif_rows.append(
                    {
                        "group": g,
                        "time": round(float(aj_g["time"][idx_t]), 5),
                        "cif_primary": round(float(aj_g["cif1"][idx_t]), 5),
                        "cif_primary_ci_low": round(float(aj_g["cif1_ci_low"][idx_t]), 5),
                        "cif_primary_ci_high": round(float(aj_g["cif1_ci_high"][idx_t]), 5),
                        "cif_competing": round(float(aj_g["cif2"][idx_t]), 5),
                        "cif_competing_ci_low": round(float(aj_g["cif2_ci_low"][idx_t]), 5),
                        "cif_competing_ci_high": round(float(aj_g["cif2_ci_high"][idx_t]), 5),
                        "n_risk": int(aj_g["n_risk"][idx_t]),
                    }
                )
    else:
        levels = ["Overall"]
        group_aj["Overall"] = aj_overall
        for idx_t in range(len(aj_overall["time"])):
            cif_rows.append(
                {
                    "group": "Overall",
                    "time": round(float(aj_overall["time"][idx_t]), 5),
                    "cif_primary": round(float(aj_overall["cif1"][idx_t]), 5),
                    "cif_primary_ci_low": round(float(aj_overall["cif1_ci_low"][idx_t]), 5),
                    "cif_primary_ci_high": round(float(aj_overall["cif1_ci_high"][idx_t]), 5),
                    "cif_competing": round(float(aj_overall["cif2"][idx_t]), 5),
                    "cif_competing_ci_low": round(float(aj_overall["cif2_ci_low"][idx_t]), 5),
                    "cif_competing_ci_high": round(float(aj_overall["cif2_ci_high"][idx_t]), 5),
                    "naive_1_minus_km": round(float(naive_arr[idx_t]), 5),
                    "n_risk": int(aj_overall["n_risk"][idx_t]),
                }
            )

    # Landmarks
    landmark_rows: list[dict[str, Any]] = []
    if args.landmarks:
        for hz in sorted(set(args.landmarks)):
            if hz <= 0:
                continue
            for g in levels:
                c1, l1, h1, c2, l2, h2, nr = cif_at(group_aj[g], float(hz))
                landmark_rows.append(
                    {
                        "horizon": float(hz),
                        "group": g,
                        "cif_primary": round(c1, 4),
                        "cif_primary_ci_low": round(l1, 4),
                        "cif_primary_ci_high": round(h1, 4),
                        "cif_competing": round(c2, 4),
                        "cif_competing_ci_low": round(l2, 4),
                        "cif_competing_ci_high": round(h2, 4),
                        "n_at_risk": nr,
                    }
                )

    cif_df = pd.DataFrame(cif_rows)
    cif_path = os.path.join(args.output_dir, "cif_estimates.tsv")
    cif_df.to_csv(cif_path, sep="\t", index=False)

    # Fine-Gray regression when covariates or group_col provided
    fg_covs = covs if covs else ([args.group_col] if args.group_col else [])
    fg_model = None
    fg_path = None
    fg_terms: list[dict[str, Any]] = []
    artifacts: list[str] = []

    if fg_covs and n_ev1 >= 4:
        for c in fg_covs:
            if c not in sub.columns:
                fail(f"covariate '{c}' not found in clinical table: {list(sub.columns)}")
        enc_df, ref_notes = encode_covariates(sub, fg_covs)
        valid_fg = enc_df.notna().all(axis=1)
        n_drop_fg = int((~valid_fg).sum())
        if n_drop_fg > 0:
            warnings.append(f"Fine-Gray regression excluded {n_drop_fg} sample(s) with missing covariates")
        if int(valid_fg.sum()) >= max(6, enc_df.shape[1] + 2) and int((sub.loc[valid_fg, "_status"] == 1).sum()) >= 2:
            fg_fit = fit_fine_gray(
                enc_df.loc[valid_fg].to_numpy(dtype=float, copy=True),
                sub.loc[valid_fg, "_time"].to_numpy(dtype=float, copy=True),
                sub.loc[valid_fg, "_status"].to_numpy(dtype=int, copy=True),
            )
            padj = bh_fdr(fg_fit["pvalue"])
            for idx_t, term_name in enumerate(enc_df.columns):
                fg_terms.append(
                    {
                        "term": str(term_name),
                        "beta": round(float(fg_fit["beta"][idx_t]), 5),
                        "se": round(float(fg_fit["se"][idx_t]), 5),
                        "shr": round(float(fg_fit["shr"][idx_t]), 4),
                        "ci_low": round(float(fg_fit["ci_low"][idx_t]), 4),
                        "ci_high": round(float(fg_fit["ci_high"][idx_t]), 4),
                        "z": round(float(fg_fit["z"][idx_t]), 4),
                        "pvalue": float(fg_fit["pvalue"][idx_t]),
                        "padj": float(padj[idx_t]) if np.isfinite(padj[idx_t]) else None,
                    }
                )
            fg_df = pd.DataFrame(fg_terms)
            fg_path = os.path.join(args.output_dir, "fine_gray_shr.tsv")
            fg_df.to_csv(fg_path, sep="\t", index=False)
            epv_fg = round(fg_fit["n_primary_events"] / float(max(1, fg_fit["n_parameters"])), 2)
            if epv_fg < 10.0:
                warnings.append(
                    f"Fine-Gray EPV warning: {fg_fit['n_primary_events']} primary events across "
                    f"{fg_fit['n_parameters']} parameters (EPV = {epv_fg} < 10)"
                )
            fg_model = {
                "n_samples": fg_fit["n_samples"],
                "n_primary_events": fg_fit["n_primary_events"],
                "n_competing_events": fg_fit["n_competing_events"],
                "n_parameters": fg_fit["n_parameters"],
                "events_per_parameter": epv_fg,
                "converged": fg_fit["converged"],
                "reference_levels": ref_notes,
                "terms": fg_terms,
            }

    # Publication-ready 2-panel figure
    unit_lbl = args.time_label or infer_time_unit(t_arr, args.time_col).capitalize()
    fig, (ax_cif, ax_right) = plt.subplots(1, 2, figsize=(12.2, 4.9))

    if args.group_col:
        for idx_g, g in enumerate(levels):
            color = PALETTE[idx_g % len(PALETTE)]
            aj_g = group_aj[g]
            ax_cif.step(aj_g["time"], aj_g["cif1"], where="post", color=color, linewidth=2.0, label=f"{g} Primary CIF (n={aj_g['n']})")
            ax_cif.fill_between(aj_g["time"], aj_g["cif1_ci_low"], aj_g["cif1_ci_high"], step="post", color=color, alpha=0.14)
            if aj_g["n_competing"] > 0:
                ax_cif.step(aj_g["time"], aj_g["cif2"], where="post", color=color, linewidth=1.3, linestyle=":", alpha=0.85, label=f"{g} Competing CIF")
        ax_cif.step(aj_overall["time"], naive_arr, where="post", color="#475569", linewidth=1.3, linestyle="--", label="Overall Naive 1-KM")
        ttl = f"Cumulative Incidence ({args.group_col})"
        if grays and grays["pvalue"] is not None:
            ttl += f" | Gray p = {grays['pvalue']:.3e}"
        ax_cif.set_title(ttl, fontsize=10)
    else:
        ax_cif.step(aj_overall["time"], aj_overall["cif1"], where="post", color="#E11D48", linewidth=2.2, label=f"Primary CIF (n={n_ev1})")
        ax_cif.fill_between(aj_overall["time"], aj_overall["cif1_ci_low"], aj_overall["cif1_ci_high"], step="post", color="#E11D48", alpha=0.16)
        ax_cif.step(aj_overall["time"], aj_overall["cif2"], where="post", color="#2563EB", linewidth=2.0, linestyle="-.", label=f"Competing CIF (n={n_ev2})")
        ax_cif.fill_between(aj_overall["time"], aj_overall["cif2_ci_low"], aj_overall["cif2_ci_high"], step="post", color="#2563EB", alpha=0.12)
        ax_cif.step(aj_overall["time"], naive_arr, where="post", color="#475569", linewidth=1.4, linestyle="--", label="Naive 1 - KM (biased)")
        ax_cif.set_title("Aalen-Johansen CIF vs Naive 1 - KM", fontsize=10)

    ax_cif.set_ylim(-0.02, 1.02)
    ax_cif.set_xlabel(f"Follow-up Time ({unit_lbl})")
    ax_cif.set_ylabel("Cumulative Incidence Probability")
    ax_cif.legend(frameon=False, fontsize=8, loc="best")
    ax_cif.spines["top"].set_visible(False)
    ax_cif.spines["right"].set_visible(False)

    if fg_terms:
        plot_rows = fg_terms[::-1]
        y_pos = np.arange(len(plot_rows))
        l2_shr = np.array([np.log2(max(r["shr"], 1e-6)) for r in plot_rows])
        l2_lo = np.array([np.log2(max(r["ci_low"], 1e-6)) for r in plot_rows])
        l2_hi = np.array([np.log2(max(r["ci_high"], 1e-6)) for r in plot_rows])
        ax_right.axvline(0.0, color="#475569", linestyle="--", linewidth=1.0)
        for i, r in enumerate(plot_rows):
            col = "#E11D48" if (r["pvalue"] < 0.05 and r["shr"] > 1.0) else ("#2563EB" if (r["pvalue"] < 0.05 and r["shr"] < 1.0) else "#64748B")
            ax_right.plot([l2_lo[i], l2_hi[i]], [y_pos[i], y_pos[i]], color=col, linewidth=2.0)
            ax_right.scatter([l2_shr[i]], [y_pos[i]], color=col, s=46, zorder=4)
        xlim_m = max(1.2, float(np.max(np.abs(np.concatenate([l2_lo, l2_hi])))) * 1.35)
        ax_right.set_xlim(-xlim_m, xlim_m * 1.6)
        for i, r in enumerate(plot_rows):
            ax_right.text(xlim_m * 0.96, y_pos[i], f"sHR {r['shr']:.2f} ({r['ci_low']:.2f}-{r['ci_high']:.2f}), p={r['pvalue']:.2e}", va="center", ha="left", fontsize=7.8)
        ax_right.set_yticks(y_pos)
        ax_right.set_yticklabels([r["term"] for r in plot_rows], fontsize=8.2)
        ax_right.set_xlabel("log2 Subdistribution Hazard Ratio (95% CI)")
        ax_right.set_title("Fine-Gray Subdistribution HR", fontsize=10)
    else:
        tm = aj_overall["time"]
        f1_a = aj_overall["cif1"]
        f2_a = aj_overall["cif2"]
        s_a = np.maximum(0.0, 1.0 - f1_a - f2_a)
        ax_right.fill_between(tm, 0.0, f1_a, step="post", color="#E11D48", alpha=0.75, label="Primary Event")
        ax_right.fill_between(tm, f1_a, f1_a + f2_a, step="post", color="#2563EB", alpha=0.65, label="Competing Event")
        ax_right.fill_between(tm, f1_a + f2_a, f1_a + f2_a + s_a, step="post", color="#E2E8F0", alpha=0.85, label="Event-Free")
        ax_right.set_ylim(0.0, 1.0)
        ax_right.set_xlabel(f"Follow-up Time ({unit_lbl})")
        ax_right.set_ylabel("State Probability")
        ax_right.set_title("Stacked Cumulative Incidence", fontsize=10)
        ax_right.legend(frameon=False, fontsize=8, loc="best")

    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "competing_risk_cif.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {
        "time_col": args.time_col,
        "status_col": args.status_col,
        "primary_code": str(args.primary_code),
        "group_col": args.group_col,
        "covariates": fg_covs,
    }
    artifacts.append(
        write_artifact(fig_path, "clin_competing_risk", "Competing risks Aalen-Johansen CIF & Fine-Gray plot", inputs=(args.clinical,), parameters=params_meta)
    )
    artifacts.append(
        write_artifact(cif_path, "clin_competing_risk", "Aalen-Johansen CIF estimates", inputs=(args.clinical,), parameters=params_meta, shape=cif_df.shape)
    )
    if fg_path:
        artifacts.append(
            write_artifact(fg_path, "clin_competing_risk", "Fine-Gray subdistribution hazard ratios", inputs=(args.clinical,), parameters=params_meta, shape=(len(fg_terms), 9))
        )

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "time_col": args.time_col,
            "status_col": args.status_col,
            "primary_code": str(args.primary_code),
            "status_mapping": {"counts": raw_counts, "roles": role_map},
            "n_samples": int(len(sub)),
            "n_primary_events": n_ev1,
            "n_competing_events": n_ev2,
            "n_censored": n_cens,
            "overall_cif": {
                "final_cif_primary": round(float(aj_overall["cif1"][-1]), 4),
                "final_cif_competing": round(float(aj_overall["cif2"][-1]), 4),
                "final_event_free": round(float(aj_overall["surv"][-1]), 4),
            },
            "naive_1_minus_km_overestimation": bias_summary,
            "group_col": args.group_col,
            "groups": group_summaries,
            "grays_test": grays,
            "landmark_cif": landmark_rows,
            "fine_gray_model": fg_model,
            "cif_table": os.path.abspath(cif_path),
            "fine_gray_table": os.path.abspath(fg_path) if fg_path else None,
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
