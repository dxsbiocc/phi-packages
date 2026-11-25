"""Shared biostatistical & survival modeling primitives for clinical-survival (clin_*).

Includes:
- Event polarity normalization (0/1, Dead/Alive, 1:DECEASED/0:LIVING, True/False)
- Kaplan-Meier product-limit estimator with log-log Greenwood 95% CIs & Brookmeyer-Crowley median CIs
- Exact K-group hypergeometric Log-rank test
- Breslow partial-likelihood Cox proportional-hazards regression (damped Newton-Raphson, optional L2 ridge)
- Harrell's Concordance Index (C-index) & Grambsch-Therneau Schoenfeld residual PH test
- IPCW (Inverse Probability of Censoring Weighting) time-dependent ROC & AUC
- Contract-1.0.0 .phi-artifact.json descriptor writer
"""

from __future__ import annotations

import datetime
import hashlib
import json
import math
import os
import struct
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats


EVENT_TRUE_TOKENS = {
    "1", "1.0", "true", "t", "yes", "y", "dead", "deceased", "death",
    "event", "progressed", "recurrence", "relapse", "1:deceased", "1:progressed", "1:event",
}
EVENT_FALSE_TOKENS = {
    "0", "0.0", "false", "f", "no", "n", "alive", "living", "censored",
    "censor", "diseasefree", "disease-free", "0:living", "0:diseasefree", "0:censored",
}


def emit(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)


def fail(message: str, code: int = 1) -> None:
    print(json.dumps({"error": message}, indent=2, ensure_ascii=False))
    sys.exit(code)


def read_table(path: str) -> pd.DataFrame:
    if not os.path.isfile(path):
        fail(f"file not found: {path}")
    sep = "," if path.lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, low_memory=False)


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    n = int(valid.sum())
    if n == 0:
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1, dtype=float)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(n, dtype=float)
    res[order] = np.clip(adj, 0.0, 1.0)
    out[valid] = res
    return out


def parse_event_series(series: pd.Series, col_name: str = "event") -> tuple[pd.Series, dict[str, int]]:
    """Map clinical event encodings to {0, 1, NaN}."""
    mapped = []
    counts: dict[str, int] = {}
    unknown: set[str] = set()
    for val in series:
        if pd.isna(val):
            mapped.append(np.nan)
            continue
        raw = str(val).strip()
        if not raw or raw.lower() in ("na", "nan", "null", "none", "."):
            mapped.append(np.nan)
            continue
        low = raw.lower()
        counts[raw] = counts.get(raw, 0) + 1
        if low in EVENT_TRUE_TOKENS:
            mapped.append(1.0)
        elif low in EVENT_FALSE_TOKENS:
            mapped.append(0.0)
        else:
            unknown.add(raw)
            mapped.append(np.nan)
    if unknown:
        fail(
            f"unrecognized event values in '{col_name}': {sorted(unknown)[:8]}; "
            "expected 0/1, Alive/Dead, 0:LIVING/1:DECEASED, or Censored/Event"
        )
    return pd.Series(mapped, index=series.index, dtype=float), counts


def infer_time_unit(times: np.ndarray, col_name: str) -> str:
    low = col_name.lower()
    if "day" in low:
        return "days"
    if "month" in low or "_mo" in low or low.endswith("mos"):
        return "months"
    if "year" in low or "_yr" in low:
        return "years"
    finite = times[np.isfinite(times) & (times > 0)]
    if len(finite) == 0:
        return "unknown"
    q90 = float(np.percentile(finite, 90))
    if q90 > 240:
        return "days"
    if q90 > 20:
        return "months"
    return "years"


def join_molecular_features(
    clin: pd.DataFrame,
    sample_col: str,
    matrix_path: str | None,
    needed_features: list[str],
) -> pd.DataFrame:
    if not matrix_path:
        return clin
    mat_df = read_table(matrix_path)
    if mat_df.shape[1] < 2:
        fail(f"molecular matrix has < 2 columns: {matrix_path}")
    feat_col = mat_df.columns[0]
    mat_df[feat_col] = mat_df[feat_col].astype(str).str.strip()
    mat_df = mat_df.drop_duplicates(subset=[feat_col], keep="first").set_index(feat_col)
    mat_df.columns = [str(c).strip() for c in mat_df.columns]

    out = clin.copy()
    s_ids = out[sample_col].astype(str).str.strip()
    for feat in needed_features:
        if feat in out.columns:
            continue
        if feat not in mat_df.index:
            fail(f"feature '{feat}' not found in clinical columns or molecular matrix index")
        f_series = pd.to_numeric(mat_df.loc[feat], errors="coerce")
        f_map = f_series.to_dict()
        out[feat] = [f_map.get(s, np.nan) for s in s_ids]
    return out


def kaplan_meier(times: np.ndarray, events: np.ndarray) -> dict[str, Any]:
    """Product-limit KM curve with log-log Greenwood 95% CIs and Brookmeyer-Crowley median CIs."""
    t = np.asarray(times, dtype=float)
    e = np.asarray(events, dtype=int)
    order = np.argsort(t)
    t = t[order]
    e = e[order]
    n_total = len(t)
    n_events = int(e.sum())

    uniq_t = np.unique(t)
    curve_t = [0.0]
    curve_s = [1.0]
    curve_lo = [1.0]
    curve_hi = [1.0]
    curve_at_risk = [n_total]
    curve_events = [0]
    curve_censored = [0]

    surv = 1.0
    greenwood_sum = 0.0
    for ut in uniq_t:
        at_risk = int(np.sum(t >= ut))
        d_i = int(np.sum((t == ut) & (e == 1)))
        c_i = int(np.sum((t == ut) & (e == 0)))
        if d_i > 0 and at_risk > 0:
            surv *= 1.0 - (d_i / float(at_risk))
            if at_risk > d_i:
                greenwood_sum += d_i / float(at_risk * (at_risk - d_i))
            else:
                greenwood_sum = float("inf")

        if surv <= 0.0:
            lo, hi = 0.0, 0.0
        elif surv >= 1.0 or greenwood_sum == 0.0:
            lo, hi = 1.0, 1.0
        elif not math.isfinite(greenwood_sum):
            lo, hi = 0.0, 1.0
        else:
            log_s = math.log(max(surv, 1e-15))
            se_theta = math.sqrt(greenwood_sum) / abs(log_s)
            c_fac = math.exp(1.959964 * se_theta)
            lo = float(surv ** c_fac)
            hi = float(surv ** (1.0 / c_fac))

        curve_t.append(float(ut))
        curve_s.append(float(surv))
        curve_lo.append(float(np.clip(lo, 0.0, 1.0)))
        curve_hi.append(float(np.clip(hi, 0.0, 1.0)))
        curve_at_risk.append(at_risk)
        curve_events.append(d_i)
        curve_censored.append(c_i)

    arr_t = np.array(curve_t, dtype=float)
    arr_s = np.array(curve_s, dtype=float)
    arr_lo = np.array(curve_lo, dtype=float)
    arr_hi = np.array(curve_hi, dtype=float)

    med = float(arr_t[np.argmax(arr_s <= 0.5)]) if np.any(arr_s <= 0.5) else None
    med_lo = float(arr_t[np.argmax(arr_lo <= 0.5)]) if np.any(arr_lo <= 0.5) else None
    med_hi = float(arr_t[np.argmax(arr_hi <= 0.5)]) if np.any(arr_hi <= 0.5) else None

    return {
        "n": n_total,
        "events": n_events,
        "time": arr_t,
        "surv": arr_s,
        "ci_low": arr_lo,
        "ci_high": arr_hi,
        "n_risk": np.array(curve_at_risk, dtype=int),
        "n_event": np.array(curve_events, dtype=int),
        "n_censor": np.array(curve_censored, dtype=int),
        "median": round(med, 4) if med is not None else None,
        "median_ci_low": round(med_lo, 4) if med_lo is not None else None,
        "median_ci_high": round(med_hi, 4) if med_hi is not None else None,
    }


def km_surv_at(km: dict[str, Any], horizon: float) -> tuple[float, float, float, int]:
    """Evaluate step-function S(horizon), 95% CI, and number at risk just before horizon."""
    t = km["time"]
    idx = int(np.searchsorted(t, horizon, side="right") - 1)
    idx = max(0, min(idx, len(t) - 1))
    # Number at risk at horizon: check first time >= horizon
    r_idx = int(np.searchsorted(t[1:], horizon, side="left")) + 1
    n_risk = int(km["n_risk"][r_idx]) if r_idx < len(km["n_risk"]) else 0
    return float(km["surv"][idx]), float(km["ci_low"][idx]), float(km["ci_high"][idx]), n_risk


def logrank_test(times: np.ndarray, events: np.ndarray, groups: np.ndarray) -> dict[str, Any]:
    """Exact K-group Mantel-Haenszel / hypergeometric log-rank test."""
    t = np.asarray(times, dtype=float)
    e = np.asarray(events, dtype=int)
    g = np.asarray(groups, dtype=str)
    levels = sorted(set(g.tolist()))
    k = len(levels)
    if k < 2:
        return {"chisq": None, "df": None, "pvalue": None, "levels": levels}

    g_idx = np.array([levels.index(val) for val in g], dtype=int)
    event_times = np.unique(t[e == 1])
    u_vec = np.zeros(k, dtype=float)
    v_mat = np.zeros((k, k), dtype=float)

    for et in event_times:
        at_risk_mask = t >= et
        n_j = int(at_risk_mask.sum())
        if n_j < 2:
            continue
        event_mask = (t == et) & (e == 1)
        d_j = int(event_mask.sum())
        if d_j == 0:
            continue
        n_jg = np.bincount(g_idx[at_risk_mask], minlength=k).astype(float)
        d_jg = np.bincount(g_idx[event_mask], minlength=k).astype(float)
        p_j = n_jg / float(n_j)
        u_vec += d_jg - d_j * p_j
        factor = (d_j * (n_j - d_j)) / float(n_j - 1)
        v_mat += factor * (np.diag(p_j) - np.outer(p_j, p_j))

    u_sub = u_vec[:-1]
    v_sub = v_mat[:-1, :-1]
    chisq = float(u_sub.T @ np.linalg.pinv(v_sub) @ u_sub)
    chisq = max(0.0, chisq)
    df = k - 1
    pval = float(stats.chi2.sf(chisq, df=df))
    return {"chisq": round(chisq, 4), "df": df, "pvalue": pval, "levels": levels}


def harrell_c_index(times: np.ndarray, events: np.ndarray, risk_scores: np.ndarray) -> float:
    """Harrell's concordance index (higher risk_score -> shorter survival / earlier event)."""
    t = np.asarray(times, dtype=float)
    e = np.asarray(events, dtype=int)
    r = np.asarray(risk_scores, dtype=float)
    n = len(t)
    concordant = 0.0
    permissible = 0.0
    for i in range(n):
        if e[i] != 1:
            continue
        # j survived strictly longer than i's event time, or censored at same time
        mask = (t > t[i]) | ((t == t[i]) & (e == 0))
        if not np.any(mask):
            continue
        r_j = r[mask]
        permissible += float(len(r_j))
        concordant += float(np.sum(r[i] > r_j)) + 0.5 * float(np.sum(r[i] == r_j))
    if permissible == 0.0:
        return 0.5
    return float(concordant / permissible)


def fit_cox_breslow(
    X_raw: np.ndarray,
    times: np.ndarray,
    events: np.ndarray,
    *,
    l2_penalty: float = 0.0,
    max_iter: int = 60,
    tol: float = 1e-7,
) -> dict[str, Any]:
    """Fit Cox PH via Breslow partial likelihood with damped Newton-Raphson on standardized X,
    returning coefficients, SEs, covariance, C-index, and Schoenfeld PH diagnostics on original scale.
    """
    X = np.asarray(X_raw, dtype=float)
    t = np.asarray(times, dtype=float)
    e = np.asarray(events, dtype=int)
    n, p = X.shape
    n_events = int(e.sum())
    if n_events < 2:
        fail("Cox regression requires >= 2 observed events")

    means = X.mean(axis=0)
    sds = X.std(axis=0, ddof=0)
    sds_safe = np.where(sds > 1e-12, sds, 1.0)
    Z = (X - means) / sds_safe

    uniq_event_t = np.unique(t[e == 1])
    risk_sets = [np.flatnonzero(t >= ut) for ut in uniq_event_t]
    event_sets = [np.flatnonzero((t == ut) & (e == 1)) for ut in uniq_event_t]

    def neg_loglik_grad_hess(beta_z: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        eta = np.clip(Z @ beta_z, -50.0, 50.0)
        nll = 0.5 * l2_penalty * float(beta_z @ beta_z)
        grad = l2_penalty * beta_z.copy()
        hess = l2_penalty * np.eye(p, dtype=float)

        for r_idx, d_idx in zip(risk_sets, event_sets):
            d_m = len(d_idx)
            eta_r = eta[r_idx]
            shift = float(np.max(eta_r))
            w = np.exp(eta_r - shift)
            s0 = float(np.sum(w))
            Z_r = Z[r_idx]
            w_norm = w / s0
            x_bar = w_norm @ Z_r  # (p,)

            # Breslow negative log-likelihood contribution
            nll -= float(np.sum(eta[d_idx])) - d_m * (shift + math.log(s0))
            grad -= np.sum(Z[d_idx], axis=0) - d_m * x_bar
            # Weighted covariance matrix sum(w_norm * Z_r^T Z_r) - outer(x_bar, x_bar)
            Z_centered = Z_r - x_bar[None, :]
            v_bar = (Z_centered * w_norm[:, None]).T @ Z_centered
            hess += d_m * v_bar

        return nll, grad, hess

    beta_z = np.zeros(p, dtype=float)
    converged = False
    nll_cur, grad_cur, hess_cur = neg_loglik_grad_hess(beta_z)

    for _ in range(max_iter):
        # Stabilize Hessian if near-singular
        reg_hess = hess_cur + 1e-9 * np.eye(p)
        try:
            step = np.linalg.solve(reg_hess, grad_cur)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(reg_hess) @ grad_cur

        # Backtracking line search
        alpha = 1.0
        accepted = False
        for _ls in range(20):
            cand = beta_z - alpha * step
            if np.max(np.abs(cand)) <= 25.0:
                nll_new, grad_new, hess_new = neg_loglik_grad_hess(cand)
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

    cov_z = np.linalg.pinv(hess_cur + 1e-10 * np.eye(p))
    inv_sd = 1.0 / sds_safe
    beta = beta_z * inv_sd
    cov = cov_z * np.outer(inv_sd, inv_sd)
    se = np.sqrt(np.maximum(np.diag(cov), 1e-18))
    z_stat = beta / se
    pvals = 2.0 * stats.norm.sf(np.abs(z_stat))
    lin_pred = (X - means) @ beta
    c_idx = harrell_c_index(t, e, lin_pred)

    # Grambsch-Therneau Schoenfeld residual test against rank-transformed failure times
    eta_hat = np.clip(Z @ beta_z, -50.0, 50.0)
    sch_res_list = []
    sch_times = []
    for ut, r_idx, d_idx in zip(uniq_event_t, risk_sets, event_sets):
        eta_r = eta_hat[r_idx]
        w = np.exp(eta_r - float(np.max(eta_r)))
        w_norm = w / float(np.sum(w))
        x_bar = w_norm @ X[r_idx]
        for i in d_idx:
            sch_res_list.append(X[i] - x_bar)
            sch_times.append(float(ut))

    ph_pvalues = [None] * p
    global_ph_p = None
    if len(sch_res_list) >= 5:
        R_mat = np.vstack(sch_res_list)  # (d, p)
        g_t = stats.rankdata(sch_times).astype(float)
        g_centered = g_t - g_t.mean()
        g_ss = float(np.sum(g_centered ** 2))
        d_ev = len(sch_times)
        if g_ss > 0:
            # Scaled Schoenfeld test per variable and global: (g^T R)_j^2 / (d * V_jj * sum(g^2))
            gt_r = g_centered @ R_mat  # (p,)
            ph_pvalues = []
            for j in range(p):
                denom = d_ev * cov[j, j] * g_ss
                if denom > 1e-18:
                    chi_j = float((gt_r[j] ** 2) / denom)
                    ph_pvalues.append(float(stats.chi2.sf(chi_j, df=1)))
                else:
                    ph_pvalues.append(None)
            global_mat = d_ev * g_ss * cov
            global_chi = float(gt_r.T @ np.linalg.pinv(global_mat) @ gt_r)
            global_ph_p = float(stats.chi2.sf(max(0.0, global_chi), df=p))

    return {
        "beta": beta,
        "se": se,
        "z": z_stat,
        "pvalue": pvals,
        "hr": np.exp(np.clip(beta, -20.0, 20.0)),
        "ci_low": np.exp(np.clip(beta - 1.959964 * se, -20.0, 20.0)),
        "ci_high": np.exp(np.clip(beta + 1.959964 * se, -20.0, 20.0)),
        "c_index": round(c_idx, 4),
        "converged": converged,
        "n_samples": n,
        "n_events": n_events,
        "n_parameters": p,
        "events_per_parameter": round(n_events / float(max(1, p)), 2),
        "ph_pvalues": ph_pvalues,
        "global_ph_pvalue": global_ph_p,
        "means": means,
    }


def ipcw_time_dependent_roc(
    times: np.ndarray, events: np.ndarray, scores: np.ndarray, horizon: float
) -> dict[str, Any]:
    """Uno / Chiang cumulative/dynamic time-dependent ROC with KM censoring IPCW weights."""
    t = np.asarray(times, dtype=float)
    e = np.asarray(events, dtype=int)
    s = np.asarray(scores, dtype=float)

    # Estimate censoring survival distribution G(u) via KM where censoring (1 - e) is the event
    km_cens = kaplan_meier(t, 1 - e)
    g_times = km_cens["time"]
    g_surv = np.maximum(km_cens["surv"], 0.02)  # stabilize tail weights at 0.02 floor

    def g_at(u: float, left_limit: bool = False) -> float:
        side = "left" if left_limit else "right"
        idx = int(np.searchsorted(g_times, u, side=side) - 1)
        idx = max(0, min(idx, len(g_surv) - 1))
        return float(g_surv[idx])

    case_mask = (t <= horizon) & (e == 1)
    ctrl_mask = t > horizon
    n_cases = int(case_mask.sum())
    n_ctrls = int(ctrl_mask.sum())
    if n_cases < 3 or n_ctrls < 3:
        return {
            "horizon": horizon,
            "n_cases": n_cases,
            "n_controls": n_ctrls,
            "skipped": f"insufficient cases ({n_cases}) or controls ({n_ctrls}) at horizon {horizon}; need >= 3 each",
        }

    w_case = np.array([1.0 / g_at(ti, left_limit=True) for ti in t[case_mask]], dtype=float)
    w_ctrl = np.full(n_ctrls, 1.0 / g_at(horizon, left_limit=False), dtype=float)
    s_case = s[case_mask]
    s_ctrl = s[ctrl_mask]

    thresholds = np.unique(np.concatenate([[-np.inf], np.sort(s)[::-1], [np.inf]]))
    sum_wc = float(w_case.sum())
    sum_w0 = float(w_ctrl.sum())
    tpr = []
    fpr = []
    for thr in thresholds:
        tpr.append(float(w_case[s_case > thr].sum() / sum_wc))
        fpr.append(float(w_ctrl[s_ctrl > thr].sum() / sum_w0))

    fpr_arr = np.array(fpr, dtype=float)
    tpr_arr = np.array(tpr, dtype=float)
    order = np.argsort(fpr_arr, kind="mergesort")
    fpr_arr = fpr_arr[order]
    tpr_arr = tpr_arr[order]
    auc = float(np.trapezoid(tpr_arr, fpr_arr) if hasattr(np, "trapezoid") else np.trapz(tpr_arr, fpr_arr))

    return {
        "horizon": float(horizon),
        "n_cases": n_cases,
        "n_controls": n_ctrls,
        "auc": round(auc, 4),
        "fpr": fpr_arr,
        "tpr": tpr_arr,
    }


def write_artifact(
    path: str | os.PathLike[str],
    tool: str,
    title: str,
    *,
    skill: str = "survival-workflow",
    inputs: tuple[str | None, ...] | list[str | None] = (),
    parameters: dict[str, Any] | None = None,
    shape: tuple[int, int] | None = None,
) -> str:
    p = Path(path)
    ext = p.suffix.lower()
    figures = {".png": "image/png", ".pdf": "application/pdf", ".svg": "image/svg+xml"}
    tables = {".tsv": "text/tab-separated-values", ".csv": "text/csv"}
    if ext in figures:
        kind, media = "figure", figures[ext]
    elif ext in tables:
        kind, media = "table", tables[ext]
    else:
        raise ValueError(f"unsupported artifact extension: {ext}")

    descriptor: dict[str, Any] = {
        "contractVersion": "1.0.0",
        "kind": kind,
        "file": p.name,
        "mediaType": media,
        "title": title[:200],
        "provenance": {
            "createdAt": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "tool": tool,
            "skill": skill,
        },
    }
    if kind == "figure":
        block: dict[str, Any] = {"format": ext[1:]}
        if ext == ".png":
            with open(p, "rb") as fh:
                hdr = fh.read(24)
            if hdr[:8] == b"\x89PNG\r\n\x1a\n":
                w, h = struct.unpack(">II", hdr[16:24])
                block["widthPx"], block["heightPx"] = int(w), int(h)
        descriptor["figure"] = block
    elif shape is not None:
        descriptor["table"] = {"rows": int(shape[0]), "columns": int(shape[1])}

    cwd = Path.cwd().resolve()
    prov_in = []
    for raw in inputs:
        if not raw:
            continue
        ip = Path(raw)
        try:
            rel = ip.resolve().relative_to(cwd)
        except ValueError:
            continue
        entry: dict[str, Any] = {"path": str(rel)}
        if ip.is_file():
            entry["sha256"] = hashlib.sha256(ip.read_bytes()).hexdigest()
        prov_in.append(entry)
    if prov_in:
        descriptor["provenance"]["inputs"] = prov_in
    if parameters:
        descriptor["provenance"]["parameters"] = parameters

    sidecar = p.with_name(p.name + ".phi-artifact.json")
    tmp = sidecar.with_name(sidecar.name + ".tmp")
    tmp.write_text(json.dumps(descriptor, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, sidecar)
    return str(p)
