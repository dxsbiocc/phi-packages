"""Shared statistical, regression, and artifact helpers for prediction-causal-workflow (clin_*)."""

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


POS_TOKENS = {
    "1", "1.0", "true", "t", "yes", "y", "pos", "positive", "case",
    "disease", "tumor", "malignant", "dead", "deceased", "event", "treated",
    "treatment", "exposed", "1:deceased", "1:event", "1:case",
}
NEG_TOKENS = {
    "0", "0.0", "false", "f", "no", "n", "neg", "negative", "control",
    "ctrl", "normal", "healthy", "benign", "alive", "living", "censored",
    "untreated", "unexposed", "0:living", "0:censored", "0:control",
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


def wilson_ci(k: int, n: int, z: float = 1.959964) -> tuple[float, float, float]:
    if n <= 0:
        return float("nan"), float("nan"), float("nan")
    p = k / float(n)
    denom = 1.0 + (z * z) / n
    center = (p + (z * z) / (2.0 * n)) / denom
    half = (z * math.sqrt((p * (1.0 - p) + (z * z) / (4.0 * n)) / n)) / denom
    return p, max(0.0, center - half), min(1.0, center + half)


def parse_binary_series(
    series: pd.Series,
    col_name: str,
    positive_label: str | None = None,
) -> tuple[pd.Series, dict[str, str]]:
    """Map binary outcome/treatment column to {0.0, 1.0, NaN}."""
    cleaned = series.where(series.notna(), np.nan)
    raw_strs = []
    valid_idx = []
    for idx, val in cleaned.items():
        if pd.isna(val):
            continue
        s = str(val).strip()
        if not s or s.lower() in ("na", "nan", "null", "none", "."):
            continue
        raw_strs.append(s)
        valid_idx.append(idx)

    distinct = sorted(set(raw_strs))
    if len(distinct) != 2 and positive_label is None:
        fail(f"column '{col_name}' must have exactly 2 distinct non-missing values (found {len(distinct)}: {distinct[:6]})")
    if len(distinct) < 2:
        fail(f"column '{col_name}' has < 2 distinct non-missing values")

    out = pd.Series(np.nan, index=series.index, dtype=float)
    if positive_label is not None:
        pos_str = str(positive_label).strip()
        pos_low = pos_str.lower()
        pos_matches = [v for v in distinct if v == pos_str or v.lower() == pos_low]
        if not pos_matches:
            fail(f"positive_label '{positive_label}' not found in '{col_name}' values: {distinct}")
        pos_set = set(pos_matches)
        for idx, s in zip(valid_idx, raw_strs):
            out.loc[idx] = 1.0 if s in pos_set else 0.0
        neg_labels = [v for v in distinct if v not in pos_set]
        return out, {"positive": pos_matches[0], "negative": ", ".join(neg_labels)}

    # Check standard tokens
    pos_found = [v for v in distinct if v.lower() in POS_TOKENS]
    neg_found = [v for v in distinct if v.lower() in NEG_TOKENS]
    if len(pos_found) == 1 and len(neg_found) == 1:
        pos_val, neg_val = pos_found[0], neg_found[0]
    else:
        neg_val, pos_val = distinct[0], distinct[1]

    for idx, s in zip(valid_idx, raw_strs):
        out.loc[idx] = 1.0 if s == pos_val else 0.0
    return out, {"positive": pos_val, "negative": neg_val}


def encode_covariates(
    df: pd.DataFrame, covariates: list[str]
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Encode continuous & categorical covariates (one-hot vs first reference).

    Returns (encoded_df, reference_notes, is_continuous_flags_names).
    """
    parts = []
    notes = []
    cont_cols = []
    for cov in covariates:
        s = df[cov]
        num = pd.to_numeric(s, errors="coerce")
        is_cont = num.notna().sum() >= 0.85 * max(1, s.notna().sum()) and num.dropna().nunique() > 2
        if is_cont:
            col_s = num.astype(float).rename(cov)
            parts.append(col_s)
            cont_cols.append(cov)
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
    return pd.concat(parts, axis=1), notes, cont_cols


def fit_logistic_weighted(
    X_no_int: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray | None = None,
    l2_penalty: float = 0.0,
    max_iter: int = 80,
    tol: float = 1e-8,
) -> dict[str, Any]:
    """Fit weighted logistic regression with intercept via damped Newton-Raphson + sandwich SE."""
    X_in = np.asarray(X_no_int, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    n, p_in = X_in.shape
    w = np.ones(n, dtype=float) if weights is None else np.asarray(weights, dtype=float)

    # Standardize design columns internally for numerical stability, keeping intercept unpenalized
    means = X_in.mean(axis=0)
    sds = X_in.std(axis=0, ddof=0)
    sds_safe = np.where(sds > 1e-12, sds, 1.0)
    Z_in = (X_in - means) / sds_safe
    Z = np.column_stack([np.ones(n, dtype=float), Z_in])
    p = p_in + 1

    pen_diag = np.zeros(p, dtype=float)
    pen_diag[1:] = l2_penalty

    def nll_grad_hess(beta_z: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        eta = np.clip(Z @ beta_z, -30.0, 30.0)
        prob = 1.0 / (1.0 + np.exp(-eta))
        # stable log-loss: log(1+exp(eta)) - y*eta
        loss_vec = np.logaddexp(0.0, eta) - y_arr * eta
        nll = float(np.sum(w * loss_vec)) + 0.5 * float(np.sum(pen_diag * (beta_z ** 2)))
        resid = w * (prob - y_arr)
        grad = Z.T @ resid + pen_diag * beta_z
        wt = w * prob * (1.0 - prob)
        hess = (Z * wt[:, None]).T @ Z + np.diag(pen_diag)
        return nll, grad, hess

    beta_z = np.zeros(p, dtype=float)
    p0 = float(np.clip(np.sum(w * y_arr) / max(1e-12, float(np.sum(w))), 1e-4, 1.0 - 1e-4))
    beta_z[0] = math.log(p0 / (1.0 - p0))
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
        for _ls in range(25):
            cand = beta_z - alpha * step
            if np.max(np.abs(cand)) <= 30.0:
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

    eta_hat = np.clip(Z @ beta_z, -30.0, 30.0)
    prob_hat = 1.0 / (1.0 + np.exp(-eta_hat))
    score_i = Z * (w * (y_arr - prob_hat))[:, None]
    inv_h = np.linalg.pinv(hess_cur + 1e-10 * np.eye(p))
    cov_model_z = inv_h
    cov_rob_z = inv_h @ (score_i.T @ score_i) @ inv_h

    # Transform back to original X scale
    # beta_0 = beta_z0 - sum(beta_zj * mean_j / sd_j); beta_j = beta_zj / sd_j
    T_mat = np.zeros((p, p), dtype=float)
    T_mat[0, 0] = 1.0
    T_mat[0, 1:] = -means / sds_safe
    for j in range(p_in):
        T_mat[j + 1, j + 1] = 1.0 / sds_safe[j]

    beta = T_mat @ beta_z
    cov_model = T_mat @ cov_model_z @ T_mat.T
    cov_rob = T_mat @ cov_rob_z @ T_mat.T
    return {
        "beta": beta,
        "cov_model": cov_model,
        "cov_robust": cov_rob,
        "se_model": np.sqrt(np.maximum(np.diag(cov_model), 1e-18)),
        "se_robust": np.sqrt(np.maximum(np.diag(cov_rob), 1e-18)),
        "probs": prob_hat,
        "converged": converged,
    }


def fit_cox_weighted(
    X_raw: np.ndarray,
    times: np.ndarray,
    events: np.ndarray,
    weights: np.ndarray | None = None,
    max_iter: int = 60,
    tol: float = 1e-7,
) -> dict[str, Any]:
    """Fit weighted Breslow Cox PH regression with model & robust sandwich covariance."""
    X = np.asarray(X_raw, dtype=float)
    t = np.asarray(times, dtype=float)
    e = np.asarray(events, dtype=int)
    n, p = X.shape
    w = np.ones(n, dtype=float) if weights is None else np.asarray(weights, dtype=float)
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

    def nll_grad_hess(beta_z: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        eta = np.clip(Z @ beta_z, -50.0, 50.0)
        nll = 0.0
        grad = np.zeros(p, dtype=float)
        hess = np.zeros((p, p), dtype=float)
        for r_idx, d_idx in zip(risk_sets, event_sets):
            w_d = float(np.sum(w[d_idx]))
            eta_r = eta[r_idx]
            shift = float(np.max(eta_r))
            ew = w[r_idx] * np.exp(eta_r - shift)
            s0 = float(np.sum(ew))
            if s0 <= 1e-15:
                continue
            w_norm = ew / s0
            Z_r = Z[r_idx]
            x_bar = w_norm @ Z_r
            nll -= float(np.sum(w[d_idx] * eta[d_idx])) - w_d * (shift + math.log(s0))
            grad -= (w[d_idx, None] * Z[d_idx]).sum(axis=0) - w_d * x_bar
            Z_c = Z_r - x_bar[None, :]
            hess += w_d * ((Z_c * w_norm[:, None]).T @ Z_c)
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

    eta_hat = np.clip(Z @ beta_z, -50.0, 50.0)
    u_subj = np.zeros((n, p), dtype=float)
    for r_idx, d_idx in zip(risk_sets, event_sets):
        w_d = float(np.sum(w[d_idx]))
        ew = w[r_idx] * np.exp(eta_hat[r_idx] - float(np.max(eta_hat[r_idx])))
        s0 = float(np.sum(ew))
        if s0 <= 1e-15:
            continue
        w_norm = ew / s0
        Z_r = Z[r_idx]
        x_bar = w_norm @ Z_r
        u_subj[d_idx] += w[d_idx, None] * (Z[d_idx] - x_bar[None, :])
        u_subj[r_idx] -= (w_d * w_norm)[:, None] * (Z_r - x_bar[None, :])

    inv_h = np.linalg.pinv(hess_cur + 1e-10 * np.eye(p))
    cov_model_z = inv_h
    cov_rob_z = inv_h @ (u_subj.T @ u_subj) @ inv_h

    inv_sd = 1.0 / sds_safe
    scale_mat = np.outer(inv_sd, inv_sd)
    beta = beta_z * inv_sd
    cov_model = cov_model_z * scale_mat
    cov_rob = cov_rob_z * scale_mat
    return {
        "beta": beta,
        "cov_model": cov_model,
        "cov_robust": cov_rob,
        "se_model": np.sqrt(np.maximum(np.diag(cov_model), 1e-18)),
        "se_robust": np.sqrt(np.maximum(np.diag(cov_rob), 1e-18)),
        "converged": converged,
        "n_samples": n,
        "n_events": n_events,
    }


def write_artifact(
    path: str | os.PathLike[str],
    tool: str,
    title: str,
    *,
    skill: str = "prediction-causal-workflow",
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
