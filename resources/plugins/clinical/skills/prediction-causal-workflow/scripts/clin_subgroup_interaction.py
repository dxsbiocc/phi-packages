#!/usr/bin/env python3
"""clin_subgroup_interaction: Subgroup forest plot analysis with formal
multiplicative interaction testing (P for interaction via LRT / Wald chi-square)
and BH-FDR multiplicity control for Cox PH, Logistic, and Linear models.
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
    bh_fdr,
    emit,
    encode_covariates,
    fail,
    fit_cox_weighted,
    fit_logistic_weighted,
    parse_binary_series,
    read_table,
    write_artifact,
)


def _fit_ols_nll(X_no_int: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    n = len(y)
    X = np.column_stack([np.ones(n, dtype=float), X_no_int])
    p = X.shape[1]
    inv_xtx = np.linalg.pinv(X.T @ X + 1e-12 * np.eye(p))
    beta = inv_xtx @ (X.T @ y)
    resid = y - X @ beta
    rss = max(1e-18, float(np.sum(resid ** 2)))
    df_resid = max(1, n - p)
    sigma2 = rss / float(df_resid)
    cov = sigma2 * inv_xtx
    # Gaussian NLL = 0.5 * n * log(rss / n) (+ const)
    nll = 0.5 * float(n) * math.log(rss / float(n))
    return {"beta": beta, "cov": cov, "se": np.sqrt(np.maximum(np.diag(cov), 1e-18)), "nll": nll}


def _fit_single_model(
    mtype: str,
    X_design: np.ndarray,
    y_arr: np.ndarray,
    tm_arr: np.ndarray | None,
) -> dict[str, Any]:
    """Fit Cox / Logistic / Linear model and return exposure effect at column 0 of X_design."""
    if mtype == "cox":
        assert tm_arr is not None
        fit = fit_cox_weighted(X_design, tm_arr, y_arr.astype(int))
        b0 = float(fit["beta"][0])
        se0 = float(fit["se_model"][0])
        z = b0 / max(se0, 1e-12)
        p = float(2.0 * stats.norm.sf(abs(z)))
        return {
            "beta": b0,
            "se": se0,
            "estimate": float(math.exp(np.clip(b0, -20.0, 20.0))),
            "ci_low": float(math.exp(np.clip(b0 - 1.959964 * se0, -20.0, 20.0))),
            "ci_high": float(math.exp(np.clip(b0 + 1.959964 * se0, -20.0, 20.0))),
            "pvalue": p,
            "nll": float(fit["nll"]),
            "cov_model": fit["cov_model"],
            "beta_vec": fit["beta"],
        }
    if mtype == "logistic":
        fit = fit_logistic_weighted(X_design, y_arr, l2_penalty=0.0)
        b0 = float(fit["beta"][1])
        se0 = float(fit["se_model"][1])
        z = b0 / max(se0, 1e-12)
        p = float(2.0 * stats.norm.sf(abs(z)))
        return {
            "beta": b0,
            "se": se0,
            "estimate": float(math.exp(np.clip(b0, -20.0, 20.0))),
            "ci_low": float(math.exp(np.clip(b0 - 1.959964 * se0, -20.0, 20.0))),
            "ci_high": float(math.exp(np.clip(b0 + 1.959964 * se0, -20.0, 20.0))),
            "pvalue": p,
            "nll": float(fit["nll"]),
            "cov_model": fit["cov_model"][1:, 1:],
            "beta_vec": fit["beta"][1:],
        }
    fit = _fit_ols_nll(X_design, y_arr)
    b0 = float(fit["beta"][1])
    se0 = float(fit["se"][1])
    df_r = max(1, len(y_arr) - (X_design.shape[1] + 1))
    t_stat = b0 / max(se0, 1e-12)
    p = float(2.0 * stats.t.sf(abs(t_stat), df=df_r))
    t_crit = float(stats.t.ppf(0.975, df=df_r))
    return {
        "beta": b0,
        "se": se0,
        "estimate": b0,
        "ci_low": b0 - t_crit * se0,
        "ci_high": b0 + t_crit * se0,
        "pvalue": p,
        "nll": float(fit["nll"]),
        "cov_model": fit["cov"][1:, 1:],
        "beta_vec": fit["beta"][1:],
    }


def _discretize_subgroup(series: pd.Series, col_name: str) -> tuple[pd.Series, bool, float | None]:
    """Return (string-level series with NaN, dichotomized_at_median, median_cutoff)."""
    cleaned = series.where(series.notna() & (series.astype(str).str.strip() != ""), np.nan)
    num = pd.to_numeric(cleaned, errors="coerce")
    n_nonmiss = int(cleaned.notna().sum())
    if n_nonmiss < 4:
        fail(f"subgroup '{col_name}' has < 4 non-missing values")
    # Continuous if >=85% numeric and >5 unique values
    if num.notna().sum() >= 0.85 * n_nonmiss and num.dropna().nunique() > 5:
        med = float(num.dropna().median())
        below = num < med
        if int(below.sum()) == 0 or int((num >= med).sum()) == 0:
            below = num <= med
            lbl_lo = f"<= {med:.3g}"
            lbl_hi = f"> {med:.3g}"
        else:
            lbl_lo = f"< {med:.3g}"
            lbl_hi = f">= {med:.3g}"
        out = pd.Series(np.nan, index=series.index, dtype=object)
        valid_m = num.notna()
        out.loc[valid_m] = np.where(below.loc[valid_m], lbl_lo, lbl_hi)
        return out, True, med

    out = cleaned.where(cleaned.isna(), cleaned.astype(str).str.strip())
    return out, False, None


def _drop_collinear_cols(X: np.ndarray) -> np.ndarray:
    keep = []
    for j in range(X.shape[1]):
        if float(np.std(X[:, j])) > 1e-12:
            keep.append(j)
    if not keep:
        return np.empty((X.shape[0], 0), dtype=float)
    return X[:, keep]


def _render_subgroup_forest(
    overall: dict[str, Any],
    subgroups: list[dict[str, Any]],
    effect_measure: str,
    title: str,
    out_path: str,
) -> None:
    is_ratio = effect_measure in ("HR", "OR")
    rows: list[dict[str, Any]] = [
        {
            "label": "Overall",
            "is_header": False,
            "is_overall": True,
            "n": overall["n"],
            "n_events": overall["n_events"],
            "est": overall["estimate"],
            "lo": overall["ci_low"],
            "hi": overall["ci_high"],
            "p_int": None,
        }
    ]
    for sg in subgroups:
        p_int = sg["p_interaction"]
        rows.append(
            {
                "label": sg["subgroup_var"] + (" (median split)" if sg["dichotomized_at_median"] else ""),
                "is_header": True,
                "is_overall": False,
                "n": None,
                "n_events": None,
                "est": None,
                "lo": None,
                "hi": None,
                "p_int": p_int,
            }
        )
        for lv in sg["levels"]:
            rows.append(
                {
                    "label": f"    {lv['level']}",
                    "is_header": False,
                    "is_overall": False,
                    "n": lv["n"],
                    "n_events": lv["n_events"],
                    "est": lv["estimate"],
                    "lo": lv["ci_low"],
                    "hi": lv["ci_high"],
                    "p_int": None,
                }
            )

    n_rows = len(rows)
    fig_h = max(4.2, 0.36 * n_rows + 1.8)
    fig = plt.figure(figsize=(11.5, fig_h))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.55, 1.35, 1.35], wspace=0.04)
    ax_left = fig.add_subplot(gs[0, 0])
    ax_mid = fig.add_subplot(gs[0, 1], sharey=ax_left)
    ax_right = fig.add_subplot(gs[0, 2], sharey=ax_left)

    y_pos = np.arange(n_rows)[::-1]
    ax_left.set_ylim(-0.8, n_rows + 0.35)

    # Header row at y = n_rows - 0.1
    hdr_y = n_rows - 0.15
    ax_left.text(0.0, hdr_y, "Subgroup", fontweight="bold", fontsize=9.0, va="center", ha="left")
    ax_left.text(0.96, hdr_y, "N (Events)" if overall["n_events"] is not None else "N", fontweight="bold", fontsize=9.0, va="center", ha="right")
    ax_mid.set_title(f"{effect_measure} (95% CI)", fontweight="bold", fontsize=9.5, pad=8)
    ax_right.text(0.04, hdr_y, f"{effect_measure} (95% CI)", fontweight="bold", fontsize=9.0, va="center", ha="left")
    ax_right.text(0.98, hdr_y, "P for interaction", fontweight="bold", fontsize=9.0, va="center", ha="right")

    # Determine x-axis limits for forest panel
    valid_ests = [r for r in rows if r["est"] is not None and r["lo"] is not None and r["hi"] is not None]
    if is_ratio:
        null_v = 1.0
        lo_vals = [max(0.05, r["lo"]) for r in valid_ests]
        hi_vals = [min(20.0, r["hi"]) for r in valid_ests]
        xmin = min(0.5, min(lo_vals) * 0.8) if lo_vals else 0.25
        xmax = max(2.0, max(hi_vals) * 1.2) if hi_vals else 4.0
        ax_mid.set_xscale("log")
        ax_mid.set_xlim(xmin, xmax)
    else:
        null_v = 0.0
        lo_vals = [r["lo"] for r in valid_ests]
        hi_vals = [r["hi"] for r in valid_ests]
        span = max(0.5, (max(hi_vals) - min(lo_vals)) if lo_vals else 1.0)
        xmin = min(-0.2 * span, (min(lo_vals) - 0.15 * span) if lo_vals else -1.0)
        xmax = max(0.2 * span, (max(hi_vals) + 0.15 * span) if hi_vals else 1.0)
        ax_mid.set_xlim(xmin, xmax)

    ax_mid.axvline(null_v, color="#475569", linestyle="--", linewidth=1.1, zorder=1)

    for idx, r in enumerate(rows):
        yp = y_pos[idx]
        if idx % 2 == 1 and not r["is_header"]:
            for ax in (ax_left, ax_mid, ax_right):
                ax.axhspan(yp - 0.45, yp + 0.45, color="#F8FAFC", zorder=0)

        fw = "bold" if (r["is_header"] or r["is_overall"]) else "normal"
        ax_left.text(0.0, yp, r["label"], fontsize=8.5, fontweight=fw, va="center", ha="left", color="#0F172A")
        if r["n"] is not None:
            cnt_txt = f"{r['n']} ({r['n_events']})" if r["n_events"] is not None else f"{r['n']}"
            ax_left.text(0.96, yp, cnt_txt, fontsize=8.3, fontweight=fw, va="center", ha="right", color="#334155")

        if r["est"] is not None and r["lo"] is not None and r["hi"] is not None:
            est_v, lo_v, hi_v = r["est"], r["lo"], r["hi"]
            lo_clip = max(xmin, lo_v)
            hi_clip = min(xmax, hi_v)
            if r["is_overall"]:
                # Diamond for overall effect
                if xmin <= est_v <= xmax:
                    ax_mid.fill(
                        [lo_clip, est_v, hi_clip, est_v],
                        [yp, yp + 0.25, yp, yp - 0.25],
                        color="#E11D48",
                        alpha=0.88,
                        zorder=4,
                    )
            else:
                ax_mid.plot([lo_clip, hi_clip], [yp, yp], color="#1E293B", linewidth=1.6, zorder=3)
                if xmin <= est_v <= xmax:
                    ax_mid.scatter([est_v], [yp], color="#2563EB", marker="s", s=38, zorder=4)
            ci_txt = f"{est_v:.2f} ({lo_v:.2f} to {hi_v:.2f})"
            ax_right.text(0.04, yp, ci_txt, fontsize=8.3, fontweight=fw, va="center", ha="left", color="#0F172A")
        elif not r["is_header"]:
            ax_right.text(0.04, yp, "NE (sparse)", fontsize=8.0, style="italic", va="center", ha="left", color="#64748B")

        if r["p_int"] is not None:
            p_str = "<0.001" if r["p_int"] < 0.001 else f"{r['p_int']:.3f}"
            p_col = "#E11D48" if r["p_int"] < 0.05 else "#0F172A"
            p_fw = "bold" if r["p_int"] < 0.05 else "normal"
            ax_right.text(0.98, yp, p_str, fontsize=8.5, fontweight=p_fw, va="center", ha="right", color=p_col)

    for ax in (ax_left, ax_right):
        ax.set_xlim(0.0, 1.0)
        ax.axis("off")
    ax_mid.set_yticks([])
    ax_mid.spines["top"].set_visible(False)
    ax_mid.spines["right"].set_visible(False)
    ax_mid.spines["left"].set_visible(False)
    ax_mid.set_xlabel(f"{effect_measure} (95% CI)", fontsize=9.0)

    fig.suptitle(title, fontsize=10.8, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.10)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--exposure_col", required=True)
    parser.add_argument("--subgroups", action="extend", nargs="+", default=[], required=True)
    parser.add_argument("--outcome_col", required=True)
    parser.add_argument("--model_type", choices=["auto", "cox", "logistic", "linear"], default="auto")
    parser.add_argument("--time_col", default=None)
    parser.add_argument("--adjust_covariates", action="extend", nargs="+", default=[])
    parser.add_argument("--exposed_label", default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.clinical)
    for req in (args.exposure_col, args.outcome_col):
        if req not in df.columns:
            fail(f"column '{req}' not found in clinical table: {list(df.columns)}")

    sub_cols = []
    for item in args.subgroups:
        for s in str(item).split(","):
            s = s.strip()
            if s and s not in sub_cols:
                sub_cols.append(s)
    if not sub_cols:
        fail("at least one subgroup variable is required")
    for s in sub_cols:
        if s not in df.columns:
            fail(f"subgroup column '{s}' not found in clinical table")

    adj_cols = []
    for item in args.adjust_covariates:
        for c in str(item).split(","):
            c = c.strip()
            if c and c != args.exposure_col and c not in adj_cols:
                adj_cols.append(c)
    for c in adj_cols:
        if c not in df.columns:
            fail(f"adjust_covariate '{c}' not found in clinical table")

    # Determine model_type
    mtype = args.model_type
    if mtype == "auto":
        if args.time_col:
            mtype = "cox"
        else:
            s_out = df[args.outcome_col].dropna()
            num_out = pd.to_numeric(s_out, errors="coerce")
            if s_out.astype(str).str.strip().nunique() == 2:
                mtype = "logistic"
            elif num_out.notna().sum() >= 0.85 * max(1, len(s_out)):
                mtype = "linear"
            else:
                fail(f"cannot infer model_type for '{args.outcome_col}'; pass --model_type")

    tm_s = None
    if mtype == "cox":
        if not args.time_col or args.time_col not in df.columns:
            fail("time_col is required and must exist in clinical table when model_type='cox'")
        tm_s = pd.to_numeric(df[args.time_col], errors="coerce")
        y_s, _ = parse_binary_series(df[args.outcome_col], args.outcome_col)
        valid_base = tm_s.notna() & (tm_s > 0) & y_s.notna()
        effect_measure = "HR"
    elif mtype == "logistic":
        y_s, _ = parse_binary_series(df[args.outcome_col], args.outcome_col)
        valid_base = y_s.notna()
        effect_measure = "OR"
    else:
        y_s = pd.to_numeric(df[args.outcome_col], errors="coerce")
        valid_base = y_s.notna()
        effect_measure = "beta"

    # Exposure parsing (binary vs continuous)
    exp_cleaned = df[args.exposure_col].where(
        df[args.exposure_col].notna() & (df[args.exposure_col].astype(str).str.strip() != ""),
        np.nan,
    )
    exp_distinct = exp_cleaned.dropna().astype(str).str.strip().nunique()
    if exp_distinct == 2 or args.exposed_label is not None:
        exp_s, exp_map = parse_binary_series(df[args.exposure_col], args.exposure_col, args.exposed_label)
        exp_type = "binary"
        exposed_lbl = exp_map["positive"]
        ref_lbl = exp_map["negative"]
    else:
        exp_s = pd.to_numeric(df[args.exposure_col], errors="coerce")
        if exp_s.dropna().nunique() < 2:
            fail(f"exposure_col '{args.exposure_col}' has < 2 distinct numeric values")
        exp_type = "continuous"
        exposed_lbl = None
        ref_lbl = None

    valid_base = valid_base & exp_s.notna()
    adj_df = None
    if adj_cols:
        adj_df, _, _ = encode_covariates(df, adj_cols)
        valid_base = valid_base & adj_df.notna().all(axis=1)

    warnings: list[str] = []
    n_drop_base = int((~valid_base).sum())
    if n_drop_base > 0:
        warnings.append(f"excluded {n_drop_base} row(s) with missing exposure, outcome, time, or adjustment covariates")

    n_base = int(valid_base.sum())
    if n_base < 10:
        fail(f"subgroup analysis requires >= 10 valid base observations (got {n_base})")

    # Fit overall model
    t_base = exp_s.loc[valid_base].to_numpy(dtype=float, copy=True)
    y_base = y_s.loc[valid_base].to_numpy(dtype=float, copy=True)
    tm_base = tm_s.loc[valid_base].to_numpy(dtype=float, copy=True) if tm_s is not None else None
    if adj_df is not None:
        X_adj_base = _drop_collinear_cols(adj_df.loc[valid_base].to_numpy(dtype=float, copy=True))
        X_overall = np.column_stack([t_base, X_adj_base]) if X_adj_base.shape[1] > 0 else t_base[:, None]
    else:
        X_overall = t_base[:, None]

    ov_fit = _fit_single_model(mtype, X_overall, y_base, tm_base)
    ov_events = int(np.sum(y_base == 1)) if mtype in ("cox", "logistic") else None
    overall_out = {
        "n": n_base,
        "n_events": ov_events,
        "estimate": round(ov_fit["estimate"], 4),
        "ci_low": round(ov_fit["ci_low"], 4),
        "ci_high": round(ov_fit["ci_high"], 4),
        "se": round(ov_fit["se"], 5),
        "pvalue": ov_fit["pvalue"],
    }

    # Evaluate each subgroup variable
    subgroup_results: list[dict[str, Any]] = []
    raw_p_ints: list[float] = []

    for sg_col in sub_cols:
        sg_ser, is_med, med_cut = _discretize_subgroup(df[sg_col], sg_col)
        if is_med:
            warnings.append(f"continuous subgroup '{sg_col}' auto-dichotomized at median ({med_cut:.4g})")

        valid_sg = valid_base & sg_ser.notna()
        n_sg = int(valid_sg.sum())
        levels = sorted(set(str(v) for v in sg_ser.loc[valid_sg]))
        if len(levels) < 2:
            warnings.append(f"subgroup '{sg_col}' has < 2 non-missing levels after base filtering; skipped")
            continue

        # Exclude sg_col from adjustment covariates if present
        sg_adj_cols = [c for c in adj_cols if c != sg_col]
        sg_adj_df = None
        if sg_adj_cols:
            sg_adj_df, _, _ = encode_covariates(df, sg_adj_cols)

        # Stratified estimates per level
        level_rows: list[dict[str, Any]] = []
        for lvl in levels:
            m_lvl = valid_sg & (sg_ser == lvl)
            n_l = int(m_lvl.sum())
            y_l = y_s.loc[m_lvl].to_numpy(dtype=float, copy=True)
            t_l = exp_s.loc[m_lvl].to_numpy(dtype=float, copy=True)
            tm_l = tm_s.loc[m_lvl].to_numpy(dtype=float, copy=True) if tm_s is not None else None
            ev_l = int(np.sum(y_l == 1)) if mtype in ("cox", "logistic") else None

            sparse = (
                n_l < 6
                or float(np.std(t_l)) < 1e-12
                or (mtype == "cox" and (ev_l is None or ev_l < 2))
                or (mtype == "logistic" and (ev_l is None or ev_l < 2 or (n_l - ev_l) < 2))
            )
            if sparse:
                warnings.append(f"stratum '{sg_col}={lvl}' is too sparse (n={n_l}, events={ev_l}); effect omitted")
                level_rows.append(
                    {
                        "level": lvl,
                        "n": n_l,
                        "n_events": ev_l,
                        "estimate": None,
                        "ci_low": None,
                        "ci_high": None,
                        "pvalue": None,
                    }
                )
                continue

            if sg_adj_df is not None:
                Xa_l = _drop_collinear_cols(sg_adj_df.loc[m_lvl].to_numpy(dtype=float, copy=True))
                X_l = np.column_stack([t_l, Xa_l]) if Xa_l.shape[1] > 0 else t_l[:, None]
            else:
                X_l = t_l[:, None]

            l_fit = _fit_single_model(mtype, X_l, y_l, tm_l)
            level_rows.append(
                {
                    "level": lvl,
                    "n": n_l,
                    "n_events": ev_l,
                    "estimate": round(l_fit["estimate"], 4),
                    "ci_low": round(l_fit["ci_low"], 4),
                    "ci_high": round(l_fit["ci_high"], 4),
                    "pvalue": l_fit["pvalue"],
                }
            )

        # Interaction test: Y ~ T + S + (T x S) + adj vs Y ~ T + S + adj
        t_sg = exp_s.loc[valid_sg].to_numpy(dtype=float, copy=True)
        y_sg = y_s.loc[valid_sg].to_numpy(dtype=float, copy=True)
        tm_sg = tm_s.loc[valid_sg].to_numpy(dtype=float, copy=True) if tm_s is not None else None
        sg_vals = sg_ser.loc[valid_sg].astype(str).to_numpy(copy=True)
        # L-1 dummy columns for subgroup levels[1:]
        S_dum = np.column_stack([(sg_vals == lvl).astype(float) for lvl in levels[1:]])
        Int_cols = S_dum * t_sg[:, None]
        df_int = S_dum.shape[1]

        if sg_adj_df is not None:
            Xa_sg = _drop_collinear_cols(sg_adj_df.loc[valid_sg].to_numpy(dtype=float, copy=True))
            X_red = np.column_stack([t_sg, S_dum, Xa_sg]) if Xa_sg.shape[1] > 0 else np.column_stack([t_sg, S_dum])
            X_full = np.column_stack([t_sg, S_dum, Int_cols, Xa_sg]) if Xa_sg.shape[1] > 0 else np.column_stack([t_sg, S_dum, Int_cols])
        else:
            X_red = np.column_stack([t_sg, S_dum])
            X_full = np.column_stack([t_sg, S_dum, Int_cols])

        red_fit = _fit_single_model(mtype, X_red, y_sg, tm_sg)
        full_fit = _fit_single_model(mtype, X_full, y_sg, tm_sg)

        lrt_chi2 = 2.0 * (red_fit["nll"] - full_fit["nll"])
        if math.isfinite(lrt_chi2) and lrt_chi2 >= 0.0:
            chisq_int = float(lrt_chi2)
        else:
            # Joint Wald chi-square fallback on the L-1 interaction terms
            int_slice = slice(1 + df_int, 1 + 2 * df_int)
            b_int = full_fit["beta_vec"][int_slice]
            c_int = full_fit["cov_model"][int_slice, int_slice]
            chisq_int = float(max(0.0, b_int.T @ np.linalg.pinv(c_int + 1e-12 * np.eye(df_int)) @ b_int))

        p_int = float(stats.chi2.sf(chisq_int, df=df_int))
        raw_p_ints.append(p_int)

        subgroup_results.append(
            {
                "subgroup_var": sg_col,
                "dichotomized_at_median": is_med,
                "median_cutoff": round(med_cut, 5) if med_cut is not None else None,
                "n_complete": n_sg,
                "chisq_interaction": round(chisq_int, 4),
                "df_interaction": df_int,
                "p_interaction": p_int,
                "padj_interaction": None,
                "levels": level_rows,
            }
        )

    if not subgroup_results:
        fail("no valid subgroup variables had >= 2 levels")

    if len(raw_p_ints) > 1:
        padj = bh_fdr(np.array(raw_p_ints, dtype=float))
        for sg_obj, qv in zip(subgroup_results, padj):
            sg_obj["padj_interaction"] = float(qv) if math.isfinite(qv) else None
        warnings.append(
            f"tested {len(subgroup_results)} subgroup interactions; inspect padj_interaction (BH-FDR) to guard against subgroup multiplicity"
        )
    else:
        subgroup_results[0]["padj_interaction"] = subgroup_results[0]["p_interaction"]

    # Export TSV summary table
    tsv_rows: list[dict[str, Any]] = [
        {
            "row_type": "overall",
            "subgroup_var": "Overall",
            "level": "All",
            "n": overall_out["n"],
            "n_events": overall_out["n_events"],
            "effect_measure": effect_measure,
            "estimate": overall_out["estimate"],
            "ci_low": overall_out["ci_low"],
            "ci_high": overall_out["ci_high"],
            "pvalue": overall_out["pvalue"],
            "p_interaction": None,
            "padj_interaction": None,
        }
    ]
    for sg in subgroup_results:
        for idx_l, lv in enumerate(sg["levels"]):
            tsv_rows.append(
                {
                    "row_type": "subgroup_level",
                    "subgroup_var": sg["subgroup_var"],
                    "level": lv["level"],
                    "n": lv["n"],
                    "n_events": lv["n_events"],
                    "effect_measure": effect_measure,
                    "estimate": lv["estimate"],
                    "ci_low": lv["ci_low"],
                    "ci_high": lv["ci_high"],
                    "pvalue": lv["pvalue"],
                    "p_interaction": sg["p_interaction"] if idx_l == 0 else None,
                    "padj_interaction": sg["padj_interaction"] if idx_l == 0 else None,
                }
            )
    sum_df = pd.DataFrame(tsv_rows)
    sum_path = os.path.join(args.output_dir, "subgroup_interaction.tsv")
    sum_df.to_csv(sum_path, sep="\t", index=False)

    # Render Forest Plot
    fig_path = os.path.join(args.output_dir, "subgroup_forest_plot.png")
    exp_desc = f"{args.exposure_col} ({exposed_lbl} vs {ref_lbl})" if exp_type == "binary" else f"{args.exposure_col} (continuous)"
    _render_subgroup_forest(
        overall_out,
        subgroup_results,
        effect_measure,
        f"Subgroup Forest Plot & Interaction Tests ({mtype.capitalize()}: {exp_desc} -> {args.outcome_col})",
        fig_path,
    )

    params_meta = {
        "exposure_col": args.exposure_col,
        "outcome_col": args.outcome_col,
        "model_type": mtype,
        "subgroups": [s["subgroup_var"] for s in subgroup_results],
        "adjust_covariates": adj_cols,
    }
    artifacts = [
        write_artifact(fig_path, "clin_subgroup_interaction", "Subgroup interaction forest plot", inputs=(args.clinical,), parameters=params_meta),
        write_artifact(sum_path, "clin_subgroup_interaction", "Subgroup stratified effects and P for interaction", inputs=(args.clinical,), parameters=params_meta, shape=sum_df.shape),
    ]

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "exposure_col": args.exposure_col,
            "exposure_type": exp_type,
            "exposed_label": exposed_lbl,
            "reference_label": ref_lbl,
            "outcome_col": args.outcome_col,
            "model_type": mtype,
            "effect_measure": effect_measure,
            "n_samples": n_base,
            "n_events": ov_events,
            "adjust_covariates": adj_cols,
            "overall": overall_out,
            "subgroups": subgroup_results,
            "summary_table": os.path.abspath(sum_path),
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
