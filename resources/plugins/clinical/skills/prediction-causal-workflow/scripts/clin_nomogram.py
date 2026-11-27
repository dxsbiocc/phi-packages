#!/usr/bin/env python3
"""clin_nomogram: Multivariable Cox PH (with Breslow baseline survival at
landmark horizons) or Logistic clinical prediction nomogram builder: constructs a
0-100 point scorecard, exact Total Points -> Risk / Survival lookup table, and
publication-ready alignment nomogram chart.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pred_common import (
    emit,
    fail,
    fit_cox_weighted,
    fit_logistic_weighted,
    parse_binary_series,
    read_table,
    write_artifact,
)


def _infer_time_unit(times: np.ndarray, col_name: str) -> str:
    low = col_name.lower()
    if "day" in low:
        return "days"
    if "month" in low or "_mo" in low or low.endswith("mos"):
        return "months"
    if "year" in low or "_yr" in low:
        return "years"
    finite = times[np.isfinite(times) & (times > 0)]
    if len(finite) == 0:
        return "time"
    q90 = float(np.percentile(finite, 90))
    if q90 > 240:
        return "days"
    if q90 > 20:
        return "months"
    return "years"


def _harrell_c_index(times: np.ndarray, events: np.ndarray, risk_scores: np.ndarray) -> float:
    n = len(times)
    conc = 0.0
    perm = 0.0
    for i in range(n):
        if events[i] != 1:
            continue
        mask = (times > times[i]) | ((times == times[i]) & (events == 0))
        if not np.any(mask):
            continue
        rj = risk_scores[mask]
        perm += float(len(rj))
        conc += float(np.sum(risk_scores[i] > rj)) + 0.5 * float(np.sum(risk_scores[i] == rj))
    return float(conc / perm) if perm > 0.0 else 0.5


def _breslow_baseline_survival(
    times: np.ndarray,
    events: np.ndarray,
    eta_centered: np.ndarray,
    horizons: list[float],
) -> dict[float, float]:
    """Compute Breslow cumulative baseline hazard H0(t) and S0(t) = exp(-H0(t)) at landmark horizons."""
    uniq_et = np.unique(times[events == 1])
    exp_eta = np.exp(np.clip(eta_centered, -30.0, 30.0))
    jumps = []
    for ut in uniq_et:
        d_k = float(np.sum((times == ut) & (events == 1)))
        denom = float(np.sum(exp_eta[times >= ut]))
        if denom > 1e-15:
            jumps.append((float(ut), d_k / denom))

    out: dict[float, float] = {}
    for h in horizons:
        h0 = sum(dh for t_k, dh in jumps if t_k <= h)
        out[float(h)] = float(np.clip(math.exp(-h0), 1e-6, 1.0))
    return out


def _encode_predictors_structured(
    df: pd.DataFrame, pred_cols: list[str]
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """Encode predictors while preserving grouping per original variable for nomogram construction."""
    parts: list[pd.Series] = []
    specs: list[dict[str, Any]] = []
    col_idx = 0

    for pred in pred_cols:
        s = df[pred]
        num = pd.to_numeric(s, errors="coerce")
        is_cont = num.notna().sum() >= 0.85 * max(1, s.notna().sum()) and num.dropna().nunique() > 2
        if is_cont:
            cs = num.astype(float).rename(pred)
            parts.append(cs)
            specs.append(
                {
                    "predictor": pred,
                    "type": "continuous",
                    "term_names": [pred],
                    "design_indices": [col_idx],
                }
            )
            col_idx += 1
        else:
            cats = s.where(s.notna() & (s.astype(str).str.strip() != ""), np.nan).astype(object)
            levels = sorted(set(str(x).strip() for x in cats.dropna()))
            if len(levels) < 2:
                fail(f"categorical predictor '{pred}' has < 2 non-missing levels")
            ref = levels[0]
            t_names = []
            d_indices = []
            for lvl in levels[1:]:
                tname = f"{pred} [{lvl} vs {ref}]"
                ind = pd.Series(
                    np.where(cats.isna(), np.nan, (cats.astype(str).str.strip() == lvl).astype(float)),
                    index=df.index,
                    name=tname,
                )
                parts.append(ind)
                t_names.append(tname)
                d_indices.append(col_idx)
                col_idx += 1
            specs.append(
                {
                    "predictor": pred,
                    "type": "categorical",
                    "levels": levels,
                    "ref": ref,
                    "term_names": t_names,
                    "design_indices": d_indices,
                }
            )
    return pd.concat(parts, axis=1), specs


def _nice_ticks(vmin: float, vmax: float, n_bins: int = 5) -> np.ndarray:
    loc = MaxNLocator(nbins=n_bins, steps=[1, 2, 2.5, 5, 10])
    ticks = loc.tick_values(vmin, vmax)
    ticks = ticks[(ticks >= vmin - 1e-9) & (ticks <= vmax + 1e-9)]
    if len(ticks) < 3:
        ticks = np.linspace(vmin, vmax, 5)
    return np.asarray(ticks, dtype=float)


def _render_nomogram(
    predictor_items: list[dict[str, Any]],
    tp_max: float,
    lp_zero: float,
    pts_per_lp: float,
    mtype: str,
    landmarks: list[float],
    s0_map: dict[float, float],
    time_label: str | None,
    title: str,
    out_path: str,
) -> None:
    n_pred = len(predictor_items)
    n_risk_axes = 1 if mtype == "logistic" else len(landmarks)
    # Rows from top to bottom: Points, each predictor, separator gap, Total Points, Linear Predictor, Risk/Surv axes
    total_rows = 1 + n_pred + 3 + n_risk_axes
    fig_h = max(5.2, 0.52 * total_rows + 1.6)
    fig, ax = plt.subplots(figsize=(11.5, fig_h))

    y_cur = float(total_rows)
    x_left_label = -22.0
    ax.set_xlim(x_left_label - 3.0, 106.0)
    ax.set_ylim(-0.6, total_rows + 0.8)

    def _draw_ruler(
        y_val: float,
        label: str,
        x0: float,
        x1: float,
        tick_xs: list[float],
        tick_lbls: list[str],
        *,
        color: str = "#1E293B",
        label_color: str = "#0F172A",
        lw: float = 1.4,
        bold: bool = False,
        alternate_stagger: bool = False,
    ) -> None:
        ax.text(
            x_left_label,
            y_val,
            label,
            va="center",
            ha="left",
            fontsize=8.8,
            fontweight="bold" if bold else "normal",
            color=label_color,
        )
        ax.plot([x0, x1], [y_val, y_val], color=color, linewidth=lw)
        for i_t, (tx, tl) in enumerate(zip(tick_xs, tick_lbls)):
            up = True
            if alternate_stagger and len(tick_xs) >= 2 and i_t % 2 == 1:
                up = False
            dy_tick = 0.12 if up else -0.12
            dy_txt = 0.18 if up else -0.22
            va = "bottom" if up else "top"
            ax.plot([tx, tx], [y_val, y_val + dy_tick], color=color, linewidth=1.0)
            ax.text(tx, y_val + dy_txt, tl, va=va, ha="center", fontsize=7.4, color=color)

    # 1) Top Points ruler (0..100)
    pts_ticks = list(range(0, 101, 10))
    _draw_ruler(
        y_cur,
        "Points",
        0.0,
        100.0,
        [float(v) for v in pts_ticks],
        [str(v) for v in pts_ticks],
        color="#2563EB",
        label_color="#1E3A8A",
        lw=1.8,
        bold=True,
    )
    y_cur -= 1.0

    # 2) Predictor axes
    for p_item in predictor_items:
        p_ticks = p_item["ticks"]
        xs = [float(t["points"]) for t in p_ticks]
        lbls = [str(t["label"]) for t in p_ticks]
        x_hi = max(xs) if xs else 0.0
        stagger = p_item["type"] == "categorical" and len(xs) >= 2
        _draw_ruler(
            y_cur,
            p_item["predictor"],
            0.0,
            max(0.5, x_hi),
            xs,
            lbls,
            color="#334155",
            lw=1.3,
            alternate_stagger=stagger,
        )
        y_cur -= 1.0

    # Divider line
    ax.axhline(y_cur + 0.35, xmin=0.02, xmax=0.98, color="#CBD5E1", linestyle="--", linewidth=0.9)
    y_cur -= 0.45

    # 3) Total Points ruler mapped onto [0, 100] canvas
    tp_nice = _nice_ticks(0.0, max(1.0, tp_max), n_bins=7)
    tp_xs = [float(tp / max(tp_max, 1e-9) * 100.0) for tp in tp_nice]
    _draw_ruler(
        y_cur,
        "Total Points",
        0.0,
        100.0,
        tp_xs,
        [f"{tp:g}" for tp in tp_nice],
        color="#0F172A",
        lw=1.7,
        bold=True,
    )
    y_cur -= 1.0

    # 4) Linear Predictor ruler aligned directly below Total Points
    lp_min = lp_zero
    lp_max = lp_zero + tp_max / pts_per_lp
    lp_nice = _nice_ticks(lp_min, lp_max, n_bins=6)
    lp_xs = [float(((lp - lp_zero) * pts_per_lp) / max(tp_max, 1e-9) * 100.0) for lp in lp_nice]
    _draw_ruler(
        y_cur,
        "Linear Predictor",
        0.0,
        100.0,
        lp_xs,
        [f"{lp:.2g}" for lp in lp_nice],
        color="#475569",
        lw=1.3,
    )
    y_cur -= 1.0

    # 5) Predicted probability / landmark survival rulers aligned to Total Points [0, 100]
    prob_grid = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    if mtype == "logistic":
        p_xs = []
        p_lbls = []
        for pr in prob_grid:
            lp_val = math.log(pr / (1.0 - pr))
            tp_val = (lp_val - lp_zero) * pts_per_lp
            if 0.0 <= tp_val <= tp_max:
                p_xs.append(float(tp_val / max(tp_max, 1e-9) * 100.0))
                p_lbls.append(f"{pr:g}")
        x0_p = min(p_xs) if p_xs else 0.0
        x1_p = max(p_xs) if p_xs else 100.0
        _draw_ruler(
            y_cur,
            "Event Probability",
            x0_p,
            x1_p,
            p_xs,
            p_lbls,
            color="#E11D48",
            label_color="#9F1239",
            lw=1.8,
            bold=True,
        )
    else:
        t_unit = time_label or "time"
        for h in landmarks:
            s0_h = s0_map[float(h)]
            h0_h = -math.log(max(s0_h, 1e-9))
            s_xs = []
            s_lbls = []
            for sv in prob_grid:
                # sv = exp(-H0 * exp(lp)) -> exp(lp) = -ln(sv) / H0
                lp_val = math.log(max(1e-9, -math.log(sv) / max(h0_h, 1e-12)))
                tp_val = (lp_val - lp_zero) * pts_per_lp
                if 0.0 <= tp_val <= tp_max:
                    s_xs.append(float(tp_val / max(tp_max, 1e-9) * 100.0))
                    s_lbls.append(f"{sv:g}")
            x0_s = min(s_xs) if s_xs else 0.0
            x1_s = max(s_xs) if s_xs else 100.0
            _draw_ruler(
                y_cur,
                f"{h:g}-{t_unit} Survival",
                x0_s,
                x1_s,
                s_xs,
                s_lbls,
                color="#E11D48",
                label_color="#9F1239",
                lw=1.7,
                bold=True,
            )
            y_cur -= 1.0

    ax.axis("off")
    ax.set_title(title, fontsize=10.5, fontweight="bold", pad=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--predictors", action="extend", nargs="+", default=[], required=True)
    parser.add_argument("--outcome_col", required=True)
    parser.add_argument("--model_type", choices=["auto", "cox", "logistic"], default="auto")
    parser.add_argument("--time_col", default=None)
    parser.add_argument("--landmarks", action="extend", nargs="+", default=[])
    parser.add_argument("--time_label", default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.clinical)
    if args.outcome_col not in df.columns:
        fail(f"outcome_col '{args.outcome_col}' not found in clinical table")

    pred_cols: list[str] = []
    for item in args.predictors:
        for p in str(item).split(","):
            p = p.strip()
            if p and p not in pred_cols:
                pred_cols.append(p)
    if len(pred_cols) < 2 or len(pred_cols) > 12:
        fail(f"nomogram requires 2 to 12 predictors (got {len(pred_cols)})")
    for p in pred_cols:
        if p not in df.columns:
            fail(f"predictor '{p}' not found in clinical table")

    mtype = args.model_type
    if mtype == "auto":
        mtype = "cox" if args.time_col else "logistic"

    tm_s = None
    if mtype == "cox":
        if not args.time_col or args.time_col not in df.columns:
            fail("time_col is required and must exist in clinical table when model_type='cox'")
        tm_s = pd.to_numeric(df[args.time_col], errors="coerce")
        y_s, _ = parse_binary_series(df[args.outcome_col], args.outcome_col)
        valid = tm_s.notna() & (tm_s > 0) & y_s.notna()
    else:
        y_s, _ = parse_binary_series(df[args.outcome_col], args.outcome_col)
        valid = y_s.notna()

    X_df, pred_specs = _encode_predictors_structured(df, pred_cols)
    valid = valid & X_df.notna().all(axis=1)

    warnings: list[str] = []
    n_drop = int((~valid).sum())
    if n_drop > 0:
        warnings.append(f"excluded {n_drop} record(s) with missing outcome, time, or predictors")

    X = X_df.loc[valid].to_numpy(dtype=float, copy=True)
    y = y_s.loc[valid].to_numpy(dtype=float, copy=True)
    n_samples = len(y)
    n_events = int(np.sum(y == 1))
    n_params = X.shape[1]
    if n_samples < 20 or n_events < 6:
        fail(f"nomogram requires >= 20 samples and >= 6 events (got n={n_samples}, events={n_events})")

    epv = round(n_events / float(max(1, n_params)), 2)
    if epv < 10.0:
        warnings.append(f"low Events-Per-Variable (EPV={epv} < 10); apparent discrimination may be optimistic — validate externally")

    landmarks: list[float] = []
    s0_map: dict[float, float] = {}
    t_label: str | None = None

    if mtype == "cox":
        assert tm_s is not None
        tm = tm_s.loc[valid].to_numpy(dtype=float, copy=True)
        t_label = args.time_label or _infer_time_unit(tm, args.time_col)
        fit = fit_cox_weighted(X, tm, y.astype(int))
        beta = fit["beta"]
        se = fit["se_model"]
        mean_x = X.mean(axis=0)
        eta_bar = float(mean_x @ beta)
        eta_centered = X @ beta - eta_bar
        disc_metric = "c_index"
        disc_val = _harrell_c_index(tm, y.astype(int), eta_centered)

        # Parse or infer landmarks
        raw_lm = []
        for item in args.landmarks:
            for tok in str(item).split(","):
                tok = tok.strip()
                if tok:
                    raw_lm.append(float(tok))
        max_t = float(np.max(tm))
        if raw_lm:
            landmarks = [h for h in raw_lm if 0.0 < h <= max_t]
            dropped_lm = [h for h in raw_lm if h not in landmarks]
            if dropped_lm:
                warnings.append(f"landmarks {dropped_lm} outside observed follow-up (0, {max_t:.2f}] omitted")
        if not landmarks:
            q_lms = np.unique(np.round(np.percentile(tm, [25.0, 50.0, 75.0]), 1))
            landmarks = [float(h) for h in q_lms if h > 0.0]

        s0_map = _breslow_baseline_survival(tm, y.astype(int), eta_centered, landmarks)
        effect_ratio_name = "HR"
    else:
        fit = fit_logistic_weighted(X, y, l2_penalty=0.0)
        if not fit["converged"]:
            fit = fit_logistic_weighted(X, y, l2_penalty=0.05)
            warnings.append("unpenalized logistic regression hit near-separation; refitted with tiny L2 ridge (0.05)")
        beta_0 = float(fit["beta"][0])
        beta = fit["beta"][1:]
        se = fit["se_model"][1:]
        disc_metric = "roc_auc"
        disc_val = float(roc_auc_score(y, fit["probs"]))
        effect_ratio_name = "OR"

    # Model terms summary
    term_names = list(X_df.columns)
    model_terms = []
    for idx_t, tname in enumerate(term_names):
        b_j = float(beta[idx_t])
        se_j = float(se[idx_t])
        z_j = b_j / max(se_j, 1e-12)
        p_j = float(2.0 * stats.norm.sf(abs(z_j)))
        model_terms.append(
            {
                "term": tname,
                "beta": round(b_j, 5),
                "se": round(se_j, 5),
                "effect_measure": effect_ratio_name,
                "effect_ratio": round(float(math.exp(np.clip(b_j, -20.0, 20.0))), 4),
                "ci_low": round(float(math.exp(np.clip(b_j - 1.959964 * se_j, -20.0, 20.0))), 4),
                "ci_high": round(float(math.exp(np.clip(b_j + 1.959964 * se_j, -20.0, 20.0))), 4),
                "pvalue": p_j,
            }
        )

    # Compute linear-predictor spans per predictor to scale max span -> 100 points
    spans = []
    c_mins = []
    c_maxs = []
    cont_bounds: dict[str, tuple[float, float, np.ndarray]] = {}

    for spec in pred_specs:
        pname = spec["predictor"]
        d_idxs = spec["design_indices"]
        if spec["type"] == "continuous":
            b_j = float(beta[d_idxs[0]])
            col_vals = X[:, d_idxs[0]]
            vmin = float(np.min(col_vals))
            vmax = float(np.max(col_vals))
            ticks_v = _nice_ticks(vmin, vmax, n_bins=5)
            # Ensure span covers plotted tick range and observed range cleanly
            r_lo = min(vmin, float(ticks_v[0]))
            r_hi = max(vmax, float(ticks_v[-1]))
            cont_bounds[pname] = (r_lo, r_hi, ticks_v)
            c_lo = b_j * r_lo
            c_hi = b_j * r_hi
            c_min = min(c_lo, c_hi)
            c_max = max(c_lo, c_hi)
        else:
            contribs = [0.0] + [float(beta[i]) for i in d_idxs]
            c_min = min(contribs)
            c_max = max(contribs)
        spans.append(max(0.0, c_max - c_min))
        c_mins.append(c_min)
        c_maxs.append(c_max)

    max_span = max(spans)
    if max_span <= 1e-12:
        fail("all predictor effects are zero; cannot construct nomogram scale")

    pts_per_lp = 100.0 / max_span
    c_min_total = sum(c_mins)
    tp_max = sum((cmx - cmn) * pts_per_lp for cmn, cmx in zip(c_mins, c_maxs))
    lp_zero = (c_min_total - eta_bar) if mtype == "cox" else (beta_0 + c_min_total)

    # Build predictor point tables
    predictor_items: list[dict[str, Any]] = []
    scorecard_rows: list[dict[str, Any]] = []
    for spec, c_min in zip(pred_specs, c_mins):
        pname = spec["predictor"]
        d_idxs = spec["design_indices"]
        ticks_out = []
        if spec["type"] == "continuous":
            b_j = float(beta[d_idxs[0]])
            _, _, ticks_v = cont_bounds[pname]
            for tv in ticks_v:
                c_val = b_j * float(tv)
                pts = max(0.0, (c_val - c_min) * pts_per_lp)
                ticks_out.append(
                    {
                        "label": f"{tv:g}",
                        "value": round(float(tv), 5),
                        "beta_contribution": round(c_val, 5),
                        "points": round(pts, 2),
                    }
                )
                scorecard_rows.append(
                    {
                        "predictor": pname,
                        "type": "continuous",
                        "value_or_level": round(float(tv), 5),
                        "beta": round(b_j, 5),
                        "points_0_to_100": round(pts, 2),
                    }
                )
        else:
            lvls = spec["levels"]
            betas_lvl = [0.0] + [float(beta[i]) for i in d_idxs]
            for lv, b_val in zip(lvls, betas_lvl):
                pts = max(0.0, (b_val - c_min) * pts_per_lp)
                ticks_out.append(
                    {
                        "label": lv,
                        "value": lv,
                        "beta_contribution": round(b_val, 5),
                        "points": round(pts, 2),
                    }
                )
                scorecard_rows.append(
                    {
                        "predictor": pname,
                        "type": "categorical",
                        "value_or_level": lv,
                        "beta": round(b_val, 5),
                        "points_0_to_100": round(pts, 2),
                    }
                )
        max_pts_j = max(t["points"] for t in ticks_out)
        predictor_items.append(
            {
                "predictor": pname,
                "type": spec["type"],
                "max_points": round(max_pts_j, 2),
                "ticks": ticks_out,
            }
        )

    # Total Points -> Risk / Landmark Survival lookup table
    tp_grid = np.linspace(0.0, tp_max, 51)
    lookup_rows = []
    for tp in tp_grid:
        lp = lp_zero + float(tp) / pts_per_lp
        row: dict[str, Any] = {
            "total_points": round(float(tp), 2),
            "linear_predictor": round(lp, 5),
        }
        if mtype == "logistic":
            prob = 1.0 / (1.0 + math.exp(-float(np.clip(lp, -30.0, 30.0))))
            row["predicted_event_prob"] = round(prob, 5)
        else:
            exp_lp = math.exp(float(np.clip(lp, -20.0, 20.0)))
            for h in landmarks:
                s0_h = s0_map[float(h)]
                surv_h = float(np.clip(s0_h ** exp_lp, 0.0, 1.0))
                row[f"surv_at_{h:g}"] = round(surv_h, 5)
        lookup_rows.append(row)

    scorecard_df = pd.DataFrame(scorecard_rows)
    lookup_df = pd.DataFrame(lookup_rows)
    scorecard_path = os.path.join(args.output_dir, "nomogram_points_scorecard.tsv")
    lookup_path = os.path.join(args.output_dir, "nomogram_total_points_lookup.tsv")
    scorecard_df.to_csv(scorecard_path, sep="\t", index=False)
    lookup_df.to_csv(lookup_path, sep="\t", index=False)

    fig_path = os.path.join(args.output_dir, "clinical_nomogram.png")
    _render_nomogram(
        predictor_items,
        tp_max,
        lp_zero,
        pts_per_lp,
        mtype,
        landmarks,
        s0_map,
        t_label,
        f"Clinical Prediction Nomogram ({mtype.capitalize()}: {args.outcome_col}, {disc_metric.upper()}={disc_val:.3f})",
        fig_path,
    )

    params_meta = {
        "outcome_col": args.outcome_col,
        "model_type": mtype,
        "predictors": pred_cols,
        "landmarks": landmarks,
    }
    artifacts = [
        write_artifact(fig_path, "clin_nomogram", f"Clinical nomogram ({mtype})", inputs=(args.clinical,), parameters=params_meta),
        write_artifact(scorecard_path, "clin_nomogram", "Nomogram 0-100 points scorecard", inputs=(args.clinical,), parameters=params_meta, shape=scorecard_df.shape),
        write_artifact(lookup_path, "clin_nomogram", "Nomogram Total Points risk/survival lookup", inputs=(args.clinical,), parameters=params_meta, shape=lookup_df.shape),
    ]

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "outcome_col": args.outcome_col,
            "model_type": mtype,
            "time_col": args.time_col if mtype == "cox" else None,
            "time_label": t_label,
            "n_samples": n_samples,
            "n_events": n_events,
            "events_per_variable": epv,
            "discrimination_metric": disc_metric,
            "apparent_discrimination": round(disc_val, 4),
            "total_points_max": round(tp_max, 2),
            "lp_at_zero_points": round(lp_zero, 5),
            "points_per_unit_lp": round(pts_per_lp, 5),
            "landmarks": landmarks,
            "baseline_survival_at_landmarks": {f"{k:g}": round(v, 5) for k, v in s0_map.items()},
            "predictors": predictor_items,
            "model_terms": model_terms,
            "scorecard_table": os.path.abspath(scorecard_path),
            "lookup_table": os.path.abspath(lookup_path),
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
