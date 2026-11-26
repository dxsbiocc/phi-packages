#!/usr/bin/env python3
"""clin_calibration_dca: TRIPOD calibration diagnostics (Brier, scaled Brier IPA,
E/O ratio, calibration-in-the-large, calibration slope, Hosmer-Lemeshow test) and
Vickers Decision Curve Analysis (DCA) net benefit & interventions avoided.
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
from scipy import optimize, stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pred_common import (
    emit,
    encode_covariates,
    fail,
    fit_logistic_weighted,
    parse_binary_series,
    read_table,
    wilson_ci,
    write_artifact,
)


def fit_citl_offset(y: np.ndarray, logit_p: np.ndarray) -> tuple[float, float, float]:
    """Fit logistic calibration-in-the-large: logit(P(Y=1)) = alpha + 1 * logit_p."""
    def nll(alpha: float) -> float:
        eta = np.clip(alpha + logit_p, -30.0, 30.0)
        return float(np.sum(np.logaddexp(0.0, eta) - y * eta))

    res = optimize.minimize_scalar(nll, bounds=(-10.0, 10.0), method="bounded")
    alpha_hat = float(res.x)
    eta_hat = np.clip(alpha_hat + logit_p, -30.0, 30.0)
    p_tilde = 1.0 / (1.0 + np.exp(-eta_hat))
    info = float(np.sum(p_tilde * (1.0 - p_tilde)))
    se = 1.0 / math.sqrt(max(info, 1e-12))
    return alpha_hat, alpha_hat - 1.959964 * se, alpha_hat + 1.959964 * se


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--outcome_col", required=True)
    parser.add_argument("--prob_col", default=None)
    parser.add_argument("--features", action="extend", nargs="+", default=[])
    parser.add_argument("--positive_label", default=None)
    parser.add_argument("--n_bins", type=int, default=10)
    parser.add_argument("--n_folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    feats = [f.strip() for item in args.features for f in str(item).split(",") if f.strip()]
    if bool(args.prob_col) == bool(feats):
        fail("specify exactly one of --prob_col (precomputed probabilities) or --features (OOF CV logistic model)")

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.clinical)
    if args.outcome_col not in df.columns:
        fail(f"outcome_col '{args.outcome_col}' not found in columns: {list(df.columns)}")

    y_ser, _ = parse_binary_series(df[args.outcome_col], args.outcome_col, args.positive_label)
    warnings: list[str] = []

    if args.prob_col:
        if args.prob_col not in df.columns:
            fail(f"prob_col '{args.prob_col}' not found in columns")
        p_ser = pd.to_numeric(df[args.prob_col], errors="coerce")
        valid = y_ser.notna() & p_ser.notna()
        y = y_ser.loc[valid].to_numpy(dtype=int, copy=True)
        p_hat = p_ser.loc[valid].to_numpy(dtype=float, copy=True)
        if np.any((p_hat < 0.0) | (p_hat > 1.0)):
            fail("prob_col must contain predicted probabilities strictly within [0, 1]")
        eval_mode = "precomputed_prob"
    else:
        for f in feats:
            if f not in df.columns:
                fail(f"feature '{f}' not found in columns: {list(df.columns)}")
        enc_df, _, _ = encode_covariates(df, feats)
        valid = y_ser.notna() & enc_df.notna().all(axis=1)
        y = y_ser.loc[valid].to_numpy(dtype=int, copy=True)
        X = enc_df.loc[valid].to_numpy(dtype=float, copy=True)
        n_pos = int(np.sum(y == 1))
        n_neg = int(np.sum(y == 0))
        k_folds = min(args.n_folds, n_pos, n_neg)
        if k_folds < 3:
            fail(f"need >= 3 positives and >= 3 negatives for stratified CV (got pos={n_pos}, neg={n_neg})")
        skf = StratifiedKFold(n_splits=k_folds, shuffle=True, random_state=args.seed)
        p_hat = np.zeros(len(y), dtype=float)
        for tr_idx, te_idx in skf.split(X, y):
            mu = X[tr_idx].mean(axis=0)
            sd = X[tr_idx].std(axis=0, ddof=0)
            sd = np.where(sd > 1e-12, sd, 1.0)
            X_tr = (X[tr_idx] - mu) / sd
            X_te = (X[te_idx] - mu) / sd
            clf = LogisticRegression(C=1.0, max_iter=1000)
            clf.fit(X_tr, y[tr_idx])
            p_hat[te_idx] = clf.predict_proba(X_te)[:, 1]
        eval_mode = "oof_cv_logistic"

    n_dropped = int((~valid).sum())
    if n_dropped > 0:
        warnings.append(f"excluded {n_dropped} row(s) with missing outcome/predictors")

    n = len(y)
    n_pos = int(np.sum(y == 1))
    if n < 15 or n_pos < 3 or (n - n_pos) < 3:
        fail(f"insufficient samples for calibration/DCA (n={n}, pos={n_pos})")

    prevalence = n_pos / float(n)
    p_clip = np.clip(p_hat, 1e-6, 1.0 - 1e-6)
    logit_p = np.log(p_clip / (1.0 - p_clip))

    roc_auc = float(roc_auc_score(y, p_clip))
    brier = float(np.mean((y - p_hat) ** 2))
    brier_null = prevalence * (1.0 - prevalence)
    ipa = 1.0 - (brier / brier_null) if brier_null > 1e-12 else 0.0
    eo_ratio = float(np.sum(p_hat) / max(1, n_pos))

    citl, citl_lo, citl_hi = fit_citl_offset(y, logit_p)
    cal_fit = fit_logistic_weighted(logit_p[:, None], y, l2_penalty=0.0)
    cal_intercept = float(cal_fit["beta"][0])
    cal_slope = float(cal_fit["beta"][1])
    slope_se = float(cal_fit["se_model"][1])
    slope_lo = cal_slope - 1.959964 * slope_se
    slope_hi = cal_slope + 1.959964 * slope_se

    # Calibration bins & Hosmer-Lemeshow goodness-of-fit test
    n_bins_req = max(4, min(args.n_bins, n // 5))
    try:
        bins = pd.qcut(p_clip, q=n_bins_req, labels=False, duplicates="drop")
    except ValueError:
        bins = pd.cut(p_clip, bins=min(n_bins_req, 5), labels=False)

    bin_rows: list[dict[str, Any]] = []
    hl_stat = 0.0
    for b_id in sorted(np.unique(bins)):
        m = bins == b_id
        n_b = int(np.sum(m))
        obs_cnt = int(np.sum(y[m]))
        exp_cnt = float(np.sum(p_clip[m]))
        p_mean = float(np.mean(p_clip[m]))
        obs_prop, obs_lo, obs_hi = wilson_ci(obs_cnt, n_b)
        denom_hl = n_b * p_mean * (1.0 - p_mean)
        if denom_hl > 1e-9:
            hl_stat += ((obs_cnt - exp_cnt) ** 2) / denom_hl
        bin_rows.append(
            {
                "bin": int(b_id) + 1,
                "n": n_b,
                "observed_events": obs_cnt,
                "expected_events": round(exp_cnt, 3),
                "mean_predicted_prob": round(p_mean, 5),
                "observed_proportion": round(obs_prop, 5),
                "observed_ci_low": round(obs_lo, 5),
                "observed_ci_high": round(obs_hi, 5),
            }
        )

    n_bins_used = len(bin_rows)
    hl_df = max(1, n_bins_used - 2)
    hl_p = float(stats.chi2.sf(hl_stat, df=hl_df))

    # Decision Curve Analysis (Vickers & Elkin 2006)
    thresholds = np.round(np.arange(0.01, 1.00, 0.01), 2)
    dca_rows: list[dict[str, Any]] = []
    useful_pts: list[float] = []
    for pt in thresholds:
        odds = pt / (1.0 - pt)
        pred_pos = p_hat >= pt
        tp = float(np.sum(pred_pos & (y == 1)))
        fp = float(np.sum(pred_pos & (y == 0)))
        nb_model = (tp / n) - (fp / n) * odds
        nb_all = prevalence - (1.0 - prevalence) * odds
        nb_none = 0.0
        interventions_avoided = ((nb_model - nb_all) / odds) * 100.0 if odds > 1e-9 else 0.0
        beats_both = nb_model > max(nb_all, nb_none) + 1e-6
        if beats_both:
            useful_pts.append(float(pt))
        dca_rows.append(
            {
                "threshold_prob": float(pt),
                "net_benefit_model": round(nb_model, 5),
                "net_benefit_treat_all": round(nb_all, 5),
                "net_benefit_treat_none": 0.0,
                "net_interventions_avoided_per_100": round(interventions_avoided, 3),
                "beats_all_and_none": bool(beats_both),
            }
        )

    # Longest contiguous useful threshold range
    useful_range = None
    if useful_pts:
        # group consecutive 0.01 steps
        best_seg = [useful_pts[0], useful_pts[0]]
        cur_start = useful_pts[0]
        prev = useful_pts[0]
        for pt in useful_pts[1:]:
            if pt - prev <= 0.015:
                prev = pt
            else:
                if (prev - cur_start) > (best_seg[1] - best_seg[0]):
                    best_seg = [cur_start, prev]
                cur_start = pt
                prev = pt
        if (prev - cur_start) > (best_seg[1] - best_seg[0]):
            best_seg = [cur_start, prev]
        useful_range = [round(best_seg[0], 2), round(best_seg[1], 2)]

    max_nb = max(r["net_benefit_model"] for r in dca_rows)
    max_avoid = max(r["net_interventions_avoided_per_100"] for r in dca_rows)

    cal_df = pd.DataFrame(bin_rows)
    dca_df = pd.DataFrame(dca_rows)
    cal_path = os.path.join(args.output_dir, "calibration_bins.tsv")
    dca_path = os.path.join(args.output_dir, "dca_curve.tsv")
    cal_df.to_csv(cal_path, sep="\t", index=False)
    dca_df.to_csv(dca_path, sep="\t", index=False)

    # 2-panel Publication Figure: Calibration + DCA
    fig, (ax_cal, ax_dca) = plt.subplots(1, 2, figsize=(11.8, 5.0))

    ax_cal.plot([0, 1], [0, 1], color="#64748B", linestyle="--", linewidth=1.2, label="Ideal (45 deg)")
    p_grid = np.linspace(0.01, 0.99, 150)
    logit_grid = np.log(p_grid / (1.0 - p_grid))
    smooth_cal = 1.0 / (1.0 + np.exp(-(cal_intercept + cal_slope * logit_grid)))
    ax_cal.plot(p_grid, smooth_cal, color="#2563EB", linewidth=2.0, label=f"Logistic cal (slope={cal_slope:.2f}, CITL={citl:.2f})")
    ax_cal.errorbar(
        cal_df["mean_predicted_prob"],
        cal_df["observed_proportion"],
        yerr=[
            np.maximum(0.0, cal_df["observed_proportion"] - cal_df["observed_ci_low"]),
            np.maximum(0.0, cal_df["observed_ci_high"] - cal_df["observed_proportion"]),
        ],
        fmt="o",
        color="#E11D48",
        ecolor="#E11D48",
        elinewidth=1.2,
        capsize=3,
        markersize=5,
        label=f"Decile bins (HL p={hl_p:.3f})",
    )
    # Spike histogram along bottom
    hist_vals, hist_edges = np.histogram(p_hat, bins=25, range=(0.0, 1.0))
    hist_scaled = 0.12 * (hist_vals / max(1, int(hist_vals.max())))
    ax_cal.bar(hist_edges[:-1], hist_scaled, width=0.036, align="edge", color="#94A3B8", alpha=0.45)
    ax_cal.set_xlim(-0.02, 1.02)
    ax_cal.set_ylim(-0.02, 1.02)
    ax_cal.set_xlabel("Predicted Probability")
    ax_cal.set_ylabel("Observed Proportion")
    ax_cal.set_title(f"TRIPOD Calibration ({eval_mode}, Brier={brier:.3f}, IPA={ipa:.2f})", fontsize=9.8)
    ax_cal.legend(frameon=False, loc="upper left", fontsize=8)
    ax_cal.spines["top"].set_visible(False)
    ax_cal.spines["right"].set_visible(False)

    pts = dca_df["threshold_prob"].to_numpy()
    nb_m = dca_df["net_benefit_model"].to_numpy()
    nb_a = dca_df["net_benefit_treat_all"].to_numpy()
    ax_dca.plot(pts, nb_m, color="#2563EB", linewidth=2.2, label="Prediction Model")
    ax_dca.plot(pts, nb_a, color="#E11D48", linestyle="--", linewidth=1.5, label="Treat All")
    ax_dca.axhline(0.0, color="#475569", linestyle=":", linewidth=1.3, label="Treat None")
    y_min_dca = max(-0.25, -0.3 * prevalence)
    y_max_dca = max(prevalence * 1.15, float(np.max(nb_m)) * 1.15, 0.05)
    ax_dca.set_xlim(0.0, 1.0)
    ax_dca.set_ylim(y_min_dca, y_max_dca)
    ax_dca.set_xlabel("Threshold Probability (p_t)")
    ax_dca.set_ylabel("Net Benefit")
    rng_txt = f"useful {useful_range[0]:.2f}-{useful_range[1]:.2f}" if useful_range else "no positive range"
    ax_dca.set_title(f"Decision Curve Analysis ({rng_txt})", fontsize=9.8)
    ax_dca.legend(frameon=False, loc="upper right", fontsize=8)
    ax_dca.spines["top"].set_visible(False)
    ax_dca.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "calibration_dca.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {"outcome_col": args.outcome_col, "evaluation_mode": eval_mode, "n_bins": n_bins_used}
    artifacts = [
        write_artifact(fig_path, "clin_calibration_dca", "TRIPOD Calibration & Decision Curve Analysis", inputs=(args.clinical,), parameters=params_meta),
        write_artifact(cal_path, "clin_calibration_dca", "Calibration quantile bins table", inputs=(args.clinical,), parameters=params_meta, shape=cal_df.shape),
        write_artifact(dca_path, "clin_calibration_dca", "Decision Curve Analysis net benefit table", inputs=(args.clinical,), parameters=params_meta, shape=dca_df.shape),
    ]

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "outcome_col": args.outcome_col,
            "evaluation_mode": eval_mode,
            "n_samples": n,
            "n_positive": n_pos,
            "prevalence": round(prevalence, 4),
            "roc_auc": round(roc_auc, 5),
            "calibration": {
                "brier_score": round(brier, 5),
                "brier_null": round(brier_null, 5),
                "scaled_brier_ipa": round(ipa, 5),
                "expected_observed_ratio": round(eo_ratio, 4),
                "citl": round(citl, 4),
                "citl_ci_low": round(citl_lo, 4),
                "citl_ci_high": round(citl_hi, 4),
                "calibration_intercept": round(cal_intercept, 4),
                "calibration_slope": round(cal_slope, 4),
                "slope_ci_low": round(slope_lo, 4),
                "slope_ci_high": round(slope_hi, 4),
                "hosmer_lemeshow_chisq": round(hl_stat, 4),
                "hosmer_lemeshow_df": hl_df,
                "hosmer_lemeshow_pvalue": round(hl_p, 5),
                "n_bins_used": n_bins_used,
            },
            "dca": {
                "useful_threshold_range": useful_range,
                "max_net_benefit": round(max_nb, 5),
                "max_interventions_avoided_per_100": round(max_avoid, 3),
            },
            "calibration_table": os.path.abspath(cal_path),
            "dca_table": os.path.abspath(dca_path),
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
