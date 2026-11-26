#!/usr/bin/env python3
"""clin_diagnostic_roc: Empirical ROC/PR curves, exact Sun & Xu / DeLong (1988)
variance and paired DeLong Z-tests across biomarkers, inverse-direction detection,
and 2x2 diagnostic performance metrics with Wilson 95% CIs.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_curve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pred_common import (
    bh_fdr,
    emit,
    fail,
    parse_binary_series,
    read_table,
    wilson_ci,
    write_artifact,
)

PALETTE = ["#2563EB", "#E11D48", "#10B981", "#F59E0B", "#8B5CF6", "#06B6D4"]


def _delong_components(y: np.ndarray, scores: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Exact Sun & Xu (2014) midrank computation of DeLong structural components V10, V01."""
    pos_scores = scores[y == 1]
    neg_scores = scores[y == 0]
    m = len(pos_scores)
    n = len(neg_scores)
    tx = stats.rankdata(pos_scores, method="average")
    tz = stats.rankdata(neg_scores, method="average")
    txz = stats.rankdata(np.concatenate([pos_scores, neg_scores]), method="average")
    v10 = (txz[:m] - tx) / float(n)
    v01 = 1.0 - (txz[m:] - tz) / float(m)
    auc = float(np.mean(v10))
    return auc, v10, v01


def delong_covariance(y: np.ndarray, score_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (auc_vec, cov_mat) for K predictors evaluated on the same y."""
    m = int(np.sum(y == 1))
    n = int(np.sum(y == 0))
    k = score_matrix.shape[1]
    aucs = np.zeros(k, dtype=float)
    v10_mat = np.zeros((m, k), dtype=float)
    v01_mat = np.zeros((n, k), dtype=float)
    for j in range(k):
        auc_j, v10_j, v01_j = _delong_components(y, score_matrix[:, j])
        aucs[j] = auc_j
        v10_mat[:, j] = v10_j
        v01_mat[:, j] = v01_j

    s10 = np.cov(v10_mat, rowvar=False, ddof=1) if k > 1 else np.array([[float(np.var(v10_mat[:, 0], ddof=1))]])
    s01 = np.cov(v01_mat, rowvar=False, ddof=1) if k > 1 else np.array([[float(np.var(v01_mat[:, 0], ddof=1))]])
    s10 = np.atleast_2d(np.asarray(s10, dtype=float))
    s01 = np.atleast_2d(np.asarray(s01, dtype=float))
    cov = s10 / float(m) + s01 / float(n)
    return aucs, cov


def bootstrap_pr_auc(y: np.ndarray, scores: np.ndarray, n_boot: int, seed: int) -> tuple[float, float, float]:
    ap = float(average_precision_score(y, scores))
    rng = np.random.default_rng(seed)
    pos_idx = np.flatnonzero(y == 1)
    neg_idx = np.flatnonzero(y == 0)
    boots = []
    for _ in range(n_boot):
        b_pos = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        b_neg = rng.choice(neg_idx, size=len(neg_idx), replace=True)
        idx = np.concatenate([b_pos, b_neg])
        boots.append(float(average_precision_score(y[idx], scores[idx])))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return ap, float(lo), float(hi)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--outcome_col", required=True)
    parser.add_argument("--predictors", action="extend", nargs="+", default=[], required=True)
    parser.add_argument("--positive_label", default=None)
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument("--n_bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.clinical)
    if args.outcome_col not in df.columns:
        fail(f"outcome_col '{args.outcome_col}' not found in columns: {list(df.columns)}")

    pred_cols = []
    for item in args.predictors:
        for p in str(item).split(","):
            p = p.strip()
            if p and p not in pred_cols:
                pred_cols.append(p)
    if not pred_cols:
        fail("at least one predictor column is required")
    for p in pred_cols:
        if p not in df.columns:
            fail(f"predictor '{p}' not found in columns: {list(df.columns)}")

    y_ser, outcome_map = parse_binary_series(df[args.outcome_col], args.outcome_col, args.positive_label)
    num_df = pd.DataFrame({p: pd.to_numeric(df[p], errors="coerce") for p in pred_cols}, index=df.index)
    valid = y_ser.notna() & num_df.notna().all(axis=1)

    warnings: list[str] = []
    n_dropped = int((~valid).sum())
    if n_dropped > 0:
        warnings.append(f"complete-case diagnostic ROC excluded {n_dropped} row(s) with missing outcome or predictors")

    y = y_ser.loc[valid].to_numpy(dtype=int, copy=True)
    X_raw = num_df.loc[valid].to_numpy(dtype=float, copy=True)
    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    if n_pos < 3 or n_neg < 3:
        fail(f"need >= 3 positive and >= 3 negative samples (got pos={n_pos}, neg={n_neg})")

    # Direction detection per predictor
    X_aligned = X_raw.copy()
    raw_aucs = []
    directions = []
    for j, col in enumerate(pred_cols):
        r_auc, _, _ = _delong_components(y, X_raw[:, j])
        raw_aucs.append(r_auc)
        if r_auc < 0.5:
            X_aligned[:, j] = -X_raw[:, j]
            directions.append("negative")
            warnings.append(
                f"predictor '{col}' showed inverse association (raw AUC={r_auc:.3f} < 0.5); "
                "direction aligned for discrimination & paired DeLong comparison"
            )
        else:
            directions.append("positive")

    aucs_aligned, cov_mat = delong_covariance(y, X_aligned)
    prevalence = n_pos / float(len(y))

    if args.threshold is None:
        warnings.append(
            "operating threshold selected via Youden's J on the evaluation cohort (apparent resubstitution estimate); "
            "validate externally before clinical deployment"
        )

    pred_rows: list[dict[str, Any]] = []
    roc_curves: list[tuple[str, np.ndarray, np.ndarray, float]] = []
    pr_curves: list[tuple[str, np.ndarray, np.ndarray, float]] = []

    for j, col in enumerate(pred_cols):
        auc_j = float(aucs_aligned[j])
        se_j = float(np.sqrt(max(cov_mat[j, j], 1e-18)))
        ci_lo = float(np.clip(auc_j - 1.959964 * se_j, 0.0, 1.0))
        ci_hi = float(np.clip(auc_j + 1.959964 * se_j, 0.0, 1.0))
        z_vs_50 = (auc_j - 0.5) / se_j
        p_vs_50 = float(2.0 * stats.norm.sf(abs(z_vs_50)))

        ap, ap_lo, ap_hi = bootstrap_pr_auc(y, X_aligned[:, j], args.n_bootstrap, args.seed + j)

        fpr, tpr, roc_thr = roc_curve(y, X_aligned[:, j])
        prec, rec, _ = precision_recall_curve(y, X_aligned[:, j])
        roc_curves.append((col, fpr, tpr, auc_j))
        pr_curves.append((col, rec, prec, ap))

        raw_col = X_raw[:, j]
        if args.threshold is not None:
            thr_orig = float(args.threshold)
            pred_pos = raw_col <= thr_orig if directions[j] == "negative" else raw_col >= thr_orig
            thr_mode = "prespecified"
        else:
            youden = tpr - fpr
            best_idx = int(np.argmax(youden))
            thr_aligned = float(roc_thr[best_idx])
            thr_orig = -thr_aligned if directions[j] == "negative" else thr_aligned
            pred_pos = X_aligned[:, j] >= thr_aligned
            thr_mode = "youden_j"

        tp = int(np.sum(pred_pos & (y == 1)))
        fp = int(np.sum(pred_pos & (y == 0)))
        tn = int(np.sum((~pred_pos) & (y == 0)))
        fn = int(np.sum((~pred_pos) & (y == 1)))

        sens, sens_lo, sens_hi = wilson_ci(tp, tp + fn)
        spec, spec_lo, spec_hi = wilson_ci(tn, tn + fp)
        ppv, ppv_lo, ppv_hi = wilson_ci(tp, tp + fp)
        npv, npv_lo, npv_hi = wilson_ci(tn, tn + fn)
        acc, acc_lo, acc_hi = wilson_ci(tp + tn, len(y))
        lr_pos = (sens / (1.0 - spec)) if spec < 1.0 - 1e-12 else None
        lr_neg = ((1.0 - sens) / spec) if spec > 1e-12 else None

        pred_rows.append(
            {
                "predictor": col,
                "direction": directions[j],
                "direction_inverted": directions[j] == "negative",
                "raw_auc": round(float(raw_aucs[j]), 5),
                "auc": round(auc_j, 5),
                "auc_se_delong": round(se_j, 5),
                "auc_ci_low": round(ci_lo, 5),
                "auc_ci_high": round(ci_hi, 5),
                "auc_pvalue_vs_0_5": p_vs_50,
                "pr_auc": round(ap, 5),
                "pr_auc_ci_low": round(ap_lo, 5),
                "pr_auc_ci_high": round(ap_hi, 5),
                "threshold_mode": thr_mode,
                "threshold": round(thr_orig, 5),
                "youden_j": round(sens + spec - 1.0, 5),
                "tp": tp,
                "fp": fp,
                "tn": tn,
                "fn": fn,
                "sensitivity": round(sens, 4),
                "sensitivity_ci_low": round(sens_lo, 4),
                "sensitivity_ci_high": round(sens_hi, 4),
                "specificity": round(spec, 4),
                "specificity_ci_low": round(spec_lo, 4),
                "specificity_ci_high": round(spec_hi, 4),
                "ppv": round(ppv, 4) if np.isfinite(ppv) else None,
                "ppv_ci_low": round(ppv_lo, 4) if np.isfinite(ppv_lo) else None,
                "ppv_ci_high": round(ppv_hi, 4) if np.isfinite(ppv_hi) else None,
                "npv": round(npv, 4) if np.isfinite(npv) else None,
                "npv_ci_low": round(npv_lo, 4) if np.isfinite(npv_lo) else None,
                "npv_ci_high": round(npv_hi, 4) if np.isfinite(npv_hi) else None,
                "accuracy": round(acc, 4),
                "accuracy_ci_low": round(acc_lo, 4),
                "accuracy_ci_high": round(acc_hi, 4),
                "lr_plus": round(lr_pos, 4) if lr_pos is not None else None,
                "lr_minus": round(lr_neg, 4) if lr_neg is not None else None,
            }
        )

    # Paired DeLong comparisons vs reference (first predictor)
    delong_rows: list[dict[str, Any]] = []
    ref_name = pred_cols[0]
    if len(pred_cols) >= 2:
        pvals_dl = []
        for j in range(1, len(pred_cols)):
            delta = float(aucs_aligned[j] - aucs_aligned[0])
            var_d = float(cov_mat[j, j] + cov_mat[0, 0] - 2.0 * cov_mat[0, j])
            se_d = float(np.sqrt(max(var_d, 1e-18)))
            z_d = delta / se_d if var_d > 1e-18 else 0.0
            p_d = float(2.0 * stats.norm.sf(abs(z_d))) if var_d > 1e-18 else 1.0
            pvals_dl.append(p_d)
            delong_rows.append(
                {
                    "reference": ref_name,
                    "comparator": pred_cols[j],
                    "auc_reference": round(float(aucs_aligned[0]), 5),
                    "auc_comparator": round(float(aucs_aligned[j]), 5),
                    "delta_auc": round(delta, 5),
                    "se_delta": round(se_d, 5),
                    "ci_low": round(delta - 1.959964 * se_d, 5),
                    "ci_high": round(delta + 1.959964 * se_d, 5),
                    "z": round(z_d, 4),
                    "pvalue": p_d,
                }
            )
        fdr_dl = bh_fdr(np.array(pvals_dl, dtype=float))
        for idx_r, f_val in enumerate(fdr_dl):
            delong_rows[idx_r]["padj"] = float(f_val) if np.isfinite(f_val) else None

    metrics_df = pd.DataFrame(pred_rows)
    metrics_path = os.path.join(args.output_dir, "diagnostic_metrics.tsv")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)

    delong_path = None
    if delong_rows:
        delong_df = pd.DataFrame(delong_rows)
        delong_path = os.path.join(args.output_dir, "delong_comparisons.tsv")
        delong_df.to_csv(delong_path, sep="\t", index=False)

    # 2-panel ROC + PR figure
    fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(11.6, 5.0))
    ax_roc.plot([0, 1], [0, 1], color="#94A3B8", linestyle="--", linewidth=1.0, label="Chance (0.50)")
    for idx_c, (col, fpr, tpr, auc_v) in enumerate(roc_curves):
        c = PALETTE[idx_c % len(PALETTE)]
        r = pred_rows[idx_c]
        ax_roc.plot(fpr, tpr, color=c, linewidth=2.0, label=f"{col} (AUC={auc_v:.3f} [{r['auc_ci_low']:.3f}-{r['auc_ci_high']:.3f}])")
    ax_roc.set_xlim(-0.02, 1.02)
    ax_roc.set_ylim(-0.02, 1.02)
    ax_roc.set_xlabel("1 - Specificity (False Positive Rate)")
    ax_roc.set_ylabel("Sensitivity (True Positive Rate)")
    ax_roc.set_title(f"Diagnostic ROC Curves (n={len(y)}, pos={n_pos})", fontsize=10)
    ax_roc.legend(frameon=False, loc="lower right", fontsize=8)
    ax_roc.spines["top"].set_visible(False)
    ax_roc.spines["right"].set_visible(False)

    ax_pr.axhline(prevalence, color="#94A3B8", linestyle="--", linewidth=1.0, label=f"Baseline ({prevalence:.2f})")
    for idx_c, (col, rec, prec, ap_v) in enumerate(pr_curves):
        c = PALETTE[idx_c % len(PALETTE)]
        ax_pr.plot(rec, prec, color=c, linewidth=2.0, label=f"{col} (AP={ap_v:.3f})")
    ax_pr.set_xlim(-0.02, 1.02)
    ax_pr.set_ylim(-0.02, 1.02)
    ax_pr.set_xlabel("Recall (Sensitivity)")
    ax_pr.set_ylabel("Precision (PPV)")
    ax_pr.set_title("Precision-Recall Curves", fontsize=10)
    ax_pr.legend(frameon=False, loc="best", fontsize=8)
    ax_pr.spines["top"].set_visible(False)
    ax_pr.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "diagnostic_roc_pr.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {"outcome_col": args.outcome_col, "predictors": pred_cols, "threshold": args.threshold}
    artifacts = [
        write_artifact(fig_path, "clin_diagnostic_roc", "Diagnostic ROC & Precision-Recall curves", inputs=(args.clinical,), parameters=params_meta),
        write_artifact(metrics_path, "clin_diagnostic_roc", "Diagnostic ROC/PR & 2x2 Wilson metrics", inputs=(args.clinical,), parameters=params_meta, shape=metrics_df.shape),
    ]
    if delong_path:
        artifacts.append(
            write_artifact(delong_path, "clin_diagnostic_roc", "Paired DeLong ROC AUC comparisons", inputs=(args.clinical,), parameters=params_meta, shape=(len(delong_rows), 11))
        )

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "outcome_col": args.outcome_col,
            "outcome_mapping": outcome_map,
            "n_samples": int(len(y)),
            "n_positive": n_pos,
            "n_negative": n_neg,
            "prevalence": round(prevalence, 4),
            "reference_predictor": ref_name,
            "predictors": pred_rows,
            "delong_comparisons": delong_rows,
            "metrics_table": os.path.abspath(metrics_path),
            "delong_table": os.path.abspath(delong_path) if delong_path else None,
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
