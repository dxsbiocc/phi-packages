#!/usr/bin/env python3
"""mb_drift_correct: QC-RLSC (Quality Control-based Robust LOESS Signal Correction)."""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    BLANK_PATTERN,
    emit,
    fail,
    read_matrix,
    read_tsv,
    setup_matplotlib,
    write_artifact,
)

SKILL = "lcms-workflow"


def _resolve_col(cols: list[str], candidates: list[str], label: str, *, required: bool = True) -> str | None:
    lower = {str(c).strip().lower(): str(c) for c in cols}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    if required:
        fail(f"sample_sheet missing '{label}' column (tried {candidates}); found {list(cols)}")
    return None


def _loess_1d_fit_predict(x_fit: np.ndarray, y_fit: np.ndarray, x_eval: np.ndarray, span: float) -> np.ndarray:
    """Pure-NumPy tricube locally weighted linear regression (LOESS) with 2 robust bisquare IRLS steps."""
    n = len(x_fit)
    k = max(3, min(n, int(math_ceil(span * n))))
    weights_rob = np.ones(n, dtype=float)
    y_hat_qc = np.empty(n, dtype=float)

    def _predict_at(xq: np.ndarray, rob_w: np.ndarray) -> np.ndarray:
        out = np.empty(len(xq), dtype=float)
        for j, x0 in enumerate(xq):
            dists = np.abs(x_fit - x0)
            # Distance to the k-th nearest QC point
            d_k = float(np.partition(dists, k - 1)[k - 1])
            scale = max(d_k * 1.0001, 1e-9)
            u = np.clip(dists / scale, 0.0, 1.0)
            w = ((1.0 - u ** 3) ** 3) * rob_w
            if np.sum(w > 0) < 2:
                out[j] = float(np.median(y_fit))
                continue
            dx = x_fit - x0
            sw = float(np.sum(w))
            sx = float(np.sum(w * dx))
            sy = float(np.sum(w * y_fit))
            sxx = float(np.sum(w * dx * dx))
            sxy = float(np.sum(w * dx * y_fit))
            denom = sw * sxx - sx * sx
            if abs(denom) > 1e-12:
                intercept = (sy * sxx - sx * sxy) / denom
            else:
                intercept = sy / max(sw, 1e-12)
            out[j] = intercept
        return out

    for _ in range(2):
        y_hat_qc = _predict_at(x_fit, weights_rob)
        resid = y_fit - y_hat_qc
        mad = float(np.median(np.abs(resid)))
        if mad <= 1e-12:
            break
        u_r = np.abs(resid) / (6.0 * mad)
        weights_rob = np.where(u_r < 1.0, (1.0 - u_r ** 2) ** 2, 0.0)

    y_hat_qc = _predict_at(x_fit, weights_rob)
    y_eval = _predict_at(x_eval, weights_rob)
    lo_clamp = float(np.min(y_hat_qc))
    hi_clamp = float(np.max(y_hat_qc))
    return np.clip(y_eval, lo_clamp, hi_clamp)


def math_ceil(val: float) -> int:
    return int(np.ceil(val))


def _pca_2d_svd(mat_feat_by_sample: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Center/scale rows and return 2D sample scores + variance ratios via pure NumPy SVD."""
    x = mat_feat_by_sample.copy()
    row_med = np.nanmedian(x, axis=1)
    row_med = np.where(np.isfinite(row_med), row_med, 0.0)
    nan_idx = np.where(~np.isfinite(x))
    if nan_idx[0].size > 0:
        x[nan_idx] = row_med[nan_idx[0]]
    x = np.log2(np.clip(x, 0.0, None) + 1.0)
    mu = np.mean(x, axis=1, keepdims=True)
    sd = np.std(x, axis=1, ddof=1, keepdims=True)
    keep = (sd[:, 0] > 1e-12)
    if int(keep.sum()) < 2:
        return np.zeros((x.shape[1], 2), dtype=float), np.array([0.0, 0.0])
    xz = ((x[keep] - mu[keep]) / sd[keep]).T  # samples x features
    _, s, _ = np.linalg.svd(xz, full_matrices=False)
    u, _, _ = np.linalg.svd(xz, full_matrices=False)
    scores = u[:, :2] * s[:2]
    if scores.shape[1] < 2:
        scores = np.column_stack([scores, np.zeros(scores.shape[0])])
    tot_var = float(np.sum(s ** 2))
    vr = (s[:2] ** 2) / tot_var if tot_var > 1e-12 else np.array([0.0, 0.0])
    if len(vr) < 2:
        vr = np.array([float(vr[0]), 0.0])
    return scores, vr


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--sample_sheet", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--qc_label", default="QC")
    ap.add_argument("--span", type=float, default=0.75)
    ap.add_argument("--rsd_cutoff", type=float, default=30.0)
    args = ap.parse_args()

    if not (0.25 <= float(args.span) <= 1.0):
        fail(f"span must be in [0.25, 1.0] (got {args.span})")
    if float(args.rsd_cutoff) <= 0.0:
        fail("rsd_cutoff must be > 0")

    meta, X = read_matrix(args.input)
    sheet = read_tsv(args.sample_sheet)
    s_col = _resolve_col(list(sheet.columns), ["sample", "sample_id", "name", "filename"], "sample")
    ord_col = _resolve_col(list(sheet.columns), ["injection_order", "run_order", "order", "inj_order", "seq"], "injection_order")
    type_col = _resolve_col(list(sheet.columns), ["sample_type", "type", "group", "class", "condition"], "sample_type", required=False)

    sheet = sheet.copy()
    sheet["_sample"] = sheet[s_col].astype(str).str.strip()
    sheet["_order"] = pd.to_numeric(sheet[ord_col], errors="coerce")
    sheet = sheet.dropna(subset=["_sample", "_order"]).drop_duplicates(subset=["_sample"], keep="first")

    shared = [c for c in X.columns if str(c).strip() in set(sheet["_sample"])]
    if len(shared) < 6:
        fail(f"Found only {len(shared)} sample(s) matching sample_sheet in intensity matrix (need >= 6).")

    order_map = dict(zip(sheet["_sample"], sheet["_order"]))
    type_map: dict[str, str] = {}
    qc_tok = str(args.qc_label).strip().lower()
    if type_col is not None:
        for s_name, t_raw in zip(sheet["_sample"], sheet[type_col].astype(str).str.strip()):
            type_map[s_name] = t_raw
    else:
        for s_name in shared:
            if BLANK_PATTERN.search(s_name):
                type_map[s_name] = "blank"
            elif re.search(r"(qc|pool)", s_name, re.I):
                type_map[s_name] = qc_tok
            else:
                type_map[s_name] = "sample"

    qc_cols = [
        s for s in shared
        if not BLANK_PATTERN.search(s)
        and (type_map.get(s, "").strip().lower() == qc_tok or qc_tok in type_map.get(s, "").strip().lower())
    ]
    if not qc_cols:
        qc_cols = [s for s in shared if re.search(r"(qc|pool)", s, re.I) and not BLANK_PATTERN.search(s)]
    if len(qc_cols) < 4:
        fail(f"QC-RLSC drift correction requires >= 4 pooled QC samples across injection order (found {len(qc_cols)}).")

    bio_cols = [s for s in shared if s not in qc_cols and not BLANK_PATTERN.search(s)]
    if not bio_cols:
        fail("No biological samples found alongside QC samples.")

    # Align matrix columns to shared samples
    X_sub = X[shared].copy()
    x_orders = np.array([float(order_map[s]) for s in shared], dtype=float)
    qc_mask = np.array([s in set(qc_cols) for s in shared], dtype=bool)
    qc_orders = x_orders[qc_mask]
    sort_qc = np.argsort(qc_orders)
    qc_orders_sorted = qc_orders[sort_qc]

    raw_mat = X_sub.to_numpy(dtype=float, copy=True)
    corr_mat = raw_mat.copy()

    rsd_before_list: list[float] = []
    rsd_after_list: list[float] = []
    drift_span_pct_list: list[float] = []
    corrected_flags: list[bool] = []
    warnings: list[str] = []

    for r_i in range(raw_mat.shape[0]):
        row_vals = raw_mat[r_i]
        qc_vals_unsorted = row_vals[qc_mask]
        qc_vals = qc_vals_unsorted[sort_qc]

        fin_qc = np.isfinite(qc_vals) & (qc_vals > 0)
        if int(fin_qc.sum()) < 4:
            rsd_before = float(
                100.0 * np.std(qc_vals[fin_qc], ddof=1) / max(1e-12, np.mean(qc_vals[fin_qc]))
            ) if int(fin_qc.sum()) >= 2 else float("nan")
            rsd_before_list.append(rsd_before)
            rsd_after_list.append(rsd_before)
            drift_span_pct_list.append(0.0)
            corrected_flags.append(False)
            continue

        x_q = qc_orders_sorted[fin_qc]
        y_q = qc_vals[fin_qc]
        qc_med = float(np.median(y_q))
        qc_mean_b = float(np.mean(y_q))
        rsd_b = 100.0 * float(np.std(y_q, ddof=1)) / max(1e-12, qc_mean_b)

        y_fit_all = _loess_1d_fit_predict(x_q, y_q, x_orders, float(args.span))
        floor_y = max(1e-12, 0.10 * qc_med)
        y_fit_safe = np.maximum(y_fit_all, floor_y)

        fin_all = np.isfinite(row_vals) & (row_vals > 0)
        corr_row = row_vals.copy()
        corr_row[fin_all] = row_vals[fin_all] * (qc_med / y_fit_safe[fin_all])
        corr_mat[r_i] = corr_row

        qc_corr = corr_row[qc_mask]
        fin_qc_a = np.isfinite(qc_corr) & (qc_corr > 0)
        rsd_a = 100.0 * float(np.std(qc_corr[fin_qc_a], ddof=1)) / max(1e-12, float(np.mean(qc_corr[fin_qc_a])))
        drift_span = 100.0 * (float(np.max(y_fit_safe)) - float(np.min(y_fit_safe))) / max(1e-12, qc_med)

        rsd_before_list.append(rsd_b)
        rsd_after_list.append(rsd_a)
        drift_span_pct_list.append(drift_span)
        corrected_flags.append(True)

    rsd_b_arr = np.asarray(rsd_before_list, dtype=float)
    rsd_a_arr = np.asarray(rsd_after_list, dtype=float)
    fin_pair = np.isfinite(rsd_b_arr) & np.isfinite(rsd_a_arr)
    med_b = float(np.median(rsd_b_arr[fin_pair])) if fin_pair.any() else 0.0
    med_a = float(np.median(rsd_a_arr[fin_pair])) if fin_pair.any() else 0.0

    cut = float(args.rsd_cutoff)
    n_pass_b = int(np.sum(np.isfinite(rsd_b_arr) & (rsd_b_arr <= cut)))
    n_pass_a = int(np.sum(np.isfinite(rsd_a_arr) & (rsd_a_arr <= cut)))
    n_rescued = int(np.sum(np.isfinite(rsd_b_arr) & (rsd_b_arr > cut) & np.isfinite(rsd_a_arr) & (rsd_a_arr <= cut)))

    if int(sum(corrected_flags)) < raw_mat.shape[0]:
        n_skip = raw_mat.shape[0] - int(sum(corrected_flags))
        warnings.append(f"{n_skip} feature(s) had <4 non-zero finite QC measurements and were kept uncorrected.")
    if med_a > med_b and fin_pair.any():
        warnings.append("Median QC RSD did not decrease after LOESS drift correction; inspect QC injection stability or increase --span.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    corr_df = pd.DataFrame(np.round(corr_mat, 6), index=X_sub.index, columns=shared)
    out_mat_df = pd.concat([meta.loc[corr_df.index], corr_df], axis=1) if not meta.empty else corr_df
    out_mat_df.index.name = "feature"
    mat_path = out_dir / "drift_corrected_matrix.tsv"
    out_mat_df.to_csv(mat_path, sep="\t")

    metrics_df = pd.DataFrame(
        {
            "feature": X_sub.index.astype(str),
            "corrected": corrected_flags,
            "qc_rsd_before_pct": np.round(rsd_b_arr, 3),
            "qc_rsd_after_pct": np.round(rsd_a_arr, 3),
            "drift_span_pct": np.round(np.asarray(drift_span_pct_list, dtype=float), 3),
            "rescued_below_cutoff": np.isfinite(rsd_b_arr) & (rsd_b_arr > cut) & np.isfinite(rsd_a_arr) & (rsd_a_arr <= cut),
        }
    )
    met_path = out_dir / "feature_drift_qc_metrics.tsv"
    metrics_df.to_csv(met_path, sep="\t", index=False)

    # 2-panel figure: QC RSD distribution before vs after + PCA clustering tightening
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.6), dpi=160)

    if fin_pair.any():
        ax1.scatter(rsd_b_arr[fin_pair], rsd_a_arr[fin_pair], s=22, alpha=0.75, color="#0f766e")
        lim = max(50.0, float(np.percentile(np.concatenate([rsd_b_arr[fin_pair], rsd_a_arr[fin_pair]]), 98)) * 1.08)
        ax1.plot([0, lim], [0, lim], "--", color="#64748b", lw=1.1, label="y = x (no change)")
        ax1.axhline(cut, color="#dc2626", ls=":", lw=1.1, label=f"Cutoff = {cut:g}%")
        ax1.axvline(cut, color="#dc2626", ls=":", lw=1.1)
        ax1.set_xlim(0, lim)
        ax1.set_ylim(0, lim)
    ax1.set_xlabel("QC RSD Before Correction (%)")
    ax1.set_ylabel("QC RSD After QC-RLSC (%)")
    ax1.set_title(f"QC Precision Rescue ({n_rescued} features rescued <= {cut:g}%)")
    ax1.legend(frameon=False, fontsize=8.2)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    sc_after, vr_a = _pca_2d_svd(corr_mat)
    sc_before, _ = _pca_2d_svd(raw_mat)
    ax2.scatter(sc_before[qc_mask, 0], sc_before[qc_mask, 1], c="#f59e0b", marker="x", s=45, label="QC (Raw)", alpha=0.85)
    ax2.scatter(sc_after[~qc_mask, 0], sc_after[~qc_mask, 1], c="#94a3b8", marker="o", s=30, label="Bio (Corrected)", alpha=0.65)
    ax2.scatter(sc_after[qc_mask, 0], sc_after[qc_mask, 1], c="#2563eb", marker="^", s=55, label="QC (Corrected)", alpha=0.95)
    ax2.set_xlabel(f"PC1 ({vr_a[0] * 100:.1f}% var)")
    ax2.set_ylabel(f"PC2 ({vr_a[1] * 100:.1f}% var)")
    ax2.set_title("PCA QC Tightening (Raw vs QC-RLSC)")
    ax2.legend(frameon=False, fontsize=8.0)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = out_dir / "qc_drift_pca_rsd_comparison.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {"qc_label": args.qc_label, "span": float(args.span), "rsd_cutoff": cut}
    inputs_all = [args.input, args.sample_sheet]
    artifacts = [
        write_artifact(mat_path, "mb_drift_correct", "QC-RLSC drift-corrected intensity matrix", skill=SKILL, inputs=inputs_all, parameters=params_meta, shape=out_mat_df.shape),
        write_artifact(met_path, "mb_drift_correct", "Per-feature QC RSD drift metrics", skill=SKILL, inputs=inputs_all, parameters=params_meta, shape=metrics_df.shape),
        write_artifact(fig_path, "mb_drift_correct", "QC RSD rescue and PCA tightening", skill=SKILL, inputs=inputs_all, parameters=params_meta),
    ]

    emit(
        {
            "status": "ok",
            "tool": "mb_drift_correct",
            "input": str(args.input),
            "sample_sheet": str(args.sample_sheet),
            "n_features": int(raw_mat.shape[0]),
            "n_samples_total": int(len(shared)),
            "n_qc_samples": int(len(qc_cols)),
            "n_biological_samples": int(len(bio_cols)),
            "span": float(args.span),
            "rsd_cutoff": cut,
            "median_qc_rsd_before": round(med_b, 3),
            "median_qc_rsd_after": round(med_a, 3),
            "n_features_below_cutoff_before": n_pass_b,
            "n_features_below_cutoff_after": n_pass_a,
            "n_features_rescued": n_rescued,
            "outputs": {
                "drift_corrected_matrix": str(mat_path),
                "feature_qc_metrics": str(met_path),
                "drift_qc_plot": str(fig_path),
            },
            "scientificContext": {
                str(mat_path): {
                    "valueSemantics": "normalized",
                    "sampleIds": shared,
                }
            },
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
