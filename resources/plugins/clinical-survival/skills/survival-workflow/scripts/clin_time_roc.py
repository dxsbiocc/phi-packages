#!/usr/bin/env python3
"""clin_time_roc: IPCW time-dependent ROC & landmark AUC evaluation for precomputed
risk scores or leak-free stratified K-fold L2-penalized multivariable Cox signatures.
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    fit_cox_breslow,
    harrell_c_index,
    ipcw_time_dependent_roc,
    join_molecular_features,
    logrank_test,
    parse_event_series,
    read_table,
    write_artifact,
)

PALETTE = ["#2563EB", "#E11D48", "#10B981", "#F59E0B", "#8B5CF6"]


def stratified_event_folds(events: np.ndarray, n_folds: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    folds = np.zeros(len(events), dtype=int)
    for ev_val in (0, 1):
        idx = np.flatnonzero(events == ev_val)
        rng.shuffle(idx)
        for k, i in enumerate(idx):
            folds[i] = k % n_folds
    return folds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--time_col", required=True)
    parser.add_argument("--event_col", required=True)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--score_col", default=None)
    parser.add_argument("--feature_cols", action="extend", nargs="+", default=[])
    parser.add_argument("--matrix", default=None)
    parser.add_argument("--landmarks", action="extend", nargs="+", type=float, default=[], required=True)
    parser.add_argument("--n_folds", type=int, default=5)
    parser.add_argument("--l2_penalty", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    feats = [f.strip() for item in args.feature_cols for f in str(item).split(",") if f.strip()]
    if bool(args.score_col) == bool(feats):
        fail("specify exactly one of --score_col (precomputed score) or --feature_cols (L2-Cox CV signature)")

    os.makedirs(args.output_dir, exist_ok=True)
    clin = read_table(args.clinical)
    sample_col = args.sample_col or clin.columns[0]
    needed = [args.score_col] if args.score_col else feats
    if args.matrix:
        clin = join_molecular_features(clin, sample_col, args.matrix, needed)
    for col in needed:
        if col not in clin.columns:
            fail(f"column/feature '{col}' not found in clinical table or --matrix")

    t_s = pd.to_numeric(clin[args.time_col], errors="coerce")
    e_s, _ = parse_event_series(clin[args.event_col], args.event_col)
    num_df = pd.DataFrame({c: pd.to_numeric(clin[c], errors="coerce") for c in needed}, index=clin.index)

    valid = t_s.notna() & (t_s > 0) & e_s.notna() & num_df.notna().all(axis=1)
    sub_ids = clin.loc[valid, sample_col].astype(str).to_numpy()
    times = t_s.loc[valid].to_numpy(dtype=float)
    events = e_s.loc[valid].to_numpy(dtype=int)
    X_all = num_df.loc[valid].to_numpy(dtype=float)

    n_samples = len(times)
    n_events = int(events.sum())
    if n_samples < 10 or n_events < 4:
        fail(f"need >= 10 complete samples and >= 4 events for time-dependent ROC (got n={n_samples}, events={n_events})")

    warnings: list[str] = []
    feature_weights = None

    if args.score_col:
        eval_mode = "precomputed_score"
        eval_scores = X_all[:, 0]
    else:
        eval_mode = "oof_cv_l2_cox"
        k_folds = min(args.n_folds, n_events)
        if k_folds < 3:
            fail(f"need >= 3 events for stratified CV L2-Cox (got {n_events})")
        folds = stratified_event_folds(events, k_folds, args.seed)
        eval_scores = np.zeros(n_samples, dtype=float)

        for f_id in range(k_folds):
            tr = folds != f_id
            te = folds == f_id
            mu = X_all[tr].mean(axis=0)
            sd = X_all[tr].std(axis=0, ddof=0)
            sd = np.where(sd > 1e-12, sd, 1.0)
            X_tr_std = (X_all[tr] - mu) / sd
            X_te_std = (X_all[te] - mu) / sd
            fit_f = fit_cox_breslow(X_tr_std, times[tr], events[tr], l2_penalty=max(args.l2_penalty, 1e-4))
            eval_scores[te] = X_te_std @ fit_f["beta"]

        # Fit full-cohort L2-Cox strictly to report interpretable final feature weights (never used for AUC!)
        mu_full = X_all.mean(axis=0)
        sd_full = np.where(X_all.std(axis=0, ddof=0) > 1e-12, X_all.std(axis=0, ddof=0), 1.0)
        fit_full = fit_cox_breslow((X_all - mu_full) / sd_full, times, events, l2_penalty=max(args.l2_penalty, 1e-4))
        feature_weights = {feats[i]: round(float(fit_full["beta"][i]), 5) for i in range(len(feats))}

    c_idx = round(harrell_c_index(times, events, eval_scores), 4)
    med_score = float(np.median(eval_scores))
    risk_grp = np.where(eval_scores >= med_score, "High_Risk", "Low_Risk")
    lr_oof = logrank_test(times, events, risk_grp) if len(set(risk_grp)) == 2 else {"pvalue": None}

    scores_df = pd.DataFrame(
        {
            sample_col: sub_ids,
            "time": times,
            "event": events,
            "risk_score": np.round(eval_scores, 5),
            "risk_group_median": risk_grp,
        }
    )
    scores_path = os.path.join(args.output_dir, "risk_scores.tsv")
    scores_df.to_csv(scores_path, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    ax.plot([0, 1], [0, 1], color="#94A3B8", linestyle="--", linewidth=1.0, label="Chance (AUC = 0.50)")

    landmark_results = []
    auc_rows = []
    plotted = 0
    for idx_h, hz in enumerate(sorted(set(args.landmarks))):
        roc = ipcw_time_dependent_roc(times, events, eval_scores, float(hz))
        if "skipped" in roc:
            warnings.append(roc["skipped"])
            landmark_results.append({k: v for k, v in roc.items() if k not in ("fpr", "tpr")})
            continue
        color = PALETTE[plotted % len(PALETTE)]
        plotted += 1
        ax.plot(roc["fpr"], roc["tpr"], color=color, linewidth=2.0, label=f"t = {hz:g} (AUC = {roc['auc']:.3f})")
        entry = {
            "horizon": roc["horizon"],
            "n_cases": roc["n_cases"],
            "n_controls": roc["n_controls"],
            "auc": roc["auc"],
        }
        landmark_results.append(entry)
        auc_rows.append(entry)

    if plotted == 0:
        ax.text(0.5, 0.5, "No landmark had >= 3 cases and >= 3 controls", ha="center", va="center")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("1 - Specificity (IPCW False Positive Rate)")
    ax.set_ylabel("Sensitivity (IPCW True Positive Rate)")
    mode_tag = "Out-of-Fold CV L2-Cox" if eval_mode == "oof_cv_l2_cox" else args.score_col
    ax.set_title(f"Time-Dependent ROC ({mode_tag}, C-index = {c_idx:.3f})", fontsize=10)
    ax.legend(frameon=False, loc="lower right", fontsize=8.5)
    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "time-dependent-roc.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    auc_df = pd.DataFrame(auc_rows if auc_rows else [{"horizon": np.nan, "n_cases": 0, "n_controls": 0, "auc": np.nan}])
    auc_path = os.path.join(args.output_dir, "landmark_auc.tsv")
    auc_df.to_csv(auc_path, sep="\t", index=False)

    params_meta = {
        "evaluation_mode": eval_mode,
        "landmarks": sorted(set(args.landmarks)),
        "l2_penalty": args.l2_penalty if eval_mode == "oof_cv_l2_cox" else None,
    }
    in_files = (args.clinical, args.matrix)
    artifacts = [
        write_artifact(fig_path, "clin_time_roc", f"IPCW Time-Dependent ROC ({mode_tag})", inputs=in_files, parameters=params_meta),
        write_artifact(auc_path, "clin_time_roc", "Time-dependent ROC landmark AUCs", inputs=in_files, parameters=params_meta, shape=auc_df.shape),
        write_artifact(scores_path, "clin_time_roc", "Per-sample prognostic risk scores", inputs=in_files, parameters=params_meta, shape=scores_df.shape),
    ]

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "evaluation_mode": eval_mode,
            "n_samples": n_samples,
            "n_events": n_events,
            "c_index": c_idx,
            "oof_km_logrank_pvalue": lr_oof["pvalue"],
            "feature_weights": feature_weights,
            "landmarks": landmark_results,
            "scores_table": os.path.abspath(scores_path),
            "auc_table": os.path.abspath(auc_path),
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
