#!/usr/bin/env python3
"""mg_classify: binary case-control biomarker classification using leak-free
nested/stratified k-fold cross-validation, in-fold mRMR or Mann-Whitney
feature selection, and L2-regularized logistic regression (pure numpy/scipy).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import mannwhitneyu

from _common import emit, fail, read_matrix, read_metadata, setup_matplotlib, write_artifact


def stratified_folds(y, k, seed):
    rng = np.random.default_rng(seed)
    folds = [[] for _ in range(k)]
    for cls in (0, 1):
        idx = np.where(y == cls)[0]
        rng.shuffle(idx)
        for i, sample_idx in enumerate(idx):
            folds[i % k].append(int(sample_idx))
    return [np.array(sorted(f), dtype=int) for f in folds]


def select_features_mrmr(x_tr, y_tr, top_k):
    """Greedy mRMR (F-statistic relevance minus mean |r| redundancy)."""
    n_feat = x_tr.shape[1]
    k = min(top_k, n_feat)
    # Relevance: point-biserial |r| (equivalent to sqrt(F) up to monotone transform)
    y_c = y_tr - y_tr.mean()
    y_norm = np.linalg.norm(y_c)
    x_c = x_tr - x_tr.mean(axis=0, keepdims=True)
    x_norm = np.linalg.norm(x_c, axis=0)
    valid = x_norm > 1e-12
    rel = np.zeros(n_feat, dtype=float)
    if y_norm > 0:
        rel[valid] = np.abs((x_c[:, valid].T @ y_c) / (x_norm[valid] * y_norm))

    # Pre-screen top 200 by relevance for speed on wide matrices
    cand_limit = min(n_feat, max(200, k * 4))
    candidates = list(np.argsort(rel)[::-1][:cand_limit])
    selected = []
    if not candidates:
        return selected
    first = candidates.pop(0)
    selected.append(int(first))

    # Cache normalized columns for redundancy
    x_unit = np.zeros_like(x_c)
    x_unit[:, valid] = x_c[:, valid] / x_norm[valid]

    while len(selected) < k and candidates:
        best_score = -1e9
        best_pos = 0
        sel_mat = x_unit[:, selected]  # n_tr x len(selected)
        for pos, c in enumerate(candidates):
            red = float(np.abs(sel_mat.T @ x_unit[:, c]).mean())
            score = rel[c] - red
            if score > best_score:
                best_score = score
                best_pos = pos
        selected.append(int(candidates.pop(best_pos)))
    return selected


def select_features_mannwhitney(x_tr, y_tr, top_k):
    n_feat = x_tr.shape[1]
    k = min(top_k, n_feat)
    pvals = np.ones(n_feat, dtype=float)
    x0 = x_tr[y_tr == 0]
    x1 = x_tr[y_tr == 1]
    for j in range(n_feat):
        if np.std(x_tr[:, j]) <= 1e-12:
            continue
        try:
            pvals[j] = mannwhitneyu(x1[:, j], x0[:, j], alternative="two-sided").pvalue
        except ValueError:
            pvals[j] = 1.0
    return [int(i) for i in np.argsort(pvals)[:k]]


def fit_l2_logistic(x_tr, y_tr, l2_alpha=1.0):
    n, p = x_tr.shape
    w0 = np.zeros(p + 1, dtype=float)

    def loss_and_grad(params):
        w = params[:p]
        b = params[p]
        z = x_tr @ w + b
        prob = expit(z)
        # Stable log-loss
        loss = float(np.mean(np.logaddexp(0.0, z) - y_tr * z)) + 0.5 * l2_alpha * float(np.dot(w, w)) / n
        diff = (prob - y_tr) / n
        grad_w = x_tr.T @ diff + (l2_alpha / n) * w
        grad_b = float(diff.sum())
        return loss, np.concatenate([grad_w, [grad_b]])

    res = minimize(loss_and_grad, w0, jac=True, method="L-BFGS-B")
    return res.x[:p], float(res.x[p])


def roc_auc_and_youden(y_true, y_score):
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    order = np.argsort(y_score)[::-1]
    y_sorted = y_true[order]
    s_sorted = y_score[order]
    n_pos = int((y_true == 1).sum())
    n_neg = int((y_true == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return None

    tps = np.cumsum(y_sorted == 1)
    fps = np.cumsum(y_sorted == 0)
    tpr = np.concatenate([[0.0], tps / n_pos])
    fpr = np.concatenate([[0.0], fps / n_neg])
    thresholds = np.concatenate([[s_sorted[0] + 1e-6], s_sorted])
    auroc = float(np.trapezoid(tpr, fpr))

    youden = tpr - fpr
    best_i = int(np.argmax(youden))
    return {
        "auroc": auroc,
        "youden_threshold": float(thresholds[best_i]),
        "sensitivity": float(tpr[best_i]),
        "specificity": float(1.0 - fpr[best_i]),
        "fpr": fpr,
        "tpr": tpr,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--group_col", default="group")
    ap.add_argument("--case", required=True)
    ap.add_argument("--control", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--top_features", type=int, default=15)
    ap.add_argument("--selector", choices=["mrmr", "mannwhitney"], default="mrmr")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--l2_alpha", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    matrix = read_matrix(args.input)  # features x samples
    samples = matrix.columns.tolist()
    meta = read_metadata(args.metadata, samples)
    if args.group_col not in meta.columns:
        fail(f"metadata has no column {args.group_col}; columns: {meta.columns.tolist()}")

    groups = meta[args.group_col].astype(str)
    mask = groups.isin([args.case, args.control])
    kept_samples = [s for s, m in zip(samples, mask) if m]
    if len(kept_samples) < 6:
        fail(f"classification needs at least 6 samples across case/control (got {len(kept_samples)})")

    y = np.array([1 if groups.loc[s] == args.case else 0 for s in kept_samples], dtype=int)
    n_case = int((y == 1).sum())
    n_ctrl = int((y == 0).sum())
    k = min(args.folds, n_case, n_ctrl)
    if k < 2:
        fail(f"need at least 2 samples per class for cross-validation (case={n_case}, control={n_ctrl})")

    sub_df = matrix[kept_samples]
    # Drop zero-variance features globally (unsupervised)
    sub_df = sub_df.loc[sub_df.std(axis=1) > 0]
    feat_names = sub_df.index.tolist()
    if len(feat_names) == 0:
        fail("no variable features in matrix")

    x_all = np.log1p(np.maximum(sub_df.to_numpy(dtype=float).T, 0.0))  # samples x features
    folds = stratified_folds(y, k, args.seed)
    oof_prob = np.zeros(len(y), dtype=float)
    fold_ids = np.zeros(len(y), dtype=int)
    fold_aurocs = []
    coef_records = {f: [] for f in feat_names}

    for fold_idx, test_idx in enumerate(folds):
        train_mask = np.ones(len(y), dtype=bool)
        train_mask[test_idx] = False
        train_idx = np.where(train_mask)[0]

        x_tr_raw, y_tr = x_all[train_idx], y[train_idx]
        x_te_raw, y_te = x_all[test_idx], y[test_idx]

        # Standardize strictly inside training fold
        mu = x_tr_raw.mean(axis=0)
        sd = x_tr_raw.std(axis=0)
        sd[sd < 1e-12] = 1.0
        x_tr = (x_tr_raw - mu) / sd
        x_te = (x_te_raw - mu) / sd

        if args.selector == "mrmr":
            sel = select_features_mrmr(x_tr, y_tr, args.top_features)
        else:
            sel = select_features_mannwhitney(x_tr, y_tr, args.top_features)

        w, b = fit_l2_logistic(x_tr[:, sel], y_tr, l2_alpha=args.l2_alpha)
        probs = expit(x_te[:, sel] @ w + b)
        oof_prob[test_idx] = probs
        fold_ids[test_idx] = fold_idx + 1

        fold_eval = roc_auc_and_youden(y_te, probs)
        if fold_eval is not None:
            fold_aurocs.append(fold_eval["auroc"])

        for j, weight in zip(sel, w):
            coef_records[feat_names[j]].append(float(weight))

    pooled = roc_auc_and_youden(y, oof_prob)
    if pooled is None:
        fail("unable to compute ROC AUC on out-of-fold predictions")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    params = {
        "case": args.case,
        "control": args.control,
        "selector": args.selector,
        "top_features": args.top_features,
        "folds": k,
        "l2_alpha": args.l2_alpha,
        "seed": args.seed,
    }
    inputs = [args.input, args.metadata]

    pred_df = pd.DataFrame(
        {
            "sample": kept_samples,
            "group": [groups.loc[s] for s in kept_samples],
            "y_true": y,
            "fold": fold_ids,
            "oof_prob": np.round(oof_prob, 6),
            "oof_pred_youden": (oof_prob >= pooled["youden_threshold"]).astype(int),
        }
    )
    pred_path = outdir / "oof-predictions.tsv"
    pred_df.to_csv(pred_path, sep="\t", index=False)

    imp_rows = []
    for f, ws in coef_records.items():
        if not ws:
            continue
        imp_rows.append(
            {
                "feature": f,
                "folds_selected": len(ws),
                "selection_frequency": round(len(ws) / k, 4),
                "mean_coef": round(float(np.mean(ws)), 4),
                "abs_mean_coef": round(float(abs(np.mean(ws))), 4),
            }
        )
    imp_rows.sort(key=lambda r: (-r["folds_selected"], -r["abs_mean_coef"], r["feature"]))
    imp_df = pd.DataFrame(imp_rows)
    imp_path = outdir / "feature-importance.tsv"
    imp_df.to_csv(imp_path, sep="\t", index=False)

    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.plot(
        pooled["fpr"],
        pooled["tpr"],
        color="#0EA5E9",
        linewidth=2.4,
        label=f"Pooled OOF ROC (AUROC = {pooled['auroc']:.3f})",
    )
    ax.plot([0, 1], [0, 1], linestyle="--", color="#94A3B8", linewidth=1.2)
    ax.scatter(
        [1.0 - pooled["specificity"]],
        [pooled["sensitivity"]],
        color="#E11D48",
        s=55,
        zorder=3,
        label=f"Youden (Sens={pooled['sensitivity']:.2f}, Spec={pooled['specificity']:.2f})",
    )
    ax.set_xlabel("False Positive Rate (1 - Specificity)")
    ax.set_ylabel("True Positive Rate (Sensitivity)")
    ax.set_title(f"{args.case} vs {args.control} ({k}-fold CV, {args.selector})")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    roc_path = outdir / "roc-curve.png"
    fig.savefig(roc_path, dpi=200)
    plt.close(fig)

    artifacts = [
        write_artifact(pred_path, "mg_classify", "Out-of-fold classification predictions",
                       skill="profiling-workflow", inputs=inputs, parameters=params, shape=pred_df.shape),
        write_artifact(imp_path, "mg_classify", "Cross-validated biomarker feature importance",
                       skill="profiling-workflow", inputs=inputs, parameters=params, shape=imp_df.shape),
        write_artifact(roc_path, "mg_classify", "Pooled out-of-fold ROC curve",
                       skill="profiling-workflow", inputs=inputs, parameters=params),
    ]

    emit(
        {
            "tool": "mg_classify",
            "case": args.case,
            "control": args.control,
            "n_case": n_case,
            "n_control": n_ctrl,
            "folds": int(k),
            "selector": args.selector,
            "top_features": int(args.top_features),
            "oof_auroc": round(pooled["auroc"], 4),
            "fold_auroc_mean": round(float(np.mean(fold_aurocs)), 4) if fold_aurocs else None,
            "fold_auroc_sd": round(float(np.std(fold_aurocs)), 4) if fold_aurocs else None,
            "youden_threshold": round(pooled["youden_threshold"], 4),
            "sensitivity": round(pooled["sensitivity"], 4),
            "specificity": round(pooled["specificity"], 4),
            "top_biomarkers": imp_rows[:10],
            "outputs": {
                "predictions": str(pred_path),
                "feature_importance": str(imp_path),
                "roc_figure": str(roc_path),
            },
            "artifacts": artifacts,
            "caveat": "Feature selection and scaling are performed strictly inside each training fold "
                      "to prevent leakage, but small cohorts still carry high variance; validate "
                      "diagnostic performance on an independent external cohort before clinical claims.",
        }
    )


if __name__ == "__main__":
    main()
