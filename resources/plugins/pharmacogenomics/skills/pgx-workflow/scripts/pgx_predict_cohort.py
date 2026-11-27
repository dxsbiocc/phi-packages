#!/usr/bin/env python3
"""pgx_predict_cohort: cross-cohort drug sensitivity imputation using RidgeCV (oncoPredict/Geeleher-style)."""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, read_table, setup_matplotlib, write_artifact


def _load_response_matrix(path: str):
    """Load training drug response as samples x drugs DataFrame (supports wide or long format)."""
    import pandas as pd

    raw = read_table(path)
    if raw.empty:
        fail(f"Training drug response table is empty: {path}")
    cols_lower = {c.lower().strip(): c for c in raw.columns}
    s_col = cols_lower.get("sample") or cols_lower.get("cell_line") or cols_lower.get("depmap_id")
    d_col = cols_lower.get("drug") or cols_lower.get("compound") or cols_lower.get("drug_name")
    m_col = (
        cols_lower.get("response")
        or cols_lower.get("ln_ic50")
        or cols_lower.get("ic50")
        or cols_lower.get("auc")
        or cols_lower.get("aac")
        or cols_lower.get("ic50_censored")
    )
    if s_col and d_col and m_col and len(raw.columns) <= 8:
        sub = raw[[s_col, d_col, m_col]].copy()
        sub[m_col] = pd.to_numeric(sub[m_col], errors="coerce")
        sub = sub.dropna(subset=[s_col, d_col, m_col])
        mat = sub.pivot_table(index=s_col, columns=d_col, values=m_col, aggfunc="mean")
        mat.index = mat.index.astype(str).str.strip()
        mat.columns = mat.columns.astype(str).str.strip()
        return mat

    idx_col = raw.columns[0]
    wide = raw.set_index(idx_col).apply(pd.to_numeric, errors="coerce")
    wide.index = wide.index.astype(str).str.strip()
    wide.columns = wide.columns.astype(str).str.strip()
    return wide


def _load_expr_raw(path: str):
    import pandas as pd

    raw = read_table(path)
    if raw.empty or raw.shape[1] < 2:
        fail(f"Expression matrix must contain ID column and >= 1 sample column: {path}")
    idx_col = raw.columns[0]
    mat = raw.set_index(idx_col).apply(pd.to_numeric, errors="coerce")
    mat.index = mat.index.astype(str).str.strip()
    mat.columns = mat.columns.astype(str).str.strip()
    if mat.isna().all().any():
        bad = mat.columns[mat.isna().all()].tolist()
        fail(f"Expression matrix {path} contains all-NaN columns: {bad}")
    return mat


def _to_log_scale(mat):
    import numpy as np

    X = mat.to_numpy(dtype=float, copy=True)
    finite = X[np.isfinite(X)]
    if finite.size == 0:
        return mat
    if (finite < 0).any():
        return mat
    if float(np.max(finite)) <= 35.0 and float(np.percentile(finite, 99)) <= 25.0:
        return mat
    # Convert raw counts or TPM/FPKM to log2(x + 1)
    return mat.__class__(np.log2(np.clip(X, 0.0, None) + 1.0), index=mat.index, columns=mat.columns)


def main() -> None:
    ap = argparse.ArgumentParser(description="Cross-cohort pharmacogenomic RidgeCV drug sensitivity prediction")
    ap.add_argument("--train_expr", required=True, help="Reference cell-line expression matrix TSV/CSV")
    ap.add_argument("--train_response", required=True, help="Reference cell-line drug response matrix or long table TSV/CSV")
    ap.add_argument("--clinical_expr", required=True, help="Target clinical cohort expression matrix TSV/CSV")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--min_cv_r", type=float, default=0.25, help="Minimum 5-fold CV Pearson r to mark model reliable (default 0.25)")
    ap.add_argument(
        "--response_polarity",
        choices=["lower_is_sensitive", "higher_is_sensitive"],
        default="lower_is_sensitive",
        help="Drug response polarity (lower_is_sensitive for IC50/AUC; higher_is_sensitive for AAC)",
    )
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from scipy import stats
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import KFold

    resp_df = _load_response_matrix(args.train_response)
    train_raw = _load_expr_raw(args.train_expr)
    clin_raw = _load_expr_raw(args.clinical_expr)

    # Auto-orient train_expr (genes x cell_lines) and resp_df (cell_lines x drugs)
    train_ids = set(train_raw.index) | set(train_raw.columns)
    if len(set(resp_df.columns) & train_ids) > len(set(resp_df.index) & train_ids):
        resp_df = resp_df.T

    resp_lines = set(resp_df.index)
    if len(set(train_raw.index) & resp_lines) > len(set(train_raw.columns) & resp_lines):
        train_expr = train_raw.T
    else:
        train_expr = train_raw

    # Auto-orient clinical_expr (genes x patients) by matching gene IDs with train_expr.index
    train_genes_exact = set(train_expr.index)
    train_genes_upper = {g.upper() for g in train_expr.index}
    idx_ov = max(len(set(clin_raw.index) & train_genes_exact), len({g.upper() for g in clin_raw.index} & train_genes_upper))
    col_ov = max(len(set(clin_raw.columns) & train_genes_exact), len({g.upper() for g in clin_raw.columns} & train_genes_upper))
    clin_expr = clin_raw.T if col_ov > idx_ov else clin_raw

    shared_lines = [s for s in resp_df.index if s in train_expr.columns]
    if len(shared_lines) < 5:
        fail(
            f"Too few overlapping cell lines ({len(shared_lines)}) between train_expr and train_response (need >= 5)."
        )

    warnings: list[str] = []
    if len(shared_lines) < 15:
        warnings.append(
            f"Training cohort has only {len(shared_lines)} overlapping cell lines; RidgeCV generalizability may be limited."
        )

    train_sub = train_expr[shared_lines].copy()
    resp_sub = resp_df.loc[shared_lines].copy()

    # Deduplicate/align shared genes (case-insensitive fallback if exact overlap < 20)
    exact_genes = [g for g in train_sub.index if g in clin_expr.index]
    if len(exact_genes) < 20:
        train_sub.index = pd.Index([g.upper() for g in train_sub.index])
        clin_expr = clin_expr.copy()
        clin_expr.index = pd.Index([g.upper() for g in clin_expr.index])
        train_sub = train_sub[~train_sub.index.duplicated(keep="first")]
        clin_expr = clin_expr[~clin_expr.index.duplicated(keep="first")]
        shared_genes = [g for g in train_sub.index if g in clin_expr.index]
    else:
        train_sub = train_sub[~train_sub.index.duplicated(keep="first")]
        clin_expr = clin_expr[~clin_expr.index.duplicated(keep="first")]
        shared_genes = [g for g in train_sub.index if g in clin_expr.index]

    if len(shared_genes) < 20:
        fail(
            f"Only {len(shared_genes)} shared genes between train_expr and clinical_expr (need >= 20)."
        )

    train_aligned = _to_log_scale(train_sub.loc[shared_genes]).fillna(0.0)
    clin_aligned = _to_log_scale(clin_expr.loc[shared_genes]).fillna(0.0)

    # Filter non-zero-variance genes across training (and clinical if n >= 3)
    tr_var = train_aligned.var(axis=1).to_numpy(dtype=float, copy=True)
    if clin_aligned.shape[1] >= 3:
        cl_var = clin_aligned.var(axis=1).to_numpy(dtype=float, copy=True)
        keep_mask = (tr_var > 1e-8) & (cl_var > 1e-8)
    else:
        keep_mask = tr_var > 1e-8

    kept_genes = [g for g, k in zip(shared_genes, keep_mask) if k]
    if len(kept_genes) < 20:
        fail(
            f"Only {len(kept_genes)} shared non-zero-variance genes remained after filtering (need >= 20)."
        )

    if len(kept_genes) > 5000:
        top_idx = np.argsort(-tr_var[keep_mask])[:5000]
        kept_genes = [kept_genes[i] for i in top_idx]
        warnings.append("Selected top 5,000 shared genes by training variance for RidgeCV modeling.")

    X_tr_raw = train_aligned.loc[kept_genes].to_numpy(dtype=float, copy=True).T  # (n_lines, n_genes)
    X_cl_raw = clin_aligned.loc[kept_genes].to_numpy(dtype=float, copy=True).T  # (n_patients, n_genes)
    patients = clin_aligned.columns.tolist()

    # Cross-platform gene-wise empirical Bayes / Z-score harmonization
    tr_mu = X_tr_raw.mean(axis=0, keepdims=True)
    tr_sd = X_tr_raw.std(axis=0, ddof=1, keepdims=True)
    tr_sd = np.where(tr_sd > 1e-12, tr_sd, 1.0)
    Z_train = (X_tr_raw - tr_mu) / tr_sd

    n_clin = X_cl_raw.shape[0]
    if n_clin >= 5:
        cl_mu = X_cl_raw.mean(axis=0, keepdims=True)
        cl_sd = X_cl_raw.std(axis=0, ddof=1, keepdims=True)
        cl_sd = np.where(cl_sd > 1e-12, cl_sd, 1.0)
        Z_clin = (X_cl_raw - cl_mu) / cl_sd
    elif n_clin >= 2:
        # Empirical Bayes shrinkage toward training location/scale for tiny clinical cohorts
        w_eb = float(n_clin) / float(n_clin + 5)
        cl_mu = w_eb * X_cl_raw.mean(axis=0, keepdims=True) + (1.0 - w_eb) * tr_mu
        cl_sd = w_eb * X_cl_raw.std(axis=0, ddof=1, keepdims=True) + (1.0 - w_eb) * tr_sd
        cl_sd = np.where(cl_sd > 1e-12, cl_sd, 1.0)
        Z_clin = (X_cl_raw - cl_mu) / cl_sd
        warnings.append(
            f"Clinical cohort has only {n_clin} samples; applied empirical-Bayes shrinkage toward training gene moments."
        )
    else:
        Z_clin = (X_cl_raw - tr_mu) / tr_sd
        warnings.append("Single clinical patient provided; standardized against training cohort moments.")

    alphas = np.logspace(-2, 5, 25)
    cv_metrics: list[dict] = []
    pred_dict: dict[str, np.ndarray] = {}

    for drug in resp_sub.columns:
        y_series = resp_sub[drug]
        valid = y_series.notna().to_numpy(dtype=bool, copy=True)
        n_valid = int(valid.sum())
        if n_valid < 5:
            warnings.append(f"Skipped drug '{drug}': fewer than 5 non-missing training cell lines.")
            continue

        X_d = Z_train[valid]
        y_d = y_series.to_numpy(dtype=float, copy=True)[valid]
        if float(np.std(y_d)) <= 1e-12:
            warnings.append(f"Skipped drug '{drug}': zero response variance in training cohort.")
            continue

        n_splits = min(5, n_valid)
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        oof = np.zeros(n_valid, dtype=float)
        for tr_idx, va_idx in kf.split(X_d):
            fold_model = RidgeCV(alphas=alphas)
            fold_model.fit(X_d[tr_idx], y_d[tr_idx])
            oof[va_idx] = fold_model.predict(X_d[va_idx])

        if float(np.std(oof)) > 1e-12:
            pr = stats.pearsonr(y_d, oof)
            cv_r = float(pr.statistic)
            cv_p = float(pr.pvalue) if math.isfinite(float(pr.pvalue)) else 1.0
            sp = stats.spearmanr(y_d, oof)
            cv_rho = float(sp.statistic) if math.isfinite(float(sp.statistic)) else 0.0
        else:
            cv_r, cv_p, cv_rho = 0.0, 1.0, 0.0
        cv_rmse = float(np.sqrt(np.mean((y_d - oof) ** 2)))

        full_model = RidgeCV(alphas=alphas)
        full_model.fit(X_d, y_d)
        pred_clin = full_model.predict(Z_clin)
        pred_dict[str(drug)] = pred_clin

        reliable = bool(cv_r >= float(args.min_cv_r))
        low_warn = not reliable
        if low_warn:
            warnings.append(
                f"Drug '{drug}' has low CV predictability (5-fold r={cv_r:.3f} < {args.min_cv_r:.2f}); "
                "interpret predicted cohort sensitivity with caution."
            )

        cv_metrics.append(
            {
                "drug": str(drug),
                "n_train_lines": n_valid,
                "n_features": int(len(kept_genes)),
                "cv_folds": int(n_splits),
                "cv_pearson_r": round(cv_r, 4),
                "cv_pearson_p": float(cv_p),
                "cv_spearman_rho": round(cv_rho, 4),
                "cv_rmse": round(cv_rmse, 4),
                "selected_alpha": float(full_model.alpha_),
                "reliable": reliable,
                "low_predictability_warning": low_warn,
                "response_polarity": args.response_polarity,
            }
        )

    if not pred_dict:
        fail("No drug models could be trained (check sample overlap and response values).")

    pred_df = pd.DataFrame(pred_dict, index=patients)
    pred_df.index.name = "sample"
    metrics_df = pd.DataFrame(cv_metrics).sort_values("cv_pearson_r", ascending=False).reset_index(drop=True)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pred_path = out_dir / "cohort_predicted_response.tsv"
    metrics_path = out_dir / "drug_model_cv_metrics.tsv"
    fig_path = out_dir / "cohort_drug_prediction_overview.png"

    pred_df.round(5).to_csv(pred_path, sep="\t")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)

    inputs_list = [args.train_expr, args.train_response, args.clinical_expr]
    params = {"min_cv_r": float(args.min_cv_r), "response_polarity": args.response_polarity}
    write_artifact(
        pred_path,
        "pgx_predict_cohort",
        "Predicted Clinical Cohort Drug Response Matrix (RidgeCV)",
        inputs=inputs_list,
        parameters=params,
        shape=pred_df.shape,
        scientific_context={
            "dataType": "drug_response",
            "valueScale": "predicted_response",
            "responsePolarity": args.response_polarity,
            "higherIsSensitive": args.response_polarity == "higher_is_sensitive",
        },
    )
    write_artifact(
        metrics_path,
        "pgx_predict_cohort",
        "5-Fold Cross-Validation Performance of Cell-Line RidgeCV Drug Models",
        inputs=inputs_list,
        parameters=params,
        shape=metrics_df.shape,
    )

    # Render 2-panel CV performance + cohort predicted sensitivity heatmap
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.0), dpi=180, gridspec_kw={"width_ratios": [1.0, 1.35]})

    y_pos = np.arange(len(metrics_df))
    r_vals = metrics_df["cv_pearson_r"].to_numpy(dtype=float, copy=True)
    bar_cols = ["#059669" if ok else "#f59e0b" for ok in metrics_df["reliable"]]
    ax1.barh(y_pos, r_vals, color=bar_cols, alpha=0.88)
    ax1.axvline(float(args.min_cv_r), color="#dc2626", ls="--", lw=1.1, label=f"Reliable cutoff ({args.min_cv_r:.2f})")
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(metrics_df["drug"].tolist(), fontsize=8.5)
    ax1.invert_yaxis()
    ax1.set_xlabel("5-Fold CV Pearson r")
    ax1.set_title("Training Cell-Line Model Predictability")
    ax1.legend(loc="lower right", frameon=False, fontsize=8)

    # Orient heatmap color so positive z-score = greater predicted sensitivity
    ordered_drugs = metrics_df["drug"].tolist()[: min(20, len(metrics_df))]
    mat = pred_df[ordered_drugs].to_numpy(dtype=float, copy=True)
    if args.response_polarity == "lower_is_sensitive":
        mat = -mat
    if mat.shape[0] >= 2:
        col_mu = mat.mean(axis=0, keepdims=True)
        col_sd = mat.std(axis=0, ddof=1, keepdims=True)
        col_sd = np.where(col_sd > 1e-12, col_sd, 1.0)
        z_mat = (mat - col_mu) / col_sd
    else:
        z_mat = mat - mat.mean(axis=1, keepdims=True)
    show_n = min(40, z_mat.shape[0])
    z_plot = z_mat[:show_n, :].T
    vmax = max(1.5, float(np.nanpercentile(np.abs(z_plot), 98)))
    im = ax2.imshow(z_plot, cmap="RdYlGn", vmin=-vmax, vmax=vmax, aspect="auto")
    ax2.set_yticks(range(len(ordered_drugs)))
    ax2.set_yticklabels(ordered_drugs, fontsize=8.5)
    if show_n <= 20:
        ax2.set_xticks(range(show_n))
        ax2.set_xticklabels(patients[:show_n], rotation=45, ha="right", fontsize=7.5)
    else:
        ax2.set_xticks([])
        ax2.set_xlabel(f"Clinical Patients (showing {show_n}/{len(patients)})")
    fig.colorbar(im, ax=ax2, fraction=0.04, pad=0.03, label="Predicted Sensitivity Z (+ = Sensitive)")
    ax2.set_title("Predicted Patient Cohort Drug Sensitivity")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "pgx_predict_cohort",
        "Cross-Cohort Pharmacogenomic RidgeCV Performance & Patient Sensitivity Map",
        inputs=inputs_list,
        parameters=params,
    )

    n_rel = int(metrics_df["reliable"].sum())
    emit(
        {
            "status": "ok",
            "tool": "pgx_predict_cohort",
            "response_polarity": args.response_polarity,
            "min_cv_r": float(args.min_cv_r),
            "n_train_samples": int(len(shared_lines)),
            "n_clinical_samples": int(len(patients)),
            "n_shared_genes": int(len(kept_genes)),
            "n_drugs_modeled": int(len(metrics_df)),
            "n_drugs_reliable": n_rel,
            "n_drugs_low_predictability": int(len(metrics_df) - n_rel),
            "drug_metrics": metrics_df.to_dict(orient="records"),
            "outputs": {
                "predicted_response": str(pred_path),
                "drug_cv_metrics": str(metrics_path),
                "overview_plot": str(fig_path),
            },
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
