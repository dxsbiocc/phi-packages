#!/usr/bin/env python3
"""prot_qc_impute: contaminant filtering, log2 transform, normalization, imputation & QC."""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    clean_sample_name,
    contaminant_mask,
    detect_format_and_columns,
    emit,
    fail,
    infer_scale,
    load_metadata,
    read_table,
    write_artifact,
)


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Quantile normalization preserving NaN positions."""
    vals = df.to_numpy(dtype=float)
    out = np.full_like(vals, np.nan)
    n_rows, n_cols = vals.shape
    ranks = np.empty_like(vals)
    sorted_cols = []
    for j in range(n_cols):
        col = vals[:, j]
        valid = np.isfinite(col)
        if valid.sum() == 0:
            sorted_cols.append(np.full(n_rows, np.nan))
            continue
        s = np.sort(col[valid])
        # Interpolate to common length n_rows
        interp = np.interp(np.linspace(0, 1, n_rows), np.linspace(0, 1, len(s)), s)
        sorted_cols.append(interp)
        order = np.argsort(col[valid])
        r = np.empty(int(valid.sum()), dtype=int)
        r[order] = np.round(np.linspace(0, n_rows - 1, int(valid.sum()))).astype(int)
        ranks[valid, j] = r
    ref = np.nanmean(np.column_stack(sorted_cols), axis=1)
    for j in range(n_cols):
        valid = np.isfinite(vals[:, j])
        out[valid, j] = ref[ranks[valid, j].astype(int)]
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def minprob_impute(df: pd.DataFrame, downshift: float, width: float, seed: int) -> pd.DataFrame:
    """Perseus-style left-censored Gaussian MNAR imputation per sample column."""
    rng = np.random.default_rng(seed)
    out = df.copy()
    for col in out.columns:
        vals = out[col].to_numpy(dtype=float).copy()
        miss = ~np.isfinite(vals)
        if not np.any(miss):
            continue
        obs = vals[~miss]
        mu = float(np.mean(obs)) if len(obs) else 0.0
        sd = float(np.std(obs, ddof=1)) if len(obs) > 1 else 1.0
        sd = max(sd, 1e-3)
        imp_mu = mu - downshift * sd
        imp_sd = max(width * sd, 1e-3)
        vals[miss] = rng.normal(loc=imp_mu, scale=imp_sd, size=int(miss.sum()))
        out[col] = vals
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--id_col", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default="condition")
    parser.add_argument("--min_valid_frac", type=float, default=0.5)
    parser.add_argument("--filter_mode", choices=["in_at_least_one_group", "overall"], default="in_at_least_one_group")
    parser.add_argument("--normalization", choices=["median", "quantile", "none"], default="median")
    parser.add_argument("--imputation", choices=["minprob", "knn", "none"], default="minprob")
    parser.add_argument("--downshift", type=float, default=1.8)
    parser.add_argument("--width", type=float, default=0.3)
    parser.add_argument("--knn_k", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.input)
    det = detect_format_and_columns(df, args.id_col)
    id_col = det["id_col"]
    sample_cols = det["sample_cols"]
    if len(sample_cols) < 2:
        fail("fewer than 2 quantitative sample columns detected")

    n_raw = len(df)
    cmask = contaminant_mask(df, id_col)
    df = df.loc[~cmask].copy()
    n_after_contam = len(df)

    # Deduplicate IDs by keeping row with highest mean intensity
    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    mat.columns = [clean_sample_name(c) for c in sample_cols]
    ids = df[id_col].astype(str).str.strip()
    scale_before = infer_scale(mat)
    log2_applied = False
    if scale_before["scale"] == "raw_intensity":
        mat = mat.replace(0, np.nan)
        mat = np.log2(mat.where(mat > 0))
        log2_applied = True

    mat.insert(0, "protein_id", ids.to_numpy())
    mat["_mean"] = mat.iloc[:, 1:].mean(axis=1, skipna=True)
    mat = mat.sort_values("_mean", ascending=False).drop_duplicates("protein_id").drop(columns=["_mean"])
    mat = mat.set_index("protein_id")

    # Group-aware or overall valid-value filtering
    meta_groups = None
    used_filter_mode = args.filter_mode
    if args.metadata:
        meta, scol = load_metadata(args.metadata, args.sample_col)
        if args.group_col in meta.columns:
            meta_map = dict(zip(meta[scol], meta[args.group_col].astype(str)))
            meta_groups = {s: meta_map[s] for s in mat.columns if s in meta_map}
    if used_filter_mode == "in_at_least_one_group" and meta_groups:
        keep = np.zeros(len(mat), dtype=bool)
        for grp in sorted(set(meta_groups.values())):
            cols_g = [s for s, g in meta_groups.items() if g == grp]
            if not cols_g:
                continue
            valid_frac = mat[cols_g].notna().mean(axis=1).to_numpy()
            keep |= valid_frac >= args.min_valid_frac
    else:
        used_filter_mode = "overall"
        keep = (mat.notna().mean(axis=1).to_numpy()) >= args.min_valid_frac

    mat_filt = mat.loc[keep].copy()
    n_after_filter = len(mat_filt)
    if n_after_filter < 3:
        fail(f"only {n_after_filter} proteins passed valid-value filtering; lower --min_valid_frac")

    missing_before_impute = int(mat_filt.isna().sum().sum())
    missing_frac_before_impute = round(missing_before_impute / max(1, mat_filt.size), 4)

    # Normalization on log2 scale
    if args.normalization == "median":
        sample_medians = mat_filt.median(axis=0, skipna=True)
        global_med = float(sample_medians.median())
        mat_norm = mat_filt.sub(sample_medians, axis=1).add(global_med)
    elif args.normalization == "quantile":
        mat_norm = quantile_normalize(mat_filt)
    else:
        mat_norm = mat_filt.copy()

    # Imputation
    if args.imputation == "minprob":
        mat_proc = minprob_impute(mat_norm, args.downshift, args.width, args.seed)
    elif args.imputation == "knn":
        imputer = KNNImputer(n_neighbors=min(args.knn_k, max(1, mat_norm.shape[1] - 1)))
        mat_proc = pd.DataFrame(
            imputer.fit_transform(mat_norm.to_numpy(dtype=float)),
            index=mat_norm.index,
            columns=mat_norm.columns,
        )
    else:
        mat_proc = mat_norm.copy()

    # Sample QC summary
    qc_rows = []
    for s in mat_proc.columns:
        raw_col = mat_filt[s]
        proc_col = mat_proc[s]
        qc_rows.append(
            {
                "sample": s,
                "group": meta_groups.get(s, "NA") if meta_groups else "NA",
                "detected_after_filter": int(raw_col.notna().sum()),
                "missing_fraction": round(float(raw_col.isna().mean()), 4),
                "median_log2_before_norm": round(float(raw_col.median(skipna=True)), 4),
                "median_log2_after_norm": round(float(proc_col.median(skipna=True)), 4),
                "iqr_after_norm": round(float(proc_col.quantile(0.75) - proc_col.quantile(0.25)), 4),
            }
        )
    qc_df = pd.DataFrame(qc_rows)
    qc_path = os.path.join(args.output_dir, "sample_qc.tsv")
    qc_df.to_csv(qc_path, sep="\t", index=False)

    proc_path = os.path.join(args.output_dir, "processed_matrix.tsv")
    mat_proc.reset_index().to_csv(proc_path, sep="\t", index=False)

    figures: list[str] = []

    # Figure 1: Per-sample missingness barplot
    fig, ax = plt.subplots(figsize=(max(5.5, len(mat_proc.columns) * 0.55), 3.8))
    ax.bar(qc_df["sample"], qc_df["missing_fraction"] * 100.0, color="#38BDF8")
    ax.set_ylabel("Missing values (%)")
    ax.set_title("Per-sample missingness (post-filtering, pre-imputation)")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    p_miss = os.path.join(args.output_dir, "qc-missingness.png")
    fig.savefig(p_miss, dpi=150)
    plt.close(fig)
    figures.append(p_miss)

    # Figure 2: Normalized intensity boxplot
    fig, ax = plt.subplots(figsize=(max(5.5, len(mat_proc.columns) * 0.55), 4.0))
    ax.boxplot(
        [mat_proc[c].dropna().to_numpy() for c in mat_proc.columns],
        tick_labels=list(mat_proc.columns),
        showfliers=False,
    )
    ax.set_ylabel("Normalized log2 intensity")
    ax.set_title(f"Intensity distributions ({args.normalization} norm, {args.imputation} impute)")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    p_box = os.path.join(args.output_dir, "qc-intensity-boxplot.png")
    fig.savefig(p_box, dpi=150)
    plt.close(fig)
    figures.append(p_box)

    # Figure 3: PCA on complete features
    pca_info = None
    complete_mat = mat_proc.dropna(axis=0)
    if complete_mat.shape[0] >= 3 and complete_mat.shape[1] >= 2:
        X = complete_mat.T.to_numpy(dtype=float)
        X = X - X.mean(axis=0, keepdims=True)
        n_comp = min(2, X.shape[0], X.shape[1])
        pca = PCA(n_components=n_comp, random_state=args.seed)
        pcs = pca.fit_transform(X)
        var_exp = [round(float(v), 4) for v in pca.explained_variance_ratio_]
        pca_info = {"n_features_used": int(complete_mat.shape[0]), "explained_variance_ratio": var_exp}

        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        if meta_groups:
            uniq = sorted(set(meta_groups.get(s, "NA") for s in mat_proc.columns))
            cmap = plt.cm.tab10(np.linspace(0, 1, max(1, len(uniq))))
            for idx, grp in enumerate(uniq):
                mask = np.array([meta_groups.get(s, "NA") == grp for s in mat_proc.columns])
                ax.scatter(
                    pcs[mask, 0],
                    pcs[mask, 1] if n_comp > 1 else np.zeros(mask.sum()),
                    label=grp,
                    color=cmap[idx],
                    s=55,
                )
            ax.legend(frameon=False)
        else:
            ax.scatter(pcs[:, 0], pcs[:, 1] if n_comp > 1 else np.zeros(len(pcs)), s=55, color="#818CF8")
        for i, s in enumerate(mat_proc.columns):
            ax.annotate(s, (pcs[i, 0], pcs[i, 1] if n_comp > 1 else 0.0), fontsize=7, alpha=0.8)
        ax.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)")
        if n_comp > 1:
            ax.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)")
        ax.set_title("Proteomics sample PCA")
        fig.tight_layout()
        p_pca = os.path.join(args.output_dir, "qc-pca.png")
        fig.savefig(p_pca, dpi=150)
        plt.close(fig)
        figures.append(p_pca)

    params_meta = {
        "min_valid_frac": args.min_valid_frac,
        "filter_mode": used_filter_mode,
        "normalization": args.normalization,
        "imputation": args.imputation,
    }
    in_files = (args.input, args.metadata)
    fig_titles = {
        "qc-missingness.png": "Per-sample missingness barplot",
        "qc-intensity-boxplot.png": f"Normalized log2 intensity distributions ({args.normalization}, {args.imputation})",
        "qc-pca.png": "Proteomics sample PCA",
    }
    artifacts = [
        write_artifact(
            f,
            "prot_qc_impute",
            fig_titles.get(os.path.basename(f), os.path.basename(f)),
            inputs=in_files,
            parameters=params_meta,
        )
        for f in figures
    ]
    proc_out_df = mat_proc.reset_index()
    artifacts.extend(
        [
            write_artifact(
                proc_path,
                "prot_qc_impute",
                f"Processed log2 protein matrix ({mat_proc.shape[0]} proteins x {mat_proc.shape[1]} samples)",
                inputs=in_files,
                parameters=params_meta,
                shape=proc_out_df.shape,
            ),
            write_artifact(
                qc_path,
                "prot_qc_impute",
                "Proteomics per-sample QC summary",
                inputs=in_files,
                parameters=params_meta,
                shape=qc_df.shape,
            ),
        ]
    )
    proc_abs = os.path.abspath(proc_path)

    emit(
        {
            "input": os.path.abspath(args.input),
            "processed_matrix": proc_abs,
            "sample_qc_table": os.path.abspath(qc_path),
            "detected_format": det["format"],
            "n_features_raw": n_raw,
            "n_contaminants_removed": n_raw - n_after_contam,
            "n_features_kept": n_after_filter,
            "n_samples": int(mat_proc.shape[1]),
            "log2_transformed": log2_applied,
            "filter_mode": used_filter_mode,
            "min_valid_frac": args.min_valid_frac,
            "normalization": args.normalization,
            "imputation": args.imputation,
            "missing_fraction_before_impute": missing_frac_before_impute,
            "pca": pca_info,
            "figures": [os.path.abspath(f) for f in figures],
            "scientificContext": {
                proc_abs: {
                    "valueSemantics": "log_transformed",
                    "sampleIds": [str(s) for s in mat_proc.columns],
                }
            },
            "artifacts": artifacts,
        }
    )


if __name__ == "__main__":
    main()
