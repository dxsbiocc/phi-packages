#!/usr/bin/env python3
"""Multi-omics drug sensitivity biomarker discovery (MAF / binary / continuous features)."""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import bh_fdr, bool_arg, emit, fail, read_table, setup_matplotlib, write_artifact

SILENT_MAF_CLASSES = {
    "Silent",
    "Intron",
    "3'UTR",
    "5'UTR",
    "3'Flank",
    "5'Flank",
    "IGR",
    "RNA",
    "synonymous_variant",
}


def _load_drug_matrix(path: str, sample_col: str, drug_col: str, metric_col: str):
    """Return samples x drugs numeric DataFrame."""
    import pandas as pd

    raw = read_table(path)
    cols_lower = {c.lower(): c for c in raw.columns}

    s_col = sample_col or cols_lower.get("sample") or cols_lower.get("cell_line") or cols_lower.get("depmap_id")
    d_col = drug_col or cols_lower.get("drug") or cols_lower.get("compound") or cols_lower.get("drug_name")
    m_col = (
        metric_col
        or cols_lower.get("ln_ic50")
        or cols_lower.get("ln_ic50_censored")
        or cols_lower.get("ic50_censored")
        or cols_lower.get("ic50")
        or cols_lower.get("auc")
        or cols_lower.get("aac")
        or cols_lower.get("response")
    )

    if s_col and d_col and m_col and s_col in raw.columns and d_col in raw.columns and m_col in raw.columns:
        sub = raw[[s_col, d_col, m_col]].copy()
        sub[m_col] = pd.to_numeric(sub[m_col], errors="coerce")
        sub = sub.dropna(subset=[s_col, d_col, m_col])
        mat = sub.pivot_table(index=s_col, columns=d_col, values=m_col, aggfunc="mean")
        mat.index = mat.index.astype(str)
        mat.columns = mat.columns.astype(str)
        return mat

    # Wide matrix: first column is index (samples x drugs)
    idx_col = raw.columns[0]
    wide = raw.set_index(idx_col)
    wide = wide.apply(pd.to_numeric, errors="coerce")
    wide.index = wide.index.astype(str)
    wide.columns = wide.columns.astype(str)
    return wide


def _load_features(path: str, requested_type: str, drug_samples: set[str]):
    """Return (detected_type, features_x_samples DataFrame)."""
    import numpy as np
    import pandas as pd

    raw = read_table(path)
    cols_lower = {c.lower(): c for c in raw.columns}

    is_maf = requested_type == "maf" or (
        requested_type == "auto"
        and "hugo_symbol" in cols_lower
        and ("tumor_sample_barcode" in cols_lower or "sample" in cols_lower)
    )
    if is_maf:
        g_col = cols_lower.get("hugo_symbol")
        b_col = cols_lower.get("tumor_sample_barcode") or cols_lower.get("sample")
        vc_col = cols_lower.get("variant_classification")
        if not g_col or not b_col:
            fail("MAF feature input requires Hugo_Symbol and Tumor_Sample_Barcode (or sample) columns.")
        maf = raw.copy()
        if vc_col:
            maf = maf[~maf[vc_col].astype(str).isin(SILENT_MAF_CLASSES)]
        maf["_gene"] = maf[g_col].astype(str)
        maf["_sample"] = maf[b_col].astype(str)
        all_maf_samples = sorted(set(raw[b_col].dropna().astype(str).tolist()))
        ctab = pd.crosstab(maf["_gene"], maf["_sample"]).clip(upper=1).astype(float)
        for s in all_maf_samples:
            if s not in ctab.columns:
                ctab[s] = 0.0
        ctab = ctab[all_maf_samples]
        return "maf", ctab

    idx_col = raw.columns[0]
    mat = raw.set_index(idx_col).apply(pd.to_numeric, errors="coerce")
    mat.index = mat.index.astype(str)
    mat.columns = mat.columns.astype(str)

    # Orient so columns = samples overlapping drug_samples
    col_overlap = len(set(mat.columns) & drug_samples)
    idx_overlap = len(set(mat.index) & drug_samples)
    if idx_overlap > col_overlap:
        mat = mat.T

    vals = mat.to_numpy(dtype=float)
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        fail("Feature matrix contains no finite numeric values.")

    if requested_type == "binary":
        feat_type = "binary"
        mat = (mat > 0).astype(float).where(mat.notna())
    elif requested_type == "continuous":
        feat_type = "continuous"
    else:
        uniq = np.unique(np.round(finite, 6))
        feat_type = "binary" if set(uniq.tolist()).issubset({0.0, 1.0}) else "continuous"

    return feat_type, mat


def main() -> None:
    ap = argparse.ArgumentParser(description="Multi-omics drug sensitivity biomarker screening")
    ap.add_argument("--drug_response", required=True, help="Drug sensitivity matrix or long table")
    ap.add_argument("--features", required=True, help="MAF file or binary/continuous feature matrix")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--feature_type", choices=["auto", "maf", "binary", "continuous"], default="auto")
    ap.add_argument("--metric_col", default="", help="Metric column in long drug_response table")
    ap.add_argument("--sample_col", default="", help="Sample column in long drug_response table")
    ap.add_argument("--drug_col", default="", help="Drug column in long drug_response table")
    ap.add_argument(
        "--higher_is_sensitive",
        nargs="?",
        const=True,
        type=bool_arg,
        default=False,
        help="True when higher drug metric means greater sensitivity (e.g. AAC); false for IC50/LN_IC50/AUC",
    )
    ap.add_argument("--min_altered", type=int, default=3, help="Minimum altered and WT samples per binary test")
    ap.add_argument("--fdr", type=float, default=0.05, help="BH FDR cutoff (default 0.05)")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from scipy import stats

    drug_mat = _load_drug_matrix(args.drug_response, args.sample_col, args.drug_col, args.metric_col)
    feat_type, feat_mat = _load_features(args.features, args.feature_type, set(drug_mat.index))

    # Check if wide drug_mat was actually drugs x samples
    overlap_direct = [s for s in drug_mat.index if s in feat_mat.columns]
    overlap_transposed = [s for s in drug_mat.columns if s in feat_mat.columns]
    if len(overlap_transposed) > len(overlap_direct) and len(overlap_direct) < 4:
        drug_mat = drug_mat.T
        overlap_samples = overlap_transposed
    else:
        overlap_samples = overlap_direct

    if len(overlap_samples) < 5:
        fail(
            f"Too few overlapping samples ({len(overlap_samples)}) between drug_response and features (need >= 5)."
        )

    drug_sub = drug_mat.loc[overlap_samples]
    feat_sub = feat_mat[overlap_samples]
    polarity_sign = 1.0 if args.higher_is_sensitive else -1.0

    warnings: list[str] = []
    rows: list[dict] = []

    if feat_type in ("maf", "binary"):
        min_alt = max(2, int(args.min_altered))
        for drug in drug_sub.columns:
            y_s = drug_sub[drug]
            valid_s = y_s.notna()
            if int(valid_s.sum()) < 2 * min_alt:
                continue
            y_vals = y_s[valid_s].to_numpy(dtype=float)
            f_mat_v = feat_sub.loc[:, valid_s.index[valid_s]].to_numpy(dtype=float)
            for f_idx, feat_name in enumerate(feat_sub.index):
                x = f_mat_v[f_idx]
                alt_mask = x == 1.0
                wt_mask = x == 0.0
                n_alt = int(alt_mask.sum())
                n_wt = int(wt_mask.sum())
                if n_alt < min_alt or n_wt < min_alt:
                    continue
                y_alt = y_vals[alt_mask]
                y_wt = y_vals[wt_mask]
                mean_alt = float(np.mean(y_alt))
                mean_wt = float(np.mean(y_wt))
                delta = mean_alt - mean_wt
                var_alt = float(np.var(y_alt, ddof=1)) if n_alt > 1 else 0.0
                var_wt = float(np.var(y_wt, ddof=1)) if n_wt > 1 else 0.0
                pooled_sd = math.sqrt(((n_alt - 1) * var_alt + (n_wt - 1) * var_wt) / max(1, n_alt + n_wt - 2))
                cohens_d = delta / pooled_sd if pooled_sd > 1e-12 else 0.0
                try:
                    p_welch = float(stats.ttest_ind(y_alt, y_wt, equal_var=False).pvalue)
                except Exception:
                    p_welch = 1.0
                try:
                    p_mwu = float(stats.mannwhitneyu(y_alt, y_wt, alternative="two-sided").pvalue)
                except Exception:
                    p_mwu = 1.0
                if not math.isfinite(p_welch):
                    p_welch = 1.0
                sens_eff = polarity_sign * (-cohens_d if polarity_sign < 0 else cohens_d)
                # Wait: polarity_sign is -1 when lower=sensitive -> delta < 0 means sensitivity -> sens_eff = -cohens_d > 0!
                sens_eff = polarity_sign * cohens_d
                rows.append(
                    {
                        "drug": str(drug),
                        "feature": str(feat_name),
                        "feature_type": feat_type,
                        "n_samples": n_alt + n_wt,
                        "n_altered": n_alt,
                        "n_wt": n_wt,
                        "mean_altered": round(mean_alt, 4),
                        "mean_wt": round(mean_wt, 4),
                        "delta_metric": round(delta, 4),
                        "effect_size": round(cohens_d, 4),
                        "sensitivity_effect": round(sens_eff, 4),
                        "direction": "sensitivity" if sens_eff > 0 else ("resistance" if sens_eff < 0 else "neutral"),
                        "pvalue": p_welch,
                        "mannwhitney_pvalue": p_mwu,
                    }
                )
    else:
        # Pre-filter top 2000 variable continuous features if huge matrix
        if feat_sub.shape[0] > 2000:
            variances = feat_sub.var(axis=1)
            feat_sub = feat_sub.loc[variances.nlargest(2000).index]
            warnings.append("Filtered continuous features to top 2,000 most variable rows for biomarker screening.")

        for drug in drug_sub.columns:
            y_s = drug_sub[drug]
            for feat_name in feat_sub.index:
                x_s = feat_sub.loc[feat_name]
                valid = y_s.notna() & x_s.notna()
                n_eff = int(valid.sum())
                if n_eff < 6:
                    continue
                yv = y_s[valid].to_numpy(dtype=float)
                xv = x_s[valid].to_numpy(dtype=float)
                if np.std(yv) < 1e-12 or np.std(xv) < 1e-12:
                    continue
                sp = stats.spearmanr(xv, yv)
                rho = float(sp.statistic)
                pval = float(sp.pvalue) if math.isfinite(float(sp.pvalue)) else 1.0
                sens_eff = polarity_sign * rho
                rows.append(
                    {
                        "drug": str(drug),
                        "feature": str(feat_name),
                        "feature_type": "continuous",
                        "n_samples": n_eff,
                        "spearman_rho": round(rho, 4),
                        "effect_size": round(rho, 4),
                        "sensitivity_effect": round(sens_eff, 4),
                        "direction": "sensitivity" if sens_eff > 0 else ("resistance" if sens_eff < 0 else "neutral"),
                        "pvalue": pval,
                    }
                )

    if not rows:
        fail(
            "No testable (drug, feature) pairs passed sample-count and alteration-frequency filters "
            f"(overlap={len(overlap_samples)}, min_altered={args.min_altered})."
        )

    res_df = pd.DataFrame(rows)
    res_df["padj"] = bh_fdr(res_df["pvalue"].to_numpy(dtype=float))
    res_df = res_df.sort_values(["pvalue", "sensitivity_effect"], ascending=[True, False]).reset_index(drop=True)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    table_path = out_dir / "pgx_biomarkers.tsv"
    fig_path = out_dir / "pgx-biomarkers.png"
    res_df.to_csv(table_path, sep="\t", index=False)

    write_artifact(
        table_path,
        "pgx_biomarker",
        "Pharmacogenomic Sensitivity & Resistance Biomarker Associations",
        inputs=[args.drug_response, args.features],
        parameters={
            "feature_type": feat_type,
            "higher_is_sensitive": bool(args.higher_is_sensitive),
            "fdr": args.fdr,
        },
        shape=(int(len(res_df)), int(res_df.shape[1])),
        scientific_context={
            "dataType": "de_stats",
            "effectField": "sensitivity_effect",
            "fdrField": "padj",
            "higherIsSensitive": bool(args.higher_is_sensitive),
        },
    )

    # Render 2-panel volcano + top hit figure
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 5.0), dpi=180)

    eff = res_df["sensitivity_effect"].to_numpy(dtype=float)
    neglogq = -np.log10(np.clip(res_df["padj"].to_numpy(dtype=float), 1e-300, 1.0))
    sig_sens = (res_df["padj"] < args.fdr) & (res_df["sensitivity_effect"] > 0)
    sig_res = (res_df["padj"] < args.fdr) & (res_df["sensitivity_effect"] < 0)
    ns_mask = ~(sig_sens | sig_res)

    ax1.scatter(eff[ns_mask], neglogq[ns_mask], s=18, color="#9ca3af", alpha=0.6, label="NS")
    if sig_sens.any():
        ax1.scatter(eff[sig_sens], neglogq[sig_sens], s=28, color="#059669", alpha=0.88, label="Sensitivity (FDR)")
    if sig_res.any():
        ax1.scatter(eff[sig_res], neglogq[sig_res], s=28, color="#dc2626", alpha=0.88, label="Resistance (FDR)")
    ax1.axhline(-math.log10(args.fdr), color="#6b7280", ls="--", lw=1.0)
    ax1.axvline(0.0, color="#d1d5db", lw=0.9)
    for _, r in res_df.head(6).iterrows():
        yv = -math.log10(max(float(r["padj"]), 1e-300))
        ax1.annotate(
            f"{r['feature']} -> {r['drug']}",
            (float(r["sensitivity_effect"]), yv),
            fontsize=7.2,
            xytext=(4, 3),
            textcoords="offset points",
        )
    ax1.set_xlabel("Standardized Sensitivity Effect (+ Sensitive / - Resistant)")
    ax1.set_ylabel("-log10(BH FDR)")
    ax1.set_title("Pharmacogenomic Biomarker Volcano")
    ax1.legend(loc="best", frameon=False, fontsize=8)

    top = res_df.iloc[0]
    top_drug = str(top["drug"])
    top_feat = str(top["feature"])
    y_top = drug_sub[top_drug]
    x_top = feat_sub.loc[top_feat]
    valid_top = y_top.notna() & x_top.notna()
    if feat_type in ("maf", "binary"):
        wt_v = y_top[valid_top & (x_top == 0.0)].to_numpy(dtype=float)
        alt_v = y_top[valid_top & (x_top == 1.0)].to_numpy(dtype=float)
        bp = ax2.boxplot([wt_v, alt_v], tick_labels=[f"WT (n={len(wt_v)})", f"{top_feat} Alt (n={len(alt_v)})"], patch_artist=True)
        bp["boxes"][0].set_facecolor("#cbd5e1")
        bp["boxes"][1].set_facecolor("#6ee7b7" if top["sensitivity_effect"] > 0 else "#fca5a5")
        ax2.set_ylabel(f"{top_drug} Response Metric")
        ax2.set_title(f"Top Hit: {top_feat} vs {top_drug} (FDR={top['padj']:.2g})")
    else:
        xv = x_top[valid_top].to_numpy(dtype=float)
        yv = y_top[valid_top].to_numpy(dtype=float)
        ax2.scatter(xv, yv, s=28, color="#2563eb", alpha=0.8)
        ax2.set_xlabel(f"{top_feat} Feature Level")
        ax2.set_ylabel(f"{top_drug} Response Metric")
        ax2.set_title(f"Top Hit: {top_feat} vs {top_drug} (rho={top['spearman_rho']:.2f}, FDR={top['padj']:.2g})")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "pgx_biomarker",
        "Pharmacogenomic Biomarker Volcano & Top Association Plot",
        inputs=[args.drug_response, args.features],
        parameters={"top_drug": top_drug, "top_feature": top_feat},
    )

    n_sig = int((res_df["padj"] < args.fdr).sum())
    emit(
        {
            "status": "ok",
            "tool": "pgx_biomarker",
            "detected_feature_type": feat_type,
            "higher_is_sensitive": bool(args.higher_is_sensitive),
            "n_samples_overlap": int(len(overlap_samples)),
            "n_drugs_tested": int(res_df["drug"].nunique()),
            "n_features_tested": int(res_df["feature"].nunique()),
            "n_tests_total": int(len(res_df)),
            "n_significant_fdr": n_sig,
            "top_biomarkers": res_df.head(10).to_dict(orient="records"),
            "outputs": {
                "biomarkers_table": str(table_path),
                "biomarkers_plot": str(fig_path),
            },
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
