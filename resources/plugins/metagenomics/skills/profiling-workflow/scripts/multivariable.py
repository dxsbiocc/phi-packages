#!/usr/bin/env python3
"""mg_multivariable: compositionally-aware multivariable differential abundance
via Centered Log-Ratio (CLR) transformation + OLS linear regression with
covariate adjustment and Benjamini-Hochberg FDR control (pure numpy/scipy).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from _common import bh_fdr, emit, fail, read_matrix, read_metadata, setup_matplotlib, write_artifact


def clr_transform(df, pseudocount=None):
    """Apply per-sample Centered Log-Ratio transform on features x samples DataFrame."""
    x = df.to_numpy(dtype=float)
    if (x < 0).any():
        fail("CLR transform requires non-negative abundances")
    if pseudocount is None:
        pos = x[x > 0]
        if len(pos) == 0:
            fail("matrix contains no positive abundances")
        eps = float(pos.min()) * 0.5
    else:
        eps = float(pseudocount)
        if eps <= 0:
            fail("pseudocount must be > 0")
    log_x = np.log(x + eps)
    # Subtract per-sample (column) mean log abundance
    clr = log_x - log_x.mean(axis=0, keepdims=True)
    return clr, eps


def build_design_matrix(meta, exposure, reference, covariates):
    cols = [exposure] + [c for c in covariates if c != exposure]
    for c in cols:
        if c not in meta.columns:
            fail(f"metadata missing column {c}; available: {meta.columns.tolist()}")

    sub = meta[cols].copy()
    valid_mask = pd.Series(True, index=sub.index)
    for c in cols:
        valid_mask &= sub[c].notna() & (sub[c].astype(str).str.strip() != "")
    sub = sub.loc[valid_mask]
    if len(sub) < 5:
        fail(f"fewer than 5 complete samples after dropping missing metadata ({len(sub)} left)")

    parts = [pd.Series(1.0, index=sub.index, name="intercept")]
    exposure_terms = []

    for col in cols:
        s = sub[col]
        num = pd.to_numeric(s, errors="coerce")
        is_numeric = num.notna().all() and num.nunique() > 2 and (col != exposure or reference is None)
        if is_numeric:
            series = num.astype(float)
            series.name = col
            parts.append(series)
            if col == exposure:
                exposure_terms.append(col)
        else:
            cats = s.astype(str).str.strip()
            levels = sorted(cats.unique().tolist())
            if len(levels) < 2:
                fail(f"column {col} has fewer than 2 unique levels: {levels}")
            if col == exposure and reference is not None:
                if reference not in levels:
                    fail(f"reference {reference!r} not found in {col} levels: {levels}")
                ref_lvl = reference
            else:
                ref_lvl = levels[0]
            for lvl in levels:
                if lvl == ref_lvl:
                    continue
                term = f"{col}[{lvl}_vs_{ref_lvl}]"
                parts.append(pd.Series((cats == lvl).astype(float), index=sub.index, name=term))
                if col == exposure:
                    exposure_terms.append(term)

    design = pd.concat(parts, axis=1)
    n, p = design.shape
    if n <= p + 1:
        fail(f"insufficient residual degrees of freedom: {n} samples for {p} design parameters")
    return design, exposure_terms, sub.index.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--exposure", default="group")
    ap.add_argument("--reference", default=None)
    ap.add_argument("--covariates", nargs="+", default=[])
    ap.add_argument("--min_prevalence", type=float, default=0.1)
    ap.add_argument("--pseudocount", type=float, default=None)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    matrix = read_matrix(args.input)  # features x samples
    samples = matrix.columns.tolist()
    meta = read_metadata(args.metadata, samples)

    cov_list = []
    for item in args.covariates:
        for part in str(item).split(","):
            part = part.strip()
            if part and part not in cov_list:
                cov_list.append(part)

    design, exposure_terms, kept_samples = build_design_matrix(
        meta, args.exposure, args.reference, cov_list
    )
    sub_mat = matrix[kept_samples]
    prev = (sub_mat > 0).sum(axis=1) / float(len(kept_samples))
    # Compute CLR on full non-zero feature space before prevalence filtering
    non_empty = sub_mat.loc[sub_mat.sum(axis=1) > 0]
    if len(non_empty) < 2:
        fail("fewer than 2 non-zero features in matrix")
    clr_all, eps = clr_transform(non_empty, pseudocount=args.pseudocount)
    clr_df = pd.DataFrame(clr_all, index=non_empty.index, columns=kept_samples)

    keep_feats = [f for f in non_empty.index if prev.loc[f] >= args.min_prevalence]
    if not keep_feats:
        fail(f"no features passed min_prevalence={args.min_prevalence}")
    clr_test = clr_df.loc[keep_feats]

    x_mat = design.to_numpy(dtype=float)  # n x p
    n, p = x_mat.shape
    df_resid = n - p
    xtx_inv = np.linalg.pinv(x_mat.T @ x_mat)
    term_names = design.columns.tolist()
    exp_indices = [term_names.index(t) for t in exposure_terms]

    y_all = clr_test.to_numpy(dtype=float).T  # n x m
    betas = xtx_inv @ x_mat.T @ y_all         # p x m
    resid = y_all - x_mat @ betas             # n x m
    sigma2 = (resid * resid).sum(axis=0) / float(df_resid)  # (m,)

    rows = []
    for idx_t, t_name in zip(exp_indices, exposure_terms):
        se_vec = np.sqrt(np.maximum(sigma2 * xtx_inv[idx_t, idx_t], 1e-30))
        b_vec = betas[idx_t]
        t_vec = b_vec / se_vec
        p_vec = 2.0 * student_t.sf(np.abs(t_vec), df=df_resid)
        q_vec = bh_fdr(p_vec)
        for j, feat in enumerate(keep_feats):
            rows.append(
                {
                    "feature": feat,
                    "term": t_name,
                    "beta_clr": round(float(b_vec[j]), 4),
                    "se": round(float(se_vec[j]), 4),
                    "t_stat": round(float(t_vec[j]), 4),
                    "p_value": float(p_vec[j]),
                    "q_value": float(q_vec[j]),
                    "prevalence": round(float(prev.loc[feat]), 4),
                    "significant": bool(q_vec[j] <= args.fdr),
                }
            )

    res_df = pd.DataFrame(rows).sort_values(["q_value", "p_value", "feature"]).reset_index(drop=True)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = [args.input, args.metadata]
    params = {
        "exposure": args.exposure,
        "reference": args.reference,
        "covariates": cov_list,
        "min_prevalence": args.min_prevalence,
        "pseudocount": eps,
        "fdr": args.fdr,
    }

    tsv_path = outdir / "multivariable-clr.tsv"
    res_df.to_csv(tsv_path, sep="\t", index=False)
    artifacts = [
        write_artifact(
            tsv_path,
            "mg_multivariable",
            "CLR multivariable regression differential abundance table",
            skill="profiling-workflow",
            inputs=inputs,
            parameters=params,
            shape=res_df.shape,
        )
    ]

    # Volcano plot for primary exposure term
    primary_term = exposure_terms[0]
    sub_plot = res_df.loc[res_df["term"] == primary_term]
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 5))
    sig = sub_plot["significant"].to_numpy()
    y_vals = -np.log10(np.maximum(sub_plot["p_value"].to_numpy(), 1e-300))
    x_vals = sub_plot["beta_clr"].to_numpy()
    ax.scatter(x_vals[~sig], y_vals[~sig], color="#94A3B8", s=24, alpha=0.6, label="NS")
    ax.scatter(x_vals[sig], y_vals[sig], color="#E11D48", s=34, alpha=0.85, label=f"q <= {args.fdr}")
    for _, r in sub_plot.loc[sub_plot["significant"]].head(10).iterrows():
        lbl = str(r["feature"])
        if len(lbl) > 24:
            lbl = lbl[:21] + "..."
        ax.annotate(lbl, (r["beta_clr"], -np.log10(max(r["p_value"], 1e-300))), fontsize=7)
    ax.axvline(0.0, color="#475569", linestyle="--", linewidth=1)
    ax.set_xlabel(f"CLR regression coefficient ({primary_term})")
    ax.set_ylabel("-log10(p-value)")
    cov_note = f" | adj: {', '.join(cov_list)}" if cov_list else ""
    ax.set_title(f"CLR multivariable model{cov_note}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig_path = outdir / "multivariable-volcano.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    artifacts.append(
        write_artifact(
            fig_path,
            "mg_multivariable",
            "CLR multivariable differential abundance volcano plot",
            skill="profiling-workflow",
            inputs=inputs,
            parameters=params,
        )
    )

    n_sig = int(res_df["significant"].sum())
    top_hits = res_df.head(10).to_dict(orient="records")
    emit(
        {
            "tool": "mg_multivariable",
            "exposure": args.exposure,
            "exposure_terms": exposure_terms,
            "covariates": cov_list,
            "n_samples": int(n),
            "df_residual": int(df_resid),
            "n_features_tested": int(len(keep_feats)),
            "n_significant": n_sig,
            "fdr": float(args.fdr),
            "pseudocount": float(eps),
            "top_hits": top_hits,
            "outputs": {
                "table": str(tsv_path),
                "volcano_figure": str(fig_path),
            },
            "scientificContext": {
                str(tsv_path): {
                    "valueSemantics": "de_stats",
                }
            },
            "artifacts": artifacts,
            "caveat": "CLR + OLS linear regression adjusts for compositional closure and listed "
                      "covariates under Gaussian residual assumptions; it is a lightweight Python "
                      "equivalent rather than a zero-inflated hurdle model (MaAsLin 3 / ANCOM-BC).",
        }
    )


if __name__ == "__main__":
    main()
