#!/usr/bin/env python3
"""clin_cox: Breslow partial-likelihood univariable & multivariable Cox PH regression,
Harrell's C-index, Schoenfeld PH assumption diagnostics, and HR forest plot.
"""

from __future__ import annotations

import argparse
import os
import sys

import os
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    bh_fdr,
    emit,
    fail,
    fit_cox_breslow,
    join_molecular_features,
    parse_event_series,
    read_table,
    write_artifact,
)


def encode_covariates(
    df: pd.DataFrame, covariates: list[str], standardize_cont: bool
) -> tuple[pd.DataFrame, list[str]]:
    """Encode continuous variables (optionally z-scored) and categorical variables (one-hot vs first reference)."""
    parts = []
    notes = []
    for cov in covariates:
        s = df[cov]
        num = pd.to_numeric(s, errors="coerce")
        is_cont = num.notna().sum() >= 0.85 * max(1, s.notna().sum()) and num.dropna().nunique() > 2
        if is_cont:
            col_s = num.astype(float)
            lbl = cov
            if standardize_cont:
                sd = float(col_s.std(ddof=1))
                if sd > 1e-12:
                    col_s = (col_s - float(col_s.mean())) / sd
                    lbl = f"{cov} (per 1 SD)"
            parts.append(col_s.rename(lbl))
        else:
            cats = s.where(s.notna() & (s.astype(str).str.strip() != ""), np.nan).astype(object)
            valid_levels = sorted(set(str(x).strip() for x in cats.dropna()))
            if len(valid_levels) < 2:
                fail(f"categorical covariate '{cov}' has < 2 non-missing levels")
            ref = valid_levels[0]
            notes.append(f"{cov}: reference='{ref}'")
            for lvl in valid_levels[1:]:
                ind = pd.Series(
                    np.where(cats.isna(), np.nan, (cats.astype(str).str.strip() == lvl).astype(float)),
                    index=df.index,
                    name=f"{cov} [{lvl} vs {ref}]",
                )
                parts.append(ind)
    return pd.concat(parts, axis=1), notes


def render_forest_plot(rows: list[dict], title: str, out_path: str) -> None:
    plot_rows = rows[::-1]
    fig, ax = plt.subplots(figsize=(7.6, max(3.2, 0.48 * len(plot_rows) + 1.4)))
    y_pos = np.arange(len(plot_rows))

    log2_hr = np.array([np.log2(max(r["hr"], 1e-6)) for r in plot_rows])
    log2_lo = np.array([np.log2(max(r["ci_low"], 1e-6)) for r in plot_rows])
    log2_hi = np.array([np.log2(max(r["ci_high"], 1e-6)) for r in plot_rows])

    ax.axvline(0.0, color="#475569", linestyle="--", linewidth=1.0)
    for i, r in enumerate(plot_rows):
        color = "#E11D48" if (r["pvalue"] < 0.05 and r["hr"] > 1.0) else ("#2563EB" if (r["pvalue"] < 0.05 and r["hr"] < 1.0) else "#64748B")
        ax.plot([log2_lo[i], log2_hi[i]], [y_pos[i], y_pos[i]], color=color, linewidth=2.0)
        ax.scatter([log2_hr[i]], [y_pos[i]], color=color, s=48, zorder=4)

    xlim_max = max(1.5, float(np.max(np.abs(np.concatenate([log2_lo, log2_hi])))) * 1.35)
    ax.set_xlim(-xlim_max, xlim_max * 1.55)
    for i, r in enumerate(plot_rows):
        txt = f"HR {r['hr']:.2f} ({r['ci_low']:.2f}-{r['ci_high']:.2f}), p={r['pvalue']:.2e}"
        ax.text(xlim_max * 0.98, y_pos[i], txt, va="center", ha="left", fontsize=8, color="#0F172A")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([r["term"] for r in plot_rows], fontsize=8.5)
    ax.set_xlabel("log2 Hazard Ratio (95% CI)")
    ax.set_title(title, fontsize=10.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--time_col", required=True)
    parser.add_argument("--event_col", required=True)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--covariates", action="extend", nargs="+", default=[], required=True)
    parser.add_argument("--matrix", default=None)
    parser.add_argument("--mode", choices=["both", "multivariable", "univariable"], default="both")
    parser.add_argument(
        "--standardize_continuous",
        nargs="?",
        const=True,
        type=lambda x: str(x).lower() not in ("false", "0", "no"),
        default=False,
    )
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    covs = []
    for item in args.covariates:
        for c in str(item).split(","):
            c = c.strip()
            if c and c not in covs:
                covs.append(c)

    clin = read_table(args.clinical)
    sample_col = args.sample_col or clin.columns[0]
    if args.matrix:
        clin = join_molecular_features(clin, sample_col, args.matrix, covs)
    for c in covs:
        if c not in clin.columns:
            fail(f"covariate '{c}' not found in clinical table or --matrix")

    t_s = pd.to_numeric(clin[args.time_col], errors="coerce")
    e_s, _ = parse_event_series(clin[args.event_col], args.event_col)
    base_valid = t_s.notna() & (t_s > 0) & e_s.notna()
    sub = clin.loc[base_valid].copy()
    sub["_time"] = t_s.loc[base_valid].to_numpy(dtype=float)
    sub["_event"] = e_s.loc[base_valid].to_numpy(dtype=int)

    warnings: list[str] = []
    artifacts: list[str] = []
    figures: list[str] = []
    in_files = (args.clinical, args.matrix)
    params_meta = {"mode": args.mode, "covariates": covs, "standardize_continuous": args.standardize_continuous}

    uni_summary = None
    if args.mode in ("both", "univariable"):
        uni_rows = []
        for cov in covs:
            enc_df, _ = encode_covariates(sub, [cov], args.standardize_continuous)
            valid_c = enc_df.notna().all(axis=1)
            if int(valid_c.sum()) < 5 or int(sub.loc[valid_c, "_event"].sum()) < 2:
                continue
            fit = fit_cox_breslow(
                enc_df.loc[valid_c].to_numpy(dtype=float),
                sub.loc[valid_c, "_time"].to_numpy(dtype=float),
                sub.loc[valid_c, "_event"].to_numpy(dtype=int),
            )
            for idx_t, term_name in enumerate(enc_df.columns):
                uni_rows.append(
                    {
                        "variable": cov,
                        "term": str(term_name),
                        "n": fit["n_samples"],
                        "events": fit["n_events"],
                        "beta": round(float(fit["beta"][idx_t]), 5),
                        "se": round(float(fit["se"][idx_t]), 5),
                        "hr": round(float(fit["hr"][idx_t]), 4),
                        "ci_low": round(float(fit["ci_low"][idx_t]), 4),
                        "ci_high": round(float(fit["ci_high"][idx_t]), 4),
                        "z": round(float(fit["z"][idx_t]), 4),
                        "pvalue": float(fit["pvalue"][idx_t]),
                        "c_index": fit["c_index"],
                        "ph_pvalue": fit["ph_pvalues"][idx_t],
                    }
                )
        if uni_rows:
            uni_df = pd.DataFrame(uni_rows)
            uni_df["padj"] = bh_fdr(uni_df["pvalue"].to_numpy(dtype=float))
            uni_df = uni_df.sort_values(["pvalue", "term"]).reset_index(drop=True)
            uni_path = os.path.join(args.output_dir, "cox_univariable.tsv")
            uni_df.to_csv(uni_path, sep="\t", index=False)
            artifacts.append(
                write_artifact(uni_path, "clin_cox", "Univariable Cox PH screening table", inputs=in_files, parameters=params_meta, shape=uni_df.shape)
            )
            uni_summary = {
                "table": os.path.abspath(uni_path),
                "n_terms_tested": int(len(uni_df)),
                "n_significant_fdr05": int((uni_df["padj"] <= 0.05).sum()),
                "terms": uni_df.to_dict(orient="records"),
            }

    multi_summary = None
    if args.mode in ("both", "multivariable"):
        enc_all, ref_notes = encode_covariates(sub, covs, args.standardize_continuous)
        valid_m = enc_all.notna().all(axis=1)
        n_dropped_na = int((~valid_m).sum())
        if n_dropped_na > 0:
            warnings.append(f"complete-case multivariable Cox excluded {n_dropped_na} sample(s) with missing covariates")
        if int(valid_m.sum()) < max(6, enc_all.shape[1] + 2):
            fail("insufficient complete-case samples for multivariable Cox regression")

        fit_m = fit_cox_breslow(
            enc_all.loc[valid_m].to_numpy(dtype=float),
            sub.loc[valid_m, "_time"].to_numpy(dtype=float),
            sub.loc[valid_m, "_event"].to_numpy(dtype=int),
        )
        if fit_m["events_per_parameter"] < 10.0:
            warnings.append(
                f"EPV warning: {fit_m['n_events']} events across {fit_m['n_parameters']} parameters "
                f"(EPV = {fit_m['events_per_parameter']} < 10); hazard ratio CIs may be unstable — "
                "consider fewer covariates or L2-penalized clin_time_roc"
            )
        if not fit_m["converged"]:
            warnings.append("Newton-Raphson reached max_iter before strict convergence; check quasi-complete separation")

        m_rows = []
        ph_violated = []
        for idx_t, term_name in enumerate(enc_all.columns):
            php = fit_m["ph_pvalues"][idx_t]
            if php is not None and php < 0.05:
                ph_violated.append(str(term_name))
            m_rows.append(
                {
                    "term": str(term_name),
                    "beta": round(float(fit_m["beta"][idx_t]), 5),
                    "se": round(float(fit_m["se"][idx_t]), 5),
                    "hr": round(float(fit_m["hr"][idx_t]), 4),
                    "ci_low": round(float(fit_m["ci_low"][idx_t]), 4),
                    "ci_high": round(float(fit_m["ci_high"][idx_t]), 4),
                    "z": round(float(fit_m["z"][idx_t]), 4),
                    "pvalue": float(fit_m["pvalue"][idx_t]),
                    "ph_pvalue": round(float(php), 4) if php is not None else None,
                }
            )
        if ph_violated:
            warnings.append(
                f"Schoenfeld PH assumption test p < 0.05 for {ph_violated}: hazard ratio varies over follow-up time"
            )

        m_df = pd.DataFrame(m_rows)
        m_path = os.path.join(args.output_dir, "cox_multivariable.tsv")
        m_df.to_csv(m_path, sep="\t", index=False)

        forest_path = os.path.join(args.output_dir, "cox-forest.png")
        render_forest_plot(
            m_rows,
            f"Multivariable Cox PH (n={fit_m['n_samples']}, events={fit_m['n_events']}, C-index={fit_m['c_index']:.3f})",
            forest_path,
        )
        figures.append(forest_path)
        artifacts.insert(
            0,
            write_artifact(forest_path, "clin_cox", "Multivariable Cox Hazard Ratio forest plot", inputs=in_files, parameters=params_meta),
        )
        artifacts.append(
            write_artifact(m_path, "clin_cox", "Multivariable Cox PH regression table", inputs=in_files, parameters=params_meta, shape=m_df.shape)
        )

        multi_summary = {
            "table": os.path.abspath(m_path),
            "n_samples": fit_m["n_samples"],
            "n_events": fit_m["n_events"],
            "n_parameters": fit_m["n_parameters"],
            "events_per_parameter": fit_m["events_per_parameter"],
            "c_index": fit_m["c_index"],
            "converged": fit_m["converged"],
            "global_ph_pvalue": round(float(fit_m["global_ph_pvalue"]), 4) if fit_m["global_ph_pvalue"] is not None else None,
            "ph_violations": ph_violated,
            "reference_levels": ref_notes,
            "terms": m_rows,
        }
    elif uni_summary and uni_summary["terms"]:
        forest_path = os.path.join(args.output_dir, "cox-forest.png")
        render_forest_plot(uni_summary["terms"][:20], "Univariable Cox PH Forest Plot", forest_path)
        figures.append(forest_path)
        artifacts.insert(
            0,
            write_artifact(forest_path, "clin_cox", "Univariable Cox Hazard Ratio forest plot", inputs=in_files, parameters=params_meta),
        )

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "mode": args.mode,
            "univariable": uni_summary,
            "multivariable": multi_summary,
            "figures": [os.path.abspath(f) for f in figures],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
