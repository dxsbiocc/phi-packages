#!/usr/bin/env python3
"""clin_table1: stratified Baseline Characteristics Table (Table 1) with Shapiro-Wilk
normality routing, Welch/ANOVA/Mann-Whitney/Kruskal-Wallis, Fisher/chi2, and SMD.
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
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, read_table, write_artifact


def fmt_cont(arr: np.ndarray, normal: bool) -> str:
    if len(arr) == 0:
        return "NA"
    if normal:
        sd = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
        return f"{float(np.mean(arr)):.2f} ± {sd:.2f}"
    q25, med, q75 = np.percentile(arr, [25, 50, 75])
    return f"{float(med):.2f} [{float(q25):.2f}, {float(q75):.2f}]"


def is_normal_group(arr: np.ndarray) -> bool:
    if len(arr) < 3 or len(arr) > 5000:
        return len(arr) > 30
    if float(np.std(arr)) < 1e-12:
        return False
    _, p = stats.shapiro(arr)
    return bool(p >= 0.05)


def compute_cont_smd(g1: np.ndarray, g2: np.ndarray) -> float | None:
    if len(g1) < 2 or len(g2) < 2:
        return None
    v1 = float(np.var(g1, ddof=1))
    v2 = float(np.var(g2, ddof=1))
    denom = np.sqrt(0.5 * (v1 + v2))
    if denom < 1e-12:
        return 0.0
    return round(float(abs(np.mean(g1) - np.mean(g2)) / denom), 4)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--variables", action="extend", nargs="+", default=[], required=True)
    parser.add_argument("--group_col", default=None)
    parser.add_argument("--categorical_cols", action="extend", nargs="+", default=[])
    parser.add_argument("--continuous_cols", action="extend", nargs="+", default=[])
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.clinical)

    var_list = []
    for item in args.variables:
        for v in str(item).split(","):
            v = v.strip()
            if v and v not in var_list:
                var_list.append(v)
    for v in var_list:
        if v not in df.columns:
            fail(f"variable '{v}' not found in clinical table columns: {list(df.columns)}")

    cat_force = {c.strip() for item in args.categorical_cols for c in str(item).split(",") if c.strip()}
    cont_force = {c.strip() for item in args.continuous_cols for c in str(item).split(",") if c.strip()}

    if args.group_col:
        if args.group_col not in df.columns:
            fail(f"group_col '{args.group_col}' not found in clinical columns")
        df = df[df[args.group_col].notna() & (df[args.group_col].astype(str).str.strip() != "")].copy()
        df["_grp"] = df[args.group_col].astype(str).str.strip()
        groups = sorted(df["_grp"].unique().tolist())
        if len(groups) < 2:
            fail(f"group_col '{args.group_col}' has < 2 non-empty groups")
    else:
        groups = []

    n_samples = int(len(df))
    group_counts = {g: int((df["_grp"] == g).sum()) for g in groups} if groups else {"Overall": n_samples}

    tsv_rows: list[dict[str, Any]] = []
    var_summaries: list[dict[str, Any]] = []

    # Header row with sample counts
    header_entry: dict[str, Any] = {"Variable": "N", "Level_or_Metric": "", "Overall": str(n_samples)}
    for g in groups:
        header_entry[g] = str(group_counts[g])
    header_entry["Test"] = ""
    header_entry["P_value"] = ""
    header_entry["SMD"] = ""
    tsv_rows.append(header_entry)

    for var in var_list:
        if var == args.group_col:
            continue
        s = df[var]
        n_miss = int(s.isna().sum() + (s.astype(str).str.strip() == "").sum() - s.isna().sum())
        num = pd.to_numeric(s, errors="coerce")
        is_cont = False
        if var in cont_force:
            is_cont = True
        elif var in cat_force:
            is_cont = False
        else:
            valid_num = num.dropna()
            if len(valid_num) >= 0.8 * max(1, s.notna().sum()) and valid_num.nunique() > 5:
                is_cont = True

        if is_cont:
            overall_arr = num.dropna().to_numpy(dtype=float)
            grp_arrays = [num[df["_grp"] == g].dropna().to_numpy(dtype=float) for g in groups] if groups else [overall_arr]
            normal = all(is_normal_group(a) for a in grp_arrays if len(a) >= 3)
            metric_lbl = "Mean ± SD" if normal else "Median [IQR]"

            pval = None
            test_name = "descriptive"
            smd_val = None
            if len(groups) == 2 and all(len(a) >= 2 for a in grp_arrays):
                smd_val = compute_cont_smd(grp_arrays[0], grp_arrays[1])
                if normal:
                    _, pval = stats.ttest_ind(grp_arrays[0], grp_arrays[1], equal_var=False)
                    test_name = "Welch t-test"
                else:
                    _, pval = stats.mannwhitneyu(grp_arrays[0], grp_arrays[1], alternative="two-sided")
                    test_name = "Mann-Whitney U"
            elif len(groups) > 2 and sum(len(a) >= 2 for a in grp_arrays) >= 2:
                nonempty = [a for a in grp_arrays if len(a) >= 2]
                if normal:
                    _, pval = stats.f_oneway(*nonempty)
                    test_name = "One-way ANOVA"
                else:
                    _, pval = stats.kruskal(*nonempty)
                    test_name = "Kruskal-Wallis"

            row_d: dict[str, Any] = {
                "Variable": var,
                "Level_or_Metric": metric_lbl,
                "Overall": fmt_cont(overall_arr, normal),
            }
            by_grp = {}
            for g, a in zip(groups, grp_arrays):
                txt = fmt_cont(a, normal)
                row_d[g] = txt
                by_grp[g] = txt
            row_d["Test"] = test_name if groups else ""
            row_d["P_value"] = f"{float(pval):.4g}" if pval is not None and np.isfinite(pval) else ""
            row_d["SMD"] = f"{smd_val:.3f}" if smd_val is not None else ""
            tsv_rows.append(row_d)

            var_summaries.append(
                {
                    "variable": var,
                    "type": "continuous",
                    "normality_assumed": normal,
                    "metric": metric_lbl,
                    "missing": int(num.isna().sum()),
                    "overall": row_d["Overall"],
                    "summary_by_group": by_grp,
                    "test": test_name if groups else None,
                    "pvalue": float(pval) if pval is not None and np.isfinite(pval) else None,
                    "smd": smd_val,
                }
            )
        else:
            clean_cat = s.where(s.notna(), np.nan)
            valid_mask = clean_cat.notna() & (clean_cat.astype(str).str.strip() != "")
            cats = clean_cat[valid_mask].astype(str).str.strip()
            levels = sorted(cats.unique().tolist())

            pval = None
            test_name = "descriptive"
            smd_val = None
            if groups and len(levels) >= 2:
                ct = pd.crosstab(cats, df.loc[valid_mask, "_grp"]).reindex(index=levels, columns=groups, fill_value=0)
                obs = ct.to_numpy(dtype=int)
                if obs.sum() > 0 and (obs.sum(axis=0) > 0).sum() >= 2 and (obs.sum(axis=1) > 0).sum() >= 2:
                    chi2, p_chi, _, expected = stats.chi2_contingency(obs, correction=False)
                    if obs.shape == (2, 2) and np.any(expected < 5.0):
                        _, pval = stats.fisher_exact(obs)
                        test_name = "Fisher's exact"
                    else:
                        pval = float(p_chi)
                        test_name = "Chi-square"
                    if len(groups) == 2 and len(levels) == 2:
                        p1 = obs[1, 0] / max(1, obs[:, 0].sum())
                        p2 = obs[1, 1] / max(1, obs[:, 1].sum())
                        denom = np.sqrt(0.5 * (p1 * (1 - p1) + p2 * (1 - p2)))
                        smd_val = round(float(abs(p1 - p2) / denom), 4) if denom > 1e-12 else 0.0

            level_summaries: dict[str, Any] = {}
            n_val_all = max(1, int(valid_mask.sum()))
            for idx_l, lvl in enumerate(levels):
                n_l = int((cats == lvl).sum())
                ov_txt = f"{n_l} ({100.0 * n_l / n_val_all:.1f}%)"
                row_d = {
                    "Variable": var if idx_l == 0 else "",
                    "Level_or_Metric": f"{lvl}, n (%)",
                    "Overall": ov_txt,
                }
                lvl_grp = {}
                for g in groups:
                    sub_g = clean_cat[valid_mask & (df["_grp"] == g)].astype(str).str.strip()
                    n_lg = int((sub_g == lvl).sum())
                    denom_g = max(1, int(len(sub_g)))
                    g_txt = f"{n_lg} ({100.0 * n_lg / denom_g:.1f}%)"
                    row_d[g] = g_txt
                    lvl_grp[g] = g_txt
                row_d["Test"] = (test_name if groups else "") if idx_l == 0 else ""
                row_d["P_value"] = (f"{float(pval):.4g}" if pval is not None and np.isfinite(pval) else "") if idx_l == 0 else ""
                row_d["SMD"] = (f"{smd_val:.3f}" if smd_val is not None else "") if idx_l == 0 else ""
                tsv_rows.append(row_d)
                level_summaries[lvl] = {"overall": ov_txt, "by_group": lvl_grp}

            var_summaries.append(
                {
                    "variable": var,
                    "type": "categorical",
                    "levels": levels,
                    "missing": n_samples - int(valid_mask.sum()),
                    "summary_by_level": level_summaries,
                    "test": test_name if groups else None,
                    "pvalue": float(pval) if pval is not None and np.isfinite(pval) else None,
                    "smd": smd_val,
                }
            )

    t1_df = pd.DataFrame(tsv_rows)
    t1_path = os.path.join(args.output_dir, "table1.tsv")
    t1_df.to_csv(t1_path, sep="\t", index=False)

    # Overview plot: SMD / -log10(p) balance diagnostic if stratified, else variable completeness
    fig, ax = plt.subplots(figsize=(6.4, max(3.6, len(var_summaries) * 0.42)))
    if groups and any(v["pvalue"] is not None for v in var_summaries):
        names = [v["variable"] for v in var_summaries[::-1]]
        nlp = [
            -np.log10(max(float(v["pvalue"]), 1e-20)) if v["pvalue"] is not None else 0.0
            for v in var_summaries[::-1]
        ]
        colors = ["#F43F5E" if (v["pvalue"] is not None and v["pvalue"] < 0.05) else "#38BDF8" for v in var_summaries[::-1]]
        ax.barh(names, nlp, color=colors)
        ax.axvline(-np.log10(0.05), color="#475569", linestyle="--", linewidth=1.0, label="p = 0.05")
        ax.set_xlabel("-log10(p-value) across strata")
        ax.set_title(f"Table 1 Baseline Balance ({args.group_col})")
        ax.legend(frameon=False)
    else:
        names = [v["variable"] for v in var_summaries[::-1]]
        comp = [100.0 * (1.0 - v["missing"] / max(1, n_samples)) for v in var_summaries[::-1]]
        ax.barh(names, comp, color="#38BDF8")
        ax.set_xlim(0, 105)
        ax.set_xlabel("Completeness (%)")
        ax.set_title("Table 1 Baseline Variable Completeness")
    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "table1-overview.png")
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)

    params_meta = {"group_col": args.group_col, "variables": var_list}
    artifacts = [
        write_artifact(t1_path, "clin_table1", f"Baseline Characteristics Table 1 (N={n_samples})", inputs=(args.clinical,), parameters=params_meta, shape=t1_df.shape),
        write_artifact(fig_path, "clin_table1", "Table 1 baseline balance / completeness overview", inputs=(args.clinical,), parameters=params_meta),
    ]

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "table": os.path.abspath(t1_path),
            "group_col": args.group_col,
            "groups": group_counts,
            "n_samples": n_samples,
            "variables": var_summaries,
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
        }
    )


if __name__ == "__main__":
    main()
