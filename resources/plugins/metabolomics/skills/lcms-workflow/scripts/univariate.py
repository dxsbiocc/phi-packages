#!/usr/bin/env python3
"""mb_univariate: per-feature group tests + BH FDR + volcano plot."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from _common import emit, fail, read_matrix, read_tsv, setup_matplotlib, write_artifact

SKILL = "lcms-workflow"


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--group_col", default="group")
    ap.add_argument("--case")
    ap.add_argument("--control")
    ap.add_argument("--test", choices=["welch", "wilcoxon", "kruskal"],
                    default="wilcoxon")
    ap.add_argument("--logged", action="store_true")
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--fc_threshold", type=float, default=1.0)
    args = ap.parse_args()

    meta, X = read_matrix(args.input)
    md = read_tsv(args.metadata)
    md = md.set_index(md.columns[0])
    md.index = md.index.astype(str)
    if args.group_col not in md.columns:
        fail(f"group column not in metadata: {args.group_col}")
    missing = [c for c in X.columns if str(c) not in md.index]
    if missing:
        fail(f"samples missing from metadata: {missing}")
    labels = pd.Series([str(md.loc[str(c), args.group_col]) for c in X.columns],
                       index=X.columns)
    two_group = args.test in ("welch", "wilcoxon")
    if two_group:
        if not args.case or not args.control:
            fail(f"test={args.test} requires --case and --control")
        if args.case == args.control:
            fail("case and control must differ")
        keep = labels.isin([args.case, args.control])
        if not keep.any():
            fail("no samples belong to the case/control groups")
        X = X.loc[:, keep]
        labels = labels[keep]
    levels = sorted(labels.unique())
    if len(levels) < 2:
        fail("need at least two groups")

    rows = []
    for feat in X.index:
        vals = X.loc[feat]
        groups = [vals[labels == lev].dropna().to_numpy(dtype=float) for lev in levels]
        if any(len(g) < 2 for g in groups):
            rows.append({"feature": feat, "statistic": np.nan, "pvalue": np.nan,
                         "note": "fewer than 2 observed values in a group"})
            continue
        if two_group:
            case = vals[labels == args.case].dropna().to_numpy(dtype=float)
            control = vals[labels == args.control].dropna().to_numpy(dtype=float)
            if np.all(case == case[0]) and np.all(control == control[0]):
                stat, p = 0.0, 1.0
            elif args.test == "welch":
                stat, p = stats.ttest_ind(case, control, equal_var=False)
            else:
                stat, p = stats.mannwhitneyu(case, control, alternative="two-sided")
            mc, m0 = float(np.mean(case)), float(np.mean(control))
            if args.logged:
                log2fc = mc - m0
                fc = 2.0 ** log2fc
            elif m0 > 0:
                fc = mc / m0
                log2fc = float(np.log2(fc)) if fc > 0 else np.nan
            elif mc > 0:
                fc, log2fc = np.inf, np.inf
            else:
                fc, log2fc = 1.0, 0.0
            rows.append({"feature": feat, "statistic": float(stat),
                         "pvalue": float(p), "mean_case": mc, "mean_control": m0,
                         "fold_change": float(fc), "log2_fold_change": float(log2fc)})
        else:
            if all(np.all(g == g[0]) for g in groups):
                stat, p = 0.0, 1.0
            else:
                stat, p = stats.kruskal(*groups)
            row = {"feature": feat, "statistic": float(stat), "pvalue": float(p)}
            for lev, g in zip(levels, groups):
                row[f"mean_{lev}"] = float(np.mean(g))
            rows.append(row)

    table = pd.DataFrame(rows).set_index("feature")
    tested = table["pvalue"].notna()
    table["fdr"] = np.nan
    if tested.any():
        table.loc[tested, "fdr"] = bh_fdr(table.loc[tested, "pvalue"].to_numpy())
    if two_group:
        table["significant"] = (
            (table["fdr"] < args.fdr)
            & (table["log2_fold_change"].abs() >= args.fc_threshold))
    else:
        table["significant"] = table["fdr"] < args.fdr
    table = table.sort_values("pvalue", na_position="last")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    table_path = out / "univariate-table.tsv"
    table.to_csv(table_path, sep="\t")

    artifacts = [
        write_artifact(table_path, "mb_univariate", "Univariate test table",
                       skill=SKILL, inputs=[args.input, args.metadata],
                       parameters={"test": args.test, "logged": args.logged,
                                   "fdr": args.fdr, "fc_threshold": args.fc_threshold},
                       shape=table.shape),
    ]
    fig_path = None
    if two_group:
        plt = setup_matplotlib()
        fig, ax = plt.subplots(figsize=(6, 5))
        x = table["log2_fold_change"].to_numpy(dtype=float)
        pv = table["pvalue"].to_numpy(dtype=float)
        with np.errstate(divide="ignore"):
            y = -np.log10(np.where(pv > 0, pv, np.nan))
        sig = table["significant"].to_numpy()
        ok = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[ok & ~sig], y[ok & ~sig], s=12, color="#9aa0a6", label="ns")
        ax.scatter(x[ok & sig], y[ok & sig], s=16, color="#d62728", label="significant")
        ax.axvline(args.fc_threshold, color="grey", ls="--", lw=0.8)
        ax.axvline(-args.fc_threshold, color="grey", ls="--", lw=0.8)
        ax.axhline(-np.log10(args.fdr), color="grey", ls=":", lw=0.8)
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10 p")
        ax.set_title(f"{args.case} vs {args.control} ({args.test})")
        ax.legend()
        fig.tight_layout()
        fig_path = out / "volcano.png"
        fig.savefig(fig_path, dpi=150)
        plt.close(fig)
        artifacts.append(
            write_artifact(fig_path, "mb_univariate", "Volcano plot",
                           skill=SKILL, inputs=[args.input, args.metadata],
                           parameters={"test": args.test, "fdr": args.fdr,
                                       "fc_threshold": args.fc_threshold}))

    top = table.head(10)
    emit({
        "input": args.input,
        "test": args.test,
        "case": args.case,
        "control": args.control,
        "groups": {lev: int((labels == lev).sum()) for lev in levels},
        "logged": args.logged,
        "n_features": int(table.shape[0]),
        "n_tested": int(tested.sum()),
        "n_significant": int(table["significant"].sum()),
        "fdr_threshold": args.fdr,
        "fc_threshold": args.fc_threshold,
        "top": [{"feature": str(i),
                 "pvalue": None if pd.isna(r["pvalue"]) else float(r["pvalue"]),
                 "fdr": None if pd.isna(r["fdr"]) else float(r["fdr"]),
                 "log2_fold_change": (None if "log2_fold_change" not in r
                                      or pd.isna(r["log2_fold_change"])
                                      else float(r["log2_fold_change"]))}
                for i, r in top.iterrows()],
        "output_dir": str(out),
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()