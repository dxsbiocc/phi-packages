#!/usr/bin/env python3
"""clin_km: Kaplan-Meier product-limit survival curves, log-log Greenwood CIs,
exact log-rank test, landmark rates, and aligned Number-at-Risk subplot.
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
    infer_time_unit,
    join_molecular_features,
    kaplan_meier,
    km_surv_at,
    logrank_test,
    parse_event_series,
    read_table,
    write_artifact,
)

PALETTE = ["#2563EB", "#E11D48", "#10B981", "#F59E0B", "#8B5CF6", "#06B6D4"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--time_col", required=True)
    parser.add_argument("--event_col", required=True)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default=None)
    parser.add_argument("--biomarker_col", default=None)
    parser.add_argument("--matrix", default=None)
    parser.add_argument("--cutoff", default="median")
    parser.add_argument("--landmarks", action="extend", nargs="+", type=float, default=[])
    parser.add_argument("--time_label", default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    if args.group_col and args.biomarker_col:
        fail("specify either --group_col or --biomarker_col, not both")

    os.makedirs(args.output_dir, exist_ok=True)
    clin = read_table(args.clinical)
    sample_col = args.sample_col or clin.columns[0]
    if args.biomarker_col and args.matrix:
        clin = join_molecular_features(clin, sample_col, args.matrix, [args.biomarker_col])

    if args.time_col not in clin.columns or args.event_col not in clin.columns:
        fail(f"time_col '{args.time_col}' or event_col '{args.event_col}' missing from clinical table")

    t_series = pd.to_numeric(clin[args.time_col], errors="coerce")
    e_series, _ = parse_event_series(clin[args.event_col], args.event_col)
    valid = t_series.notna() & (t_series > 0) & e_series.notna()

    warnings: list[str] = []
    n_nonpos = int((t_series.notna() & (t_series <= 0)).sum())
    if n_nonpos > 0:
        warnings.append(f"excluded {n_nonpos} record(s) with non-positive survival time")

    sub = clin.loc[valid].copy()
    sub["_time"] = t_series.loc[valid].to_numpy(dtype=float)
    sub["_event"] = e_series.loc[valid].to_numpy(dtype=int)

    cutoff_info = None
    if args.biomarker_col:
        if args.biomarker_col not in sub.columns:
            fail(f"biomarker_col '{args.biomarker_col}' not found in clinical table or --matrix")
        b_vals = pd.to_numeric(sub[args.biomarker_col], errors="coerce")
        sub = sub.loc[b_vals.notna()].copy()
        bv = b_vals.loc[sub.index].to_numpy(dtype=float)
        if len(bv) < 6:
            fail("fewer than 6 valid observations with non-missing biomarker")
        rule = str(args.cutoff).strip().lower()
        if rule == "median":
            thr = float(np.median(bv))
            sub["_stratum"] = np.where(bv >= thr, f"{args.biomarker_col} High (>= {thr:.3g})", f"{args.biomarker_col} Low (< {thr:.3g})")
            cutoff_info = {"rule": "median", "thresholds": [round(thr, 5)]}
        elif rule == "tertile":
            q1, q2 = np.percentile(bv, [33.3333, 66.6667])
            sub["_stratum"] = np.where(bv < q1, "Low (T1)", np.where(bv < q2, "Mid (T2)", "High (T3)"))
            cutoff_info = {"rule": "tertile", "thresholds": [round(float(q1), 5), round(float(q2), 5)]}
        elif rule == "quartile":
            q25, q75 = np.percentile(bv, [25.0, 75.0])
            sub["_stratum"] = np.where(bv <= q25, "Bottom Quartile (Q1)", np.where(bv >= q75, "Top Quartile (Q4)", "Middle (Q2-Q3)"))
            cutoff_info = {"rule": "quartile", "thresholds": [round(float(q25), 5), round(float(q75), 5)]}
        else:
            try:
                thr = float(rule)
            except ValueError:
                fail(f"unsupported --cutoff '{args.cutoff}'; use median, tertile, quartile, or a numeric value")
            sub["_stratum"] = np.where(bv >= thr, f"{args.biomarker_col} >= {thr:g}", f"{args.biomarker_col} < {thr:g}")
            cutoff_info = {"rule": "explicit", "thresholds": [thr]}
        strat_label = args.biomarker_col
    elif args.group_col:
        if args.group_col not in sub.columns:
            fail(f"group_col '{args.group_col}' not found in clinical table")
        sub = sub[sub[args.group_col].notna() & (sub[args.group_col].astype(str).str.strip() != "")].copy()
        sub["_stratum"] = sub[args.group_col].astype(str).str.strip()
        strat_label = args.group_col
    else:
        sub["_stratum"] = "Overall"
        strat_label = "Overall"

    if len(sub) < 3:
        fail("fewer than 3 valid survival records after filtering")

    levels = sorted(sub["_stratum"].unique().tolist())
    lr = logrank_test(sub["_time"].to_numpy(), sub["_event"].to_numpy(), sub["_stratum"].to_numpy()) if len(levels) >= 2 else {"chisq": None, "df": None, "pvalue": None}

    kms: dict[str, dict] = {}
    summary_rows = []
    curve_rows = []
    for g in levels:
        sg = sub[sub["_stratum"] == g]
        km = kaplan_meier(sg["_time"].to_numpy(), sg["_event"].to_numpy())
        kms[g] = km
        summary_rows.append(
            {
                "group": g,
                "n": km["n"],
                "events": km["events"],
                "median_survival": km["median"],
                "median_ci_low": km["median_ci_low"],
                "median_ci_high": km["median_ci_high"],
            }
        )
        for idx_t in range(len(km["time"])):
            curve_rows.append(
                {
                    "group": g,
                    "time": round(float(km["time"][idx_t]), 5),
                    "survival": round(float(km["surv"][idx_t]), 5),
                    "ci_low": round(float(km["ci_low"][idx_t]), 5),
                    "ci_high": round(float(km["ci_high"][idx_t]), 5),
                    "n_risk": int(km["n_risk"][idx_t]),
                    "n_event": int(km["n_event"][idx_t]),
                    "n_censor": int(km["n_censor"][idx_t]),
                }
            )

    landmark_rows = []
    for hz in sorted(set(args.landmarks)):
        if hz <= 0:
            continue
        for g in levels:
            s_val, lo_val, hi_val, n_r = km_surv_at(kms[g], hz)
            landmark_rows.append(
                {
                    "horizon": hz,
                    "group": g,
                    "survival": round(s_val, 4),
                    "ci_low": round(lo_val, 4),
                    "ci_high": round(hi_val, 4),
                    "n_at_risk": n_r,
                }
            )

    sum_df = pd.DataFrame(summary_rows)
    cur_df = pd.DataFrame(curve_rows)
    sum_path = os.path.join(args.output_dir, "km_summary.tsv")
    cur_path = os.path.join(args.output_dir, "km_curves.tsv")
    sum_df.to_csv(sum_path, sep="\t", index=False)
    cur_df.to_csv(cur_path, sep="\t", index=False)

    # Two-panel KM Curve + Number-at-Risk subplot
    t_max = float(sub["_time"].max())
    unit_lbl = args.time_label or infer_time_unit(sub["_time"].to_numpy(), args.time_col).capitalize()
    fig, (ax_km, ax_nr) = plt.subplots(
        2,
        1,
        figsize=(7.2, 4.8 + 0.35 * len(levels)),
        gridspec_kw={"height_ratios": [3.6, max(1.0, 0.38 * len(levels))]},
        sharex=True,
    )

    for idx_g, g in enumerate(levels):
        color = PALETTE[idx_g % len(PALETTE)]
        km = kms[g]
        med_txt = f"{km['median']:.2f}" if km["median"] is not None else "NR"
        ax_km.step(km["time"], km["surv"], where="post", color=color, linewidth=2.0, label=f"{g} (n={km['n']}, med={med_txt})")
        ax_km.fill_between(km["time"], km["ci_low"], km["ci_high"], step="post", color=color, alpha=0.14)
        # Censoring ticks
        cens_mask = km["n_censor"] > 0
        if np.any(cens_mask):
            ax_km.scatter(km["time"][cens_mask], km["surv"][cens_mask], marker="|", s=45, color=color, linewidths=1.3, zorder=4)

    ax_km.set_ylim(-0.03, 1.05)
    ax_km.set_xlim(0.0, t_max * 1.04)
    ax_km.set_ylabel("Survival probability")
    ax_km.axhline(0.5, color="#94A3B8", linestyle=":", linewidth=0.9)
    title = f"Kaplan-Meier Survival ({strat_label})"
    if lr["pvalue"] is not None:
        title += f"  |  Log-rank p = {lr['pvalue']:.3e}"
    ax_km.set_title(title, fontsize=10.5)
    ax_km.legend(frameon=False, fontsize=8.5, loc="best")
    ax_km.spines["top"].set_visible(False)
    ax_km.spines["right"].set_visible(False)

    # Aligned Number-at-Risk table on shared x-axis ticks
    ticks = np.linspace(0.0, t_max, 6)
    ax_km.set_xticks(ticks)
    ax_nr.set_yticks(np.arange(len(levels)))
    ax_nr.set_yticklabels([g[:22] for g in levels[::-1]], fontsize=8)
    ax_nr.set_ylim(-0.6, len(levels) - 0.4)
    for row_y, g in enumerate(levels[::-1]):
        km = kms[g]
        for tx in ticks:
            _, _, _, n_r = km_surv_at(km, float(tx))
            ax_nr.text(tx, row_y, str(n_r), ha="center", va="center", fontsize=8, color="#0F172A")
    ax_nr.set_xlabel(f"Follow-up Time ({unit_lbl})")
    ax_nr.set_title("Number at risk", loc="left", fontsize=8.5, color="#475569")
    for sp in ("top", "right", "left"):
        ax_nr.spines[sp].set_visible(False)
    ax_nr.tick_params(axis="y", length=0)

    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "km-curve.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {
        "time_col": args.time_col,
        "event_col": args.event_col,
        "stratifier": strat_label,
        "cutoff": cutoff_info,
    }
    in_files = (args.clinical, args.matrix)
    artifacts = [
        write_artifact(fig_path, "clin_km", f"Kaplan-Meier curve & risk table ({strat_label})", inputs=in_files, parameters=params_meta),
        write_artifact(sum_path, "clin_km", "Kaplan-Meier group median survival summary", inputs=in_files, parameters=params_meta, shape=sum_df.shape),
        write_artifact(cur_path, "clin_km", "Kaplan-Meier step coordinates and 95% CIs", inputs=in_files, parameters=params_meta, shape=cur_df.shape),
    ]

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "summary_table": os.path.abspath(sum_path),
            "curves_table": os.path.abspath(cur_path),
            "stratifier": strat_label,
            "cutoff_applied": cutoff_info,
            "n_samples": int(len(sub)),
            "n_events": int(sub["_event"].sum()),
            "logrank_chisq": lr["chisq"],
            "logrank_df": lr["df"],
            "logrank_pvalue": lr["pvalue"],
            "groups": summary_rows,
            "landmark_survival": landmark_rows,
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
