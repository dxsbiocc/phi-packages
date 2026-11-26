#!/usr/bin/env python3
"""clin_meta: Study-level Fixed-Effect and Random-Effects (DerSimonian-Laird,
Paule-Mandel, Hartung-Knapp-Sidik-Jonkman) meta-analysis with Q/I^2/tau^2,
95% prediction interval, Egger's asymmetry test, subgroup Q_between, Leave-One-Out
sensitivity, and publication-grade forest + funnel plots.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pred_common import emit, fail, read_table, write_artifact


def _fmt_est(val: float, ratio: bool) -> float:
    if ratio:
        return round(float(math.exp(np.clip(val, -20.0, 20.0))), 4)
    return round(float(val), 5)


def pool_studies(y: np.ndarray, se: np.ndarray, ratio: bool) -> dict[str, Any]:
    k = len(y)
    v = se ** 2
    w_fe = 1.0 / np.maximum(v, 1e-18)
    sum_w = float(np.sum(w_fe))
    mu_fe = float(np.sum(w_fe * y) / sum_w)
    se_fe = 1.0 / math.sqrt(sum_w)
    z_fe = mu_fe / se_fe
    p_fe = float(2.0 * stats.norm.sf(abs(z_fe)))

    q_stat = float(np.sum(w_fe * ((y - mu_fe) ** 2)))
    df_q = max(0, k - 1)
    p_q = float(stats.chi2.sf(q_stat, df=df_q)) if df_q > 0 else 1.0
    i2 = max(0.0, (q_stat - df_q) / q_stat) * 100.0 if (q_stat > 0 and df_q > 0) else 0.0

    # DerSimonian-Laird tau^2
    c_dl = sum_w - float(np.sum(w_fe ** 2)) / sum_w
    tau2_dl = max(0.0, (q_stat - df_q) / c_dl) if (df_q > 0 and c_dl > 1e-18) else 0.0

    # Paule-Mandel tau^2
    def pm_eq(t2: float) -> float:
        w_t = 1.0 / np.maximum(v + t2, 1e-18)
        mu_t = float(np.sum(w_t * y) / np.sum(w_t))
        return float(np.sum(w_t * ((y - mu_t) ** 2))) - float(df_q)

    if df_q == 0 or pm_eq(0.0) <= 0.0:
        tau2_pm = 0.0
    else:
        hi = max(1.0, tau2_dl * 10.0 + 1.0)
        while pm_eq(hi) > 0.0 and hi < 1e6:
            hi *= 4.0
        tau2_pm = float(optimize.brentq(pm_eq, 0.0, hi)) if pm_eq(hi) <= 0.0 else tau2_dl

    def _re_stats(tau2: float) -> tuple[float, float, float, float, np.ndarray]:
        w_re = 1.0 / np.maximum(v + tau2, 1e-18)
        sw = float(np.sum(w_re))
        mu_r = float(np.sum(w_re * y) / sw)
        se_r = 1.0 / math.sqrt(sw)
        z_r = mu_r / se_r
        p_r = float(2.0 * stats.norm.sf(abs(z_r)))
        return mu_r, se_r, p_r, sw, w_re

    mu_dl, se_dl, p_dl, sw_dl, w_dl = _re_stats(tau2_dl)
    mu_pm, se_pm, p_pm, _, _ = _re_stats(tau2_pm)

    # HKSJ (modified Knapp-Hartung with q floored at 1.0 so CI isn't artificially narrow)
    if k >= 2:
        q_hksj = float(np.sum(w_dl * ((y - mu_dl) ** 2)) / float(k - 1))
        se_hksj = math.sqrt(max(1.0, q_hksj) / sw_dl)
        t_crit = float(stats.t.ppf(0.975, df=k - 1))
        t_stat = mu_dl / se_hksj
        p_hksj = float(2.0 * stats.t.sf(abs(t_stat), df=k - 1))
        hksj_lo = mu_dl - t_crit * se_hksj
        hksj_hi = mu_dl + t_crit * se_hksj
    else:
        se_hksj, p_hksj = se_dl, p_dl
        hksj_lo = mu_dl - 1.959964 * se_dl
        hksj_hi = mu_dl + 1.959964 * se_dl

    # 95% Prediction Interval (k >= 3)
    pi95 = None
    if k >= 3:
        t_pi = float(stats.t.ppf(0.975, df=k - 2))
        pi_half = t_pi * math.sqrt(tau2_dl + se_dl ** 2)
        pi95 = {
            "ci_low": _fmt_est(mu_dl - pi_half, ratio),
            "ci_high": _fmt_est(mu_dl + pi_half, ratio),
            "log_ci_low": round(mu_dl - pi_half, 5),
            "log_ci_high": round(mu_dl + pi_half, 5),
        }

    return {
        "fixed_effect": {
            "estimate": _fmt_est(mu_fe, ratio),
            "ci_low": _fmt_est(mu_fe - 1.959964 * se_fe, ratio),
            "ci_high": _fmt_est(mu_fe + 1.959964 * se_fe, ratio),
            "log_estimate": round(mu_fe, 5),
            "se": round(se_fe, 5),
            "pvalue": p_fe,
        },
        "random_effects_dl": {
            "estimate": _fmt_est(mu_dl, ratio),
            "ci_low": _fmt_est(mu_dl - 1.959964 * se_dl, ratio),
            "ci_high": _fmt_est(mu_dl + 1.959964 * se_dl, ratio),
            "log_estimate": round(mu_dl, 5),
            "se": round(se_dl, 5),
            "pvalue": p_dl,
        },
        "random_effects_pm": {
            "estimate": _fmt_est(mu_pm, ratio),
            "ci_low": _fmt_est(mu_pm - 1.959964 * se_pm, ratio),
            "ci_high": _fmt_est(mu_pm + 1.959964 * se_pm, ratio),
            "log_estimate": round(mu_pm, 5),
            "se": round(se_pm, 5),
            "tau2": round(tau2_pm, 5),
            "pvalue": p_pm,
        },
        "hksj": {
            "estimate": _fmt_est(mu_dl, ratio),
            "ci_low": _fmt_est(hksj_lo, ratio),
            "ci_high": _fmt_est(hksj_hi, ratio),
            "se": round(se_hksj, 5),
            "pvalue": p_hksj,
        },
        "prediction_interval_95": pi95,
        "heterogeneity": {
            "q": round(q_stat, 4),
            "df": df_q,
            "q_pvalue": p_q,
            "i2_percent": round(i2, 2),
            "tau2_dl": round(tau2_dl, 5),
            "tau2_pm": round(tau2_pm, 5),
        },
        "weights_fe": w_fe / sum_w,
        "weights_re": w_dl / sw_dl,
        "_mu_dl": mu_dl,
        "_se_dl": se_dl,
        "_mu_fe": mu_fe,
        "_se_fe": se_fe,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--studies", required=True)
    parser.add_argument("--study_col", required=True)
    parser.add_argument("--effect_col", default=None)
    parser.add_argument("--se_col", default=None)
    parser.add_argument("--ci_low_col", default=None)
    parser.add_argument("--ci_high_col", default=None)
    parser.add_argument("--is_ratio_scale", action="store_true")
    parser.add_argument("--ai_col", default=None)
    parser.add_argument("--bi_col", default=None)
    parser.add_argument("--ci_col", default=None)
    parser.add_argument("--di_col", default=None)
    parser.add_argument("--measure", choices=["OR", "RR", "RD"], default="OR")
    parser.add_argument("--subgroup_col", default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.studies)
    if args.study_col not in df.columns:
        fail(f"study_col '{args.study_col}' not found in columns: {list(df.columns)}")

    has_2x2 = all(c is not None for c in (args.ai_col, args.bi_col, args.ci_col, args.di_col))
    has_eff = args.effect_col is not None
    if has_2x2 == has_eff:
        fail("specify either --effect_col (+ --se_col or --ci_low_col/--ci_high_col) OR all four 2x2 columns (ai/bi/ci/di)")

    warnings: list[str] = []
    if has_eff:
        for col in (args.effect_col, args.se_col, args.ci_low_col, args.ci_high_col):
            if col is not None and col not in df.columns:
                fail(f"column '{col}' not found in studies table")
        if args.se_col is None and not (args.ci_low_col and args.ci_high_col):
            fail("provide --se_col or both --ci_low_col and --ci_high_col")
        eff_raw = pd.to_numeric(df[args.effect_col], errors="coerce")
        if args.se_col:
            se_raw = pd.to_numeric(df[args.se_col], errors="coerce")
            valid = df[args.study_col].notna() & eff_raw.notna() & se_raw.notna() & (se_raw > 0)
            if args.is_ratio_scale:
                valid = valid & (eff_raw > 0)
            sub = df.loc[valid].copy()
            y = np.log(eff_raw.loc[valid].to_numpy(dtype=float, copy=True)) if args.is_ratio_scale else eff_raw.loc[valid].to_numpy(dtype=float, copy=True)
            se = se_raw.loc[valid].to_numpy(dtype=float, copy=True)
        else:
            lo_raw = pd.to_numeric(df[args.ci_low_col], errors="coerce")
            hi_raw = pd.to_numeric(df[args.ci_high_col], errors="coerce")
            valid = df[args.study_col].notna() & eff_raw.notna() & lo_raw.notna() & hi_raw.notna() & (hi_raw > lo_raw)
            if args.is_ratio_scale:
                valid = valid & (eff_raw > 0) & (lo_raw > 0) & (hi_raw > 0)
            sub = df.loc[valid].copy()
            if args.is_ratio_scale:
                y = np.log(eff_raw.loc[valid].to_numpy(dtype=float, copy=True))
                se = (np.log(hi_raw.loc[valid].to_numpy(dtype=float, copy=True)) - np.log(lo_raw.loc[valid].to_numpy(dtype=float, copy=True))) / (2.0 * 1.959964)
            else:
                y = eff_raw.loc[valid].to_numpy(dtype=float, copy=True)
                se = (hi_raw.loc[valid].to_numpy(dtype=float, copy=True) - lo_raw.loc[valid].to_numpy(dtype=float, copy=True)) / (2.0 * 1.959964)
        is_ratio = bool(args.is_ratio_scale)
        scale_label = "ratio" if is_ratio else "difference"
    else:
        for col in (args.ai_col, args.bi_col, args.ci_col, args.di_col):
            if col not in df.columns:
                fail(f"column '{col}' not found in studies table")
        ai_s = pd.to_numeric(df[args.ai_col], errors="coerce")
        bi_s = pd.to_numeric(df[args.bi_col], errors="coerce")
        ci_s = pd.to_numeric(df[args.ci_col], errors="coerce")
        di_s = pd.to_numeric(df[args.di_col], errors="coerce")
        valid = df[args.study_col].notna() & ai_s.notna() & bi_s.notna() & ci_s.notna() & di_s.notna()
        sub = df.loc[valid].copy()
        a = ai_s.loc[valid].to_numpy(dtype=float, copy=True)
        b = bi_s.loc[valid].to_numpy(dtype=float, copy=True)
        c = ci_s.loc[valid].to_numpy(dtype=float, copy=True)
        d = di_s.loc[valid].to_numpy(dtype=float, copy=True)

        zero_mask = (a == 0) | (b == 0) | (c == 0) | (d == 0)
        if np.any(zero_mask) and args.measure in ("OR", "RR"):
            a[zero_mask] += 0.5
            b[zero_mask] += 0.5
            c[zero_mask] += 0.5
            d[zero_mask] += 0.5
            warnings.append(f"applied Haldane-Anscombe +0.5 continuity correction to {int(np.sum(zero_mask))} study/studies with zero cells")

        if args.measure == "OR":
            y = np.log((a * d) / (b * c))
            se = np.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
            is_ratio = True
            scale_label = "OR"
        elif args.measure == "RR":
            n1 = a + b
            n0 = c + d
            y = np.log((a / n1) / (c / n0))
            se = np.sqrt(1.0 / a - 1.0 / n1 + 1.0 / c - 1.0 / n0)
            is_ratio = True
            scale_label = "RR"
        else:
            n1 = a + b
            n0 = c + d
            p1 = a / n1
            p0 = c / n0
            y = p1 - p0
            se = np.sqrt(p1 * (1.0 - p1) / n1 + p0 * (1.0 - p0) / n0)
            se = np.maximum(se, 1e-6)
            is_ratio = False
            scale_label = "RD"

    k = len(y)
    if k < 2:
        fail(f"meta-analysis requires >= 2 valid studies (got {k})")

    studies = sub[args.study_col].astype(str).str.strip().tolist()
    pooled = pool_studies(y, se, is_ratio)

    study_rows: list[dict[str, Any]] = []
    for i, sname in enumerate(studies):
        study_rows.append(
            {
                "study": sname,
                "estimate": _fmt_est(float(y[i]), is_ratio),
                "ci_low": _fmt_est(float(y[i] - 1.959964 * se[i]), is_ratio),
                "ci_high": _fmt_est(float(y[i] + 1.959964 * se[i]), is_ratio),
                "log_estimate": round(float(y[i]), 5),
                "se": round(float(se[i]), 5),
                "weight_fe_pct": round(float(pooled["weights_fe"][i] * 100.0), 2),
                "weight_re_pct": round(float(pooled["weights_re"][i] * 100.0), 2),
            }
        )

    # Egger's regression asymmetry test (k >= 3)
    eggers = None
    if k >= 3:
        prec = 1.0 / se
        snd = y / se
        slope, intercept, _, _, _ = stats.linregress(prec, snd)
        resid = snd - (intercept + slope * prec)
        s2 = float(np.sum(resid ** 2) / max(1, k - 2))
        x_bar = float(np.mean(prec))
        ss_x = float(np.sum((prec - x_bar) ** 2))
        se_int = math.sqrt(s2 * (1.0 / k + (x_bar ** 2) / max(ss_x, 1e-18)))
        t_egg = intercept / max(se_int, 1e-12)
        p_egg = float(2.0 * stats.t.sf(abs(t_egg), df=k - 2))
        eggers = {
            "intercept": round(float(intercept), 4),
            "se": round(se_int, 4),
            "t": round(float(t_egg), 4),
            "df": k - 2,
            "pvalue": round(p_egg, 5),
        }
        if k < 10:
            warnings.append(f"Egger's test has low statistical power with k={k} (< 10) studies")

    # Leave-One-Out sensitivity
    loo_rows: list[dict[str, Any]] = []
    for i, sname in enumerate(studies):
        mask = np.ones(k, dtype=bool)
        mask[i] = False
        sub_p = pool_studies(y[mask], se[mask], is_ratio)
        re = sub_p["random_effects_dl"]
        het = sub_p["heterogeneity"]
        loo_rows.append(
            {
                "omitted_study": sname,
                "pooled_re_estimate": re["estimate"],
                "ci_low": re["ci_low"],
                "ci_high": re["ci_high"],
                "pvalue": re["pvalue"],
                "i2_percent": het["i2_percent"],
                "tau2_dl": het["tau2_dl"],
            }
        )

    # Optional subgroup analysis
    subgroups_out = None
    if args.subgroup_col:
        if args.subgroup_col not in sub.columns:
            fail(f"subgroup_col '{args.subgroup_col}' not found in studies table")
        sg_vals = sub[args.subgroup_col].astype(str).str.strip()
        sg_levels = sorted(sg_vals.unique().tolist())
        sg_list = []
        sg_mus = []
        sg_ws = []
        for sg_name in sg_levels:
            m = (sg_vals == sg_name).to_numpy()
            if int(np.sum(m)) == 0:
                continue
            sp = pool_studies(y[m], se[m], is_ratio)
            sg_list.append(
                {
                    "subgroup": sg_name,
                    "k": int(np.sum(m)),
                    "random_effects_dl": sp["random_effects_dl"],
                    "heterogeneity": sp["heterogeneity"],
                }
            )
            sg_mus.append(sp["_mu_dl"])
            sg_ws.append(1.0 / max(sp["_se_dl"] ** 2, 1e-18))
        q_between = None
        p_between = None
        if len(sg_mus) >= 2:
            w_arr = np.array(sg_ws, dtype=float)
            m_arr = np.array(sg_mus, dtype=float)
            mu_bar = float(np.sum(w_arr * m_arr) / np.sum(w_arr))
            q_between = float(np.sum(w_arr * ((m_arr - mu_bar) ** 2)))
            p_between = float(stats.chi2.sf(q_between, df=len(sg_mus) - 1))
        subgroups_out = {
            "subgroup_col": args.subgroup_col,
            "q_between": round(q_between, 4) if q_between is not None else None,
            "df": max(0, len(sg_mus) - 1),
            "pvalue_between": round(p_between, 5) if p_between is not None else None,
            "groups": sg_list,
        }

    sum_df = pd.DataFrame(study_rows)
    loo_df = pd.DataFrame(loo_rows)
    sum_path = os.path.join(args.output_dir, "meta_summary.tsv")
    loo_path = os.path.join(args.output_dir, "meta_leave_one_out.tsv")
    sum_df.to_csv(sum_path, sep="\t", index=False)
    loo_df.to_csv(loo_path, sep="\t", index=False)

    # 2-panel Forest + Funnel Plot
    fig, (ax_for, ax_fun) = plt.subplots(1, 2, figsize=(13.0, max(4.5, 0.42 * k + 2.5)), gridspec_kw={"width_ratios": [1.55, 1.0]})

    # Forest plot on log/raw scale with exponentiated ticks if ratio
    null_val = 0.0
    ax_for.axvline(null_val, color="#475569", linestyle="--", linewidth=1.0)
    y_idx = np.arange(k)[::-1] + 2  # leave rows 0 and 1 for pooled RE and FE
    for idx_s in range(k):
        yi = y[idx_s]
        sei = se[idx_s]
        yp = y_idx[idx_s]
        ax_for.plot([yi - 1.959964 * sei, yi + 1.959964 * sei], [yp, yp], color="#334155", linewidth=1.5)
        sz = 35.0 + 180.0 * float(pooled["weights_re"][idx_s])
        ax_for.scatter([yi], [yp], color="#2563EB", marker="s", s=sz, zorder=4)

    # Pooled diamonds (FE at row 1, RE at row 0)
    for row_y, mu_v, se_v, col, lbl in [
        (1, pooled["_mu_fe"], pooled["_se_fe"], "#64748B", "Fixed-Effect"),
        (0, pooled["_mu_dl"], pooled["_se_dl"], "#E11D48", "Random-Effects (DL)"),
    ]:
        lo_v = mu_v - 1.959964 * se_v
        hi_v = mu_v + 1.959964 * se_v
        ax_for.fill([lo_v, mu_v, hi_v, mu_v], [row_y, row_y + 0.28, row_y, row_y - 0.28], color=col, alpha=0.85)

    all_lo = np.concatenate([y - 1.96 * se, [pooled["_mu_dl"] - 1.96 * pooled["_se_dl"]]])
    all_hi = np.concatenate([y + 1.96 * se, [pooled["_mu_dl"] + 1.96 * pooled["_se_dl"]]])
    span = max(0.5, float(np.max(all_hi) - np.min(all_lo)))
    ax_for.set_xlim(float(np.min(all_lo)) - 0.15 * span, float(np.max(all_hi)) + 0.65 * span)
    txt_x = float(np.max(all_hi)) + 0.08 * span
    for idx_s, r in enumerate(study_rows):
        ax_for.text(txt_x, y_idx[idx_s], f"{r['estimate']:.2f} [{r['ci_low']:.2f}, {r['ci_high']:.2f}] ({r['weight_re_pct']:.1f}%)", va="center", ha="left", fontsize=7.8)
    for row_y, blk in [(1, pooled["fixed_effect"]), (0, pooled["random_effects_dl"])]:
        ax_for.text(txt_x, row_y, f"{blk['estimate']:.2f} [{blk['ci_low']:.2f}, {blk['ci_high']:.2f}]", va="center", ha="left", fontsize=8.0, fontweight="bold")

    ax_for.set_yticks(list(y_idx) + [1, 0])
    ax_for.set_yticklabels(studies + ["Fixed-Effect", "Random-Effects (DL)"], fontsize=8.2)
    ax_for.set_xlabel(f"log({scale_label})" if is_ratio else scale_label)
    het = pooled["heterogeneity"]
    ax_for.set_title(f"Forest Plot (k={k}, I2={het['i2_percent']:.1f}%, tau2={het['tau2_dl']:.3f})", fontsize=10)
    ax_for.spines["top"].set_visible(False)
    ax_for.spines["right"].set_visible(False)

    # Funnel plot
    mu_re = pooled["_mu_dl"]
    se_max = float(np.max(se)) * 1.15
    se_line = np.linspace(0.0, se_max, 50)
    ax_fun.plot(mu_re - 1.959964 * se_line, se_line, color="#64748B", linestyle="--", linewidth=1.1, label="Pseudo 95% CI")
    ax_fun.plot(mu_re + 1.959964 * se_line, se_line, color="#64748B", linestyle="--", linewidth=1.1)
    ax_fun.axvline(mu_re, color="#E11D48", linewidth=1.2, label="Pooled RE")
    ax_fun.scatter(y, se, color="#2563EB", s=38, zorder=4)
    ax_fun.set_ylim(se_max, -0.02 * se_max)
    ax_fun.set_xlabel(f"log({scale_label})" if is_ratio else scale_label)
    ax_fun.set_ylabel("Standard Error")
    egg_txt = f"Egger p={eggers['pvalue']:.3f}" if eggers else "k < 3"
    ax_fun.set_title(f"Funnel Plot ({egg_txt})", fontsize=10)
    ax_fun.legend(frameon=False, loc="lower right", fontsize=8)
    ax_fun.spines["top"].set_visible(False)
    ax_fun.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "meta_forest_funnel.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {"scale": scale_label, "is_ratio_scale": is_ratio, "n_studies": k}
    artifacts = [
        write_artifact(fig_path, "clin_meta", "Meta-analysis Forest & Funnel plot", inputs=(args.studies,), parameters=params_meta),
        write_artifact(sum_path, "clin_meta", "Meta-analysis study weights & estimates", inputs=(args.studies,), parameters=params_meta, shape=sum_df.shape),
        write_artifact(loo_path, "clin_meta", "Leave-One-Out sensitivity table", inputs=(args.studies,), parameters=params_meta, shape=loo_df.shape),
    ]

    emit(
        {
            "studies": os.path.abspath(args.studies),
            "n_studies": k,
            "scale": scale_label,
            "is_ratio_scale": is_ratio,
            "fixed_effect": pooled["fixed_effect"],
            "random_effects_dl": pooled["random_effects_dl"],
            "random_effects_pm": pooled["random_effects_pm"],
            "hksj": pooled["hksj"],
            "prediction_interval_95": pooled["prediction_interval_95"],
            "heterogeneity": pooled["heterogeneity"],
            "eggers_test": eggers,
            "subgroups": subgroups_out,
            "study_results": study_rows,
            "leave_one_out": loo_rows,
            "summary_table": os.path.abspath(sum_path),
            "leave_one_out_table": os.path.abspath(loo_path),
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
