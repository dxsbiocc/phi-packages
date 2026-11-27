#!/usr/bin/env python3
"""clin_sample_size: Clinical study sample-size & statistical power calculator
supporting survival log-rank (Schoenfeld), two proportions (Fleiss/Casagrande
superiority & Farrington-Manning non-inferiority), two means (exact non-central t),
and diagnostic ROC AUC (Hanley & McNeil / Obuchowski), with power-curve export.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from typing import Any, Callable

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pred_common import emit, fail, write_artifact


def _z_alpha(alpha: float, one_sided: bool) -> float:
    return float(stats.norm.ppf(1.0 - alpha)) if one_sided else float(stats.norm.ppf(1.0 - alpha / 2.0))


def _split_arms(n_total: int, k: float) -> tuple[int, int]:
    n_ctrl = max(1, int(round(n_total / (1.0 + k))))
    n_exp = max(1, n_total - n_ctrl)
    if n_ctrl + n_exp != n_total:
        n_ctrl = max(1, n_total - n_exp)
    return n_ctrl, n_exp


def _survival_event_prob(
    hr: float,
    k: float,
    event_rate_ctrl: float | None,
    median_surv_ctrl: float | None,
    accrual_time: float | None,
    followup_time: float | None,
) -> tuple[float, float, float, str]:
    if event_rate_ctrl is not None:
        if not (0.0 < event_rate_ctrl <= 1.0):
            fail(f"event_rate_ctrl must be in (0, 1] (got {event_rate_ctrl})")
        p_c = float(event_rate_ctrl)
        p_e = 1.0 if p_c >= 1.0 else 1.0 - ((1.0 - p_c) ** hr)
        p_pool = (p_c + k * p_e) / (1.0 + k)
        return p_c, p_e, p_pool, "direct event_rate_ctrl with exponential PH scaling"

    if median_surv_ctrl is not None:
        if median_surv_ctrl <= 0.0:
            fail(f"median_surv_ctrl must be > 0 (got {median_surv_ctrl})")
        acc = 0.0 if accrual_time is None else float(accrual_time)
        fu = 0.0 if followup_time is None else float(followup_time)
        if acc < 0.0 or fu < 0.0 or (acc + fu) <= 0.0:
            fail("when using median_surv_ctrl, provide accrual_time >= 0 and followup_time >= 0 with accrual_time + followup_time > 0")
        lam_c = math.log(2.0) / float(median_surv_ctrl)
        lam_e = lam_c * hr

        def _simpson_ev(lam: float) -> float:
            if acc == 0.0:
                return 1.0 - math.exp(-lam * fu)
            s_f = math.exp(-lam * fu)
            s_m = math.exp(-lam * (fu + 0.5 * acc))
            s_t = math.exp(-lam * (fu + acc))
            return 1.0 - (s_f + 4.0 * s_m + s_t) / 6.0

        p_c = _simpson_ev(lam_c)
        p_e = _simpson_ev(lam_e)
        p_pool = (p_c + k * p_e) / (1.0 + k)
        return p_c, p_e, p_pool, "exponential survival via Schoenfeld/Simpson rule"

    fail("survival_logrank requires either --event_rate_ctrl or --median_surv_ctrl (with --accrual_time / --followup_time)")
    raise AssertionError("unreachable")


def _fm_null_probs(p_c: float, p_e: float, k: float, delta: float) -> tuple[float, float]:
    """Restricted MLE (p_c_tilde, p_e_tilde) under H0: p_e - p_c = -delta (i.e. p_e_tilde = p_c_tilde - delta)."""
    lo = max(1e-8, delta + 1e-8)
    hi = min(1.0 - 1e-8, 1.0 - 1e-8)
    if lo >= hi:
        return min(max(p_c, 1e-4), 1.0 - 1e-4), min(max(p_c - delta, 1e-4), 1.0 - 1e-4)

    def neg_ll(pc_t: float) -> float:
        pe_t = pc_t - delta
        ll1 = p_c * math.log(pc_t) + (1.0 - p_c) * math.log(1.0 - pc_t)
        ll2 = k * (p_e * math.log(pe_t) + (1.0 - p_e) * math.log(1.0 - pe_t))
        return -(ll1 + ll2)

    res = optimize.minimize_scalar(neg_ll, bounds=(lo, hi), method="bounded")
    pc_t = float(res.x)
    pe_t = pc_t - delta
    return pc_t, pe_t


def _hanley_mcneil_var(auc: float, n_pos: float, n_neg: float) -> float:
    q1 = auc / (2.0 - auc)
    q2 = (2.0 * auc * auc) / (1.0 + auc)
    num = auc * (1.0 - auc) + max(0.0, n_pos - 1.0) * (q1 - auc * auc) + max(0.0, n_neg - 1.0) * (q2 - auc * auc)
    return max(1e-18, num / max(1e-9, n_pos * n_neg))


def _build_power_fn(
    design: str,
    alpha: float,
    one_sided: bool,
    k: float,
    params: dict[str, Any],
) -> tuple[Callable[[int], float], Callable[[float], tuple[int, int, int, dict[str, Any]]]]:
    """Return (power_from_n_eval(n_eval), solve_n_eval(target_power)) for a given scenario."""
    z_a = _z_alpha(alpha, one_sided)

    if design == "survival_logrank":
        hr = float(params["hr"])
        p_pool = float(params["p_event_pooled"])
        log_hr = math.log(hr)
        fac = ((1.0 + k) ** 2) / (k * (log_hr ** 2))

        def power_fn(n_eval: int) -> float:
            d_exp = max(1e-9, float(n_eval) * p_pool)
            z_b = math.sqrt(d_exp / fac) - z_a
            return float(stats.norm.cdf(z_b))

        def solve_fn(target_pow: float) -> tuple[int, int, int, dict[str, Any]]:
            z_b = float(stats.norm.ppf(target_pow))
            d_exact = ((z_a + z_b) ** 2) * fac
            d_req = int(math.ceil(d_exact))
            n_eval = max(2, int(math.ceil(d_req / max(p_pool, 1e-9))))
            while power_fn(n_eval) < target_pow and n_eval < 10_000_000:
                n_eval += 1
            n_ctrl, n_exp = _split_arms(n_eval, k)
            return n_eval, n_ctrl, n_exp, {"required_events": d_req, "required_events_exact": round(d_exact, 3)}

        return power_fn, solve_fn

    if design == "two_proportions":
        p_c = float(params["p_ctrl"])
        p_e = float(params["p_exp"])
        ni_margin = params.get("non_inferiority_margin")

        if ni_margin is None:
            delta = abs(p_e - p_c)
            p_bar = (p_c + k * p_e) / (1.0 + k)
            v0 = (1.0 + 1.0 / k) * p_bar * (1.0 - p_bar)
            v1 = p_c * (1.0 - p_c) + p_e * (1.0 - p_e) / k

            def _n_ctrl_fleiss_cont(z_b: float) -> float:
                num = z_a * math.sqrt(v0) + z_b * math.sqrt(v1)
                if num <= 0.0:
                    return 0.5
                n_raw = (num ** 2) / (delta ** 2)
                return 0.25 * n_raw * ((1.0 + math.sqrt(1.0 + (2.0 * (1.0 + 1.0 / k)) / (n_raw * delta))) ** 2)

            def power_fn(n_eval: int) -> float:
                n_c, _ = _split_arms(n_eval, k)
                if n_c <= 1:
                    return float(alpha if one_sided else alpha / 2.0)
                if _n_ctrl_fleiss_cont(-8.0) >= n_c:
                    return 0.0
                if _n_ctrl_fleiss_cont(8.0) <= n_c:
                    return 1.0
                z_b = float(optimize.brentq(lambda zb: _n_ctrl_fleiss_cont(zb) - n_c, -8.0, 8.0))
                return float(stats.norm.cdf(z_b))

            def solve_fn(target_pow: float) -> tuple[int, int, int, dict[str, Any]]:
                z_b = float(stats.norm.ppf(target_pow))
                n_c = max(2, int(math.ceil(_n_ctrl_fleiss_cont(z_b))))
                n_e = max(2, int(math.ceil(k * n_c)))
                n_eval = n_c + n_e
                while power_fn(n_eval) < target_pow and n_eval < 10_000_000:
                    n_eval += 1
                n_c, n_e = _split_arms(n_eval, k)
                return n_eval, n_c, n_e, {}

            return power_fn, solve_fn

        margin = float(ni_margin)
        sep = (p_e - p_c) + margin
        pc_t, pe_t = _fm_null_probs(p_c, p_e, k, margin)
        sd0 = math.sqrt(pc_t * (1.0 - pc_t) + pe_t * (1.0 - pe_t) / k)
        sd1 = math.sqrt(p_c * (1.0 - p_c) + p_e * (1.0 - p_e) / k)

        def power_fn(n_eval: int) -> float:
            n_c, _ = _split_arms(n_eval, k)
            z_b = (math.sqrt(max(1.0, float(n_c))) * sep - z_a * sd0) / max(sd1, 1e-12)
            return float(stats.norm.cdf(z_b))

        def solve_fn(target_pow: float) -> tuple[int, int, int, dict[str, Any]]:
            z_b = float(stats.norm.ppf(target_pow))
            n_c_exact = ((z_a * sd0 + z_b * sd1) / sep) ** 2
            n_c = max(2, int(math.ceil(n_c_exact)))
            n_e = max(2, int(math.ceil(k * n_c)))
            n_eval = n_c + n_e
            while power_fn(n_eval) < target_pow and n_eval < 10_000_000:
                n_eval += 1
            n_c, n_e = _split_arms(n_eval, k)
            return n_eval, n_c, n_e, {"null_p_ctrl": round(pc_t, 5), "null_p_exp": round(pe_t, 5)}

        return power_fn, solve_fn

    if design == "two_means":
        diff = abs(float(params["mean_diff"]))
        sd = float(params["sd"])
        cohen_d = diff / sd

        def _t_power(n_c: int, n_e: int) -> float:
            df = max(1, n_c + n_e - 2)
            ncp = cohen_d / math.sqrt(1.0 / n_c + 1.0 / n_e)
            if one_sided:
                t_crit = float(stats.t.ppf(1.0 - alpha, df=df))
                return float(stats.nct.sf(t_crit, df=df, nc=ncp))
            t_crit = float(stats.t.ppf(1.0 - alpha / 2.0, df=df))
            return float(stats.nct.sf(t_crit, df=df, nc=ncp) + stats.nct.cdf(-t_crit, df=df, nc=ncp))

        def power_fn(n_eval: int) -> float:
            n_c, n_e = _split_arms(max(4, n_eval), k)
            return _t_power(n_c, n_e)

        def solve_fn(target_pow: float) -> tuple[int, int, int, dict[str, Any]]:
            z_b = float(stats.norm.ppf(target_pow))
            n_c_approx = ((z_a + z_b) / cohen_d) ** 2 * (1.0 + 1.0 / k)
            n_eval = max(4, int(math.floor((1.0 + k) * n_c_approx)) - 4)
            while power_fn(n_eval) < target_pow and n_eval < 10_000_000:
                # Jump quickly if far away
                step = 1 if n_eval < 200 else max(1, int(0.02 * n_eval))
                if step > 1 and power_fn(n_eval + step) < target_pow:
                    n_eval += step
                else:
                    n_eval += 1
            n_c, n_e = _split_arms(n_eval, k)
            return n_eval, n_c, n_e, {"cohens_d": round(cohen_d, 4), "df": n_c + n_e - 2}

        return power_fn, solve_fn

    if design == "diagnostic_auc":
        auc_1 = float(params["auc_target"])
        auc_0 = float(params["auc_null"])
        prev = float(params["prevalence"])
        diff = abs(auc_1 - auc_0)

        def power_fn(n_eval: int) -> float:
            n_pos = max(1.0, prev * float(n_eval))
            n_neg = max(1.0, (1.0 - prev) * float(n_eval))
            v0 = _hanley_mcneil_var(auc_0, n_pos, n_neg)
            v1 = _hanley_mcneil_var(auc_1, n_pos, n_neg)
            z_b = (diff - z_a * math.sqrt(v0)) / math.sqrt(v1)
            return float(stats.norm.cdf(z_b))

        def solve_fn(target_pow: float) -> tuple[int, int, int, dict[str, Any]]:
            z_b = float(stats.norm.ppf(target_pow))
            # Large-sample Hanley-McNeil initializer
            def _unit_v(th: float) -> float:
                q1 = th / (2.0 - th)
                q2 = (2.0 * th * th) / (1.0 + th)
                return (q1 - th * th) / (1.0 - prev) + (q2 - th * th) / prev

            n_init = ((z_a * math.sqrt(_unit_v(auc_0)) + z_b * math.sqrt(_unit_v(auc_1))) / diff) ** 2
            n_eval = max(10, int(math.floor(n_init)) - 10)
            while power_fn(n_eval) < target_pow and n_eval < 10_000_000:
                step = 1 if n_eval < 300 else max(1, int(0.02 * n_eval))
                if step > 1 and power_fn(n_eval + step) < target_pow:
                    n_eval += step
                else:
                    n_eval += 1
            n_cases = max(1, int(math.ceil(prev * n_eval)))
            n_ctrls = max(1, n_eval - n_cases)
            n_eval = n_cases + n_ctrls
            return n_eval, n_ctrls, n_cases, {}

        return power_fn, solve_fn

    fail(f"unsupported design: {design}")
    raise AssertionError("unreachable")


def _scenario_params(design: str, base: dict[str, Any], mult: float) -> dict[str, Any]:
    out = dict(base)
    if design == "survival_logrank":
        hr_new = float(base["hr"]) ** mult
        out["hr"] = hr_new
        # Recompute pooled event probability under scenario HR when median_surv_ctrl mode was used, or keep scaled PH
        if base.get("event_rate_ctrl") is not None:
            p_c = float(base["event_rate_ctrl"])
            p_e = 1.0 if p_c >= 1.0 else 1.0 - ((1.0 - p_c) ** hr_new)
            out["p_event_pooled"] = (p_c + float(base["k"]) * p_e) / (1.0 + float(base["k"]))
        else:
            _, _, p_pool, _ = _survival_event_prob(
                hr_new,
                float(base["k"]),
                None,
                base.get("median_surv_ctrl"),
                base.get("accrual_time"),
                base.get("followup_time"),
            )
            out["p_event_pooled"] = p_pool
        return out

    if design == "two_proportions":
        p_c = float(base["p_ctrl"])
        p_e = float(base["p_exp"])
        ni = base.get("non_inferiority_margin")
        if ni is None:
            p_e_new = float(np.clip(p_c + mult * (p_e - p_c), 0.005, 0.995))
            if abs(p_e_new - p_c) < 1e-4:
                p_e_new = p_e
            out["p_exp"] = p_e_new
        else:
            margin = float(ni)
            sep = (p_e - p_c) + margin
            sep_new = max(1e-3, mult * sep)
            p_e_new = float(np.clip(p_c - margin + sep_new, 0.005, 0.995))
            out["p_exp"] = p_e_new
        return out

    if design == "two_means":
        out["mean_diff"] = float(base["mean_diff"]) * mult
        return out

    if design == "diagnostic_auc":
        a0 = float(base["auc_null"])
        a1 = float(base["auc_target"])
        out["auc_target"] = float(np.clip(a0 + mult * (a1 - a0), a0 + 1e-3, 0.995))
        return out

    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", required=True, choices=["survival_logrank", "two_proportions", "two_means", "diagnostic_auc"])
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--solve_for", choices=["sample_size", "power"], default="sample_size")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--power", type=float, default=0.80)
    parser.add_argument("--n_total", type=int, default=None)
    parser.add_argument("--allocation_ratio", type=float, default=1.0)
    parser.add_argument("--one_sided", action="store_true")
    parser.add_argument("--dropout_rate", type=float, default=0.10)
    # survival_logrank
    parser.add_argument("--hr", type=float, default=None)
    parser.add_argument("--event_rate_ctrl", type=float, default=None)
    parser.add_argument("--median_surv_ctrl", type=float, default=None)
    parser.add_argument("--accrual_time", type=float, default=None)
    parser.add_argument("--followup_time", type=float, default=None)
    # two_proportions
    parser.add_argument("--p_ctrl", type=float, default=None)
    parser.add_argument("--p_exp", type=float, default=None)
    parser.add_argument("--non_inferiority_margin", type=float, default=None)
    # two_means
    parser.add_argument("--mean_diff", type=float, default=None)
    parser.add_argument("--sd", type=float, default=None)
    # diagnostic_auc
    parser.add_argument("--auc_target", type=float, default=None)
    parser.add_argument("--auc_null", type=float, default=0.50)
    parser.add_argument("--prevalence", type=float, default=0.30)
    args = parser.parse_args()

    if not (0.0 < args.alpha < 0.5):
        fail(f"alpha must be in (0, 0.5) (got {args.alpha})")
    if not (0.0 < args.power < 1.0):
        fail(f"power must be in (0, 1) (got {args.power})")
    if args.allocation_ratio <= 0.0:
        fail(f"allocation_ratio must be > 0 (got {args.allocation_ratio})")
    if not (0.0 <= args.dropout_rate < 1.0):
        fail(f"dropout_rate must be in [0, 1) (got {args.dropout_rate})")
    if args.solve_for == "power" and (args.n_total is None or args.n_total < 4):
        fail("n_total >= 4 is required when solve_for='power'")

    os.makedirs(args.output_dir, exist_ok=True)
    warnings: list[str] = []
    k = float(args.allocation_ratio)
    d = float(args.dropout_rate)

    base_params: dict[str, Any] = {"k": k}
    design_details: dict[str, Any] = {}

    if args.design == "survival_logrank":
        if args.hr is None or args.hr <= 0.0 or abs(math.log(args.hr)) < 1e-9:
            fail("survival_logrank requires --hr > 0 and != 1.0")
        p_c, p_e, p_pool, ev_note = _survival_event_prob(
            args.hr, k, args.event_rate_ctrl, args.median_surv_ctrl, args.accrual_time, args.followup_time
        )
        method = "Schoenfeld (1983) log-rank"
        base_params.update(
            {
                "hr": float(args.hr),
                "event_rate_ctrl": args.event_rate_ctrl,
                "median_surv_ctrl": args.median_surv_ctrl,
                "accrual_time": args.accrual_time,
                "followup_time": args.followup_time,
                "p_event_pooled": p_pool,
            }
        )
        design_details = {
            "hr": float(args.hr),
            "p_event_ctrl": round(p_c, 5),
            "p_event_exp": round(p_e, 5),
            "p_event_pooled": round(p_pool, 5),
            "event_prob_source": ev_note,
        }
        g1_lbl, g2_lbl = "Control", "Experimental"

    elif args.design == "two_proportions":
        if args.p_ctrl is None or args.p_exp is None:
            fail("two_proportions requires --p_ctrl and --p_exp")
        if not (0.0 < args.p_ctrl < 1.0) or not (0.0 < args.p_exp < 1.0):
            fail("p_ctrl and p_exp must be strictly in (0, 1)")
        if args.non_inferiority_margin is not None:
            margin = float(args.non_inferiority_margin)
            if margin <= 0.0 or margin >= 1.0:
                fail("non_inferiority_margin must be in (0, 1)")
            if (args.p_exp - args.p_ctrl) + margin <= 1e-9:
                fail("non-inferiority requires p_exp - p_ctrl > -non_inferiority_margin (favorable direction)")
            method = "Farrington-Manning (1990) non-inferiority"
            if not args.one_sided:
                warnings.append("non-inferiority margin evaluated at two-sided alpha (equivalent to one-sided alpha/2)")
        else:
            if abs(args.p_exp - args.p_ctrl) < 1e-9:
                fail("two_proportions superiority requires p_exp != p_ctrl")
            method = "Fleiss / Casagrande continuity-corrected normal-binomial"
        base_params.update(
            {
                "p_ctrl": float(args.p_ctrl),
                "p_exp": float(args.p_exp),
                "non_inferiority_margin": args.non_inferiority_margin,
            }
        )
        design_details = {
            "p_ctrl": float(args.p_ctrl),
            "p_exp": float(args.p_exp),
            "risk_difference": round(float(args.p_exp - args.p_ctrl), 5),
            "non_inferiority_margin": args.non_inferiority_margin,
        }
        g1_lbl, g2_lbl = "Control", "Experimental"

    elif args.design == "two_means":
        if args.mean_diff is None or args.sd is None:
            fail("two_means requires --mean_diff and --sd")
        if args.sd <= 0.0 or abs(args.mean_diff) < 1e-12:
            fail("two_means requires --sd > 0 and --mean_diff != 0")
        method = "Exact non-central t-distribution (scipy.stats.nct)"
        base_params.update({"mean_diff": float(args.mean_diff), "sd": float(args.sd)})
        design_details = {
            "mean_diff": float(args.mean_diff),
            "sd": float(args.sd),
            "cohens_d": round(abs(float(args.mean_diff)) / float(args.sd), 4),
        }
        g1_lbl, g2_lbl = "Control", "Experimental"

    else:  # diagnostic_auc
        if args.auc_target is None:
            fail("diagnostic_auc requires --auc_target")
        if not (0.0 < args.auc_null < 1.0) or not (0.0 < args.auc_target < 1.0):
            fail("auc_target and auc_null must be strictly in (0, 1)")
        if abs(args.auc_target - args.auc_null) < 1e-6:
            fail("auc_target must differ from auc_null")
        if not (0.0 < args.prevalence < 1.0):
            fail("prevalence must be strictly in (0, 1)")
        method = "Hanley & McNeil (1982) / Obuchowski AUC"
        base_params.update(
            {
                "auc_target": float(args.auc_target),
                "auc_null": float(args.auc_null),
                "prevalence": float(args.prevalence),
            }
        )
        design_details = {
            "auc_target": float(args.auc_target),
            "auc_null": float(args.auc_null),
            "prevalence": float(args.prevalence),
        }
        g1_lbl, g2_lbl = "Controls", "Cases"

    power_fn, solve_fn = _build_power_fn(args.design, args.alpha, args.one_sided, k, base_params)

    req_events: int | None = None
    if args.solve_for == "sample_size":
        n_eval, g1_eval, g2_eval, extra = solve_fn(args.power)
        achieved_pow = power_fn(n_eval)
        target_pow: float | None = float(args.power)
        n_enroll = int(math.ceil(n_eval / (1.0 - d))) if d > 0.0 else n_eval
        if args.design == "diagnostic_auc":
            g2_enr = max(1, int(math.ceil(float(args.prevalence) * n_enroll)))
            g1_enr = max(1, n_enroll - g2_enr)
            n_enroll = g1_enr + g2_enr
        else:
            g1_enr, g2_enr = _split_arms(n_enroll, k)
    else:
        n_enroll = int(args.n_total)
        n_eval = max(4, int(math.floor(n_enroll * (1.0 - d))))
        if args.design == "diagnostic_auc":
            g2_eval = max(1, int(round(float(args.prevalence) * n_eval)))
            g1_eval = max(1, n_eval - g2_eval)
            g2_enr = max(1, int(round(float(args.prevalence) * n_enroll)))
            g1_enr = max(1, n_enroll - g2_enr)
        else:
            g1_eval, g2_eval = _split_arms(n_eval, k)
            g1_enr, g2_enr = _split_arms(n_enroll, k)
        achieved_pow = power_fn(n_eval)
        target_pow = None
        extra = {}
        if args.design == "survival_logrank":
            req_events = int(round(n_eval * float(base_params["p_event_pooled"])))
            extra = {"expected_events": req_events}

    if args.design == "survival_logrank" and "required_events" in extra:
        req_events = int(extra["required_events"])
    design_details.update(extra)

    if d > 0.0:
        warnings.append(f"enrolled sample size inflated by 1 / (1 - {d:.2f}) for anticipated {d * 100:.1f}% dropout/attrition")

    # Build power curve grid across conservative (0.8x), target (1.0x), optimistic (1.2x)
    n_lo = max(12, int(round(0.25 * n_enroll)))
    n_hi = max(n_lo + 40, int(round(2.0 * n_enroll)))
    n_grid_enr = np.unique(np.round(np.linspace(n_lo, n_hi, 45)).astype(int))
    if n_enroll not in n_grid_enr:
        n_grid_enr = np.sort(np.append(n_grid_enr, n_enroll))

    scen_specs = [
        ("conservative", 0.80),
        ("target", 1.00),
        ("optimistic", 1.20),
    ]
    scen_fns = {}
    for sname, mult in scen_specs:
        sparams = _scenario_params(args.design, base_params, mult)
        pfn, _ = _build_power_fn(args.design, args.alpha, args.one_sided, k, sparams)
        scen_fns[sname] = pfn

    grid_rows = []
    for n_enr_pt in n_grid_enr:
        n_ev_pt = max(4, int(math.floor(int(n_enr_pt) * (1.0 - d))))
        row = {
            "n_total_enrolled": int(n_enr_pt),
            "n_total_evaluable": int(n_ev_pt),
            "power_conservative": round(scen_fns["conservative"](n_ev_pt), 5),
            "power_target": round(scen_fns["target"](n_ev_pt), 5),
            "power_optimistic": round(scen_fns["optimistic"](n_ev_pt), 5),
        }
        if args.design == "survival_logrank":
            row["expected_events"] = round(n_ev_pt * float(base_params["p_event_pooled"]), 2)
        grid_rows.append(row)

    grid_df = pd.DataFrame(grid_rows)
    grid_path = os.path.join(args.output_dir, "power_curve_grid.tsv")
    grid_df.to_csv(grid_path, sep="\t", index=False)

    # Power curve figure
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    n_vals = grid_df["n_total_enrolled"].to_numpy(dtype=float, copy=True)
    ax.plot(n_vals, grid_df["power_conservative"].to_numpy(dtype=float, copy=True), color="#94A3B8", linestyle="--", linewidth=1.8, label="Conservative (0.8x effect)")
    ax.plot(n_vals, grid_df["power_target"].to_numpy(dtype=float, copy=True), color="#2563EB", linewidth=2.4, label="Target effect")
    ax.plot(n_vals, grid_df["power_optimistic"].to_numpy(dtype=float, copy=True), color="#10B981", linestyle="-.", linewidth=1.8, label="Optimistic (1.2x effect)")

    ref_p = target_pow if target_pow is not None else achieved_pow
    ax.axhline(ref_p, color="#E11D48", linestyle=":", linewidth=1.2, label=f"Power = {ref_p:.2f}")
    ax.axvline(n_enroll, color="#334155", linestyle=":", linewidth=1.2, label=f"N_enroll = {n_enroll}")
    ax.scatter([n_enroll], [achieved_pow], color="#E11D48", s=55, zorder=5)

    ax.set_ylim(0.0, 1.03)
    ax.set_xlabel("Total Enrolled Sample Size (N)")
    ax.set_ylabel("Statistical Power (1 - beta)")
    side_str = "1-sided" if args.one_sided else "2-sided"
    ev_title = f" | Events = {req_events}" if req_events is not None else ""
    ax.set_title(
        f"Power vs Sample Size ({args.design}: {method})\n"
        f"alpha={args.alpha} ({side_str}), dropout={d * 100:.0f}%, N_enroll={n_enroll}, Power={achieved_pow:.3f}{ev_title}",
        fontsize=9.5,
    )
    ax.legend(frameon=False, loc="lower right", fontsize=8.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "sample_size_power_curve.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    params_meta = {
        "design": args.design,
        "solve_for": args.solve_for,
        "alpha": args.alpha,
        "one_sided": args.one_sided,
        "allocation_ratio": k,
        "dropout_rate": d,
    }
    artifacts = [
        write_artifact(fig_path, "clin_sample_size", f"Sample size & power curve ({args.design})", parameters=params_meta),
        write_artifact(grid_path, "clin_sample_size", "Power curve N vs Power grid", parameters=params_meta, shape=grid_df.shape),
    ]

    emit(
        {
            "design": args.design,
            "solve_for": args.solve_for,
            "method": method,
            "alpha": float(args.alpha),
            "one_sided": bool(args.one_sided),
            "allocation_ratio": k,
            "dropout_rate": d,
            "target_power": round(target_pow, 4) if target_pow is not None else None,
            "achieved_power": round(float(achieved_pow), 4),
            "n_total_evaluable": int(n_eval),
            "n_total_enrolled": int(n_enroll),
            "arm_breakdown": {
                "group1_label": g1_lbl,
                "group2_label": g2_lbl,
                "group1_evaluable": int(g1_eval),
                "group2_evaluable": int(g2_eval),
                "group1_enrolled": int(g1_enr),
                "group2_enrolled": int(g2_enr),
            },
            "required_events": req_events,
            "design_details": design_details,
            "power_curve_table": os.path.abspath(grid_path),
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
