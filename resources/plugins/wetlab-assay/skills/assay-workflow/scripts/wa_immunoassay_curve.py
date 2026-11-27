#!/usr/bin/env python3
"""Weighted 4PL/5PL immunoassay calibration (ELISA/Luminex/MSD/BCA) with LOD/LLOQ/ULOQ."""
from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
from scipy import optimize

from _common import (
    akaike_weights,
    compute_aicc,
    emit,
    fail,
    jacobian_se_ci,
    read_table,
    resolve_column,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Weighted 4PL/5PL immunoassay calibration (wa_immunoassay_curve)")
    ap.add_argument("--standards", required=True, help="Standards TSV/CSV with concentration and signal")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--unknowns", default=None, help="Optional unknowns TSV/CSV with sample, signal, optional dilution_factor")
    ap.add_argument("--model", default="auto", choices=["auto", "4PL", "5PL"], help="Model family (default auto via AICc)")
    ap.add_argument("--weighting", default="1_y2", choices=["1_y2", "1_y", "none"], help="Residual weight scheme (default 1_y2)")
    return ap.parse_args()


def eval_5pl(x: np.ndarray, a: float, d: float, log10_c: float, b: float, s: float) -> np.ndarray:
    c = 10.0 ** log10_c
    ratio = np.maximum(x / c, 1e-18)
    term = 1.0 + np.power(ratio, b)
    denom = np.power(np.clip(term, 1.0, 1e12), s)
    return d + (a - d) / denom


def invert_5pl(y: float, a: float, d: float, log10_c: float, b: float, s: float, increasing: bool) -> tuple[float | None, str]:
    """Analytically invert 5PL/4PL curve; return (concentration, bound_flag)."""
    c = 10.0 ** log10_c
    lo_y, hi_y = (a, d) if increasing else (d, a)
    if increasing:
        if y <= lo_y:
            return None, "below_lower_asymptote"
        if y >= hi_y:
            return None, "above_upper_asymptote"
    else:
        if y >= hi_y:
            return None, "below_lower_asymptote"
        if y <= lo_y:
            return None, "above_upper_asymptote"

    u = (a - d) / (y - d)
    if u <= 1.0:
        return None, "out_of_asymptote"
    inner = (u ** (1.0 / s)) - 1.0
    if inner <= 0.0:
        return None, "out_of_asymptote"
    x = c * (inner ** (1.0 / b))
    if not math.isfinite(x) or x <= 0.0:
        return None, "numerical_overflow"
    return float(x), "inverted"


def _compute_weights(x: np.ndarray, y: np.ndarray, scheme: str) -> np.ndarray:
    if scheme == "none":
        return np.ones_like(y, dtype=float)
    span = max(1e-6, float(np.max(y) - np.min(y)))
    floor = max(1e-4, 0.01 * span)
    # Use level-mean signal per nominal concentration so replicate outliers don't skew weights
    level_scale = np.empty_like(y, dtype=float)
    for val in np.unique(x):
        mask = x == val
        level_scale[mask] = max(floor, abs(float(np.mean(y[mask]))))
    if scheme == "1_y":
        return 1.0 / level_scale
    return 1.0 / (level_scale ** 2)


def _fit_pl_model(
    x: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    model_name: str,
    increasing: bool,
    blank_mean: float | None,
) -> dict[str, Any]:
    y_min, y_max = float(np.min(y)), float(np.max(y))
    y_span = max(1e-6, y_max - y_min)
    log_x = np.log10(x)
    min_lx, max_lx = float(np.min(log_x)), float(np.max(log_x))

    # Bounds for [A, D, log10_C, B (, S)]
    if increasing:
        a_init = blank_mean if blank_mean is not None else max(0.0, y_min - 0.02 * y_span)
        d_init = y_max + 0.05 * y_span
        a_lo, a_hi = min(-0.5 * y_span, y_min - y_span), y_max
        d_lo, d_hi = y_min, y_max + 1.5 * y_span
    else:
        a_init = blank_mean if blank_mean is not None else (y_max + 0.02 * y_span)
        d_init = max(0.0, y_min - 0.05 * y_span)
        a_lo, a_hi = y_min, y_max + 1.5 * y_span
        d_lo, d_hi = min(-0.5 * y_span, y_min - y_span), y_max

    sqrt_w = np.sqrt(weights)
    mid_lx = 0.5 * (min_lx + max_lx)
    b_lo, b_hi = 0.15, 6.0
    s_lo, s_hi = 0.1, 10.0
    lx_lo, lx_hi = min_lx - 2.0, max_lx + 2.0

    if model_name == "4PL":
        lo = np.array([a_lo, d_lo, lx_lo, b_lo], dtype=float)
        hi = np.array([a_hi, d_hi, lx_hi, b_hi], dtype=float)
        starts = [
            [a_init, d_init, mid_lx, 1.0],
            [a_init, d_init, min_lx + 0.35 * (max_lx - min_lx), 1.3],
            [a_init, d_init, min_lx + 0.65 * (max_lx - min_lx), 0.8],
        ]

        def resid(theta: np.ndarray) -> np.ndarray:
            yhat = eval_5pl(x, theta[0], theta[1], theta[2], theta[3], 1.0)
            return sqrt_w * (yhat - y)

    else:
        lo = np.array([a_lo, d_lo, lx_lo, b_lo, s_lo], dtype=float)
        hi = np.array([a_hi, d_hi, lx_hi, b_hi, s_hi], dtype=float)
        starts = [
            [a_init, d_init, mid_lx, 1.0, 1.0],
            [a_init, d_init, mid_lx - 0.3, 1.2, 2.0],
            [a_init, d_init, mid_lx + 0.3, 1.2, 0.5],
            [a_init, d_init, min_lx + 0.4 * (max_lx - min_lx), 1.5, 1.5],
        ]

        def resid(theta: np.ndarray) -> np.ndarray:
            yhat = eval_5pl(x, theta[0], theta[1], theta[2], theta[3], theta[4])
            return sqrt_w * (yhat - y)

    best_res = None
    best_cost = float("inf")
    for p0 in starts:
        p0_c = np.clip(np.asarray(p0, dtype=float), lo + 1e-6, hi - 1e-6)
        try:
            res = optimize.least_squares(resid, p0_c, bounds=(lo, hi), method="trf", max_nfev=12000)
            if res.cost < best_cost:
                best_cost = float(res.cost)
                best_res = res
        except Exception:
            continue

    if best_res is None:
        fail(f"Nonlinear regression failed to converge for {model_name}.")

    theta = best_res.x
    a = float(theta[0])
    d = float(theta[1])
    log10_c = float(theta[2])
    b = float(theta[3])
    s = float(theta[4]) if model_name == "5PL" else 1.0
    c_val = float(10.0 ** log10_c)

    y_hat = eval_5pl(x, a, d, log10_c, b, s)
    wrss = float(np.sum(weights * ((y - y_hat) ** 2)))
    rss = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(max(0.0, 1.0 - rss / ss_tot)) if ss_tot > 1e-12 else 0.0
    n_params = 4 if model_name == "4PL" else 5
    aicc = compute_aicc(wrss, len(y), n_params)
    se, ci_lo, ci_hi = jacobian_se_ci(best_res, len(y))

    return {
        "model": model_name,
        "a": a,
        "d": d,
        "c": c_val,
        "log10_c": log10_c,
        "b": b,
        "s": s,
        "wrss": wrss,
        "rss": rss,
        "r2": r2,
        "aicc": aicc,
        "n_params": n_params,
        "se": [float(v) if math.isfinite(v) else None for v in se],
        "ci_lo": [float(v) if math.isfinite(v) else None for v in ci_lo],
        "ci_hi": [float(v) if math.isfinite(v) else None for v in ci_hi],
        "y_hat": y_hat,
    }


def _plot_calibration(
    x_pos: np.ndarray,
    y_pos: np.ndarray,
    fit: dict[str, Any],
    rec_df: pd.DataFrame,
    unk_df: pd.DataFrame | None,
    lod_conc: float | None,
    lloq: float | None,
    uloq: float | None,
    out_png: Path,
) -> None:
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.8), gridspec_kw={"width_ratios": [1.3, 1.0]})

    x_lo = float(np.min(x_pos)) * 0.5
    x_hi = float(np.max(x_pos)) * 1.5
    grid_x = np.logspace(math.log10(x_lo), math.log10(x_hi), 300)
    grid_y = eval_5pl(grid_x, fit["a"], fit["d"], fit["log10_c"], fit["b"], fit["s"])

    ax1.plot(grid_x, grid_y, color="#0f766e", linewidth=2.1, label=f"{fit['model']} fit ($R^2$={fit['r2']:.4f})")
    ax1.scatter(x_pos, y_pos, color="#0f172a", s=32, zorder=4, label="Standards")

    if lloq is not None and uloq is not None and uloq > lloq:
        ax1.axvspan(lloq, uloq, color="#10b981", alpha=0.10, label=f"Quantifiable [{lloq:g}, {uloq:g}]")
        ax1.axvline(lloq, color="#059669", linestyle="--", linewidth=1.0)
        ax1.axvline(uloq, color="#059669", linestyle="--", linewidth=1.0)
    if lod_conc is not None and lod_conc > 0:
        ax1.axvline(lod_conc, color="#d97706", linestyle=":", linewidth=1.3, label=f"LOD = {lod_conc:.3g}")

    if unk_df is not None and not unk_df.empty:
        in_rng = unk_df[unk_df["range_status"] == "in_quantifiable_range"]
        if not in_rng.empty:
            ax1.scatter(
                in_rng["undiluted_conc"],
                in_rng["signal"],
                color="#0284c7",
                marker="^",
                s=42,
                zorder=5,
                label="Unknowns (in range)",
            )

    ax1.set_xscale("log")
    ax1.set_xlabel("Concentration (log scale)")
    ax1.set_ylabel("Assay Signal")
    ax1.set_title(f"A. {fit['model']} Standard Calibration Curve")
    ax1.legend(frameon=False, fontsize=8)

    # Panel B: Standard %Recovery across nominal concentrations
    nom = rec_df["nominal_conc"].to_numpy(copy=True, dtype=float)
    pct_rec = rec_df["pct_recovery"].to_numpy(copy=True, dtype=float)
    passed = rec_df["meets_fda_ema"].to_numpy(copy=True, dtype=bool)

    ax2.axhspan(80.0, 120.0, color="#10b981", alpha=0.12, label="80%-120% acceptance")
    ax2.axhline(100.0, color="#475569", linestyle="-", linewidth=0.9)
    ax2.axhline(80.0, color="#dc2626", linestyle="--", linewidth=1.0)
    ax2.axhline(120.0, color="#dc2626", linestyle="--", linewidth=1.0)

    if np.any(passed):
        ax2.scatter(nom[passed], pct_rec[passed], color="#0f766e", s=45, zorder=4, label="Pass (LLOQ-ULOQ)")
    if np.any(~passed):
        ax2.scatter(nom[~passed], pct_rec[~passed], color="#dc2626", marker="x", s=50, zorder=4, label="Fail recovery/CV")

    ax2.set_xscale("log")
    ax2.set_xlabel("Nominal Standard Concentration")
    ax2.set_ylabel("Back-Calculated %Recovery")
    ax2.set_title("B. Standard Recovery & FDA/EMA Range")
    ax2.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    std_df = read_table(args.standards)
    c_col = resolve_column(std_df, ["concentration", "conc", "nominal_conc", "std_conc"], "standards concentration")
    s_col = resolve_column(std_df, ["signal", "od", "od450", "mfi", "response", "absorbance", "rfu"], "standards signal")

    work = pd.DataFrame(
        {
            "conc": pd.to_numeric(std_df[c_col], errors="coerce"),
            "signal": pd.to_numeric(std_df[s_col], errors="coerce"),
        }
    )
    work = work[np.isfinite(work["conc"]) & np.isfinite(work["signal"]) & (work["conc"] >= 0)].copy()
    if work.empty:
        fail("No valid non-negative standard rows found.")

    blanks = work[work["conc"] == 0]
    pos_stds = work[work["conc"] > 0].sort_values("conc").copy()
    if len(pos_stds) < 4 or pos_stds["conc"].nunique() < 4:
        fail(f"At least 4 distinct positive standard concentrations are required (found {pos_stds['conc'].nunique()}).")

    x_pos = pos_stds["conc"].to_numpy(copy=True, dtype=float)
    y_pos = pos_stds["signal"].to_numpy(copy=True, dtype=float)

    low_sig = float(np.mean(y_pos[x_pos == np.min(x_pos)]))
    high_sig = float(np.mean(y_pos[x_pos == np.max(x_pos)]))
    increasing = high_sig >= low_sig

    blank_mean = float(blanks["signal"].mean()) if not blanks.empty else None
    blank_sd = float(np.std(blanks["signal"].to_numpy(copy=True, dtype=float), ddof=1)) if len(blanks) > 1 else 0.0

    weights = _compute_weights(x_pos, y_pos, args.weighting)

    candidates: list[str]
    if args.model == "auto":
        if len(pos_stds) >= 7 and pos_stds["conc"].nunique() >= 5:
            candidates = ["4PL", "5PL"]
        else:
            candidates = ["4PL"]
            warnings.append("Fewer than 7 positive standard points or <5 levels; restricting auto model to 4PL.")
    elif args.model == "5PL" and len(pos_stds) <= 5:
        fail("5PL fitting requires at least 6 positive standard observations.")
    else:
        candidates = [args.model]

    fits = [_fit_pl_model(x_pos, y_pos, weights, m, increasing, blank_mean) for m in candidates]
    aiccs = [f["aicc"] for f in fits]
    ak_weights = akaike_weights(aiccs)
    for f, w in zip(fits, ak_weights):
        f["akaike_weight"] = w

    best_fit = min(fits, key=lambda f: f["aicc"])

    # Compute LOD signal and concentration
    if blank_mean is not None and blank_sd > 0.0:
        lod_signal = blank_mean + (3.3 * blank_sd if increasing else -3.3 * blank_sd)
        lod_source = "blank_mean_plus_3.3sd"
    else:
        res_sd = float(np.std(y_pos - best_fit["y_hat"], ddof=1))
        base_sig = blank_mean if blank_mean is not None else best_fit["a"]
        lod_signal = base_sig + (3.3 * res_sd if increasing else -3.3 * res_sd)
        lod_source = "asymptote_plus_3.3_residual_sd"
        warnings.append("No replicate zero blanks with SD > 0; LOD estimated from residual SD.")

    lod_conc, _ = invert_5pl(
        lod_signal, best_fit["a"], best_fit["d"], best_fit["log10_c"], best_fit["b"], best_fit["s"], increasing
    )

    # Back-calculate every positive standard well and evaluate FDA/EMA criteria per nominal level
    back_vals = []
    for sig in y_pos:
        val, _ = invert_5pl(
            float(sig), best_fit["a"], best_fit["d"], best_fit["log10_c"], best_fit["b"], best_fit["s"], increasing
        )
        back_vals.append(val if val is not None else float("nan"))
    pos_stds["back_calc"] = back_vals

    rec_rows = []
    for nom_c, sub in pos_stds.groupby("conc", sort=True):
        bc = sub["back_calc"].to_numpy(copy=True, dtype=float)
        fin = bc[np.isfinite(bc)]
        mean_bc = float(np.mean(fin)) if fin.size else float("nan")
        sd_bc = float(np.std(fin, ddof=1)) if fin.size > 1 else 0.0
        pct_rec = (100.0 * mean_bc / float(nom_c)) if math.isfinite(mean_bc) else float("nan")
        pct_cv = (100.0 * sd_bc / mean_bc) if (math.isfinite(mean_bc) and mean_bc > 0 and fin.size > 1) else (0.0 if fin.size == 1 else float("nan"))
        meets = bool(math.isfinite(pct_rec) and 80.0 <= pct_rec <= 120.0 and math.isfinite(pct_cv) and pct_cv <= 20.0)
        rec_rows.append(
            {
                "nominal_conc": float(nom_c),
                "n_replicates": int(len(sub)),
                "mean_signal": round(float(sub["signal"].mean()), 5),
                "mean_back_calc_conc": round(mean_bc, 5) if math.isfinite(mean_bc) else None,
                "pct_recovery": round(pct_rec, 2) if math.isfinite(pct_rec) else None,
                "pct_cv": round(pct_cv, 2) if math.isfinite(pct_cv) else None,
                "meets_fda_ema": meets,
            }
        )
    rec_df = pd.DataFrame(rec_rows)

    passing = rec_df[rec_df["meets_fda_ema"]]
    if passing.empty:
        lloq, uloq = None, None
        warnings.append("No standard levels met FDA/EMA bioanalytical criteria (80-120% recovery & <=20% CV).")
    else:
        lloq = float(passing["nominal_conc"].min())
        uloq = float(passing["nominal_conc"].max())

    # Interpolate unknowns if provided
    unk_df_out: pd.DataFrame | None = None
    p_unk: str | None = None
    status_counts: dict[str, int] = {}
    inputs = [args.standards]

    if args.unknowns:
        inputs.append(args.unknowns)
        unk_raw = read_table(args.unknowns)
        u_samp_col = resolve_column(unk_raw, ["sample", "sample_id", "unknown", "id"], "unknowns sample")
        u_sig_col = resolve_column(unk_raw, ["signal", "od", "od450", "mfi", "response", "absorbance", "rfu"], "unknowns signal")
        u_dil_col = resolve_column(unk_raw, ["dilution_factor", "dilution", "df"], "dilution_factor", required=False)
        u_grp_col = resolve_column(unk_raw, ["group", "condition", "cohort"], "group", required=False)

        unk_rows = []
        for idx, row in unk_raw.iterrows():
            samp_name = str(row[u_samp_col]).strip()
            sig_val = pd.to_numeric(row[u_sig_col], errors="coerce")
            dil_val = pd.to_numeric(row[u_dil_col], errors="coerce") if u_dil_col else 1.0
            if not math.isfinite(dil_val) or dil_val <= 0:
                dil_val = 1.0
            grp_val = str(row[u_grp_col]).strip() if u_grp_col and pd.notna(row[u_grp_col]) else None

            if not math.isfinite(sig_val):
                status = "invalid_signal"
                undil = None
                rep_conc = None
                rep_label = "NA"
            else:
                undil, inv_flag = invert_5pl(
                    float(sig_val),
                    best_fit["a"],
                    best_fit["d"],
                    best_fit["log10_c"],
                    best_fit["b"],
                    best_fit["s"],
                    increasing,
                )
                is_below_lod_sig = (sig_val < lod_signal) if increasing else (sig_val > lod_signal)
                if inv_flag == "below_lower_asymptote" or is_below_lod_sig or (
                    undil is not None and lod_conc is not None and undil < lod_conc
                ):
                    status = "below_LOD"
                    rep_conc = None
                    rep_label = f"< LOD ({lod_conc:.4g})" if lod_conc is not None else "< LOD"
                elif inv_flag == "above_upper_asymptote" or (undil is not None and uloq is not None and undil > uloq):
                    status = "above_ULOQ"
                    rep_conc = None
                    rep_label = f"> ULOQ ({uloq * float(dil_val):.4g})" if uloq is not None else "> ULOQ"
                elif undil is None:
                    status = "out_of_range"
                    rep_conc = None
                    rep_label = "out_of_range"
                elif lloq is not None and undil < lloq:
                    status = "below_LLOQ"
                    rep_conc = None
                    rep_label = f"< LLOQ ({lloq * float(dil_val):.4g})"
                else:
                    status = "in_quantifiable_range"
                    rep_conc = float(undil) * float(dil_val)
                    rep_label = f"{rep_conc:.6g}"

            status_counts[status] = status_counts.get(status, 0) + 1
            rec_item: dict[str, Any] = {
                "sample": samp_name,
                "signal": round(float(sig_val), 5) if math.isfinite(sig_val) else None,
                "dilution_factor": float(dil_val),
                "undiluted_conc": round(undil, 6) if undil is not None else None,
                "reported_conc": round(rep_conc, 6) if rep_conc is not None else None,
                "reported_conc_label": rep_label,
                "range_status": status,
            }
            if grp_val is not None:
                rec_item["group"] = grp_val
            unk_rows.append(rec_item)

        unk_df_out = pd.DataFrame(unk_rows)
        n_censored = sum(v for k, v in status_counts.items() if k != "in_quantifiable_range")
        if n_censored > 0:
            warnings.append(
                f"{n_censored}/{len(unk_df_out)} unknown samples fall outside [LLOQ, ULOQ] and were censored (no silent extrapolation)."
            )

    params_meta = {"model": args.model, "selected_model": best_fit["model"], "weighting": args.weighting}
    p_rec = write_tsv(
        rec_df,
        out_dir / "standard_curve_recovery.tsv",
        tool="wa_immunoassay_curve",
        title=f"Immunoassay {best_fit['model']} Standard Curve Back-Calculated Recovery & CV",
        inputs=inputs,
        parameters=params_meta,
    )
    artifacts = [p_rec]
    if unk_df_out is not None:
        p_unk = write_tsv(
            unk_df_out,
            out_dir / "interpolated_unknowns.tsv",
            tool="wa_immunoassay_curve",
            title="Immunoassay Unknown Sample Interpolations & Range Censoring",
            inputs=inputs,
            parameters=params_meta,
        )
        artifacts.append(p_unk)

    fig_path = out_dir / "immunoassay_calibration.png"
    _plot_calibration(x_pos, y_pos, best_fit, rec_df, unk_df_out, lod_conc, lloq, uloq, fig_path)
    p_fig = write_artifact(
        fig_path,
        tool="wa_immunoassay_curve",
        title=f"Immunoassay {best_fit['model']} Calibration Curve & %Recovery Profile",
        inputs=inputs,
        parameters=params_meta,
    )
    artifacts.append(p_fig)

    model_comp = [
        {
            "model": f["model"],
            "n_params": f["n_params"],
            "wrss": round(f["wrss"], 6),
            "r2": round(f["r2"], 5),
            "aicc": round(f["aicc"], 4),
            "akaike_weight": round(f["akaike_weight"], 4),
        }
        for f in sorted(fits, key=lambda x: x["aicc"])
    ]

    emit(
        {
            "status": "ok",
            "tool": "wa_immunoassay_curve",
            "selected_model": best_fit["model"],
            "weighting": args.weighting,
            "assay_direction": "increasing" if increasing else "decreasing",
            "parameters": {
                "A": round(best_fit["a"], 6),
                "B": round(best_fit["b"], 6),
                "C_ec50": round(best_fit["c"], 6),
                "D": round(best_fit["d"], 6),
                "S_asymmetry": round(best_fit["s"], 6),
            },
            "model_comparison": model_comp,
            "r2": round(best_fit["r2"], 5),
            "aicc": round(best_fit["aicc"], 4),
            "lod_signal": round(lod_signal, 5),
            "lod_source": lod_source,
            "lod_conc": round(lod_conc, 6) if lod_conc is not None else None,
            "lloq": lloq,
            "uloq": uloq,
            "n_standards": int(len(pos_stds)),
            "n_blanks": int(len(blanks)),
            "n_unknowns": int(len(unk_df_out)) if unk_df_out is not None else 0,
            "unknown_status_counts": status_counts,
            "outputs": {
                "standard_curve_recovery": p_rec,
                "interpolated_unknowns": p_unk,
                "immunoassay_calibration": p_fig,
            },
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
