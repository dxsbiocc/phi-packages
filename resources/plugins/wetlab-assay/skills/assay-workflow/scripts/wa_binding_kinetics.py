#!/usr/bin/env python3
"""wa_binding_kinetics: Global multi-concentration 1:1 Langmuir SPR/BLI binding
kinetics (kon, koff, KD, half-life, Rmax, reduced chi2) with mass-transport &
slow-dissociation guardrails and steady-state KD_eq isotherm cross-validation.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
import pandas as pd
from scipy import optimize, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    jacobian_se_ci,
    read_table,
    resolve_column,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)

UNIT_TO_MOLAR = {"M": 1.0, "uM": 1e-6, "nM": 1e-9}


def _eval_trace(
    t: np.ndarray,
    c_m: float,
    t0: float,
    td: float,
    kon: float,
    koff: float,
    rmax: float,
) -> np.ndarray:
    k_obs = kon * c_m + koff
    r_eq = (kon * c_m * rmax) / max(k_obs, 1e-30)
    dt_assoc = np.clip(np.minimum(t, td) - t0, 0.0, None)
    r_assoc = r_eq * (1.0 - np.exp(-np.clip(k_obs * dt_assoc, 0.0, 700.0)))
    r_at_td = r_eq * (1.0 - math.exp(-min(700.0, k_obs * max(0.0, td - t0))))
    out = np.where(
        t <= td,
        r_assoc,
        r_at_td * np.exp(-np.clip(koff * (t - td), 0.0, 700.0)),
    )
    return np.where(t < t0, 0.0, out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sensorgrams", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--dissociation_start", type=float, default=None)
    ap.add_argument("--conc_unit", choices=["nM", "uM", "M"], default="nM")
    ap.add_argument("--local_rmax", action="store_true", default=False)
    args = ap.parse_args()

    df = read_table(args.sensorgrams)
    t_col = resolve_column(df, ["time", "time_s", "t", "seconds"], "time")
    r_col = resolve_column(df, ["response", "ru", "nm", "signal", "binding"], "response")
    c_col = resolve_column(df, ["conc_m", "conc_nm", "conc_um", "concentration", "conc", "analyte_conc"], "concentration")
    ph_col = resolve_column(df, ["phase", "step"], "phase", required=False)

    sub = df.copy()
    sub["_t"] = pd.to_numeric(sub[t_col], errors="coerce")
    sub["_r"] = pd.to_numeric(sub[r_col], errors="coerce")
    sub["_c_raw"] = pd.to_numeric(sub[c_col], errors="coerce")
    sub = sub.dropna(subset=["_t", "_r", "_c_raw"])
    sub = sub[sub["_c_raw"] > 0].sort_values(["_c_raw", "_t"]).copy()

    if len(sub) < 12:
        fail("Sensorgram table requires >= 12 valid positive-concentration time points.")

    # Infer unit if column name explicitly specifies conc_m / conc_nm / conc_um
    c_col_low = c_col.lower()
    if c_col_low == "conc_m":
        unit_scale = 1.0
        eff_unit = "M"
    elif c_col_low == "conc_um":
        unit_scale = 1e-6
        eff_unit = "uM"
    elif c_col_low == "conc_nm":
        unit_scale = 1e-9
        eff_unit = "nM"
    else:
        unit_scale = UNIT_TO_MOLAR[args.conc_unit]
        eff_unit = args.conc_unit

    sub["_c_m"] = sub["_c_raw"] * unit_scale

    if args.dissociation_start is not None:
        td = float(args.dissociation_start)
    elif ph_col is not None:
        ph_low = sub[ph_col].astype(str).str.strip().str.lower()
        diss_mask = ph_low.str.startswith("diss")
        if not diss_mask.any():
            fail("phase column has no 'dissoc'/'dissociation' rows; pass --dissociation_start.")
        td = float(sub.loc[diss_mask, "_t"].min())
    else:
        fail("Pass --dissociation_start (seconds) when sensorgrams table has no phase column.")

    t0 = float(sub["_t"].min())
    if td <= t0 or td >= float(sub["_t"].max()):
        fail(f"dissociation_start ({td}) must lie strictly inside time span ({t0}, {float(sub['_t'].max())}).")

    conc_levels_m = sorted(np.unique(sub["_c_m"].to_numpy(dtype=float, copy=True)))
    n_conc = len(conc_levels_m)
    if n_conc < 2:
        fail(f"Global SPR/BLI kinetic fitting requires >= 2 analyte concentrations (got {n_conc}).")

    traces = []
    for cm in conc_levels_m:
        g = sub[np.isclose(sub["_c_m"], cm, rtol=1e-9, atol=0.0)]
        traces.append(
            (
                float(cm),
                g["_t"].to_numpy(dtype=float, copy=True),
                g["_r"].to_numpy(dtype=float, copy=True),
            )
        )

    r_peak = max(float(sub["_r"].max()), 1e-6)
    c_med = float(np.median(conc_levels_m))

    # Estimate initial koff from highest-concentration dissociation slope
    _, t_hi, r_hi = traces[-1]
    d_mask = (t_hi > td) & (r_hi > 0.05 * r_peak)
    if int(d_mask.sum()) >= 3:
        slope_d = stats.linregress(t_hi[d_mask] - td, np.log(np.maximum(r_hi[d_mask], 1e-9))).slope
        koff_init = float(np.clip(-slope_d, 1e-5, 0.5))
    else:
        koff_init = 1e-3
    kon_init = float(np.clip(koff_init / max(c_med, 1e-12), 1e2, 1e7))

    rmax_mode = "local" if args.local_rmax else "global"
    if not args.local_rmax:
        # Parameters: [log10_kon, log10_koff, rmax]
        x0 = np.array([math.log10(kon_init), math.log10(koff_init), r_peak * 1.15], dtype=float)
        bounds = ([1.0, -7.0, 1e-6], [9.0, 1.0, r_peak * 25.0])

        def resid_fn(p: np.ndarray) -> np.ndarray:
            kon = 10.0 ** p[0]
            koff = 10.0 ** p[1]
            rm = p[2]
            chunks = [_eval_trace(t_i, c_i, t0, td, kon, koff, rm) - r_i for c_i, t_i, r_i in traces]
            return np.concatenate(chunks)
    else:
        x0 = np.array([math.log10(kon_init), math.log10(koff_init)] + [r_peak * 1.15] * n_conc, dtype=float)
        bounds = ([1.0, -7.0] + [1e-6] * n_conc, [9.0, 1.0] + [r_peak * 25.0] * n_conc)

        def resid_fn(p: np.ndarray) -> np.ndarray:
            kon = 10.0 ** p[0]
            koff = 10.0 ** p[1]
            chunks = [
                _eval_trace(t_i, c_i, t0, td, kon, koff, p[2 + idx]) - r_i
                for idx, (c_i, t_i, r_i) in enumerate(traces)
            ]
            return np.concatenate(chunks)

    res = optimize.least_squares(resid_fn, x0=x0, bounds=bounds, x_scale="jac", max_nfev=4000)
    n_obs = len(res.fun)
    n_params = len(res.x)
    se, ci_lo, ci_hi = jacobian_se_ci(res, n_obs)

    kon = float(10.0 ** res.x[0])
    koff = float(10.0 ** res.x[1])
    kon_ci = [float(10.0 ** ci_lo[0]), float(10.0 ** ci_hi[0])] if np.isfinite(se[0]) else [kon, kon]
    koff_ci = [float(10.0 ** ci_lo[1]), float(10.0 ** ci_hi[1])] if np.isfinite(se[1]) else [koff, koff]

    kd_m = koff / max(kon, 1e-30)
    kd_nm = kd_m * 1e9
    # Delta method in log10 space: log10(KD) = log10(koff) - log10(kon)
    try:
        df_deg = max(1, n_obs - n_params)
        s2 = float(np.sum(res.fun ** 2)) / df_deg
        cov = np.linalg.pinv(res.jac.T @ res.jac) * s2
        var_log10_kd = max(0.0, float(cov[1, 1] + cov[0, 0] - 2.0 * cov[0, 1]))
        t_crit = float(stats.t.ppf(0.975, df_deg))
        se_log10_kd = math.sqrt(var_log10_kd)
        kd_ci_nm = [
            float((10.0 ** (res.x[1] - res.x[0] - t_crit * se_log10_kd)) * 1e9),
            float((10.0 ** (res.x[1] - res.x[0] + t_crit * se_log10_kd)) * 1e9),
        ]
    except Exception:
        kd_ci_nm = [kd_nm, kd_nm]

    half_life_s = math.log(2.0) / max(koff, 1e-30)
    rss = float(np.sum(res.fun ** 2))
    red_chi2 = rss / max(1, n_obs - n_params)
    all_r = np.concatenate([r_i for _, _, r_i in traces])
    ss_tot = float(np.sum((all_r - np.mean(all_r)) ** 2))
    r2 = 1.0 - rss / ss_tot if ss_tot > 1e-18 else 1.0

    if not args.local_rmax:
        rmax_val: float | None = round(float(res.x[2]), 4)
        local_rmax_map = None
    else:
        rmax_val = None
        local_rmax_map = {
            f"{c_i / unit_scale:g}_{eff_unit}": round(float(res.x[2 + idx]), 4)
            for idx, (c_i, _, _) in enumerate(traces)
        }

    # Steady-state equilibrium isotherm cross-check across late association window
    req_obs = []
    for c_i, t_i, r_i in traces:
        assoc_m = t_i <= td
        if np.any(assoc_m):
            t_a = t_i[assoc_m]
            r_a = r_i[assoc_m]
            win_lo = t0 + 0.85 * (td - t0)
            late = t_a >= win_lo
            req_obs.append((c_i, float(np.mean(r_a[late]) if np.any(late) else r_a[-1])))

    kd_eq_m = None
    kd_eq_nm = None
    rmax_eq = None
    if len(req_obs) >= 4:
        c_eq = np.array([x[0] for x in req_obs], dtype=float)
        r_eq_y = np.array([x[1] for x in req_obs], dtype=float)
        try:
            res_eq = optimize.least_squares(
                lambda p: (c_eq * p[0]) / ((10.0 ** p[1]) + c_eq) - r_eq_y,
                x0=np.array([max(float(np.max(r_eq_y)) * 1.2, 1e-6), math.log10(max(kd_m, 1e-14))]),
                bounds=([1e-6, -15.0], [r_peak * 50.0, 0.0]),
                x_scale="jac",
            )
            rmax_eq = round(float(res_eq.x[0]), 4)
            kd_eq_m = float(10.0 ** res_eq.x[1])
            kd_eq_nm = round(kd_eq_m * 1e9, 4)
        except Exception:
            kd_eq_m = None

    warnings: list[str] = []
    if kon > 1e7:
        warnings.append(
            f"Fitted kon ({kon:.2e} M^-1 s^-1) > 1e7 M^-1 s^-1: suspect mass-transport limitation (reduce ligand immobilization density or increase flow rate)."
        )
    if koff < 1e-5:
        warnings.append(
            f"Fitted koff ({koff:.2e} s^-1) < 1e-5 s^-1: insufficient dissociation curvature; extend dissociation wash phase (>=1800 s)."
        )
    ref_rmax = rmax_val if rmax_val is not None else float(np.mean(res.x[2:]))
    if math.sqrt(max(red_chi2, 0.0)) > 0.05 * ref_rmax:
        warnings.append(
            f"Elevated fit RMSE ({math.sqrt(red_chi2):.2f} > 5% of Rmax={ref_rmax:.1f}): inspect 1:1 stoichiometry, bulk refractive index jumps, or surface heterogeneity."
        )
    if kd_eq_m is not None:
        ratio_kd = max(kd_m / max(kd_eq_m, 1e-30), kd_eq_m / max(kd_m, 1e-30))
        if ratio_kd > 3.0:
            warnings.append(
                f"Kinetic KD ({kd_nm:.2f} nM) and steady-state KD_eq ({kd_eq_nm:.2f} nM) differ by {ratio_kd:.1f}x (>3x); association traces may not have reached equilibrium."
            )

    # Build output tables
    fit_rows: list[dict[str, Any]] = []
    for idx, (c_i, t_i, r_i) in enumerate(traces):
        rm_i = float(res.x[2 + idx]) if args.local_rmax else float(res.x[2])
        pred_i = _eval_trace(t_i, c_i, t0, td, kon, koff, rm_i)
        c_disp = round(c_i / unit_scale, 6)
        for tv, rv, pv in zip(t_i, r_i, pred_i):
            fit_rows.append(
                {
                    "conc": c_disp,
                    "conc_unit": eff_unit,
                    "conc_m": c_i,
                    "time": float(tv),
                    "observed_response": round(float(rv), 5),
                    "fitted_response": round(float(pv), 5),
                    "residual": round(float(rv - pv), 5),
                }
            )

    sum_df = pd.DataFrame(
        [
            {
                "rmax_mode": rmax_mode,
                "n_concentrations": n_conc,
                "kon_m_s": kon,
                "kon_ci95_low": kon_ci[0],
                "kon_ci95_high": kon_ci[1],
                "koff_s": koff,
                "koff_ci95_low": koff_ci[0],
                "koff_ci95_high": koff_ci[1],
                "kd_m": kd_m,
                "kd_nm": round(kd_nm, 4),
                "kd_ci95_low_nm": round(kd_ci_nm[0], 4),
                "kd_ci95_high_nm": round(kd_ci_nm[1], 4),
                "half_life_s": round(half_life_s, 3),
                "rmax": rmax_val,
                "reduced_chi2": round(red_chi2, 5),
                "r2": round(r2, 5),
                "kd_eq_nm": kd_eq_nm,
                "rmax_eq": rmax_eq,
            }
        ]
    )
    fit_df = pd.DataFrame(fit_rows)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sum_path = out_dir / "binding_kinetics_summary.tsv"
    fits_path = out_dir / "sensorgram_fits.tsv"
    fig_path = out_dir / "binding_kinetics_sensorgrams.png"

    write_tsv(
        sum_df,
        sum_path,
        "wa_binding_kinetics",
        "1:1 Langmuir SPR/BLI Kinetic & Steady-State Summary",
        inputs=[args.sensorgrams],
        parameters={"rmax_mode": rmax_mode, "dissociation_start": td},
    )
    write_tsv(
        fit_df,
        fits_path,
        "wa_binding_kinetics",
        "SPR/BLI Multi-Concentration Fitted Sensorgram Traces",
        inputs=[args.sensorgrams],
        parameters={"rmax_mode": rmax_mode},
    )

    plt = setup_matplotlib()
    fig = plt.figure(figsize=(12.2, 5.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[3.0, 1.1], width_ratios=[1.45, 1.0], hspace=0.28, wspace=0.26)
    ax_sens = fig.add_subplot(gs[0, 0])
    ax_res = fig.add_subplot(gs[1, 0], sharex=ax_sens)
    ax_iso = fig.add_subplot(gs[:, 1])

    palette = ["#2563eb", "#0d9488", "#d97706", "#dc2626", "#7c3aed", "#059669", "#db2777"]
    for idx, (c_i, t_i, r_i) in enumerate(traces):
        col = palette[idx % len(palette)]
        rm_i = float(res.x[2 + idx]) if args.local_rmax else float(res.x[2])
        pred_i = _eval_trace(t_i, c_i, t0, td, kon, koff, rm_i)
        lbl = f"{c_i / unit_scale:g} {eff_unit}"
        ax_sens.scatter(t_i, r_i, s=8, color=col, alpha=0.55)
        ax_sens.plot(t_i, pred_i, lw=1.8, color=col, label=lbl)
        ax_res.plot(t_i, r_i - pred_i, lw=1.0, color=col, alpha=0.8)

    ax_sens.axvline(td, color="#475569", ls="--", lw=1.0)
    ax_sens.set_ylabel("Response (RU / nm)")
    ax_sens.set_title(
        f"1:1 Langmuir Global Fit (KD={kd_nm:.2f} nM, kon={kon:.2e} M^-1s^-1, koff={koff:.2e} s^-1)"
    )
    ax_sens.legend(frameon=False, loc="best", ncol=2, fontsize=7.5)

    ax_res.axhline(0.0, color="#0f172a", ls=":", lw=0.9)
    ax_res.axvline(td, color="#475569", ls="--", lw=1.0)
    ax_res.set_xlabel("Time (s)")
    ax_res.set_ylabel("Residual")

    c_disp_pts = np.array([x[0] / unit_scale for x in req_obs], dtype=float)
    r_disp_pts = np.array([x[1] for x in req_obs], dtype=float)
    ax_iso.scatter(c_disp_pts, r_disp_pts, color="#e11d48", s=44, zorder=4, label="Observed Req")
    if kd_eq_m is not None and rmax_eq is not None:
        c_grid_m = np.geomspace(min(conc_levels_m) * 0.4, max(conc_levels_m) * 2.0, 150)
        r_grid_eq = (c_grid_m * rmax_eq) / (kd_eq_m + c_grid_m)
        ax_iso.plot(c_grid_m / unit_scale, r_grid_eq, color="#2563eb", lw=2.0, label=f"Steady-state (KD_eq={kd_eq_nm:.2f} nM)")
        ax_iso.axvline(kd_eq_m / unit_scale, color="#64748b", ls=":", lw=1.0)
    ax_iso.set_xscale("log")
    ax_iso.set_xlabel(f"Analyte Concentration ({eff_unit}, log scale)")
    ax_iso.set_ylabel("Equilibrium Response (Req)")
    ax_iso.set_title("Steady-State Binding Isotherm")
    ax_iso.legend(frameon=False, loc="best")

    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "wa_binding_kinetics",
        "1:1 Langmuir SPR/BLI Sensorgrams & Steady-State Isotherm",
        inputs=[args.sensorgrams],
        parameters={"rmax_mode": rmax_mode, "kd_nm": round(kd_nm, 4)},
    )

    emit(
        {
            "status": "ok",
            "tool": "wa_binding_kinetics",
            "rmax_mode": rmax_mode,
            "conc_unit": eff_unit,
            "dissociation_start": td,
            "n_concentrations": n_conc,
            "kon_m_s": kon,
            "kon_ci95": kon_ci,
            "koff_s": koff,
            "koff_ci95": koff_ci,
            "kd_m": kd_m,
            "kd_nm": round(kd_nm, 4),
            "kd_ci95_nm": [round(kd_ci_nm[0], 4), round(kd_ci_nm[1], 4)],
            "half_life_s": round(half_life_s, 3),
            "rmax": rmax_val,
            "local_rmax_values": local_rmax_map,
            "reduced_chi2": round(red_chi2, 5),
            "r2": round(r2, 5),
            "kd_eq_m": kd_eq_m,
            "kd_eq_nm": kd_eq_nm,
            "rmax_eq": rmax_eq,
            "outputs": {
                "binding_kinetics_summary": str(sum_path),
                "sensorgram_fits": str(fits_path),
                "binding_kinetics_sensorgrams": str(fig_path),
            },
            "artifacts": [str(sum_path), str(fits_path), str(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
