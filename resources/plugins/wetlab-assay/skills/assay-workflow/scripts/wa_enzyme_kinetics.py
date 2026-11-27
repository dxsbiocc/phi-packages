#!/usr/bin/env python3
"""wa_enzyme_kinetics: Steady-state Michaelis-Menten & Hill kinetics, automated
progress-curve initial-velocity (v0) extraction, AICc-ranked reversible inhibitor
mechanism discrimination (Competitive, Uncompetitive, Noncompetitive, Mixed),
kcat/Km, and Cheng-Prusoff Ki conversion.
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


def _extract_v0_from_progress(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Extract linear initial velocity v0 from (time, signal, substrate_conc, [inhibitor_conc])."""
    t_col = resolve_column(df, ["time", "time_s", "time_min", "t"], "time")
    sig_col = resolve_column(df, ["signal", "absorbance", "fluorescence", "product", "od", "rfu"], "signal")
    s_col = resolve_column(df, ["substrate_conc", "substrate", "s_conc", "s"], "substrate_conc")
    i_col = resolve_column(df, ["inhibitor_conc", "inhibitor", "i_conc", "i"], "inhibitor_conc", required=False)
    rep_col = resolve_column(df, ["replicate", "rep", "well", "curve_id"], "replicate", required=False)

    sub = df.copy()
    sub["_t"] = pd.to_numeric(sub[t_col], errors="coerce")
    sub["_sig"] = pd.to_numeric(sub[sig_col], errors="coerce")
    sub["_s"] = pd.to_numeric(sub[s_col], errors="coerce")
    sub["_i"] = pd.to_numeric(sub[i_col], errors="coerce").fillna(0.0) if i_col else 0.0
    sub["_rep"] = sub[rep_col].astype(str) if rep_col else "rep1"
    sub = sub.dropna(subset=["_t", "_sig", "_s"]).sort_values("_t")

    if len(sub) < 6:
        fail("Progress-curve table requires >= 6 valid time-point measurements.")

    warnings: list[str] = []
    out_rows: list[dict[str, Any]] = []
    low_r2_cnt = 0

    for (s_val, i_val, rep_val), grp in sub.groupby(["_s", "_i", "_rep"], sort=True):
        t_arr = grp["_t"].to_numpy(dtype=float, copy=True)
        y_arr = grp["_sig"].to_numpy(dtype=float, copy=True)
        n_pts = len(t_arr)
        if n_pts < 3:
            continue
        # Search early windows [0..k] (k from 3 to max(4, 60% of trace)) preferring high R^2 and steep early slope
        max_k = max(3, min(n_pts, max(4, int(math.ceil(0.6 * n_pts)))))
        best_slope = None
        best_r2 = -1.0
        for k in range(3, max_k + 1):
            lr = stats.linregress(t_arr[:k], y_arr[:k])
            r2 = float(lr.rvalue ** 2) if math.isfinite(lr.rvalue) else 0.0
            if r2 >= 0.98 or best_slope is None or r2 > best_r2:
                best_slope = abs(float(lr.slope))
                best_r2 = r2
                if k >= 4 and r2 >= 0.99:
                    break
        if best_r2 < 0.95:
            low_r2_cnt += 1
        out_rows.append(
            {
                "substrate_conc": float(s_val),
                "inhibitor_conc": float(i_val),
                "velocity": max(0.0, float(best_slope or 0.0)),
                "initial_r2": round(best_r2, 4),
                "replicate": str(rep_val),
            }
        )

    if not out_rows:
        fail("Could not extract initial velocities (need >= 3 time points per curve).")
    if low_r2_cnt > 0:
        warnings.append(
            f"{low_r2_cnt} progress curve(s) had early linear fit R^2 < 0.95; check substrate depletion (>15%) or lag phase."
        )
    return pd.DataFrame(out_rows), warnings


def _r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 1e-18 else 1.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kinetics_data", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--from_progress_curves", action="store_true", default=False)
    ap.add_argument("--enzyme_conc", type=float, default=None)
    ap.add_argument("--ic50", type=float, default=None)
    ap.add_argument("--cheng_prusoff_substrate", type=float, default=None)
    args = ap.parse_args()

    raw_df = read_table(args.kinetics_data)
    warnings: list[str] = []

    if args.from_progress_curves:
        mode = "progress_curves"
        vel_df, w_prog = _extract_v0_from_progress(raw_df)
        warnings.extend(w_prog)
    else:
        mode = "velocity_table"
        s_col = resolve_column(raw_df, ["substrate_conc", "substrate", "s_conc", "s"], "substrate_conc")
        v_col = resolve_column(raw_df, ["velocity", "v0", "rate", "v"], "velocity")
        i_col = resolve_column(raw_df, ["inhibitor_conc", "inhibitor", "i_conc", "i"], "inhibitor_conc", required=False)
        vel_df = pd.DataFrame(
            {
                "substrate_conc": pd.to_numeric(raw_df[s_col], errors="coerce"),
                "velocity": pd.to_numeric(raw_df[v_col], errors="coerce"),
                "inhibitor_conc": pd.to_numeric(raw_df[i_col], errors="coerce").fillna(0.0) if i_col else 0.0,
            }
        ).dropna()

    vel_df = vel_df[(vel_df["substrate_conc"] > 0) & (vel_df["velocity"] >= 0) & (vel_df["inhibitor_conc"] >= 0)].copy()
    if len(vel_df) < 5:
        fail(f"Need >= 5 valid (substrate_conc > 0, velocity >= 0) observations (got {len(vel_df)}).")

    S = vel_df["substrate_conc"].to_numpy(dtype=float, copy=True)
    I = vel_df["inhibitor_conc"].to_numpy(dtype=float, copy=True)
    V = vel_df["velocity"].to_numpy(dtype=float, copy=True)
    n_obs = len(V)

    pos_i_levels = sorted(set(float(x) for x in I if x > 1e-12))
    has_inhibitor = len(pos_i_levels) >= 1

    v_max_init = max(float(np.max(V)) * 1.15, 1e-6)
    s_med = max(float(np.median(S)), 1e-9)
    i_med = max(float(np.median(pos_i_levels)) if pos_i_levels else 1.0, 1e-9)

    model_candidates: list[dict[str, Any]] = []

    if not has_inhibitor:
        # 1) Michaelis-Menten: v = Vmax * S / (Km + S)
        def mm_res(p: np.ndarray) -> np.ndarray:
            vmax, km = p
            return (vmax * S) / (km + S) - V

        res_mm = optimize.least_squares(
            mm_res,
            x0=np.array([v_max_init, s_med]),
            bounds=([1e-12, 1e-12], [v_max_init * 50.0, float(np.max(S)) * 100.0]),
        )
        se_mm, lo_mm, hi_mm = jacobian_se_ci(res_mm, n_obs)
        rss_mm = float(np.sum(res_mm.fun ** 2))
        pred_mm = V + res_mm.fun
        model_candidates.append(
            {
                "model": "Michaelis-Menten",
                "n_params": 2,
                "rss": rss_mm,
                "aicc": compute_aicc(rss_mm, n_obs, 2),
                "r2": _r2_score(V, pred_mm),
                "pred": pred_mm,
                "params": {
                    "Vmax": float(res_mm.x[0]),
                    "Vmax_se": float(se_mm[0]),
                    "Vmax_ci95": [float(max(0.0, lo_mm[0])), float(hi_mm[0])],
                    "Km": float(res_mm.x[1]),
                    "Km_se": float(se_mm[1]),
                    "Km_ci95": [float(max(0.0, lo_mm[1])), float(hi_mm[1])],
                },
            }
        )

        # 2) Hill cooperativity: v = Vmax * S^nH / (K05^nH + S^nH)
        if len(np.unique(S)) >= 4:
            def hill_res(p: np.ndarray) -> np.ndarray:
                vmax, k05, nh = p
                sn = np.power(S, nh)
                kn = math.pow(k05, nh)
                return (vmax * sn) / (kn + sn) - V

            res_h = optimize.least_squares(
                hill_res,
                x0=np.array([float(res_mm.x[0]), float(res_mm.x[1]), 1.0]),
                bounds=([1e-12, 1e-12, 0.2], [v_max_init * 50.0, float(np.max(S)) * 100.0, 8.0]),
            )
            se_h, lo_h, hi_h = jacobian_se_ci(res_h, n_obs)
            rss_h = float(np.sum(res_h.fun ** 2))
            pred_h = V + res_h.fun
            model_candidates.append(
                {
                    "model": "Hill",
                    "n_params": 3,
                    "rss": rss_h,
                    "aicc": compute_aicc(rss_h, n_obs, 3),
                    "r2": _r2_score(V, pred_h),
                    "pred": pred_h,
                    "params": {
                        "Vmax": float(res_h.x[0]),
                        "Km": float(res_h.x[1]),
                        "K05": float(res_h.x[1]),
                        "hill_n": float(res_h.x[2]),
                        "hill_n_ci95": [float(lo_h[2]), float(hi_h[2])],
                    },
                }
            )
    else:
        # Global reversible inhibition models
        specs = [
            (
                "Competitive",
                3,
                lambda p: (p[0] * S) / (p[1] * (1.0 + I / p[2]) + S),
                [v_max_init, s_med, i_med],
                ([1e-12, 1e-12, 1e-12], [v_max_init * 50.0, float(np.max(S)) * 100.0, float(np.max(I)) * 200.0]),
                ["Vmax", "Km", "Ki"],
            ),
            (
                "Uncompetitive",
                3,
                lambda p: (p[0] * S) / (p[1] + S * (1.0 + I / p[2])),
                [v_max_init, s_med, i_med],
                ([1e-12, 1e-12, 1e-12], [v_max_init * 50.0, float(np.max(S)) * 100.0, float(np.max(I)) * 200.0]),
                ["Vmax", "Km", "Ki_prime"],
            ),
            (
                "Noncompetitive",
                3,
                lambda p: ((p[0] / (1.0 + I / p[2])) * S) / (p[1] + S),
                [v_max_init, s_med, i_med],
                ([1e-12, 1e-12, 1e-12], [v_max_init * 50.0, float(np.max(S)) * 100.0, float(np.max(I)) * 200.0]),
                ["Vmax", "Km", "Ki"],
            ),
            (
                "Mixed",
                4,
                lambda p: (p[0] * S) / (p[1] * (1.0 + I / p[2]) + S * (1.0 + I / p[3])),
                [v_max_init, s_med, i_med, i_med * 1.5],
                (
                    [1e-12, 1e-12, 1e-12, 1e-12],
                    [v_max_init * 50.0, float(np.max(S)) * 100.0, float(np.max(I)) * 200.0, float(np.max(I)) * 200.0],
                ),
                ["Vmax", "Km", "Ki", "Ki_prime"],
            ),
        ]
        for m_name, n_p, fn, x0, bnds, p_names in specs:
            if n_obs <= n_p + 1:
                continue
            res = optimize.least_squares(lambda p, _f=fn: _f(p) - V, x0=np.array(x0, dtype=float), bounds=bnds)
            se, lo, hi = jacobian_se_ci(res, n_obs)
            rss = float(np.sum(res.fun ** 2))
            pred = V + res.fun
            p_dict: dict[str, Any] = {}
            for idx_p, pn in enumerate(p_names):
                p_dict[pn] = round(float(res.x[idx_p]), 6)
                p_dict[f"{pn}_se"] = round(float(se[idx_p]), 6) if np.isfinite(se[idx_p]) else None
                p_dict[f"{pn}_ci95"] = [
                    round(float(max(0.0, lo[idx_p])), 6) if np.isfinite(lo[idx_p]) else None,
                    round(float(hi[idx_p]), 6) if np.isfinite(hi[idx_p]) else None,
                ]
            model_candidates.append(
                {
                    "model": m_name,
                    "n_params": n_p,
                    "rss": rss,
                    "aicc": compute_aicc(rss, n_obs, n_p),
                    "r2": _r2_score(V, pred),
                    "pred": pred,
                    "params": p_dict,
                }
            )

    model_candidates.sort(key=lambda d: d["aicc"])
    weights = akaike_weights([m["aicc"] for m in model_candidates])
    min_aicc = model_candidates[0]["aicc"]

    comp_rows: list[dict[str, Any]] = []
    for m, w in zip(model_candidates, weights):
        m["delta_aicc"] = round(float(m["aicc"] - min_aicc), 4)
        m["akaike_weight"] = round(float(w), 4)
        comp_rows.append(
            {
                "model": m["model"],
                "n_params": m["n_params"],
                "rss": round(m["rss"], 6),
                "r2": round(m["r2"], 5),
                "aicc": round(m["aicc"], 4),
                "delta_aicc": m["delta_aicc"],
                "akaike_weight": m["akaike_weight"],
                "Vmax": m["params"].get("Vmax"),
                "Km": m["params"].get("Km"),
                "Ki": m["params"].get("Ki"),
                "Ki_prime": m["params"].get("Ki_prime"),
                "hill_n": m["params"].get("hill_n"),
            }
        )

    best = model_candidates[0]
    if len(model_candidates) >= 2 and model_candidates[1]["delta_aicc"] < 2.0:
        warnings.append(
            f"Runner-up model '{model_candidates[1]['model']}' has delta_AICc = {model_candidates[1]['delta_aicc']:.2f} (< 2.0); "
            "collect wider [S] and [I] brackets to resolve mechanism unambiguously."
        )

    vmax_best = float(best["params"]["Vmax"])
    km_best = float(best["params"]["Km"])

    kcat = None
    cat_eff = None
    if args.enzyme_conc is not None:
        if args.enzyme_conc <= 0:
            fail("--enzyme_conc must be > 0.")
        kcat = round(vmax_best / float(args.enzyme_conc), 6)
        cat_eff = round(kcat / max(km_best, 1e-18), 6)

    cp_ki = None
    if args.ic50 is not None:
        if args.cheng_prusoff_substrate is None or args.cheng_prusoff_substrate < 0 or args.ic50 <= 0:
            fail("--ic50 requires --cheng_prusoff_substrate >= 0 and --ic50 > 0.")
        # Cheng-Prusoff competitive conversion: Ki = IC50 / (1 + [S]/Km)
        cp_ki = round(float(args.ic50) / (1.0 + float(args.cheng_prusoff_substrate) / max(km_best, 1e-18)), 6)

    vel_df["fitted_velocity"] = np.round(best["pred"], 6)
    vel_df["residual"] = np.round(V - best["pred"], 6)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    comp_path = out_dir / "kinetics_model_comparison.tsv"
    fit_path = out_dir / "velocity_fit_table.tsv"
    fig_path = out_dir / "enzyme_kinetics_plots.png"

    comp_df = pd.DataFrame(comp_rows)
    write_tsv(
        comp_df,
        comp_path,
        "wa_enzyme_kinetics",
        "Enzyme Kinetics AICc Model Comparison",
        inputs=[args.kinetics_data],
        parameters={"mode": mode, "has_inhibitor": has_inhibitor},
    )
    write_tsv(
        vel_df,
        fit_path,
        "wa_enzyme_kinetics",
        f"Enzyme Velocity Observations & {best['model']} Fits",
        inputs=[args.kinetics_data],
        parameters={"best_model": best["model"]},
    )

    # 2-panel plot: Direct nonlinear curves + Diagnostic Lineweaver-Burk overlay
    plt = setup_matplotlib()
    fig, (ax_mm, ax_lb) = plt.subplots(1, 2, figsize=(11.8, 4.8))
    palette = ["#2563eb", "#0d9488", "#d97706", "#dc2626", "#7c3aed", "#059669"]
    uniq_i = sorted(np.unique(I))
    s_grid = np.linspace(0.0, float(np.max(S)) * 1.08, 200)

    def _eval_best(s_arr: np.ndarray, i_val: float) -> np.ndarray:
        p = best["params"]
        vm = float(p["Vmax"])
        km = float(p["Km"])
        mname = best["model"]
        if mname == "Hill":
            nh = float(p["hill_n"])
            return (vm * np.power(s_arr, nh)) / (math.pow(km, nh) + np.power(s_arr, nh))
        if mname == "Competitive":
            return (vm * s_arr) / (km * (1.0 + i_val / float(p["Ki"])) + s_arr)
        if mname == "Uncompetitive":
            return (vm * s_arr) / (km + s_arr * (1.0 + i_val / float(p["Ki_prime"])))
        if mname == "Noncompetitive":
            return ((vm / (1.0 + i_val / float(p["Ki"]))) * s_arr) / (km + s_arr)
        if mname == "Mixed":
            return (vm * s_arr) / (km * (1.0 + i_val / float(p["Ki"])) + s_arr * (1.0 + i_val / float(p["Ki_prime"])))
        return (vm * s_arr) / (km + s_arr)

    inv_s_max = float(np.max(1.0 / S)) * 1.15
    inv_s_grid = np.linspace(-0.35 * inv_s_max, inv_s_max, 150)

    for idx_i, iv in enumerate(uniq_i):
        col = palette[idx_i % len(palette)]
        mask = np.isclose(I, iv)
        lbl = f"[I] = {iv:g}" if has_inhibitor else "Observed"
        ax_mm.scatter(S[mask], V[mask], color=col, s=34, alpha=0.88, label=lbl)
        ax_mm.plot(s_grid, _eval_best(s_grid, float(iv)), color=col, lw=2.0)

        pos_v = mask & (V > 1e-12)
        if np.any(pos_v):
            ax_lb.scatter(1.0 / S[pos_v], 1.0 / V[pos_v], color=col, s=32, alpha=0.88, label=lbl)
            # Analytical Lineweaver-Burk line from nonlinear parameters: 1/v = (Km_app/Vmax_app)*(1/S) + 1/Vmax_app
            p = best["params"]
            vm, km = float(p["Vmax"]), float(p["Km"])
            if best["model"] == "Competitive":
                slope_lb = (km * (1.0 + iv / float(p["Ki"]))) / vm
                int_lb = 1.0 / vm
            elif best["model"] == "Uncompetitive":
                slope_lb = km / vm
                int_lb = (1.0 + iv / float(p["Ki_prime"])) / vm
            elif best["model"] == "Noncompetitive":
                fac = 1.0 + iv / float(p["Ki"])
                slope_lb = (km * fac) / vm
                int_lb = fac / vm
            elif best["model"] == "Mixed":
                slope_lb = (km * (1.0 + iv / float(p["Ki"]))) / vm
                int_lb = (1.0 + iv / float(p["Ki_prime"])) / vm
            else:
                slope_lb = km / vm
                int_lb = 1.0 / vm
            y_lb = int_lb + slope_lb * inv_s_grid
            valid_lb = y_lb >= -0.05 * float(np.max(1.0 / V[V > 1e-12]))
            ax_lb.plot(inv_s_grid[valid_lb], y_lb[valid_lb], color=col, lw=1.6, ls="--")

    ax_mm.set_xlabel("Substrate Concentration [S]")
    ax_mm.set_ylabel("Initial Velocity (v0)")
    ax_mm.set_title(f"Nonlinear Fit: {best['model']} (R2={best['r2']:.4f})")
    ax_mm.legend(frameon=False, loc="best")

    ax_lb.axhline(0.0, color="#64748b", lw=0.9)
    ax_lb.axvline(0.0, color="#64748b", lw=0.9)
    ax_lb.set_xlabel("1 / [S]")
    ax_lb.set_ylabel("1 / v0")
    ax_lb.set_title("Diagnostic Lineweaver-Burk Overlay (Nonlinear Params)")
    ax_lb.legend(frameon=False, loc="best")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "wa_enzyme_kinetics",
        f"Enzyme Kinetics Curves & Lineweaver-Burk Overlay ({best['model']})",
        inputs=[args.kinetics_data],
        parameters={"best_model": best["model"], "has_inhibitor": has_inhibitor},
    )

    emit(
        {
            "status": "ok",
            "tool": "wa_enzyme_kinetics",
            "mode": mode,
            "has_inhibitor": has_inhibitor,
            "best_model": best["model"],
            "best_params": best["params"],
            "kcat": kcat,
            "catalytic_efficiency": cat_eff,
            "cheng_prusoff_ki": cp_ki,
            "model_comparison": comp_rows,
            "outputs": {
                "kinetics_model_comparison": str(comp_path),
                "velocity_fit_table": str(fit_path),
                "enzyme_kinetics_plots": str(fig_path),
            },
            "artifacts": [str(comp_path), str(fit_path), str(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
