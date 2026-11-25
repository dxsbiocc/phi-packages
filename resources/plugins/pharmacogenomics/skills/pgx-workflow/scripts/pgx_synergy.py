#!/usr/bin/env python3
"""2D drug-combination synergy analysis across ZIP, Bliss, Loewe, and HSA models."""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    detect_viability_scale,
    emit,
    fail,
    fit_hill_4pl,
    hill_4pl_log10,
    read_table,
    setup_matplotlib,
    write_artifact,
)


def _pick_col(cols: list[str], requested: str, candidates: list[str], *, required: bool = True) -> str | None:
    if requested:
        if requested not in cols:
            fail(f"Requested column '{requested}' not found in input columns: {cols}")
        return requested
    lower_map = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand in lower_map:
            return lower_map[cand]
    if required:
        fail(f"Could not auto-detect column from candidates {candidates} in input columns: {cols}")
    return None


def _inv_hill_dose(inhib: float, fit: dict) -> float:
    """Concentration achieving fractional inhibition `inhib` = 1 - y under a 4PL Hill fit."""
    y_target = 1.0 - inhib
    e0, einf = fit["e0"], fit["einf"]
    if y_target >= e0 - 1e-6:
        return 1e-12
    if y_target <= einf + 1e-6:
        return float("inf")
    ratio = (e0 - y_target) / (y_target - einf)
    if ratio <= 0:
        return float("inf")
    ld = fit["log10_ec50"] + math.log10(ratio) / fit["hill"]
    if ld > 25.0:
        return float("inf")
    return float(10.0 ** ld)


def _loewe_expected_inhibition(da: float, db: float, fit_a: dict, fit_b: dict, hsa_exp: float) -> float:
    """Solve da / D_A(E) + db / D_B(E) = 1 via bisection on fractional inhibition E."""
    max_e_a = max(0.0, 1.0 - fit_a["einf"])
    max_e_b = max(0.0, 1.0 - fit_b["einf"])
    e_ceil = max(max_e_a, max_e_b) - 1e-4
    if e_ceil <= 0.01:
        return hsa_exp

    def iso_lhs(e_val: float) -> float:
        d_req_a = _inv_hill_dose(e_val, fit_a)
        d_req_b = _inv_hill_dose(e_val, fit_b)
        term_a = da / d_req_a if math.isfinite(d_req_a) and d_req_a > 0 else 0.0
        term_b = db / d_req_b if math.isfinite(d_req_b) and d_req_b > 0 else 0.0
        return term_a + term_b

    lo_e, hi_e = 1e-4, min(0.999, e_ceil)
    if iso_lhs(lo_e) <= 1.0:
        return lo_e
    if iso_lhs(hi_e) >= 1.0:
        return max(hi_e, hsa_exp)
    for _ in range(40):
        mid = 0.5 * (lo_e + hi_e)
        if iso_lhs(mid) > 1.0:
            lo_e = mid
        else:
            hi_e = mid
    return 0.5 * (lo_e + hi_e)


def _classify_synergy(score: float) -> str:
    if score > 10.0:
        return "synergistic"
    if score < -10.0:
        return "antagonistic"
    return "additive"


def main() -> None:
    ap = argparse.ArgumentParser(description="2D drug-combination synergy scoring (ZIP, Bliss, Loewe, HSA)")
    ap.add_argument("--input", required=True, help="Long-format 2D combination TSV/CSV")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--dose_a_col", default="", help="Column for Drug A dose")
    ap.add_argument("--dose_b_col", default="", help="Column for Drug B dose")
    ap.add_argument("--viability_col", default="", help="Column for cell viability")
    ap.add_argument("--drug_a_name", default="", help="Display name for Drug A")
    ap.add_argument("--drug_b_name", default="", help="Display name for Drug B")
    ap.add_argument("--sample_col", default="", help="Optional sample/cell-line column")
    ap.add_argument("--viability_scale", choices=["auto", "fraction", "percent"], default="auto")
    ap.add_argument("--dose_unit", default="uM", help="Concentration unit label")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    df = read_table(args.input)
    cols = list(df.columns)
    da_col = _pick_col(cols, args.dose_a_col, ["dose_a", "conc_a", "conc1", "dose1", "drug_a_conc"])
    db_col = _pick_col(cols, args.dose_b_col, ["dose_b", "conc_b", "conc2", "dose2", "drug_b_conc"])
    v_col = _pick_col(cols, args.viability_col, ["viability", "response", "cell_viability", "pct_viability"])
    samp_col = _pick_col(cols, args.sample_col, ["sample", "cell_line", "cellline", "model"], required=False)
    drug_a_col = _pick_col(cols, "", ["drug_a", "drug1", "compound_a"], required=False)
    drug_b_col = _pick_col(cols, "", ["drug_b", "drug2", "compound_b"], required=False)

    df = df.copy()
    df["_da"] = pd.to_numeric(df[da_col], errors="coerce")
    df["_db"] = pd.to_numeric(df[db_col], errors="coerce")
    df["_v"] = pd.to_numeric(df[v_col], errors="coerce")
    df = df.dropna(subset=["_da", "_db", "_v"])
    if df.empty:
        fail("No valid numeric (dose_a, dose_b, viability) rows found.")
    if (df["_da"] < 0).any() or (df["_db"] < 0).any():
        fail("Negative concentrations encountered in combination matrix.")

    scale_name, scale_factor = detect_viability_scale(df["_v"].to_numpy(), args.viability_scale)
    df["_v_frac"] = df["_v"] * scale_factor
    df["_sample"] = df[samp_col].astype(str) if samp_col else "Sample"
    df["_drug_a"] = (
        args.drug_a_name
        if args.drug_a_name
        else (df[drug_a_col].astype(str) if drug_a_col else "DrugA")
    )
    df["_drug_b"] = (
        args.drug_b_name
        if args.drug_b_name
        else (df[drug_b_col].astype(str) if drug_b_col else "DrugB")
    )

    warnings: list[str] = []
    summary_rows: list[dict] = []
    grid_rows: list[dict] = []
    first_plot_data = None

    for (samp, drug_a, drug_b), grp in df.groupby(["_sample", "_drug_a", "_drug_b"], sort=False):
        agg = grp.groupby(["_da", "_db"], as_index=False)["_v_frac"].mean()
        veh = agg[(agg["_da"] == 0.0) & (agg["_db"] == 0.0)]["_v_frac"]
        if len(veh) > 0 and float(veh.iloc[0]) > 0.1:
            agg["_v_norm"] = agg["_v_frac"] / float(veh.iloc[0])
        else:
            agg["_v_norm"] = agg["_v_frac"]
            warnings.append(
                f"({samp}: {drug_a} + {drug_b}) had no (0, 0) vehicle well; assuming input viability is already vehicle-normalized."
            )

        mono_a = agg[(agg["_da"] > 0.0) & (agg["_db"] == 0.0)].sort_values("_da")
        mono_b = agg[(agg["_db"] > 0.0) & (agg["_da"] == 0.0)].sort_values("_db")
        combo = agg[(agg["_da"] > 0.0) & (agg["_db"] > 0.0)].copy()

        if len(mono_a) < 3 or len(mono_b) < 3 or combo.empty:
            fail(
                f"Combination ({samp}: {drug_a} + {drug_b}) requires >= 3 positive monotherapy doses along both "
                "dose_b=0 and dose_a=0 margins plus interior combination wells."
            )

        # Duplicate points if exactly 3 monotherapy concentrations so 4PL has >=4 observations
        da_m = mono_a["_da"].to_numpy(dtype=float)
        va_m = mono_a["_v_norm"].to_numpy(dtype=float)
        if da_m.size == 3:
            da_m = np.repeat(da_m, 2)
            va_m = np.repeat(va_m, 2)
        fit_a = fit_hill_4pl(da_m, va_m, fix_e0=1.0)

        db_m = mono_b["_db"].to_numpy(dtype=float)
        vb_m = mono_b["_v_norm"].to_numpy(dtype=float)
        if db_m.size == 3:
            db_m = np.repeat(db_m, 2)
            vb_m = np.repeat(vb_m, 2)
        fit_b = fit_hill_4pl(db_m, vb_m, fix_e0=1.0)

        doses_a = sorted(combo["_da"].unique().tolist())
        doses_b = sorted(combo["_db"].unique().tolist())

        # Conditional row/column Hill fits for ZIP smoothing
        cond_y_ab: dict[tuple[float, float], float] = {}
        for db_val in doses_b:
            yb_hat = float(np.clip(hill_4pl_log10([math.log10(db_val)], fit_b["e0"], fit_b["einf"], fit_b["log10_ec50"], fit_b["hill"])[0], 0.02, 1.2))
            slice_a = combo[combo["_db"] == db_val].sort_values("_da")
            if len(slice_a) >= 3:
                s_d = slice_a["_da"].to_numpy(dtype=float)
                s_v = slice_a["_v_norm"].to_numpy(dtype=float)
                if s_d.size == 3:
                    s_d, s_v = np.repeat(s_d, 2), np.repeat(s_v, 2)
                c_fit = fit_hill_4pl(s_d, s_v, fix_e0=yb_hat)
                for da_val in slice_a["_da"]:
                    y_pred = float(hill_4pl_log10([math.log10(da_val)], c_fit["e0"], c_fit["einf"], c_fit["log10_ec50"], c_fit["hill"])[0])
                    cond_y_ab[(float(da_val), float(db_val))] = float(np.clip(y_pred, 0.0, 1.0))

        cond_y_ba: dict[tuple[float, float], float] = {}
        for da_val in doses_a:
            ya_hat = float(np.clip(hill_4pl_log10([math.log10(da_val)], fit_a["e0"], fit_a["einf"], fit_a["log10_ec50"], fit_a["hill"])[0], 0.02, 1.2))
            slice_b = combo[combo["_da"] == da_val].sort_values("_db")
            if len(slice_b) >= 3:
                s_d = slice_b["_db"].to_numpy(dtype=float)
                s_v = slice_b["_v_norm"].to_numpy(dtype=float)
                if s_d.size == 3:
                    s_d, s_v = np.repeat(s_d, 2), np.repeat(s_v, 2)
                c_fit = fit_hill_4pl(s_d, s_v, fix_e0=ya_hat)
                for db_val in slice_b["_db"]:
                    y_pred = float(hill_4pl_log10([math.log10(db_val)], c_fit["e0"], c_fit["einf"], c_fit["log10_ec50"], c_fit["hill"])[0])
                    cond_y_ba[(float(da_val), float(db_val))] = float(np.clip(y_pred, 0.0, 1.0))

        # Observed single-agent inhibition maps (fallback to fitted Hill if exact dose missing on margin)
        obs_a_map = {float(r["_da"]): float(1.0 - np.clip(r["_v_norm"], 0.0, 1.0)) for _, r in mono_a.iterrows()}
        obs_b_map = {float(r["_db"]): float(1.0 - np.clip(r["_v_norm"], 0.0, 1.0)) for _, r in mono_b.iterrows()}

        sub_grid = []
        for _, r in combo.iterrows():
            da_v = float(r["_da"])
            db_v = float(r["_db"])
            v_obs = float(np.clip(r["_v_norm"], 0.0, 1.0))
            e_obs = 1.0 - v_obs

            ya_fit = float(np.clip(hill_4pl_log10([math.log10(da_v)], fit_a["e0"], fit_a["einf"], fit_a["log10_ec50"], fit_a["hill"])[0], 0.0, 1.0))
            yb_fit = float(np.clip(hill_4pl_log10([math.log10(db_v)], fit_b["e0"], fit_b["einf"], fit_b["log10_ec50"], fit_b["hill"])[0], 0.0, 1.0))
            ea_fit = 1.0 - ya_fit
            eb_fit = 1.0 - yb_fit

            ea_obs = obs_a_map.get(da_v, ea_fit)
            eb_obs = obs_b_map.get(db_v, eb_fit)

            e_hsa = max(ea_obs, eb_obs)
            e_bliss = ea_obs + eb_obs - ea_obs * eb_obs
            e_loewe = _loewe_expected_inhibition(da_v, db_v, fit_a, fit_b, max(ea_fit, eb_fit))

            y1 = cond_y_ab.get((da_v, db_v), v_obs)
            y2 = cond_y_ba.get((da_v, db_v), v_obs)
            e_zip_obs = 1.0 - 0.5 * (y1 + y2)
            e_zip_exp = ea_fit + eb_fit - ea_fit * eb_fit

            cell = {
                "sample": samp,
                "drug_a": drug_a,
                "drug_b": drug_b,
                "dose_a": da_v,
                "dose_b": db_v,
                "viability_norm": round(v_obs, 4),
                "inhibition_pct": round(100.0 * e_obs, 2),
                "zip_delta": round(100.0 * (e_zip_obs - e_zip_exp), 2),
                "bliss_delta": round(100.0 * (e_obs - e_bliss), 2),
                "loewe_delta": round(100.0 * (e_obs - e_loewe), 2),
                "hsa_delta": round(100.0 * (e_obs - e_hsa), 2),
            }
            sub_grid.append(cell)
            grid_rows.append(cell)

        sub_df = pd.DataFrame(sub_grid)
        zip_piv = sub_df.pivot(index="dose_b", columns="dose_a", values="zip_delta").sort_index(ascending=True)
        inh_piv = agg.copy()
        inh_piv["inh_pct"] = 100.0 * (1.0 - np.clip(inh_piv["_v_norm"], 0.0, 1.0))
        inh_mat = inh_piv.pivot(index="_db", columns="_da", values="inh_pct").sort_index(ascending=True)

        # Peak 2x2 window mean ZIP score (most synergistic dose region)
        z_arr = zip_piv.to_numpy(dtype=float)
        if z_arr.shape[0] >= 2 and z_arr.shape[1] >= 2:
            win_scores = []
            for r_i in range(z_arr.shape[0] - 1):
                for c_i in range(z_arr.shape[1] - 1):
                    patch = z_arr[r_i : r_i + 2, c_i : c_i + 2]
                    if np.isfinite(patch).all():
                        win_scores.append(float(np.mean(patch)))
            peak_win_zip = max(win_scores) if win_scores else float(np.nanmax(z_arr))
        else:
            peak_win_zip = float(np.nanmax(z_arr))

        mean_zip = float(sub_df["zip_delta"].mean())
        mean_bliss = float(sub_df["bliss_delta"].mean())
        mean_loewe = float(sub_df["loewe_delta"].mean())
        mean_hsa = float(sub_df["hsa_delta"].mean())
        summary_rows.append(
            {
                "sample": samp,
                "drug_a": drug_a,
                "drug_b": drug_b,
                "n_combo_wells": int(len(sub_df)),
                "mean_zip": round(mean_zip, 2),
                "peak_window_zip": round(peak_win_zip, 2),
                "mean_bliss": round(mean_bliss, 2),
                "mean_loewe": round(mean_loewe, 2),
                "mean_hsa": round(mean_hsa, 2),
                "drug_a_ec50": round(fit_a["ec50"], 5),
                "drug_b_ec50": round(fit_b["ec50"], 5),
                "classification": _classify_synergy(mean_zip),
            }
        )
        if first_plot_data is None:
            first_plot_data = (samp, drug_a, drug_b, inh_mat, zip_piv, mean_zip, peak_win_zip)

    sum_df = pd.DataFrame(summary_rows)
    grid_df = pd.DataFrame(grid_rows)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sum_path = out_dir / "synergy_summary.tsv"
    grid_path = out_dir / "synergy_grid.tsv"
    fig_path = out_dir / "drug-synergy.png"

    sum_df.to_csv(sum_path, sep="\t", index=False)
    grid_df.to_csv(grid_path, sep="\t", index=False)

    write_artifact(
        sum_path,
        "pgx_synergy",
        "Drug Combination Synergy Summary (ZIP / Bliss / Loewe / HSA)",
        inputs=[args.input],
        parameters={"dose_unit": args.dose_unit},
        shape=(int(len(sum_df)), int(sum_df.shape[1])),
    )
    write_artifact(
        grid_path,
        "pgx_synergy",
        "2D Dose-Combination Inhibition & Synergy Delta Matrix",
        inputs=[args.input],
        parameters={"dose_unit": args.dose_unit},
        shape=(int(len(grid_df)), int(grid_df.shape[1])),
    )

    # Render 2-panel dose-response inhibition + ZIP delta synergy heatmap
    samp, drug_a, drug_b, inh_mat, zip_piv, mean_zip, peak_win_zip = first_plot_data
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 5.0), dpi=180)

    im1 = ax1.imshow(inh_mat.to_numpy(dtype=float), cmap="YlOrRd", vmin=0, vmax=100, origin="lower", aspect="auto")
    ax1.set_xticks(range(inh_mat.shape[1]))
    ax1.set_xticklabels([f"{c:g}" for c in inh_mat.columns], rotation=35, ha="right", fontsize=8)
    ax1.set_yticks(range(inh_mat.shape[0]))
    ax1.set_yticklabels([f"{r:g}" for r in inh_mat.index], fontsize=8)
    for i in range(inh_mat.shape[0]):
        for j in range(inh_mat.shape[1]):
            v = float(inh_mat.iloc[i, j])
            if math.isfinite(v):
                ax1.text(j, i, f"{v:.0f}%", ha="center", va="center", fontsize=7.5, color="white" if v > 60 else "#111827")
    ax1.set_xlabel(f"{drug_a} ({args.dose_unit})")
    ax1.set_ylabel(f"{drug_b} ({args.dose_unit})")
    ax1.set_title(f"{samp}: Observed Inhibition (%)")
    fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04, label="% Inhibition")

    z_mat = zip_piv.to_numpy(dtype=float)
    vmax_z = max(25.0, float(np.nanmax(np.abs(z_mat))))
    im2 = ax2.imshow(z_mat, cmap="RdBu_r", vmin=-vmax_z, vmax=vmax_z, origin="lower", aspect="auto")
    ax2.set_xticks(range(zip_piv.shape[1]))
    ax2.set_xticklabels([f"{c:g}" for c in zip_piv.columns], rotation=35, ha="right", fontsize=8)
    ax2.set_yticks(range(zip_piv.shape[0]))
    ax2.set_yticklabels([f"{r:g}" for r in zip_piv.index], fontsize=8)
    for i in range(zip_piv.shape[0]):
        for j in range(zip_piv.shape[1]):
            zv = float(zip_piv.iloc[i, j])
            if math.isfinite(zv):
                ax2.text(
                    j,
                    i,
                    f"{zv:+.1f}",
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    color="white" if abs(zv) > 0.55 * vmax_z else "#111827",
                )
    ax2.set_xlabel(f"{drug_a} ({args.dose_unit})")
    ax2.set_ylabel(f"{drug_b} ({args.dose_unit})")
    ax2.set_title(f"ZIP Synergy Landscape (Mean={mean_zip:+.1f}, Peak={peak_win_zip:+.1f})")
    fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04, label="ZIP Synergy Score (delta %)")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "pgx_synergy",
        "2D Drug-Combination Inhibition & ZIP Synergy Landscape",
        inputs=[args.input],
        parameters={"drug_a": drug_a, "drug_b": drug_b, "sample": samp},
    )

    emit(
        {
            "status": "ok",
            "tool": "pgx_synergy",
            "dose_unit": args.dose_unit,
            "detected_viability_scale": scale_name,
            "n_combinations": int(len(sum_df)),
            "combinations": sum_df.to_dict(orient="records"),
            "outputs": {
                "synergy_summary": str(sum_path),
                "synergy_grid": str(grid_path),
                "synergy_plot": str(fig_path),
            },
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
