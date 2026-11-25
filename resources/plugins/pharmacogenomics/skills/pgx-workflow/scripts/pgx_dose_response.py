#!/usr/bin/env python3
"""Fit 4-parameter log-logistic Hill dose-response curves with IC50 censoring guards."""
from __future__ import annotations

import argparse
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


def main() -> None:
    ap = argparse.ArgumentParser(description="Fit 4PL Hill dose-response curves")
    ap.add_argument("--input", required=True, help="Long-format dose-viability TSV/CSV")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--compound_col", default="", help="Compound/drug column name")
    ap.add_argument("--sample_col", default="", help="Sample/cell-line column name")
    ap.add_argument("--dose_col", default="", help="Dose/concentration column name")
    ap.add_argument("--viability_col", default="", help="Viability column name")
    ap.add_argument("--viability_scale", choices=["auto", "fraction", "percent"], default="auto")
    ap.add_argument("--dose_unit", default="uM", help="Concentration unit label (default uM)")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    df = read_table(args.input)
    cols = list(df.columns)
    dose_col = _pick_col(cols, args.dose_col, ["dose", "conc", "concentration", "drug_conc"])
    viab_col = _pick_col(cols, args.viability_col, ["viability", "cell_viability", "response", "survival", "pct_viability"])
    comp_col = _pick_col(cols, args.compound_col, ["compound", "drug", "drug_name", "agent", "inhibitor"], required=False)
    samp_col = _pick_col(cols, args.sample_col, ["sample", "cell_line", "cellline", "model", "depmap_id"], required=False)

    df = df.copy()
    df["_dose"] = pd.to_numeric(df[dose_col], errors="coerce")
    df["_viab"] = pd.to_numeric(df[viab_col], errors="coerce")
    df = df.dropna(subset=["_dose", "_viab"])
    if df.empty:
        fail("No valid numeric (dose, viability) rows found.")
    if (df["_dose"] < 0).any():
        fail("Negative drug concentrations encountered; doses must be >= 0.")

    scale_name, scale_factor = detect_viability_scale(df["_viab"].to_numpy(), args.viability_scale)
    df["_viab_frac"] = df["_viab"] * scale_factor
    df["_compound"] = df[comp_col].astype(str) if comp_col else "Compound"
    df["_sample"] = df[samp_col].astype(str) if samp_col else "Sample"

    warnings: list[str] = []
    fit_rows: list[dict] = []
    plot_series: list[dict] = []

    for (comp, samp), grp in df.groupby(["_compound", "_sample"], sort=False):
        veh = grp[grp["_dose"] == 0.0]["_viab_frac"]
        pos = grp[grp["_dose"] > 0.0].copy()
        if len(pos) < 4 or pos["_dose"].nunique() < 3:
            warnings.append(
                f"Skipped ({comp}, {samp}): fewer than 4 positive points or 3 distinct concentrations."
            )
            continue

        if len(veh) > 0 and float(veh.mean()) > 0.1:
            pos["_viab_norm"] = pos["_viab_frac"] / float(veh.mean())
            fix_e0 = 1.0
        else:
            pos["_viab_norm"] = pos["_viab_frac"]
            fix_e0 = None

        res = fit_hill_4pl(
            pos["_dose"].to_numpy(dtype=float),
            pos["_viab_norm"].to_numpy(dtype=float),
            fix_e0=fix_e0,
        )
        label = f"{comp} | {samp}" if (comp_col and samp_col) else (comp if comp_col else samp)
        ic50_display = f">{res['max_dose']:.4g}" if res["beyond_tested_range"] else f"{res['ic50_raw']:.4g}"
        row = {
            "compound": comp,
            "sample": samp,
            "dose_unit": args.dose_unit,
            "n_points": res["n_points"],
            "n_doses": res["n_doses"],
            "min_dose": round(res["min_dose"], 6),
            "max_dose": round(res["max_dose"], 6),
            "e0": round(res["e0"], 4),
            "einf": round(res["einf"], 4),
            "emax_inhibition": round(res["emax_inhibition"], 4),
            "hill_slope": round(res["hill"], 4),
            "ec50": round(res["ec50"], 6),
            "ic50": round(res["ic50_raw"], 6) if res["ic50_raw"] is not None else np.nan,
            "ic50_censored": round(res["ic50_censored"], 6),
            "ln_ic50_censored": round(res["ln_ic50_censored"], 5),
            "ic50_report": ic50_display,
            "beyond_tested_range": res["beyond_tested_range"],
            "auc": round(res["auc"], 4),
            "aac": round(res["aac"], 4),
            "r2": round(res["r2"], 4),
        }
        fit_rows.append(row)
        plot_series.append(
            {
                "label": label,
                "doses": pos["_dose"].to_numpy(dtype=float),
                "viab": pos["_viab_norm"].to_numpy(dtype=float),
                "fit": res,
                "ic50_report": ic50_display,
            }
        )

    if not fit_rows:
        fail("No curves had sufficient dose points (>= 4 points across >= 3 positive concentrations) to fit.")

    fits_df = pd.DataFrame(fit_rows)
    n_censored = int(fits_df["beyond_tested_range"].sum())
    if n_censored > 0:
        warnings.append(
            f"{n_censored}/{len(fits_df)} curve(s) did not reach 50% viability loss within the tested concentration "
            "range (beyond_tested_range=true); report IC50 as '> max_dose' and compare lines via AUC/AAC."
        )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fits_path = out_dir / "dose_response_fits.tsv"
    fig_path = out_dir / "dose-response-curves.png"
    fits_df.to_csv(fits_path, sep="\t", index=False)

    write_artifact(
        fits_path,
        "pgx_dose_response",
        "4-Parameter Hill Dose-Response Curve Parameters & Censoring Status",
        inputs=[args.input],
        parameters={"dose_unit": args.dose_unit, "viability_scale": scale_name},
        shape=(int(len(fits_df)), int(fits_df.shape[1])),
        scientific_context={
            "dataType": "drug_response",
            "valueScale": "ic50_and_auc",
            "doseUnit": args.dose_unit,
            "censoringField": "beyond_tested_range",
        },
    )

    # Plot up to 12 curves cleanly
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(8.4, 5.2), dpi=180)
    palette = [
        "#2563eb", "#dc2626", "#059669", "#d97706", "#7c3aed", "#0891b2",
        "#db2777", "#4f46e5", "#65a30d", "#ea580c", "#0d9488", "#4b5563",
    ]
    shown = plot_series[:12]
    for i, item in enumerate(shown):
        col = palette[i % len(palette)]
        d_arr = item["doses"]
        v_arr = item["viab"]
        f_res = item["fit"]
        grid_ld = np.linspace(np.log10(f_res["min_dose"]) - 0.1, np.log10(f_res["max_dose"]) + 0.1, 150)
        grid_y = hill_4pl_log10(grid_ld, f_res["e0"], f_res["einf"], f_res["log10_ec50"], f_res["hill"])
        lbl = f"{item['label']} (IC50={item['ic50_report']} {args.dose_unit}, AUC={f_res['auc']:.2f})"
        ax.scatter(d_arr, v_arr, s=22, color=col, alpha=0.75)
        ax.plot(10.0 ** grid_ld, grid_y, lw=1.8, color=col, label=lbl)
        if not f_res["beyond_tested_range"] and f_res["ic50_raw"] is not None:
            ax.scatter([f_res["ic50_raw"]], [0.50], marker="D", s=30, color=col, zorder=5)

    ax.axhline(0.50, color="#6b7280", ls="--", lw=1.0, label="50% Viability (IC50)")
    ax.set_xscale("log")
    ax.set_ylim(-0.05, 1.25)
    ax.set_xlabel(f"Drug Concentration ({args.dose_unit}, log scale)")
    ax.set_ylabel("Normalized Fractional Viability")
    ax.set_title("4-Parameter Log-Logistic Dose-Response Curves")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=7.8)
    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "pgx_dose_response",
        "4PL Hill Dose-Response Curves",
        inputs=[args.input],
        parameters={"dose_unit": args.dose_unit, "n_curves_shown": len(shown)},
    )

    emit(
        {
            "status": "ok",
            "tool": "pgx_dose_response",
            "dose_unit": args.dose_unit,
            "detected_viability_scale": scale_name,
            "n_curves_fitted": int(len(fits_df)),
            "n_beyond_tested_range": n_censored,
            "median_r2": round(float(fits_df["r2"].median()), 4),
            "curves": fits_df.head(25).replace({np.nan: None}).to_dict(orient="records"),
            "outputs": {
                "dose_response_fits": str(fits_path),
                "curves_plot": str(fig_path),
            },
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
