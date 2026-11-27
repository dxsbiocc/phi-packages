#!/usr/bin/env python3
"""pgx_hts_qc: high-throughput screening plate QC (Z'-factor, SSMD, edge effects, Tukey B-score, hit calling)."""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, read_table, setup_matplotlib, write_artifact

WELL_TYPE_ALIASES = {
    "pos_ctrl": "pos_ctrl",
    "pos": "pos_ctrl",
    "positive": "pos_ctrl",
    "posctrl": "pos_ctrl",
    "positive_control": "pos_ctrl",
    "neg_ctrl": "neg_ctrl",
    "neg": "neg_ctrl",
    "negative": "neg_ctrl",
    "negctrl": "neg_ctrl",
    "negative_control": "neg_ctrl",
    "dmso": "neg_ctrl",
    "vehicle": "neg_ctrl",
    "sample": "sample",
    "compound": "sample",
    "test": "sample",
    "cmpd": "sample",
    "blank": "blank",
    "empty": "blank",
}


def _row_letters_to_idx(lbl: str) -> int:
    s = re.sub(r"[^A-Za-z]", "", str(lbl).strip().upper())
    if not s:
        return -1
    val = 0
    for ch in s:
        val = val * 26 + (ord(ch) - ord("A") + 1)
    return val


def _idx_to_row_letters(idx: int) -> str:
    n = int(idx)
    chars = []
    while n > 0:
        n, rem = divmod(n - 1, 26)
        chars.append(chr(ord("A") + rem))
    return "".join(reversed(chars)) or "A"


def _parse_well(well_str: str) -> tuple[int, int]:
    m = re.match(r"^\s*([A-Za-z]+)\s*0*(\d+)\s*$", str(well_str))
    if not m:
        return -1, -1
    return _row_letters_to_idx(m.group(1)), int(m.group(2))


def _tukey_median_polish(mat, max_iter: int = 30, tol: float = 1e-6):
    """Two-way Tukey median polish on a 2D array with NaNs; returns (overall, row_eff, col_eff, resid)."""
    import numpy as np

    resid = np.array(mat, dtype=float, copy=True)
    n_rows, n_cols = resid.shape
    overall = 0.0
    row_eff = np.zeros(n_rows, dtype=float)
    col_eff = np.zeros(n_cols, dtype=float)

    for _ in range(max_iter):
        delta_max = 0.0
        for r in range(n_rows):
            finite = np.isfinite(resid[r, :])
            if finite.any():
                m = float(np.median(resid[r, finite]))
                resid[r, :] -= m
                row_eff[r] += m
                delta_max = max(delta_max, abs(m))
        finite_r = np.isfinite(row_eff)
        if finite_r.any():
            r_med = float(np.median(row_eff[finite_r]))
            row_eff -= r_med
            overall += r_med

        for c in range(n_cols):
            finite = np.isfinite(resid[:, c])
            if finite.any():
                m = float(np.median(resid[:, c]))
                resid[:, c] -= m
                col_eff[c] += m
                delta_max = max(delta_max, abs(m))
        finite_c = np.isfinite(col_eff)
        if finite_c.any():
            c_med = float(np.median(col_eff[finite_c]))
            col_eff -= c_med
            overall += c_med

        if delta_max < tol:
            break

    return overall, row_eff, col_eff, resid


def _pick_col(cols: list[str], candidates: list[str], *, required: bool = True) -> str | None:
    lower = {c.lower().strip(): c for c in cols}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    if required:
        fail(f"Could not find column matching {candidates} in input columns: {cols}")
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description="High-throughput screening plate QC and hit calling")
    ap.add_argument("--plates", required=True, help="Plate screening TSV/CSV")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--apply_bscore", action="store_true", help="Use Tukey median polish B-score for hit calling")
    ap.add_argument("--zscore_hit_cutoff", type=float, default=3.0, help="Robust MAD |z| cutoff for hit calling (default 3.0)")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from scipy import stats

    raw = read_table(args.plates)
    if raw.empty:
        fail("HTS plates input table is empty.")

    cols = list(raw.columns)
    plate_col = _pick_col(cols, ["plate", "plate_id", "barcode"], required=False)
    well_col = _pick_col(cols, ["well", "well_id", "position"], required=False)
    row_col = _pick_col(cols, ["row", "row_idx", "row_id"], required=False)
    col_col = _pick_col(cols, ["col", "column", "col_idx", "column_id"], required=False)
    sig_col = _pick_col(cols, ["signal", "value", "readout", "intensity", "luminescence", "viability", "raw_signal"])
    wt_col = _pick_col(cols, ["well_type", "type", "role", "control_type", "sample_type"])

    if not well_col and not (row_col and col_col):
        fail("Input must contain either 'well' (e.g. A01..P24) or both 'row' and 'col' columns.")

    df = raw.copy()
    df["_plate"] = df[plate_col].astype(str).str.strip() if plate_col else "Plate1"
    df["_signal"] = pd.to_numeric(df[sig_col], errors="coerce")
    df = df.dropna(subset=["_signal"]).copy()
    if df.empty:
        fail("No finite numeric signal values found in HTS plate file.")

    # Parse row & col indices (1-based)
    rows_idx = []
    cols_idx = []
    well_labels = []
    for _, r in df.iterrows():
        ri, ci = -1, -1
        if well_col:
            ri, ci = _parse_well(str(r[well_col]))
        if (ri <= 0 or ci <= 0) and row_col and col_col:
            rv = r[row_col]
            try:
                ri = int(float(rv))
            except Exception:
                ri = _row_letters_to_idx(str(rv))
            try:
                ci = int(float(r[col_col]))
            except Exception:
                ci = -1
        if ri <= 0 or ci <= 0:
            fail(f"Could not parse valid row/column coordinates for row: {r.to_dict()}")
        rows_idx.append(ri)
        cols_idx.append(ci)
        well_labels.append(str(r[well_col]).strip() if well_col else f"{_idx_to_row_letters(ri)}{ci:02d}")

    df["_row"] = rows_idx
    df["_col"] = cols_idx
    df["_well"] = well_labels

    norm_wt = []
    for v in df[wt_col].astype(str):
        key = v.strip().lower().replace("-", "_").replace(" ", "_")
        if key not in WELL_TYPE_ALIASES:
            fail(f"Unsupported well_type '{v}'; expected one of pos_ctrl, neg_ctrl, sample, blank.")
        norm_wt.append(WELL_TYPE_ALIASES[key])
    df["_well_type"] = norm_wt

    warnings: list[str] = []
    plate_summaries: list[dict] = []
    scored_frames: list[pd.DataFrame] = []
    plate_grids: dict[str, tuple[np.ndarray, float]] = {}

    for plate_id, p_df in df.groupby("_plate", sort=False):
        p_df = p_df.copy()
        pos_sig = p_df.loc[p_df["_well_type"] == "pos_ctrl", "_signal"].to_numpy(dtype=float, copy=True)
        neg_sig = p_df.loc[p_df["_well_type"] == "neg_ctrl", "_signal"].to_numpy(dtype=float, copy=True)
        samp_mask = (p_df["_well_type"] == "sample").to_numpy(dtype=bool, copy=True)
        samp_sig = p_df.loc[samp_mask, "_signal"].to_numpy(dtype=float, copy=True)
        n_blank = int((p_df["_well_type"] == "blank").sum())

        if pos_sig.size < 2 or neg_sig.size < 2:
            fail(f"Plate '{plate_id}' requires >= 2 pos_ctrl and >= 2 neg_ctrl wells for Z'-factor calculation.")

        mu_pos = float(np.mean(pos_sig))
        sd_pos = float(np.std(pos_sig, ddof=1))
        mu_neg = float(np.mean(neg_sig))
        sd_neg = float(np.std(neg_sig, ddof=1))
        cv_pos = 100.0 * sd_pos / max(1e-12, abs(mu_pos))
        cv_neg = 100.0 * sd_neg / max(1e-12, abs(mu_neg))

        sep = abs(mu_pos - mu_neg)
        if sep <= 1e-12:
            z_prime = -10.0
            ssmd = 0.0
        else:
            z_prime = 1.0 - 3.0 * (sd_pos + sd_neg) / sep
            denom_ssmd = math.sqrt(sd_pos ** 2 + sd_neg ** 2)
            ssmd = (mu_pos - mu_neg) / max(1e-12, denom_ssmd)

        z_grade = "excellent" if z_prime >= 0.5 else ("marginal" if z_prime >= 0.0 else "failed")
        sb_ratio = max(abs(mu_pos), abs(mu_neg)) / max(1e-9, min(abs(mu_pos), abs(mu_neg)))

        max_r = int(p_df["_row"].max())
        max_c = int(p_df["_col"].max())
        if max_r <= 8 and max_c <= 12:
            n_r, n_c, fmt = 8, 12, "96-well"
        elif max_r <= 16 and max_c <= 24:
            n_r, n_c, fmt = 16, 24, "384-well"
        else:
            n_r, n_c, fmt = max_r, max_c, f"{max_r}x{max_c}"

        # Row/col ANOVA and edge-effect diagnostics on sample wells (or neg_ctrl+sample if few samples)
        diag_df = p_df[p_df["_well_type"] == "sample"]
        if len(diag_df) < 6:
            diag_df = p_df[p_df["_well_type"].isin(["sample", "neg_ctrl"])]

        diag_r = diag_df["_row"].to_numpy(dtype=float, copy=True)
        diag_c = diag_df["_col"].to_numpy(dtype=float, copy=True)
        diag_s = diag_df["_signal"].to_numpy(dtype=float, copy=True)

        row_groups = [g["_signal"].to_numpy(dtype=float, copy=True) for _, g in diag_df.groupby("_row") if len(g) >= 2]
        col_groups = [g["_signal"].to_numpy(dtype=float, copy=True) for _, g in diag_df.groupby("_col") if len(g) >= 2]
        try:
            kw_r = float(stats.kruskal(*row_groups).pvalue) if len(row_groups) >= 2 else 1.0
            sp_r = float(stats.spearmanr(diag_r, diag_s).pvalue) if len(diag_s) >= 6 else 1.0
            row_anova_p = min(kw_r if math.isfinite(kw_r) else 1.0, sp_r if math.isfinite(sp_r) else 1.0)
        except Exception:
            row_anova_p = 1.0
        try:
            kw_c = float(stats.kruskal(*col_groups).pvalue) if len(col_groups) >= 2 else 1.0
            sp_c = float(stats.spearmanr(diag_c, diag_s).pvalue) if len(diag_s) >= 6 else 1.0
            col_anova_p = min(kw_c if math.isfinite(kw_c) else 1.0, sp_c if math.isfinite(sp_c) else 1.0)
        except Exception:
            col_anova_p = 1.0

        edge_mask = (diag_r == 1) | (diag_r == n_r) | (diag_c == 1) | (diag_c == n_c)
        if edge_mask.sum() >= 2 and (~edge_mask).sum() >= 2:
            edge_vals = diag_s[edge_mask]
            int_vals = diag_s[~edge_mask]
            try:
                mw_p = float(stats.mannwhitneyu(edge_vals, int_vals, alternative="two-sided").pvalue)
                tt_p = float(stats.ttest_ind(edge_vals, int_vals, equal_var=False).pvalue)
                edge_p = min(mw_p if math.isfinite(mw_p) else 1.0, tt_p if math.isfinite(tt_p) else 1.0)
            except Exception:
                edge_p = 1.0
            edge_bias_pct = 100.0 * (float(np.median(edge_vals)) - float(np.median(int_vals))) / max(1e-9, abs(float(np.median(int_vals))))
        else:
            edge_p = 1.0
            edge_bias_pct = 0.0

        spatial_bias = bool(row_anova_p < 0.01 or col_anova_p < 0.01 or edge_p < 0.01)

        # Build 2D grid of sample wells (controls masked NaN) for Tukey median polish B-score
        sample_grid = np.full((n_r, n_c), np.nan, dtype=float)
        full_grid = np.full((n_r, n_c), np.nan, dtype=float)
        for _, rw in p_df.iterrows():
            r_i = int(rw["_row"]) - 1
            c_i = int(rw["_col"]) - 1
            full_grid[r_i, c_i] = float(rw["_signal"])
            if rw["_well_type"] == "sample":
                sample_grid[r_i, c_i] = float(rw["_signal"])

        if np.isfinite(sample_grid).sum() >= 6:
            overall, row_eff, col_eff, _ = _tukey_median_polish(sample_grid)
        else:
            overall, row_eff, col_eff, _ = _tukey_median_polish(full_grid)

        resid_all = []
        for _, rw in p_df.iterrows():
            r_i = int(rw["_row"]) - 1
            c_i = int(rw["_col"]) - 1
            resid_all.append(float(rw["_signal"]) - (overall + row_eff[r_i] + col_eff[c_i]))
        resid_arr = np.asarray(resid_all, dtype=float)
        ref_resid = resid_arr[samp_mask] if samp_mask.sum() >= 3 else resid_arr
        med_resid = float(np.median(ref_resid))
        mad_resid = 1.4826 * float(np.median(np.abs(ref_resid - med_resid)))
        if mad_resid <= 1e-12:
            mad_resid = float(np.std(ref_resid, ddof=1)) if ref_resid.size > 1 else 1.0
        mad_resid = max(mad_resid, 1e-12)
        b_scores = (resid_arr - med_resid) / mad_resid

        # Raw MAD z_robust on sample wells
        ref_raw = samp_sig if samp_sig.size >= 3 else p_df["_signal"].to_numpy(dtype=float, copy=True)
        med_raw = float(np.median(ref_raw))
        mad_raw = 1.4826 * float(np.median(np.abs(ref_raw - med_raw)))
        if mad_raw <= 1e-12:
            mad_raw = float(np.std(ref_raw, ddof=1)) if ref_raw.size > 1 else 1.0
        mad_raw = max(mad_raw, 1e-12)
        z_raw = (p_df["_signal"].to_numpy(dtype=float, copy=True) - med_raw) / mad_raw

        z_final = b_scores if args.apply_bscore else z_raw
        if abs(mu_neg - mu_pos) > 1e-12:
            pct_inhib = 100.0 * (mu_neg - p_df["_signal"].to_numpy(dtype=float, copy=True)) / (mu_neg - mu_pos)
        else:
            pct_inhib = np.zeros(len(p_df), dtype=float)

        is_hit = samp_mask & (np.abs(z_final) >= float(args.zscore_hit_cutoff))
        n_hits = int(is_hit.sum())
        n_samp = int(samp_mask.sum())

        if z_prime < 0.0:
            warnings.append(f"Plate '{plate_id}' failed assay QC (Z' = {z_prime:.3f} < 0).")
        elif z_prime < 0.5:
            warnings.append(f"Plate '{plate_id}' has marginal assay window (0 <= Z' = {z_prime:.3f} < 0.5).")
        if spatial_bias and not args.apply_bscore:
            warnings.append(
                f"Plate '{plate_id}' exhibits spatial row/column/edge bias (row_p={row_anova_p:.2e}, "
                f"col_p={col_anova_p:.2e}, edge_p={edge_p:.2e}); consider re-running with --apply_bscore."
            )

        p_out = pd.DataFrame(
            {
                "plate": str(plate_id),
                "well": p_df["_well"].to_numpy(copy=True),
                "row": p_df["_row"].to_numpy(copy=True),
                "col": p_df["_col"].to_numpy(copy=True),
                "well_type": p_df["_well_type"].to_numpy(copy=True),
                "signal": np.round(p_df["_signal"].to_numpy(dtype=float, copy=True), 4),
                "pct_inhibition": np.round(pct_inhib, 3),
                "z_robust_raw": np.round(z_raw, 4),
                "b_score": np.round(b_scores, 4),
                "z_robust": np.round(z_final, 4),
                "is_hit": is_hit,
            }
        )
        scored_frames.append(p_out)

        # Grid for heatmap visualization (B-score if apply_bscore else % inhibition)
        vis_grid = np.full((n_r, n_c), np.nan, dtype=float)
        for _, rw in p_out.iterrows():
            vis_grid[int(rw["row"]) - 1, int(rw["col"]) - 1] = float(rw["b_score"] if args.apply_bscore else rw["pct_inhibition"])
        plate_grids[str(plate_id)] = (vis_grid, float(z_prime))

        plate_summaries.append(
            {
                "plate": str(plate_id),
                "plate_format": fmt,
                "n_wells": int(len(p_df)),
                "n_pos_ctrl": int(pos_sig.size),
                "n_neg_ctrl": int(neg_sig.size),
                "n_sample": n_samp,
                "n_blank": n_blank,
                "mean_pos": round(mu_pos, 4),
                "sd_pos": round(sd_pos, 4),
                "cv_pos_pct": round(cv_pos, 2),
                "mean_neg": round(mu_neg, 4),
                "sd_neg": round(sd_neg, 4),
                "cv_neg_pct": round(cv_neg, 2),
                "z_prime": round(float(z_prime), 4),
                "z_prime_grade": z_grade,
                "ssmd": round(float(ssmd), 4),
                "s_b_ratio": round(float(sb_ratio), 4),
                "row_anova_p": float(row_anova_p),
                "col_anova_p": float(col_anova_p),
                "edge_vs_interior_p": float(edge_p),
                "edge_bias_pct": round(float(edge_bias_pct), 2),
                "spatial_bias_detected": spatial_bias,
                "n_hits": n_hits,
                "hit_rate_pct": round(100.0 * n_hits / max(1, n_samp), 2),
            }
        )

    summary_df = pd.DataFrame(plate_summaries)
    wells_df = pd.concat(scored_frames, ignore_index=True)

    hits_df = wells_df[wells_df["is_hit"]].copy()
    if not hits_df.empty:
        hits_df["_abs_z"] = hits_df["z_robust"].abs()
        hits_df = hits_df.sort_values("_abs_z", ascending=False).drop(columns=["_abs_z"])

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_path = out_dir / "hts_plate_qc_summary.tsv"
    well_path = out_dir / "hts_well_scored.tsv"
    fig_path = out_dir / "hts_plate_heatmap_qc.png"

    summary_df.to_csv(qc_path, sep="\t", index=False)
    wells_df.to_csv(well_path, sep="\t", index=False)

    params = {"apply_bscore": bool(args.apply_bscore), "zscore_hit_cutoff": float(args.zscore_hit_cutoff)}
    write_artifact(
        qc_path,
        "pgx_hts_qc",
        "High-Throughput Screen Per-Plate QC Summary (Z'-factor, SSMD, Spatial Bias)",
        inputs=[args.plates],
        parameters=params,
        shape=summary_df.shape,
    )
    write_artifact(
        well_path,
        "pgx_hts_qc",
        "High-Throughput Screen Well-Level Normalized & B-Score Hit Calls",
        inputs=[args.plates],
        parameters=params,
        shape=wells_df.shape,
        scientific_context={
            "dataType": "drug_response",
            "valueScale": "b_score" if args.apply_bscore else "percent_inhibition",
        },
    )

    # Render spatial heatmap + Z'-factor summary figure
    plt = setup_matplotlib()
    show_plates = list(plate_grids.keys())[:3]
    n_panels = len(show_plates) + 1
    fig, axes = plt.subplots(1, n_panels, figsize=(4.6 * n_panels, 4.4), dpi=180)
    if n_panels == 1:
        axes = [axes]

    for idx, pid in enumerate(show_plates):
        ax = axes[idx]
        grid, zp = plate_grids[pid]
        cmap = "RdBu_r" if args.apply_bscore else "viridis"
        im = ax.imshow(grid, cmap=cmap, aspect="auto")
        ax.set_title(f"{pid} (Z'={zp:.2f})", fontsize=10)
        ax.set_xlabel("Column")
        ax.set_ylabel("Row")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="B-score" if args.apply_bscore else "% Inhibition")

    ax_bar = axes[-1]
    y_pos = np.arange(len(summary_df))
    z_vals = summary_df["z_prime"].to_numpy(dtype=float, copy=True)
    bar_cols = ["#059669" if z >= 0.5 else ("#f59e0b" if z >= 0.0 else "#dc2626") for z in z_vals]
    ax_bar.barh(y_pos, np.clip(z_vals, -0.5, 1.0), color=bar_cols, alpha=0.88)
    ax_bar.axvline(0.5, color="#059669", ls="--", lw=1.0, label="Excellent (0.5)")
    ax_bar.axvline(0.0, color="#dc2626", ls=":", lw=1.0, label="Fail (<0)")
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels(summary_df["plate"].tolist(), fontsize=8.5)
    ax_bar.invert_yaxis()
    ax_bar.set_xlim(-0.5, 1.05)
    ax_bar.set_xlabel("Zhang's Z'-factor")
    ax_bar.set_title("Per-Plate Assay Quality")
    ax_bar.legend(loc="lower right", frameon=False, fontsize=7.8)

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "pgx_hts_qc",
        "High-Throughput Screening Plate Spatial Heatmaps & Z'-Factor QC",
        inputs=[args.plates],
        parameters=params,
    )

    n_pass = int((summary_df["z_prime_grade"] == "excellent").sum())
    n_marg = int((summary_df["z_prime_grade"] == "marginal").sum())
    n_fail = int((summary_df["z_prime_grade"] == "failed").sum())

    emit(
        {
            "status": "ok",
            "tool": "pgx_hts_qc",
            "apply_bscore": bool(args.apply_bscore),
            "zscore_hit_cutoff": float(args.zscore_hit_cutoff),
            "n_plates": int(len(summary_df)),
            "n_wells_total": int(len(wells_df)),
            "n_passed_plates": n_pass,
            "n_marginal_plates": n_marg,
            "n_failed_plates": n_fail,
            "n_hits_total": int(len(hits_df)),
            "plate_summaries": plate_summaries,
            "top_hits": hits_df.head(15).to_dict(orient="records"),
            "outputs": {
                "plate_qc_summary": str(qc_path),
                "well_scored": str(well_path),
                "plate_heatmap_qc": str(fig_path),
            },
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
