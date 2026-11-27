#!/usr/bin/env python3
"""gwas_finemap: Multi-causal SuSiE-RSS statistical fine-mapping on summary
Z-scores and regularized LD matrix R_tilde = (1 - lambda)R + lambda*I (or exact
single-causal L=1 Wakefield ABF fallback when ld_matrix is omitted).
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
from scipy import optimize, special, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    parse_sumstats,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)


def _load_ld_matrix(path: str | Path, expected_m: int) -> np.ndarray:
    """Load M x M LD correlation matrix from whitespace, TSV, or CSV file."""
    p = Path(path)
    if not p.is_file():
        fail(f"LD matrix file not found: {path}")

    low = p.name.lower()
    sep = "," if (low.endswith(".csv") or low.endswith(".csv.gz")) else None
    try:
        raw_df = pd.read_csv(p, sep=sep if sep else r"\s+", header=None, low_memory=False)
        if raw_df.shape[1] == 1 and sep is None:
            raw_df = pd.read_csv(p, sep="\t", header=None, low_memory=False)
    except Exception as exc:
        fail(f"Failed to parse ld_matrix '{path}': {exc}")

    # Strip optional header row and/or SNP ID index column
    if raw_df.shape[0] == expected_m + 1:
        row0_num = pd.to_numeric(raw_df.iloc[0], errors="coerce")
        if row0_num.isna().sum() > expected_m // 2:
            raw_df = raw_df.iloc[1:].reset_index(drop=True)
    if raw_df.shape[1] == expected_m + 1:
        col0_num = pd.to_numeric(raw_df.iloc[:, 0], errors="coerce")
        if col0_num.isna().sum() > expected_m // 2:
            raw_df = raw_df.iloc[:, 1:].reset_index(drop=True)

    if raw_df.shape != (expected_m, expected_m):
        fail(f"LD matrix shape {raw_df.shape} does not match expected ({expected_m}, {expected_m}) from locus_sumstats")

    mat = raw_df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float, copy=True)
    if not np.all(np.isfinite(mat)):
        fail("LD matrix contains non-numeric or non-finite values")

    mat = 0.5 * (mat + mat.T)
    mat = np.clip(mat, -1.0, 1.0)
    np.fill_diagonal(mat, 1.0)
    return mat


def _regularize_ld(r_mat: np.ndarray) -> tuple[np.ndarray, float]:
    """Ensure PSD regularized LD matrix R_tilde = (1 - lam) * R + lam * I."""
    m = r_mat.shape[0]
    min_eig = float(np.min(np.linalg.eigvalsh(r_mat)))
    if min_eig >= 0.005:
        lam = 0.005
    else:
        lam = min(0.25, max(0.01, -min_eig + 0.01))
    r_reg = (1.0 - lam) * r_mat + lam * np.eye(m, dtype=float)
    np.fill_diagonal(r_reg, 1.0)
    return r_reg, float(lam)


def _wakefield_abf(z: np.ndarray, se: np.ndarray | None, n_arr: np.ndarray, coverage: float) -> dict[str, Any]:
    """Single-causal L=1 Wakefield (2009) Approximate Bayes Factor fine-mapping."""
    m = len(z)
    w = 0.04  # prior variance 0.2^2
    if se is not None and np.all(np.isfinite(se)) and np.all(se > 0):
        v_j = se * se
    else:
        v_j = 1.0 / np.maximum(n_arr, 1.0)
    r_j = w / (w + v_j)
    lbf = 0.5 * (np.log(1.0 - r_j) + (z * z) * r_j)
    log_post = lbf - special.logsumexp(lbf)
    pip = np.exp(log_post)

    order = np.argsort(-pip)
    cs_idx: list[int] = []
    cum = 0.0
    for idx in order:
        cs_idx.append(int(idx))
        cum += float(pip[idx])
        if cum >= coverage:
            break

    return {
        "pip": pip,
        "alpha": pip.reshape(1, m),
        "V": np.array([w]),
        "converged": True,
        "n_iterations": 1,
        "cs_list": [
            {
                "cs_id": 1,
                "indices": cs_idx,
                "coverage_actual": cum,
                "min_abs_corr": None,
                "mean_abs_corr": None,
                "median_abs_corr": None,
            }
        ],
    }


def _run_susie_rss(
    z: np.ndarray,
    r_reg: np.ndarray,
    r_raw: np.ndarray,
    n_med: float,
    max_causal: int,
    coverage: float,
    min_abs_corr: float,
    max_iter: int = 200,
    tol: float = 1e-4,
) -> dict[str, Any]:
    """SuSiE-RSS Iterative Bayesian Stepwise Selection on summary Z-scores and LD."""
    m = len(z)
    l_max = min(max_causal, m)
    d = float(max(n_med, 2.0))
    # PVE-scale Z adjustment (Zou et al. 2022)
    z_adj = z * np.sqrt(d / np.maximum(d - 1.0 + z * z, 1.0))

    xtx = d * r_reg
    xty = math.sqrt(d) * z_adj
    sigma2 = 1.0
    s2 = sigma2 / d
    log_pi = np.full(m, -math.log(m), dtype=float)

    alpha = np.full((l_max, m), 1.0 / m, dtype=float)
    mu = np.zeros((l_max, m), dtype=float)
    b = np.zeros((l_max, m), dtype=float)
    v_prior = np.full(l_max, 0.05, dtype=float)
    xtr = xty.copy()

    grid_lv = np.linspace(-14.0, 1.0, 16)
    converged = False
    n_iter = 0

    for it in range(1, max_iter + 1):
        n_iter = it
        alpha_prev = alpha.copy()

        for l_idx in range(l_max):
            g = xtr + xtx @ b[l_idx]
            betahat = g / d
            zhat = betahat / math.sqrt(s2)
            zhat2 = zhat * zhat

            def _neg_lbf_tot(lv: float) -> float:
                v_val = math.exp(lv)
                ratio = v_val / (s2 + v_val)
                lbf_vec = 0.5 * math.log(s2 / (s2 + v_val)) + 0.5 * zhat2 * ratio
                return -float(special.logsumexp(log_pi + lbf_vec))

            grid_vals = [_neg_lbf_tot(float(lv)) for lv in grid_lv]
            best_g = int(np.argmin(grid_vals))
            lv_lo = max(-16.0, float(grid_lv[max(0, best_g - 1)]))
            lv_hi = min(2.5, float(grid_lv[min(len(grid_lv) - 1, best_g + 1)]))
            opt = optimize.minimize_scalar(_neg_lbf_tot, bounds=(lv_lo, lv_hi), method="bounded")
            best_lbf_tot = -float(opt.fun)

            if best_lbf_tot <= 1e-6:
                v_prior[l_idx] = 0.0
                alpha[l_idx] = np.full(m, 1.0 / m, dtype=float)
                mu[l_idx] = 0.0
                b[l_idx] = 0.0
            else:
                v_opt = math.exp(float(opt.x))
                v_prior[l_idx] = v_opt
                ratio = v_opt / (s2 + v_opt)
                lbf_vec = 0.5 * math.log(s2 / (s2 + v_opt)) + 0.5 * zhat2 * ratio
                log_w = log_pi + lbf_vec
                alpha[l_idx] = np.exp(log_w - special.logsumexp(log_w))
                mu[l_idx] = ratio * betahat
                b[l_idx] = alpha[l_idx] * mu[l_idx]

            xtr = g - xtx @ b[l_idx]

        max_diff = float(np.max(np.abs(alpha - alpha_prev)))
        if it >= 5 and max_diff < tol:
            converged = True
            break

    pip = 1.0 - np.prod(1.0 - alpha, axis=0)
    pip = np.clip(pip, 0.0, 1.0)

    # Extract purity-filtered credible sets
    raw_cs: list[dict[str, Any]] = []
    for l_idx in range(l_max):
        if v_prior[l_idx] <= 1e-8:
            continue
        a_l = alpha[l_idx]
        order = np.argsort(-a_l)
        cs_idx: list[int] = []
        cum = 0.0
        for idx in order:
            cs_idx.append(int(idx))
            cum += float(a_l[idx])
            if cum >= coverage:
                break
        if len(cs_idx) == 1:
            min_r, mean_r, med_r = 1.0, 1.0, 1.0
        else:
            sub_r = np.abs(r_raw[np.ix_(cs_idx, cs_idx)])
            iu = np.triu_indices(len(cs_idx), k=1)
            pair_r = sub_r[iu]
            min_r = float(np.min(pair_r))
            mean_r = float(np.mean(pair_r))
            med_r = float(np.median(pair_r))
        if min_r >= min_abs_corr:
            lead_idx = cs_idx[0]
            raw_cs.append(
                {
                    "indices": cs_idx,
                    "coverage_actual": cum,
                    "min_abs_corr": min_r,
                    "mean_abs_corr": mean_r,
                    "median_abs_corr": med_r,
                    "lead_pip": float(pip[lead_idx]),
                }
            )

    raw_cs.sort(key=lambda c: c["lead_pip"], reverse=True)
    for k, cs in enumerate(raw_cs, start=1):
        cs["cs_id"] = k

    return {
        "pip": pip,
        "alpha": alpha,
        "V": v_prior,
        "converged": converged,
        "n_iterations": n_iter,
        "cs_list": raw_cs,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--locus_sumstats", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--ld_matrix", default=None)
    ap.add_argument("--max_causal", type=int, default=5)
    ap.add_argument("--coverage", type=float, default=0.95)
    ap.add_argument("--min_abs_corr", type=float, default=0.50)
    args = ap.parse_args()

    if not (1 <= args.max_causal <= 10):
        fail("--max_causal must be between 1 and 10")
    if not (0.5 <= args.coverage < 1.0):
        fail("--coverage must be in [0.5, 1.0)")
    if not (0.0 <= args.min_abs_corr <= 1.0):
        fail("--min_abs_corr must be in [0.0, 1.0]")

    df, _ = parse_sumstats(
        args.locus_sumstats,
        require_alleles=False,
        require_beta_se=False,
        require_n=True,
    )
    valid = np.isfinite(df["Z"].to_numpy(dtype=float, copy=True)) & np.isfinite(df["N"].to_numpy(dtype=float, copy=True)) & (df["N"].to_numpy(dtype=float, copy=True) > 1)
    if not valid.all():
        if args.ld_matrix is not None:
            fail("locus_sumstats contains non-finite Z or invalid N while ld_matrix requires exact row-order alignment")
        df = df.loc[valid].reset_index(drop=True)
    if df.empty:
        fail("No valid variants in locus_sumstats")

    m_snps = int(len(df))
    z_arr = df["Z"].to_numpy(dtype=float, copy=True)
    n_arr = df["N"].to_numpy(dtype=float, copy=True)
    se_arr = df["SE"].to_numpy(dtype=float, copy=True) if df["SE"].notna().all() else None
    n_med = float(np.median(n_arr))

    warnings: list[str] = []
    lam_reg: float | None = None
    r_raw: np.ndarray | None = None

    if args.ld_matrix:
        r_raw = _load_ld_matrix(args.ld_matrix, m_snps)
        r_reg, lam_reg = _regularize_ld(r_raw)
        res = _run_susie_rss(
            z_arr,
            r_reg,
            r_raw,
            n_med,
            args.max_causal,
            args.coverage,
            args.min_abs_corr,
        )
        method = "SuSiE-RSS"
        if not res["converged"]:
            warnings.append(f"SuSiE-RSS IBSS reached max_iter={res['n_iterations']} before strict tolerance.")
    else:
        res = _wakefield_abf(z_arr, se_arr, n_arr, args.coverage)
        method = "Wakefield_ABF_L1"
        warnings.append(
            "ld_matrix was omitted; performed single-causal (L=1) Wakefield ABF fine-mapping instead of multi-causal SuSiE-RSS."
        )

    pip = res["pip"]
    cs_list = res["cs_list"]

    # Assign cs_id per SNP (primary credible set if SNP appears in any)
    cs_id_col = np.zeros(m_snps, dtype=int)
    cs_rows: list[dict[str, Any]] = []
    snps_in_cs_set: set[int] = set()
    for cs in cs_list:
        cid = int(cs["cs_id"])
        idxs = cs["indices"]
        for i_snp in idxs:
            snps_in_cs_set.add(i_snp)
            if cs_id_col[i_snp] == 0:
                cs_id_col[i_snp] = cid
        lead_i = max(idxs, key=lambda idx: float(pip[idx]))
        cs_rows.append(
            {
                "cs_id": cid,
                "n_snps": len(idxs),
                "coverage_actual": round(float(cs["coverage_actual"]), 5),
                "min_abs_corr": round(float(cs["min_abs_corr"]), 4) if cs["min_abs_corr"] is not None else np.nan,
                "mean_abs_corr": round(float(cs["mean_abs_corr"]), 4) if cs["mean_abs_corr"] is not None else np.nan,
                "median_abs_corr": round(float(cs["median_abs_corr"]), 4) if cs["median_abs_corr"] is not None else np.nan,
                "lead_snp": str(df.loc[lead_i, "SNP"]),
                "lead_pos": int(df.loc[lead_i, "POS"]) if pd.notna(df.loc[lead_i, "POS"]) else np.nan,
                "lead_pip": round(float(pip[lead_i]), 5),
                "snps": ",".join(str(df.loc[i, "SNP"]) for i in idxs),
            }
        )

    pip_df = pd.DataFrame(
        {
            "SNP": df["SNP"],
            "CHR": df["CHR"],
            "POS": df["POS"],
            "Z": z_arr,
            "P": df["P"],
            "BETA": df["BETA"],
            "SE": df["SE"],
            "N": n_arr,
            "pip": np.round(pip, 6),
            "cs_id": cs_id_col,
        }
    )
    cs_df = pd.DataFrame(
        cs_rows,
        columns=[
            "cs_id",
            "n_snps",
            "coverage_actual",
            "min_abs_corr",
            "mean_abs_corr",
            "median_abs_corr",
            "lead_snp",
            "lead_pos",
            "lead_pip",
            "snps",
        ],
    )

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    pip_path = outdir / "finemap_pip.tsv"
    cs_path = outdir / "credible_sets.tsv"
    fig_path = outdir / "finemap_pip_locus.png"

    # LocusZoom-style 2-panel plot (-log10 P and PIP colored by credible set)
    plt = setup_matplotlib()
    fig, (ax_p, ax_pip) = plt.subplots(2, 1, figsize=(9.8, 6.2), sharex=True)

    has_pos = bool(df["POS"].notna().all())
    x_coord = df["POS"].to_numpy(dtype=float, copy=True) / 1e6 if has_pos else np.arange(1, m_snps + 1, dtype=float)
    nlp = -np.log10(np.clip(df["P"].to_numpy(dtype=float, copy=True), 1e-300, 1.0))

    palette = ["#E11D48", "#2563EB", "#059669", "#D97706", "#7C3AED", "#0891B2", "#BE185D", "#4F46E5", "#65A30D", "#EA580C"]

    bg = cs_id_col == 0
    ax_p.scatter(x_coord[bg], nlp[bg], s=22, color="#94A3B8", alpha=0.75, label="Outside CS")
    ax_pip.scatter(x_coord[bg], pip[bg], s=22, color="#94A3B8", alpha=0.75)

    for cid in range(1, len(cs_rows) + 1):
        mask = cs_id_col == cid
        col = palette[(cid - 1) % len(palette)]
        ax_p.scatter(x_coord[mask], nlp[mask], s=44, color=col, edgecolors="white", linewidths=0.6, zorder=4, label=f"CS {cid}")
        ax_pip.scatter(x_coord[mask], pip[mask], s=44, color=col, edgecolors="white", linewidths=0.6, zorder=4, label=f"CS {cid}")

    for row_cs in cs_rows:
        lead_mask = df["SNP"].astype(str) == row_cs["lead_snp"]
        if lead_mask.any():
            li = int(np.where(lead_mask.to_numpy())[0][0])
            ax_pip.annotate(
                f"{row_cs['lead_snp']} (PIP={pip[li]:.2f})",
                (x_coord[li], pip[li]),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                fontsize=7.6,
                color="#0F172A",
            )

    ax_p.axhline(-math.log10(5e-8), color="#64748B", linestyle="--", linewidth=0.9)
    ax_p.set_ylabel("−log10(P)")
    ax_p.set_title(f"Locus Fine-Mapping ({method}, {len(cs_rows)} 95% credible set(s))", fontsize=10.5)
    ax_p.legend(frameon=False, ncol=min(6, len(cs_rows) + 1), loc="upper right")

    ax_pip.set_ylabel("Posterior Inclusion Probability (PIP)")
    ax_pip.set_ylim(-0.04, 1.12)
    ax_pip.set_xlabel("Genomic Position (Mb)" if has_pos else "Variant Index")

    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    params_meta = {
        "max_causal": args.max_causal,
        "coverage": args.coverage,
        "min_abs_corr": args.min_abs_corr,
        "method": method,
        "ld_regularization_lambda": lam_reg,
    }
    in_files = [args.locus_sumstats, args.ld_matrix]
    artifacts = [
        write_tsv(pip_df, pip_path, "gwas_finemap", f"Per-SNP fine-mapping PIPs ({method})", inputs=in_files, parameters=params_meta),
        write_tsv(cs_df, cs_path, "gwas_finemap", f"Purity-filtered credible sets (K={len(cs_df)})", inputs=in_files, parameters=params_meta),
        write_artifact(fig_path, "gwas_finemap", "Locus fine-mapping PIP and -log10(P) plot", inputs=in_files, parameters=params_meta),
    ]

    top_pip_df = pip_df.sort_values("pip", ascending=False).head(15)
    emit(
        {
            "tool": "gwas_finemap",
            "method": method,
            "n_snps": m_snps,
            "sample_size": round(n_med, 1),
            "max_causal": args.max_causal,
            "coverage": args.coverage,
            "min_abs_corr": args.min_abs_corr,
            "ld_regularization_lambda": round(lam_reg, 5) if lam_reg is not None else None,
            "converged": bool(res["converged"]),
            "n_iterations": int(res["n_iterations"]),
            "n_credible_sets": int(len(cs_rows)),
            "n_snps_in_cs": int(len(snps_in_cs_set)),
            "credible_sets": cs_rows,
            "top_pip_snps": top_pip_df.to_dict(orient="records"),
            "outputs": {
                "finemap_pip": str(pip_path),
                "credible_sets": str(cs_path),
            },
            "figures": [str(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
