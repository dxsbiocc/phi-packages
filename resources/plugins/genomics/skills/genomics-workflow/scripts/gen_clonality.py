#!/usr/bin/env python3
"""gen_clonality: McGranahan/PyClone-inspired Cancer Cell Fraction (CCF) & subclonal evolution."""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    CANCER_DRIVER_GENES,
    NON_SYNONYMOUS_CLASSES,
    emit,
    extract_vaf,
    fail,
    read_allele_segments,
    read_table,
    write_artifact,
)


def load_clonality_maf(path: str) -> pd.DataFrame:
    df = read_table(path)
    required = ["Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Start_Position"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        fail(f"MAF/TSV missing required columns: {', '.join(missing)}")
    out = df.copy()
    out["Hugo_Symbol"] = out["Hugo_Symbol"].astype(str).str.strip()
    out["Tumor_Sample_Barcode"] = out["Tumor_Sample_Barcode"].astype(str).str.strip()
    out["Chromosome"] = out["Chromosome"].astype(str).str.strip().str.replace("^chr", "", regex=True)
    out["Start_Position"] = pd.to_numeric(out["Start_Position"], errors="coerce")
    if "Variant_Classification" not in out.columns:
        out["Variant_Classification"] = "Missense_Mutation"
    else:
        out["Variant_Classification"] = out["Variant_Classification"].astype(str).str.strip()
    return out.dropna(subset=["Start_Position"]).copy()


def extract_counts_and_vaf(df: pd.DataFrame, warnings: list[str]) -> pd.DataFrame:
    """Resolve t_alt_count, t_depth, and VAF per variant."""
    lower = {c.lower(): c for c in df.columns}
    alt_c = next((lower[k] for k in ("t_alt_count", "alt_count", "ad_alt", "alt") if k in lower), None)
    ref_c = next((lower[k] for k in ("t_ref_count", "ref_count", "ad_ref", "ref") if k in lower), None)
    dep_c = next((lower[k] for k in ("t_depth", "depth", "dp", "total_depth") if k in lower), None)

    vaf_series, _ = extract_vaf(df)
    alt_s = pd.to_numeric(df[alt_c], errors="coerce") if alt_c else pd.Series(np.nan, index=df.index, dtype=float)
    ref_s = pd.to_numeric(df[ref_c], errors="coerce") if ref_c else pd.Series(np.nan, index=df.index, dtype=float)
    dep_s = pd.to_numeric(df[dep_c], errors="coerce") if dep_c else pd.Series(np.nan, index=df.index, dtype=float)

    if alt_c and dep_c:
        dep_final = dep_s.fillna(alt_s + ref_s) if ref_c else dep_s
        alt_final = alt_s
        vaf_final = vaf_series if vaf_series is not None else (alt_final / dep_final.replace(0, np.nan))
    elif alt_c and ref_c:
        dep_final = alt_s + ref_s
        alt_final = alt_s
        vaf_final = vaf_series if vaf_series is not None else (alt_final / dep_final.replace(0, np.nan))
    elif vaf_series is not None and dep_c:
        dep_final = dep_s
        vaf_final = vaf_series
        alt_final = (vaf_final * dep_final).round()
    elif vaf_series is not None:
        warnings.append("no t_depth or t_alt_count/t_ref_count detected; using nominal depth=100 for Binomial/Beta CI")
        dep_final = pd.Series(100.0, index=df.index, dtype=float)
        vaf_final = vaf_series
        alt_final = (vaf_final * dep_final).round()
    else:
        fail("MAF must provide t_alt_count+t_ref_count, t_alt_count+t_depth, or vaf+t_depth")

    res = df.copy()
    res["_depth"] = pd.to_numeric(dep_final, errors="coerce")
    res["_alt"] = pd.to_numeric(alt_final, errors="coerce")
    res["_vaf"] = pd.to_numeric(vaf_final, errors="coerce").clip(0.0, 1.0)
    res = res.dropna(subset=["_depth", "_alt", "_vaf"])
    res = res[res["_depth"] > 0].copy()
    res["_depth"] = res["_depth"].round().astype(int).clip(lower=1)
    res["_alt"] = res["_alt"].round().astype(int)
    res["_alt"] = np.minimum(np.maximum(res["_alt"], 0), res["_depth"])
    return res


def load_purity_table(path: str) -> dict[str, tuple[float, float | None]]:
    tbl = read_table(path)
    lower = {c.lower(): c for c in tbl.columns}
    sample_c = next(
        (lower[k] for k in ("sample", "tumor_sample_barcode", "sample_id", "id") if k in lower),
        tbl.columns[0],
    )
    pur_c = next((lower[k] for k in ("purity", "tumor_purity", "rho", "cellularity") if k in lower), None)
    if pur_c is None:
        fail("purity_table must contain a purity column (purity, tumor_purity, rho, or cellularity)")
    plo_c = next((lower[k] for k in ("ploidy", "tumor_ploidy") if k in lower), None)

    out: dict[str, tuple[float, float | None]] = {}
    for _, row in tbl.iterrows():
        s = str(row[sample_c]).strip()
        p = pd.to_numeric(row[pur_c], errors="coerce")
        if not np.isfinite(p):
            continue
        if p > 1.0:
            p = p / 100.0
        p = float(np.clip(p, 0.05, 1.0))
        pl = None
        if plo_c and pd.notna(row[plo_c]):
            pl_val = pd.to_numeric(row[plo_c], errors="coerce")
            if np.isfinite(pl_val) and pl_val > 0:
                pl = float(pl_val)
        out[s] = (p, pl)
    return out


def lookup_total_cn(
    sample: str,
    chrom: str,
    pos: float,
    seg_by_sample: dict[str, pd.DataFrame],
    fallback_cn: float,
) -> float:
    sub = seg_by_sample.get(sample)
    if sub is None or sub.empty:
        return fallback_cn
    hit = sub[(sub["chrom"] == chrom) & (sub["start"] <= pos) & (sub["end"] >= pos)]
    if hit.empty:
        return fallback_cn
    val = float(hit.iloc[0]["total_cn"])
    return val if np.isfinite(val) and val >= 1.0 else fallback_cn


def estimate_ccf(
    k_alt: int,
    n_depth: int,
    vaf: float,
    purity: float,
    total_cn: float,
) -> tuple[int, float, float, float, float]:
    """Return (multiplicity_hat, expected_clonal_vaf, ccf, ccf_ci_low, ccf_ci_high)."""
    cn_eff = max(1.0, float(total_cn))
    max_m = max(1, int(round(cn_eff)))
    denom = purity * cn_eff + (1.0 - purity) * 2.0

    def f_exp(m: int) -> float:
        return (purity * float(m)) / denom

    f1 = f_exp(1)
    if vaf >= f1 and max_m > 1:
        best_m = 1
        best_ll = -np.inf
        for m in range(1, max_m + 1):
            p_m = min(0.999, max(1e-6, f_exp(m)))
            ll = float(stats.binom.logpmf(k_alt, n_depth, p_m))
            if ll > best_ll:
                best_ll = ll
                best_m = m
        m_hat = best_m
    else:
        m_hat = 1

    f_hat = f_exp(m_hat)
    scale = denom / (purity * float(m_hat))
    ccf = min(1.5, vaf * scale)

    # Exact Clopper-Pearson 95% binomial CI mapped onto CCF
    vaf_lo = 0.0 if k_alt <= 0 else float(stats.beta.ppf(0.025, k_alt, n_depth - k_alt + 1))
    vaf_hi = 1.0 if k_alt >= n_depth else float(stats.beta.ppf(0.975, k_alt + 1, n_depth - k_alt))
    ccf_lo = min(1.5, vaf_lo * scale)
    ccf_hi = min(1.5, vaf_hi * scale)
    return m_hat, f_hat, ccf, ccf_lo, ccf_hi


def williams_neutral_r2(vafs: np.ndarray, f_min: float = 0.12, f_max: float = 0.24) -> tuple[int, float | None]:
    """Test Williams et al. (2016) neutral tail linearity M(f) ~ 1/f over [0.12, 0.24]."""
    tail = vafs[(vafs >= f_min) & (vafs <= f_max)]
    n_tail = int(len(tail))
    if n_tail < 5:
        return n_tail, None
    unique_f = np.sort(np.unique(tail))[::-1]
    if len(unique_f) < 4:
        return n_tail, None
    cum_m = np.array([float(np.sum(tail >= f)) for f in unique_f], dtype=float)
    inv_f = 1.0 / unique_f - (1.0 / f_max)
    if np.allclose(inv_f, inv_f[0]) or np.allclose(cum_m, cum_m[0]):
        return n_tail, None
    res = stats.linregress(inv_f, cum_m)
    r2 = float(res.rvalue ** 2) if res.slope > 0 and np.isfinite(res.rvalue) else 0.0
    return n_tail, round(r2, 4)


def render_clonal_architecture(
    var_df: pd.DataFrame,
    samp_df: pd.DataFrame,
    out_path: str,
    ccf_cutoff: float,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5))
    plot_s = samp_df.head(30).copy()
    x = np.arange(len(plot_s))

    # Left panel: Per-sample clonal vs subclonal mutation burden
    clonal_c = plot_s["n_clonal"].to_numpy(dtype=float)
    sub_c = plot_s["n_subclonal"].to_numpy(dtype=float)
    axes[0].bar(x, clonal_c, color="#0284C7", label="Clonal")
    axes[0].bar(x, sub_c, bottom=clonal_c, color="#F59E0B", label="Subclonal")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(plot_s["Tumor_Sample_Barcode"], rotation=45, ha="right", fontsize=7)
    axes[0].set_ylabel("Somatic variants")
    axes[0].set_title("Per-sample Clonal vs Subclonal Burden")
    axes[0].legend(frameon=False, fontsize=8)

    # Right panel: Cohort CCF histogram split by clonality
    c_vals = var_df.loc[var_df["clonality"] == "clonal", "ccf"].to_numpy(dtype=float)
    s_vals = var_df.loc[var_df["clonality"] == "subclonal", "ccf"].to_numpy(dtype=float)
    bins = np.linspace(0.0, 1.5, 31)
    if len(c_vals):
        axes[1].hist(c_vals, bins=bins, alpha=0.75, color="#0284C7", edgecolor="#0F172A", label=f"Clonal (n={len(c_vals)})")
    if len(s_vals):
        axes[1].hist(s_vals, bins=bins, alpha=0.75, color="#F59E0B", edgecolor="#0F172A", label=f"Subclonal (n={len(s_vals)})")
    axes[1].axvline(ccf_cutoff, color="#EF4444", linestyle="--", linewidth=1.1, label=f"Cutoff ({ccf_cutoff})")
    axes[1].set_xlabel("Cancer Cell Fraction (CCF)")
    axes[1].set_ylabel("Variant count")
    axes[1].set_title("Cohort CCF Distribution")
    axes[1].legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--purity_table", default=None)
    parser.add_argument("--segments", default=None)
    parser.add_argument("--default_purity", type=float, default=0.70)
    parser.add_argument("--default_cn", type=int, default=2)
    parser.add_argument("--clonal_ccf_cutoff", type=float, default=0.85)
    args = parser.parse_args()

    if not (0.05 <= args.default_purity <= 1.0):
        fail("default_purity must be between 0.05 and 1.0")
    if args.default_cn < 1:
        fail("default_cn must be >= 1")

    os.makedirs(args.output_dir, exist_ok=True)
    warnings: list[str] = []

    raw_df = load_clonality_maf(args.maf)
    df = extract_counts_and_vaf(raw_df, warnings)
    if df.empty:
        fail("no variants with valid depth and VAF remaining")

    purity_map = load_purity_table(args.purity_table) if args.purity_table else {}
    seg_by_sample: dict[str, pd.DataFrame] = {}
    if args.segments:
        seg_df = read_allele_segments(args.segments, require_allele_cn=False)
        for s, grp in seg_df.groupby("sample"):
            seg_by_sample[str(s)] = grp

    var_rows: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        sample = str(row["Tumor_Sample_Barcode"])
        gene = str(row["Hugo_Symbol"])
        vc = str(row["Variant_Classification"])
        chrom = str(row["Chromosome"])
        pos = float(row["Start_Position"])
        k_alt = int(row["_alt"])
        n_depth = int(row["_depth"])
        vaf = float(row["_vaf"])

        pur, ploidy = purity_map.get(sample, (args.default_purity, None))
        fallback_cn = float(args.default_cn)
        tot_cn = lookup_total_cn(sample, chrom, pos, seg_by_sample, fallback_cn)
        m_hat, f_exp_val, ccf, ccf_lo, ccf_hi = estimate_ccf(k_alt, n_depth, vaf, pur, tot_cn)

        is_clonal = (ccf_hi >= 0.95) or (ccf >= args.clonal_ccf_cutoff)
        clonality = "clonal" if is_clonal else "subclonal"
        is_driver = (gene in CANCER_DRIVER_GENES) and (vc in NON_SYNONYMOUS_CLASSES)

        var_rows.append(
            {
                "Tumor_Sample_Barcode": sample,
                "Hugo_Symbol": gene,
                "Chromosome": chrom,
                "Start_Position": int(pos),
                "Variant_Classification": vc,
                "t_alt_count": k_alt,
                "t_depth": n_depth,
                "vaf": round(vaf, 4),
                "purity": round(pur, 4),
                "total_cn": round(tot_cn, 3),
                "multiplicity": m_hat,
                "expected_clonal_vaf": round(f_exp_val, 4),
                "ccf": round(ccf, 4),
                "ccf_ci_low": round(ccf_lo, 4),
                "ccf_ci_high": round(ccf_hi, 4),
                "clonality": clonality,
                "is_driver": bool(is_driver),
            }
        )

    var_df = pd.DataFrame(var_rows)

    samp_rows: list[dict[str, Any]] = []
    for sample, grp in var_df.groupby("Tumor_Sample_Barcode"):
        pur, ploidy = purity_map.get(str(sample), (args.default_purity, None))
        n_var = int(len(grp))
        n_clonal = int((grp["clonality"] == "clonal").sum())
        n_sub = int((grp["clonality"] == "subclonal").sum())
        drv = grp[grp["is_driver"]]
        clonal_drv = sorted(drv.loc[drv["clonality"] == "clonal", "Hugo_Symbol"].unique())
        sub_drv = sorted(drv.loc[drv["clonality"] == "subclonal", "Hugo_Symbol"].unique())
        n_tail, r2_neutral = williams_neutral_r2(grp["vaf"].to_numpy(dtype=float))
        samp_rows.append(
            {
                "Tumor_Sample_Barcode": str(sample),
                "purity": round(pur, 4),
                "ploidy": round(ploidy, 3) if ploidy is not None else np.nan,
                "n_variants": n_var,
                "n_clonal": n_clonal,
                "n_subclonal": n_sub,
                "clonal_fraction": round(n_clonal / n_var, 4) if n_var else 0.0,
                "median_ccf": round(float(grp["ccf"].median()), 4),
                "n_clonal_drivers": len(clonal_drv),
                "n_subclonal_drivers": len(sub_drv),
                "clonal_driver_genes": ",".join(clonal_drv),
                "subclonal_driver_genes": ",".join(sub_drv),
                "n_neutral_tail_variants": n_tail,
                "r2_neutral": r2_neutral if r2_neutral is not None else np.nan,
            }
        )
    samp_df = pd.DataFrame(samp_rows).sort_values("clonal_fraction", ascending=False)

    var_path = os.path.join(args.output_dir, "clonality_variants.tsv")
    samp_path = os.path.join(args.output_dir, "clonality_sample_summary.tsv")
    fig_path = os.path.join(args.output_dir, "ccf_clonal_architecture.png")
    var_df.to_csv(var_path, sep="\t", index=False)
    samp_df.to_csv(samp_path, sep="\t", index=False)
    render_clonal_architecture(var_df, samp_df, fig_path, args.clonal_ccf_cutoff)

    params_meta = {
        "default_purity": args.default_purity,
        "default_cn": args.default_cn,
        "clonal_ccf_cutoff": args.clonal_ccf_cutoff,
        "purity_table_provided": bool(args.purity_table),
        "segments_provided": bool(args.segments),
    }
    in_files = (args.maf, args.purity_table, args.segments)
    artifacts = [
        write_artifact(fig_path, "gen_clonality", "Cohort CCF clonal vs subclonal architecture", inputs=in_files, parameters=params_meta),
        write_artifact(var_path, "gen_clonality", "Per-variant multiplicity and Cancer Cell Fraction (CCF)", inputs=in_files, parameters=params_meta, shape=var_df.shape),
        write_artifact(samp_path, "gen_clonality", "Per-sample clonality, driver timing, and neutral tail R2", inputs=in_files, parameters=params_meta, shape=samp_df.shape),
    ]

    clean_samples = []
    for rec in samp_df.to_dict(orient="records"):
        clean_samples.append(
            {k: (None if (isinstance(v, float) and not np.isfinite(v)) else v) for k, v in rec.items()}
        )
    r2_vals = samp_df["r2_neutral"].dropna().to_numpy(dtype=float)

    emit(
        {
            "maf": os.path.abspath(args.maf),
            "variants_table": os.path.abspath(var_path),
            "sample_summary_table": os.path.abspath(samp_path),
            "n_samples": int(len(samp_df)),
            "n_variants": int(len(var_df)),
            "n_clonal_variants": int((var_df["clonality"] == "clonal").sum()),
            "n_subclonal_variants": int((var_df["clonality"] == "subclonal").sum()),
            "median_ccf": round(float(var_df["ccf"].median()), 4),
            "median_clonal_fraction": round(float(samp_df["clonal_fraction"].median()), 4),
            "median_r2_neutral": round(float(np.median(r2_vals)), 4) if len(r2_vals) else None,
            "samples": clean_samples,
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
