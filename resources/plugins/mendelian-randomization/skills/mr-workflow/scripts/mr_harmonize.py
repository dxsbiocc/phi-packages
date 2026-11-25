#!/usr/bin/env python3
"""mr_harmonize: Align exposure genetic instruments against outcome GWAS alleles,
resolve strand flips and palindromic A/T & C/G SNPs via EAF margins, filter weak
instruments (F < 10), and optionally apply Steiger directionality filtering.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    clump_by_distance,
    emit,
    fail,
    is_palindromic,
    read_summary_stats,
    strand_flip,
    write_artifact,
)


def approx_r2(beta: float, se: float, n: float) -> float | None:
    if not np.isfinite(n) or n <= 3 or se <= 0:
        return None
    f_stat = (beta / se) ** 2
    return float(f_stat / (f_stat + n - 2.0))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exposure", required=True)
    parser.add_argument("--outcome", required=True)
    parser.add_argument("--p_cutoff", type=float, default=5e-8)
    parser.add_argument("--clump_kb", type=int, default=10000)
    parser.add_argument("--min_f", type=float, default=10.0)
    parser.add_argument("--palindromic_eaf_margin", type=float, default=0.08)
    parser.add_argument(
        "--steiger_filter",
        nargs="?",
        const=True,
        type=lambda x: str(x).lower() not in ("false", "0", "no"),
        default=False,
    )
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    exp_all, _ = read_summary_stats(args.exposure)
    out_all, _ = read_summary_stats(args.outcome)

    exp_ivs = clump_by_distance(exp_all, args.p_cutoff, args.clump_kb)
    n_cand = int(len(exp_ivs))
    if n_cand == 0:
        fail(f"0 exposure variants passed p_cutoff={args.p_cutoff}")

    out_by_snp = out_all.set_index("SNP")
    matched = exp_ivs[exp_ivs["SNP"].isin(out_by_snp.index)].copy()
    n_matched = int(len(matched))
    if n_matched == 0:
        fail("0 exposure instrument SNPs overlapped with outcome GWAS SNP column")

    n_pal_drop = 0
    n_incompat = 0
    n_weak = 0
    n_steiger_drop = 0
    rows = []

    for _, erow in matched.iterrows():
        snp = str(erow["SNP"])
        orow = out_by_snp.loc[snp]
        ea_e, oa_e = str(erow["EA"]), str(erow["OA"])
        ea_o, oa_o = str(orow["EA"]), str(orow["OA"])
        b_e, se_e = float(erow["BETA"]), float(erow["SE"])
        b_o, se_o = float(orow["BETA"]), float(orow["SE"])
        eaf_e = float(erow["EAF"]) if pd.notna(erow["EAF"]) else np.nan
        eaf_o = float(orow["EAF"]) if pd.notna(orow["EAF"]) else np.nan

        pal = is_palindromic(ea_e, oa_e)
        if pal:
            # Ambiguous near 0.5 or missing EAF -> drop
            if not (np.isfinite(eaf_e) and np.isfinite(eaf_o)):
                n_pal_drop += 1
                continue
            if abs(eaf_e - 0.5) < args.palindromic_eaf_margin or abs(eaf_o - 0.5) < args.palindromic_eaf_margin:
                n_pal_drop += 1
                continue
            same_side = (eaf_e < 0.5 and eaf_o < 0.5) or (eaf_e > 0.5 and eaf_o > 0.5)
            if ea_e == ea_o and oa_e == oa_o:
                flip_sign = not same_side
            elif ea_e == oa_o and oa_e == ea_o:
                flip_sign = same_side
            else:
                n_incompat += 1
                continue
        else:
            if ea_e == ea_o and oa_e == oa_o:
                flip_sign = False
            elif ea_e == oa_o and oa_e == ea_o:
                flip_sign = True
            elif ea_e == strand_flip(ea_o) and oa_e == strand_flip(oa_o):
                flip_sign = False
            elif ea_e == strand_flip(oa_o) and oa_e == strand_flip(ea_o):
                flip_sign = True
            else:
                n_incompat += 1
                continue

        if flip_sign:
            b_o = -b_o
            if np.isfinite(eaf_o):
                eaf_o = 1.0 - eaf_o

        # Enforce positive exposure beta convention (TwoSampleMR standard so Egger intercept is invariant to coding)
        if b_e < 0:
            b_e = -b_e
            b_o = -b_o
            ea_e, oa_e = oa_e, ea_e
            if np.isfinite(eaf_e):
                eaf_e = 1.0 - eaf_e
            if np.isfinite(eaf_o):
                eaf_o = 1.0 - eaf_o

        f_stat = (b_e / se_e) ** 2
        if f_stat < args.min_f:
            n_weak += 1
            continue

        n_e = float(erow["N"]) if pd.notna(erow["N"]) else np.nan
        n_o = float(orow["N"]) if pd.notna(orow["N"]) else np.nan
        r2_e = approx_r2(b_e, se_e, n_e)
        r2_o = approx_r2(b_o, se_o, n_o)
        steiger_dir = (r2_e > r2_o) if (r2_e is not None and r2_o is not None) else None
        if args.steiger_filter and steiger_dir is False:
            n_steiger_drop += 1
            continue

        rows.append(
            {
                "SNP": snp,
                "CHR": erow["CHR"] if pd.notna(erow["CHR"]) else orow["CHR"],
                "POS": erow["POS"] if pd.notna(erow["POS"]) else orow["POS"],
                "effect_allele": ea_e,
                "other_allele": oa_e,
                "palindromic": pal,
                "beta_exp": round(b_e, 6),
                "se_exp": round(se_e, 6),
                "pval_exp": float(erow["P"]),
                "eaf_exp": round(eaf_e, 5) if np.isfinite(eaf_e) else None,
                "n_exp": int(n_e) if np.isfinite(n_e) else None,
                "beta_out": round(b_o, 6),
                "se_out": round(se_o, 6),
                "pval_out": float(orow["P"]),
                "eaf_out": round(eaf_o, 5) if np.isfinite(eaf_o) else None,
                "n_out": int(n_o) if np.isfinite(n_o) else None,
                "f_stat": round(f_stat, 3),
                "r2_exp": round(r2_e, 6) if r2_e is not None else None,
                "r2_out": round(r2_o, 6) if r2_o is not None else None,
                "steiger_pass": steiger_dir,
            }
        )

    if not rows:
        fail(
            f"0 harmonized instruments remained (matched={n_matched}, palindromic_dropped={n_pal_drop}, weak_F_dropped={n_weak})"
        )

    harm_df = pd.DataFrame(rows)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    harm_path = outdir / "harmonized_instruments.tsv"
    harm_df.to_csv(harm_path, sep="\t", index=False)

    params_meta = {
        "p_cutoff": args.p_cutoff,
        "clump_kb": args.clump_kb,
        "min_f": args.min_f,
        "palindromic_eaf_margin": args.palindromic_eaf_margin,
        "steiger_filter": args.steiger_filter,
    }
    artifacts = [
        write_artifact(
            harm_path,
            "mr_harmonize",
            f"Harmonized MR instruments ({len(harm_df)} SNPs, mean F={harm_df['f_stat'].mean():.1f})",
            inputs=(args.exposure, args.outcome),
            parameters=params_meta,
            shape=harm_df.shape,
        )
    ]

    warnings: list[str] = []
    if len(harm_df) < 3:
        warnings.append(
            f"only {len(harm_df)} harmonized SNP(s); MR-Egger, Weighted Median, and Mode estimators require >= 3 SNPs (run mr_coloc if cis-QTL)"
        )
    if args.clump_kb > 0 and harm_df["POS"].isna().all():
        warnings.append("exposure table lacked CHR/POS coordinates: ensure exposure SNPs were already LD-clumped upstream")

    emit(
        {
            "tool": "mr_harmonize",
            "n_exposure_candidates": n_cand,
            "n_matched_outcome": n_matched,
            "n_palindromic_dropped": n_pal_drop,
            "n_incompatible_dropped": n_incompat,
            "n_weak_f_dropped": n_weak,
            "n_steiger_dropped": n_steiger_drop,
            "n_harmonized_retained": int(len(harm_df)),
            "mean_f_stat": round(float(harm_df["f_stat"].mean()), 3),
            "min_f_stat": round(float(harm_df["f_stat"].min()), 3),
            "outputs": {"harmonized_table": str(harm_path)},
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
