#!/usr/bin/env python3
"""gen_inspect: validate MAF, clinical metadata, and optional CNV .seg inputs."""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    NON_SYNONYMOUS_CLASSES,
    emit,
    extract_vaf,
    load_clinical,
    read_maf,
    read_seg,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--clinical", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--vaf_col", default=None)
    parser.add_argument("--seg", default=None)
    args = parser.parse_args()

    df = read_maf(args.maf)
    n_variants = int(len(df))
    samples = sorted(df["Tumor_Sample_Barcode"].unique())
    nonsyn_mask = df["Variant_Classification"].isin(NON_SYNONYMOUS_CLASSES)
    n_nonsyn = int(nonsyn_mask.sum())
    n_silent = int((df["Variant_Classification"] == "Silent").sum())

    vc_counts = {str(k): int(v) for k, v in df["Variant_Classification"].value_counts().items()}
    vt_counts = {str(k): int(v) for k, v in df["Variant_Type"].value_counts().items()}

    vaf, vaf_source = extract_vaf(df, args.vaf_col)
    vaf_info = {"present": False, "source": None}
    if vaf is not None and vaf.notna().any():
        vaf_info = {
            "present": True,
            "source": vaf_source,
            "non_null_fraction": round(float(vaf.notna().mean()), 4),
            "median": round(float(vaf.median()), 4),
            "p05": round(float(vaf.quantile(0.05)), 4),
            "p95": round(float(vaf.quantile(0.95)), 4),
        }

    # Check if trinucleotide context column is present for SBS96
    lower_cols = {c.lower(): c for c in df.columns}
    context_col = next((lower_cols[k] for k in ("context", "triplet", "trinucleotide_context", "ref_context") if k in lower_cols), None)

    warnings: list[str] = []
    if not vaf_info["present"]:
        warnings.append("no VAF or t_alt_count/t_depth columns detected; MATH score cannot run and TMB cannot filter low-VAF artifacts")
    if context_col is None:
        warnings.append("no trinucleotide Context column detected; gen_signatures will fall back to 6-class Ti/Tv unless context is supplied")
    if len(samples) < 3:
        warnings.append(f"only {len(samples)} sample(s) in MAF; cohort interaction and differential mutation tests need larger cohorts")

    clinical_info = {"provided": False}
    if args.clinical:
        clin, scol = load_clinical(args.clinical, args.sample_col)
        clin_samples = set(clin[scol])
        maf_samples = set(samples)
        overlap = sorted(maf_samples & clin_samples)
        clinical_info = {
            "provided": True,
            "sample_col": scol,
            "n_clinical_rows": int(len(clin)),
            "n_matched_samples": len(overlap),
            "maf_only_samples": len(maf_samples - clin_samples),
            "clinical_only_samples": len(clin_samples - maf_samples),
            "columns": list(clin.columns),
        }

    seg_info = {"provided": False}
    if args.seg:
        seg_df = read_seg(args.seg)
        seg_samples = set(seg_df["Sample"])
        seg_info = {
            "provided": True,
            "n_segments": int(len(seg_df)),
            "n_samples": len(seg_samples),
            "n_matched_with_maf": len(seg_samples & set(samples)),
        }

    emit(
        {
            "maf": os.path.abspath(args.maf),
            "n_variants": n_variants,
            "n_nonsynonymous": n_nonsyn,
            "n_silent": n_silent,
            "n_samples": len(samples),
            "n_genes": int(df.loc[nonsyn_mask, "Hugo_Symbol"].nunique()),
            "variant_classifications": vc_counts,
            "variant_types": vt_counts,
            "vaf": vaf_info,
            "context_column": context_col,
            "clinical": clinical_info,
            "cnv_seg": seg_info,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
