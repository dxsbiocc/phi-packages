#!/usr/bin/env python3
"""prot_inspect: inspect MS proteomics / PTM matrices and experimental design."""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    check_collinearity,
    clean_sample_name,
    contaminant_mask,
    detect_format_and_columns,
    emit,
    infer_scale,
    load_metadata,
    read_table,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--id_col", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default="condition")
    parser.add_argument("--batch_col", default=None)
    args = parser.parse_args()

    df = read_table(args.input)
    det = detect_format_and_columns(df, args.id_col)
    id_col = det["id_col"]
    sample_cols = det["sample_cols"]

    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    # Treat exact zeros as missing in raw linear MS intensities
    scale_info = infer_scale(mat)
    if scale_info["scale"] == "raw_intensity":
        mat_na = mat.replace(0, np.nan)
    else:
        mat_na = mat

    cmask = contaminant_mask(df, id_col)
    total_cells = int(mat_na.size)
    missing_cells = int(mat_na.isna().sum().sum())
    overall_missing_frac = round(missing_cells / total_cells, 4) if total_cells else 0.0
    per_sample_missing = {
        clean_sample_name(c): round(float(mat_na[c].isna().mean()), 4) for c in sample_cols
    }
    per_sample_detected = {
        clean_sample_name(c): int(mat_na[c].notna().sum()) for c in sample_cols
    }

    warnings: list[str] = []
    if int(cmask.sum()) > 0:
        warnings.append(
            f"detected {int(cmask.sum())} contaminant/reverse/site-only rows; filter them in prot_qc_impute"
        )
    if scale_info["scale"] == "raw_intensity":
        warnings.append("values appear to be linear raw intensities; log2-transform in prot_qc_impute before DE")
    if overall_missing_frac > 0.3:
        warnings.append(
            f"high overall missingness ({overall_missing_frac:.1%}); use group-aware valid-value filtering and MinProb/kNN imputation"
        )

    design: dict = {"metadata_provided": False}
    if args.metadata:
        meta, scol = load_metadata(args.metadata, args.sample_col)
        cleaned_samples = [clean_sample_name(c) for c in sample_cols]
        meta_samples = list(meta[scol])
        overlap = [s for s in cleaned_samples if s in set(meta_samples)]
        missing_in_meta = [s for s in cleaned_samples if s not in set(meta_samples)]
        missing_in_mat = [s for s in meta_samples if s not in set(cleaned_samples)]
        group_counts = {}
        if args.group_col in meta.columns:
            sub = meta[meta[scol].isin(overlap)]
            vc = sub[args.group_col].astype(str).value_counts()
            group_counts = {str(k): int(v) for k, v in vc.items()}
            small = [g for g, cnt in group_counts.items() if cnt < 3]
            if small:
                warnings.append(f"groups with < 3 replicates: {', '.join(small)}")
        collin = check_collinearity(meta, args.group_col, args.batch_col)
        if collin.get("confounded"):
            warnings.append(
                f"batch column '{args.batch_col}' is confounded with '{args.group_col}'; do not fit additive batch adjustment"
            )
        design = {
            "metadata_provided": True,
            "sample_col": scol,
            "group_col": args.group_col,
            "matched_samples": len(overlap),
            "unmatched_matrix_samples": missing_in_meta,
            "unmatched_metadata_samples": missing_in_mat,
            "group_counts": group_counts,
            "batch_collinearity": collin,
        }

    emit(
        {
            "input": os.path.abspath(args.input),
            "detected_format": det["format"],
            "id_col": id_col,
            "n_features": int(len(df)),
            "n_samples": int(len(sample_cols)),
            "contaminant_or_decoy_rows": int(cmask.sum()),
            "scale": scale_info,
            "overall_missing_fraction": overall_missing_frac,
            "per_sample_missing_fraction": per_sample_missing,
            "per_sample_detected_features": per_sample_detected,
            "design": design,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
