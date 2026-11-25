#!/usr/bin/env python3
"""meth_inspect: inspect Bismark/MethylDackel coverage files or beta/M-value matrices."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    clean_sample_stem,
    emit,
    fail,
    is_coverage_file,
    parse_coverage_records,
    parse_cpg_coord,
    read_methylation_matrix,
    read_table,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default=None)
    args = parser.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        fail(f"input path not found: {in_path}")

    warnings: list[str] = []
    if in_path.is_dir():
        cov_files = sorted([p for p in in_path.iterdir() if p.is_file() and is_coverage_file(p)])
        if not cov_files:
            fail(f"directory '{in_path}' contains no Bismark .cov(.gz) or MethylDackel .bedGraph(.gz) files")
        samples = [clean_sample_stem(p) for p in cov_files]
        fmt, recs = parse_coverage_records(cov_files[0], max_lines=5000)
        betas = np.array([r[1] for r in recs], dtype=float) if recs else np.array([0.0])
        depths = np.array([r[2] for r in recs], dtype=float) if recs else np.array([0.0])
        input_kind = f"{fmt}_directory"
        n_cpgs = len(recs)
        coord_ok = True
        val_summary = {
            "median_beta_first_sample": round(float(np.median(betas)), 4),
            "median_depth_first_sample": round(float(np.median(depths)), 2),
        }
        if float(np.median(depths)) < 5.0:
            warnings.append("first sample has median CpG read depth < 5x; consider RRBS/targeted filtering or DMR aggregation")
    elif is_coverage_file(in_path):
        fmt, recs = parse_coverage_records(in_path, max_lines=10000)
        samples = [clean_sample_stem(in_path)]
        betas = np.array([r[1] for r in recs], dtype=float) if recs else np.array([0.0])
        depths = np.array([r[2] for r in recs], dtype=float) if recs else np.array([0.0])
        input_kind = fmt
        n_cpgs = len(recs)
        coord_ok = True
        val_summary = {
            "median_beta": round(float(np.median(betas)), 4),
            "median_depth": round(float(np.median(depths)), 2),
        }
    else:
        mat, scale = read_methylation_matrix(in_path)
        samples = list(mat.columns)
        n_cpgs = int(mat.shape[0])
        input_kind = f"{scale}_matrix"
        sample_ids = list(mat.index[: min(200, len(mat.index))])
        coord_ok = sum(1 for cid in sample_ids if parse_cpg_coord(cid) is not None) >= 0.8 * max(1, len(sample_ids))
        arr = mat.to_numpy(dtype=float)
        finite = arr[np.isfinite(arr)]
        val_summary = {
            "min": round(float(np.min(finite)), 4),
            "median": round(float(np.median(finite)), 4),
            "max": round(float(np.max(finite)), 4),
            "missing_fraction": round(float(np.isnan(arr).mean()), 4),
        }
        if not coord_ok:
            warnings.append(
                "CpG IDs are probe identifiers (e.g. cg00000029) rather than 'chr:pos'; supply --probe_coords in meth_dmr for DMR calling"
            )

    meta_overlap = None
    if args.metadata:
        meta = read_table(args.metadata)
        scol = args.sample_col or str(meta.columns[0])
        if scol not in meta.columns:
            fail(f"sample_col '{scol}' not found in metadata")
        meta_s = set(meta[scol].astype(str).str.strip())
        matched = [s for s in samples if s in meta_s]
        meta_overlap = {
            "sample_col": scol,
            "n_metadata_rows": int(len(meta)),
            "n_matched_samples": len(matched),
        }
        if len(matched) < len(samples):
            warnings.append(f"metadata matched {len(matched)}/{len(samples)} methylation samples")
        if args.group_col and args.group_col in meta.columns:
            vc = meta[meta[scol].astype(str).str.strip().isin(set(matched))][args.group_col].astype(str).value_counts()
            meta_overlap["group_counts"] = {str(k): int(v) for k, v in vc.items()}

    emit(
        {
            "tool": "meth_inspect",
            "input_kind": input_kind,
            "n_samples": len(samples),
            "samples_preview": samples[:15],
            "n_cpgs_sampled_or_total": n_cpgs,
            "coordinate_anchored_ids": coord_ok,
            "value_summary": val_summary,
            "metadata_overlap": meta_overlap,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
