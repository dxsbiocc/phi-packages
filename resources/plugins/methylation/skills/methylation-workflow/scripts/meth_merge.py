#!/usr/bin/env python3
"""meth_merge: merge per-sample Bismark .cov(.gz) or MethylDackel .bedGraph(.gz)
files into coverage-filtered beta_matrix.tsv and coverage_matrix.tsv.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    clean_sample_stem,
    emit,
    fail,
    is_coverage_file,
    parse_coverage_records,
    parse_cpg_coord,
    write_artifact,
)


def resolve_inputs(specs: list[str]) -> list[tuple[str, Path]]:
    out: list[tuple[str, Path]] = []
    for raw in specs:
        for token in str(raw).split(","):
            s = token.strip()
            if not s:
                continue
            if "=" in s and not Path(s).exists():
                lbl, p_str = s.split("=", 1)
                p = Path(p_str.strip())
                if not p.is_file():
                    fail(f"coverage file not found for '{lbl}': {p}")
                out.append((lbl.strip(), p))
            else:
                p = Path(s)
                if p.is_dir():
                    for child in sorted(p.iterdir()):
                        if child.is_file() and is_coverage_file(child):
                            out.append((clean_sample_stem(child), child))
                elif p.is_file():
                    out.append((clean_sample_stem(p), p))
                else:
                    fail(f"coverage input not found: {p}")
    if len(out) < 2:
        fail(f"meth_merge requires at least 2 sample coverage files; found {len(out)}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", action="extend", nargs="+", default=[], required=True)
    parser.add_argument("--min_cov", type=int, default=5)
    parser.add_argument("--min_samples_frac", type=float, default=0.7)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    labeled = resolve_inputs(args.inputs)
    samples = [s for s, _ in labeled]
    if len(set(samples)) != len(samples):
        fail(f"duplicate sample names in meth_merge inputs: {samples}")

    n_samples = len(samples)
    min_req = max(1, int(math.ceil(args.min_samples_frac * n_samples)))

    # Map cpg_id -> list of (sample_idx, beta, depth)
    cpg_data: dict[str, list[tuple[int, float, int]]] = {}
    for s_idx, (_, path) in enumerate(labeled):
        _, recs = parse_coverage_records(path)
        for cid, beta, depth in recs:
            if depth >= args.min_cov:
                cpg_data.setdefault(cid, []).append((s_idx, beta, depth))

    n_union = len(cpg_data)
    kept_ids = [cid for cid, entries in cpg_data.items() if len(entries) >= min_req]
    if not kept_ids:
        fail(
            f"0 CpGs passed min_cov={args.min_cov} in >= {min_req}/{n_samples} samples; "
            "lower --min_cov or --min_samples_frac"
        )

    def sort_key(cid: str) -> tuple[str, int]:
        parsed = parse_cpg_coord(cid)
        return (parsed[0], parsed[1]) if parsed else (cid, 0)

    kept_ids.sort(key=sort_key)

    beta_mat = np.full((len(kept_ids), n_samples), np.nan, dtype=float)
    cov_mat = np.zeros((len(kept_ids), n_samples), dtype=int)
    for r_i, cid in enumerate(kept_ids):
        for s_i, b_val, d_val in cpg_data[cid]:
            beta_mat[r_i, s_i] = round(b_val, 5)
            cov_mat[r_i, s_i] = d_val

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    beta_df = pd.DataFrame(beta_mat, index=kept_ids, columns=samples)
    beta_df.index.name = "cpg_id"
    beta_path = outdir / "beta_matrix.tsv"
    beta_df.to_csv(beta_path, sep="\t")

    cov_df = pd.DataFrame(cov_mat, index=kept_ids, columns=samples)
    cov_df.index.name = "cpg_id"
    cov_path = outdir / "coverage_matrix.tsv"
    cov_df.to_csv(cov_path, sep="\t")

    in_files = [p for _, p in labeled]
    params_meta = {"min_cov": args.min_cov, "min_samples_frac": args.min_samples_frac}
    artifacts = [
        write_artifact(beta_path, "meth_merge", f"Merged CpG beta matrix ({len(kept_ids)} CpGs x {n_samples} samples)", inputs=in_files, parameters=params_meta, shape=beta_df.shape),
        write_artifact(cov_path, "meth_merge", f"Merged CpG read-depth matrix ({len(kept_ids)} CpGs x {n_samples} samples)", inputs=in_files, parameters=params_meta, shape=cov_df.shape),
    ]

    warnings: list[str] = []
    if args.min_cov < 5:
        warnings.append(f"min_cov={args.min_cov} (< 5): low-depth CpG beta estimates carry high binomial sampling variance")

    emit(
        {
            "tool": "meth_merge",
            "n_samples": n_samples,
            "samples": samples,
            "min_cov": args.min_cov,
            "min_samples_required": min_req,
            "n_cpgs_union": n_union,
            "n_cpgs_retained": len(kept_ids),
            "outputs": {
                "beta_matrix": str(beta_path),
                "coverage_matrix": str(cov_path),
            },
            "scientificContext": {
                str(beta_path): {"valueSemantics": "normalized", "sampleIds": samples},
                str(cov_path): {"valueSemantics": "raw_counts", "sampleIds": samples},
            },
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
