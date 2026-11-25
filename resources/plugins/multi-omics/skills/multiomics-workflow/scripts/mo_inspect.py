#!/usr/bin/env python3
"""mo_inspect: inspect 2+ omics feature x sample matrices for sample alignment,
barcode format mismatches, missingness, value scales, and feature overlap.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, infer_matrix_scale, parse_labeled_paths, read_matrix, read_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrices", action="extend", nargs="+", default=[], required=True)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default=None)
    args = parser.parse_args()

    labeled = parse_labeled_paths(args.matrices, prefix="modality")
    if len(labeled) < 2:
        fail("mo_inspect requires at least 2 omics matrices")

    mod_records = []
    sample_sets: list[set[str]] = []
    feature_sets: list[tuple[str, set[str]]] = []
    ordered_first_samples: list[str] = []

    for idx, (lbl, path) in enumerate(labeled):
        df = read_matrix(path)
        if idx == 0:
            ordered_first_samples = list(df.columns)
        arr = df.to_numpy(dtype=float)
        miss_frac = round(float(np.isnan(arr).mean()), 4)
        zero_frac = round(float((np.nan_to_num(arr, nan=1.0) == 0.0).mean()), 4)
        scale = infer_matrix_scale(df)
        s_set = set(df.columns)
        sample_sets.append(s_set)
        feature_sets.append((lbl, set(df.index)))
        mod_records.append(
            {
                "modality": lbl,
                "path": str(path.resolve()),
                "n_features": int(df.shape[0]),
                "n_samples": int(df.shape[1]),
                "missing_fraction": miss_frac,
                "zero_fraction": zero_frac,
                "inferred_scale": scale,
                "samples_preview": list(df.columns[:5]),
            }
        )

    shared_set = set.intersection(*sample_sets)
    shared_ordered = [s for s in ordered_first_samples if s in shared_set]

    barcode_hint = None
    if len(shared_ordered) == 0:
        # Check TCGA-style 12-char patient barcode prefix match
        pref_sets = [{s[:12] for s in ss if len(s) >= 12} for ss in sample_sets]
        pref_shared = set.intersection(*pref_sets) if all(pref_sets) else set()
        if pref_shared:
            barcode_hint = (
                f"0 exact sample matches, but {len(pref_shared)} shared 12-char barcode prefixes found "
                f"(e.g. {sorted(pref_shared)[0]}): harmonize sample aliquots vs patient IDs before integration"
            )

    pairwise_feat = []
    for i in range(len(feature_sets)):
        for j in range(i + 1, len(feature_sets)):
            n1, f1 = feature_sets[i]
            n2, f2 = feature_sets[j]
            inter = len(f1 & f2)
            pairwise_feat.append(
                {
                    "pair": f"{n1}_vs_{n2}",
                    "shared_feature_ids": int(inter),
                }
            )

    warnings: list[str] = []
    if len(shared_ordered) == 0:
        warnings.append(barcode_hint or "0 shared sample IDs across all supplied matrices; check column headers")
    elif len(shared_ordered) < 6:
        warnings.append(
            f"only {len(shared_ordered)} shared samples across modalities (< 6); multivariate Procrustes/SVD will be underpowered"
        )

    for m in mod_records:
        if m["inferred_scale"] == "raw_counts":
            warnings.append(
                f"modality '{m['modality']}' appears to be raw integer counts; normalize (e.g. VST/log2CPM) or pass log1p before integration"
            )
        elif m["inferred_scale"] == "relative_compositional":
            warnings.append(
                f"modality '{m['modality']}' is compositional relative abundance; use 'clr' transformation in mo_pair_assoc / mo_factor_integration"
            )

    grp_counts = None
    if args.metadata and shared_ordered:
        meta = read_table(args.metadata)
        scol = args.sample_col or str(meta.columns[0])
        if scol not in meta.columns:
            fail(f"sample_col '{scol}' not found in metadata")
        meta[scol] = meta[scol].astype(str).str.strip()
        meta_sub = meta[meta[scol].isin(shared_set)]
        n_meta_match = int(meta_sub[scol].nunique())
        if n_meta_match < len(shared_ordered):
            warnings.append(
                f"metadata '{scol}' covers only {n_meta_match}/{len(shared_ordered)} shared omics samples"
            )
        if args.group_col:
            if args.group_col not in meta.columns:
                fail(f"group_col '{args.group_col}' not found in metadata")
            vc = meta_sub.drop_duplicates(subset=[scol])[args.group_col].astype(str).value_counts()
            grp_counts = {str(k): int(v) for k, v in vc.items()}

    emit(
        {
            "tool": "mo_inspect",
            "n_modalities": len(mod_records),
            "n_shared_samples": len(shared_ordered),
            "shared_samples_preview": shared_ordered[:15],
            "modalities": mod_records,
            "pairwise_feature_overlap": pairwise_feat,
            "group_counts_on_shared": grp_counts,
            "barcode_prefix_hint": barcode_hint,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
