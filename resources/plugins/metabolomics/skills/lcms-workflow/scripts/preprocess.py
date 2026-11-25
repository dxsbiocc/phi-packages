#!/usr/bin/env python3
"""mb_preprocess: impute -> normalize -> transform, in that fixed order."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import emit, fail, read_matrix, write_artifact

SKILL = "lcms-workflow"


def impute(X, how):
    # NaN and zero both count as missing (below detection limit).
    missing = X.isna() | (X == 0)
    if how == "none":
        return X, 0
    if how == "zero":
        filled = X.fillna(0.0)
        return filled, int(X.isna().sum().sum())
    n_missing = int(missing.sum().sum())
    obs = X.mask(missing)
    fill = obs.min(axis=1)
    if how == "half_min":
        fill = fill / 2.0
    filled = obs.T.fillna(fill).T
    return filled, int(n_missing - filled.isna().sum().sum())


def normalize(X, how, is_col):
    factors = pd.Series(1.0, index=X.columns)
    dropped = []
    if how == "none":
        return X, factors, dropped
    if how == "tic":
        factors = X.sum(axis=0)
    elif how == "median":
        factors = X.median(axis=0)
    elif how == "pqn":
        tic = X.sum(axis=0)
        ref = tic.sort_values().index[len(tic) // 2]
        ref_vals = X[ref]
        valid = ref_vals.notna() & (ref_vals > 0)
        ratios = X.loc[valid].div(ref_vals[valid], axis=0)
        factors = ratios.median(axis=0)
    elif how == "is":
        if not is_col:
            fail("normalize=is requires --is_col")
        if is_col not in X.columns:
            fail(f"is_col not found among sample columns: {is_col}")
        factors = X[is_col].copy()
        dropped = [is_col]
    if factors.isna().any() or (factors <= 0).any():
        bad = factors[factors.isna() | (factors <= 0)].index.tolist()
        fail(f"normalization factors invalid for samples: {bad}")
    X = X.div(factors, axis=1)
    if dropped:
        X = X.drop(columns=dropped)
    return X, factors, dropped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--impute", choices=["none", "half_min", "min", "zero"],
                    default="none")
    ap.add_argument("--normalize",
                    choices=["none", "tic", "median", "pqn", "is"],
                    default="none")
    ap.add_argument("--is_col")
    ap.add_argument("--transform", choices=["none", "log2", "log10"],
                    default="none")
    args = ap.parse_args()

    meta, X = read_matrix(args.input)
    if X.shape[1] == 0:
        fail("no numeric sample columns in input")
    X, n_imputed = impute(X, args.impute)
    X, factors, dropped = normalize(X, args.normalize, args.is_col)
    if args.transform != "none":
        if X.isna().any().any():
            fail("log transform needs complete data; impute missing values first")
        if (X <= 0).any().any():
            fail("log transform needs positive values; impute zeros/missing first "
                 "(impute=half_min or min)")
        X = np.log2(X) if args.transform == "log2" else np.log10(X)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df = pd.concat([meta.loc[X.index], X], axis=1)
    out_df.index.name = "feature"
    out_df.to_csv(out_path, sep="\t")

    report = pd.DataFrame({"normalization_factor": factors})
    report.index.name = "sample"
    report_path = out_path.with_name(out_path.stem + "-normalization-report.tsv")
    report.to_csv(report_path, sep="\t")

    params = {"impute": args.impute, "normalize": args.normalize,
              "transform": args.transform}
    if args.is_col:
        params["is_col"] = args.is_col
    artifacts = [
        write_artifact(out_path, "mb_preprocess", "Processed feature matrix",
                       skill=SKILL, inputs=[args.input], parameters=params,
                       shape=X.shape),
        write_artifact(report_path, "mb_preprocess",
                       "Per-sample normalization factors",
                       skill=SKILL, inputs=[args.input], parameters=params,
                       shape=report.shape),
    ]
    sem = (
        "log_transformed"
        if args.transform != "none"
        else ("normalized" if args.normalize != "none" else "raw_intensity")
    )
    emit({
        "input": args.input,
        "output": str(out_path),
        "impute": args.impute,
        "normalize": args.normalize,
        "transform": args.transform,
        "n_features": int(X.shape[0]),
        "n_samples": int(X.shape[1]),
        "n_imputed": n_imputed,
        "n_missing_remaining": int(X.isna().sum().sum()),
        "dropped_columns": dropped,
        "factors": {str(k): float(v) for k, v in factors.items()},
        "scientificContext": {
            str(out_path): {
                "valueSemantics": sem,
                "sampleIds": [str(c) for c in X.columns],
            }
        },
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()