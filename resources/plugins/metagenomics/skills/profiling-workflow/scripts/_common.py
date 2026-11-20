"""Shared helpers for the profiling-workflow script tools.

Every tool prints exactly one JSON object to stdout; diagnostics go to stderr.
"""
import json
import sys

import numpy as np
import pandas as pd


def emit(obj):
    print(json.dumps(obj, ensure_ascii=False, default=float))
    sys.exit(0)


def fail(msg):
    print(json.dumps({"error": str(msg)}, ensure_ascii=False))
    sys.exit(1)


def read_table(path):
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, index_col=0)


def read_matrix(path):
    df = read_table(path)
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.isna().any().any():
        bad = df.columns[df.isna().any()].tolist()
        fail(f"non-numeric values in matrix columns: {bad}")
    return df.fillna(0.0)


def read_metadata(path, samples):
    df = read_table(path)
    df.index = df.index.astype(str)
    missing = [s for s in samples if s not in df.index]
    if missing:
        fail(f"metadata missing samples: {missing}")
    return df.loc[samples]


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def sample_stem(path):
    from pathlib import Path

    name = Path(path).name
    for suffix in (
        ".k2report", ".kreport", ".kraken2.report", ".bracken", ".bracken.report",
        ".profile", ".metaphlan", ".mpa", ".report", ".txt", ".tsv", ".csv",
    ):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem