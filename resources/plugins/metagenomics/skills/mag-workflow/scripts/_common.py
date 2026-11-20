"""Shared helpers for the mag-workflow script tools.

Every tool prints exactly one JSON object to stdout; diagnostics go to stderr.
"""
import json
import sys

import pandas as pd


def emit(obj):
    print(json.dumps(obj, ensure_ascii=False, default=float))
    sys.exit(0)


def fail(msg):
    print(json.dumps({"error": str(msg)}, ensure_ascii=False))
    sys.exit(1)


def read_tsv(path):
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep)


def bin_stem(name):
    from pathlib import Path

    s = Path(str(name)).name
    for ext in (".fasta", ".fna", ".fa", ".fsa"):
        if s.lower().endswith(ext):
            return s[: -len(ext)]
    return s


def setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt