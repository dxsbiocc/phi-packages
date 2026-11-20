#!/usr/bin/env python3
"""Shared helpers for the scanpy-workflow script toolkit.

Rules for this toolkit (mirrors seurat-workflow/_common.R):
  * stdout carries exactly one JSON object; every progress line goes to stderr.
  * a missing optional package is reported in the JSON result, never fatal.
  * arguments arrive as `--flag value` pairs built from the declared schema;
    a flag with no value is boolean True, a repeated flag collects a list.
"""

import importlib.util
import json
import os
import sys

import numpy as np


def script_dir():
    return os.path.dirname(os.path.abspath(__file__))


def info(*parts):
    print("[scanpy] " + "".join(str(p) for p in parts), file=sys.stderr)


def fail(*parts):
    message = "".join(str(p) for p in parts)
    print("[scanpy] error: " + message, file=sys.stderr)
    print(json.dumps({"error": message}, indent=2, ensure_ascii=False))
    sys.exit(1)


def _clean(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_clean(v) for v in value.tolist()]
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_clean(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def emit(result):
    print(json.dumps(_clean(result), indent=2, ensure_ascii=False))


# --- argument parsing ---------------------------------------------------------


def parse_flags(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    out = {}
    index = 0
    while index < len(argv):
        token = argv[index]
        if not token.startswith("--"):
            index += 1
            continue
        key = token[2:].replace("-", "_")
        has_value = index + 1 < len(argv) and not argv[index + 1].startswith("--")
        if not has_value:
            out[key] = True
            index += 1
        else:
            value = argv[index + 1]
            if key not in out:
                out[key] = value
            elif isinstance(out[key], list):
                out[key].append(value)
            else:
                out[key] = [out[key], value]
            index += 2
    return out


def flag_str(flags, name, default=None):
    value = flags.get(name)
    if value is None or value is True:
        return default
    if isinstance(value, list):
        value = value[0]
    return str(value)


def flag_req(flags, name):
    value = flag_str(flags, name)
    if not value:
        fail("missing required argument --" + name.replace("_", "-"))
    return value


def flag_int(flags, name, default=None):
    value = flag_str(flags, name)
    if value is None:
        return default
    try:
        return int(value)
    except ValueError:
        fail("--" + name.replace("_", "-") + " must be an integer, got " + value)


def flag_num(flags, name, default=None):
    value = flag_str(flags, name)
    if value is None:
        return default
    try:
        return float(value)
    except ValueError:
        fail("--" + name.replace("_", "-") + " must be a number, got " + value)


def flag_bool(flags, name, default=False):
    value = flags.get(name)
    if value is None:
        return default
    if value is True:
        return True
    if isinstance(value, list):
        value = value[0]
    return str(value).lower() in ("true", "1", "yes", "on")


def flag_vec(flags, name, default=()):
    value = flags.get(name)
    if value is None or value is True:
        return list(default)
    if isinstance(value, list):
        return [str(v) for v in value]
    return [str(value)]


def flag_nums(flags, name, default=()):
    values = flag_vec(flags, name, default)
    try:
        return [float(v) for v in values]
    except ValueError:
        fail("--" + name.replace("_", "-") + " takes numbers only")


# --- packages -----------------------------------------------------------------


def optional_module(name):
    if importlib.util.find_spec(name) is not None:
        return True
    info("optional package not installed: " + name)
    return False


def require_scanpy():
    try:
        import scanpy as sc
    except ImportError:
        fail("scanpy is not installed in this environment; use the plugin:scanpy environment")
    return sc


# --- objects ------------------------------------------------------------------


def detect_kind(path):
    if os.path.isdir(path):
        markers = ("features.tsv.gz", "genes.tsv.gz", "features.tsv", "genes.tsv")
        has_markers = any(os.path.exists(os.path.join(path, m)) for m in markers)
        has_matrix = os.path.exists(os.path.join(path, "matrix.mtx.gz")) or os.path.exists(
            os.path.join(path, "matrix.mtx")
        )
        if has_markers and has_matrix:
            return "10x-matrix-dir"
        if os.path.isdir(os.path.join(path, "filtered_feature_bc_matrix")):
            return "10x-matrix-dir"
        return "unknown-directory"
    lower = path.lower()
    if lower.endswith(".h5ad"):
        return "anndata"
    if lower.endswith(".h5"):
        return "10x-h5"
    if lower.endswith(".loom"):
        return "loom"
    if lower.endswith(".mtx") or lower.endswith(".mtx.gz"):
        return "matrix-market"
    if lower.endswith((".csv", ".tsv", ".txt")):
        return "dense-table"
    return "unknown-file"


def read_object(path, var_names="gene_symbols"):
    """Load any supported input as AnnData with a counts layer."""
    sc = require_scanpy()
    if not os.path.exists(path):
        fail("input not found: " + path)
    kind = detect_kind(path)
    info("input kind: " + kind)
    if kind == "anndata":
        adata = sc.io.read_h5ad(path)
    elif kind == "10x-h5":
        adata = sc.io.read_10x_h5(path)
        adata.var_names_make_unique()
    elif kind == "10x-matrix-dir":
        directory = path
        nested = os.path.join(path, "filtered_feature_bc_matrix")
        if os.path.isdir(nested):
            directory = nested
        info("reading 10x matrix directory " + directory)
        adata = sc.io.read_10x_mtx(directory, var_names=var_names)
        adata.var_names_make_unique()
    elif kind == "loom":
        adata = sc.read_loom(path)
    elif kind == "dense-table":
        separator = "," if path.lower().endswith(".csv") else "\t"
        adata = sc.io.read_csv(path, delimiter=separator)
    elif kind == "matrix-market":
        fail(
            "--input " + path + " is a bare Matrix Market file without feature/barcode "
            "names; pass its 10x directory instead so cell and gene identifiers survive"
        )
    elif kind == "unknown-directory":
        fail(
            "--input is a directory with no 10x matrix files (matrix.mtx + features.tsv); "
            "expected a Cell Ranger filtered_feature_bc_matrix or raw_feature_bc_matrix folder"
        )
    else:
        fail(
            "unsupported --input; expected a .h5ad, a 10x matrix directory or .h5, "
            "a .loom, or a .csv/.tsv count table"
        )
    if "counts" not in adata.layers:
        info("no counts layer; snapshotting X as counts")
        adata.layers["counts"] = adata.X.copy()
    return adata


def write_object(adata, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    adata.write_h5ad(path)
    info("wrote " + path)
    return path


def ensure_dir(path, fallback):
    target = path if path else fallback
    os.makedirs(target, exist_ok=True)
    return target


# --- cluster levels -----------------------------------------------------------


def sorted_levels(values):
    uniques = sorted({str(v) for v in values})
    try:
        return sorted(uniques, key=float)
    except ValueError:
        return uniques