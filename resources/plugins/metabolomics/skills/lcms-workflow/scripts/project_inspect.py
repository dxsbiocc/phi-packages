#!/usr/bin/env python3
"""mb_inspect: identify metabolomics inputs and report what they support."""
import argparse
from pathlib import Path

import pandas as pd

from _common import (BLANK_PATTERN, MZ_ALIASES, RT_ALIASES, emit, fail,
                     qc_samples, read_tsv, split_matrix)

TABLE_EXTS = {".tsv", ".csv", ".txt"}
SKIP_SUFFIXES = (".phi-artifact.json",)


def classify_text(path):
    head = path.read_text(errors="replace")[:200000]
    upper = head.upper()
    if "BEGIN IONS" in upper:
        n = upper.count("BEGIN IONS")
        return {"kind": "mgf_spectra", "n_spectra": n}
    if ("NAME:" in upper or "NAME=" in upper) and (
            "NUM PEAKS:" in upper or "NUM PEAKS=" in upper or "PRECURSORMZ" in upper):
        n = upper.count("NUM PEAKS:") + upper.count("NUM PEAKS=")
        return {"kind": "msp_library", "n_entries": n}
    return None


def feature_format(df):
    cols = [str(c).strip().lower() for c in df.columns]
    if "alignment id" in cols and "average mz" in cols:
        return "msdial"
    if "row id" in cols or ("m/z" in cols and "retention time" in cols):
        return "mzmine"
    if any(c in MZ_ALIASES for c in cols) and any(c in RT_ALIASES for c in cols):
        return "xcms"
    if (len(cols) > 2 and cols[1] in {"label", "class", "group"}
            and any(pd.api.types.is_numeric_dtype(df[c]) for c in df.columns[2:])):
        return "metaboanalyst"
    return "generic"


def inspect_table(path):
    try:
        df = read_tsv(path)
    except Exception as exc:
        return {"kind": "unknown", "note": f"not parseable as a table: {exc}"}
    if df.shape[1] == 1:
        col = df.iloc[:, 0]
        nums = pd.to_numeric(col, errors="coerce")
        if nums.notna().all():
            return {"kind": "mass_list", "n_values": len(df)}
        return {"kind": "id_list", "n_values": len(df)}
    fmt = feature_format(df)
    if fmt == "metaboanalyst":
        numeric = [c for c in df.columns[2:]
                   if pd.api.types.is_numeric_dtype(df[c])]
        return {
            "kind": "feature_table", "format": fmt,
            "orientation": "samples-as-rows",
            "n_samples": len(df), "n_features": len(numeric),
            "note": "MetaboAnalyst layout: transpose to features-as-rows before use",
        }
    if df.shape[1] == 2:
        first = df.iloc[:, 0].astype(str)
        if first.duplicated().any():
            return {"kind": "pathway_mapping", "n_rows": len(df),
                    "columns": [str(c) for c in df.columns]}
    meta, X = split_matrix(df)
    if X.shape[1] == 0:
        lower = {str(c).strip().lower() for c in df.columns}
        if ("mass" in lower or "mz" in lower or "m/z" in lower) and (
                "compound_id" in lower or "id" in lower or "name" in lower):
            return {"kind": "compound_db", "n_rows": len(df),
                    "columns": [str(c) for c in df.columns]}
        nonnum = [c for c in df.columns[1:]
                  if not pd.api.types.is_numeric_dtype(df[c])]
        if nonnum and df.shape[1] <= 12:
            return {"kind": "metadata", "n_rows": len(df),
                    "columns": [str(c) for c in df.columns]}
        return {"kind": "unknown",
                "note": "no numeric sample columns found"}
    lower = {str(c).strip().lower() for c in df.columns}
    if ("mass" in lower or "mz" in lower) and (
            "compound_id" in lower or "name" in lower) and X.shape[1] <= 3:
        return {"kind": "compound_db", "n_rows": len(df),
                "columns": [str(c) for c in df.columns]}
    n_cells = X.size
    missing = float(X.isna().sum().sum()) / n_cells if n_cells else 0.0
    zeros = float((X == 0).sum().sum()) / n_cells if n_cells else 0.0
    rt_unit = "unknown"
    if "rt" in meta.columns:
        rt_med = float(meta["rt"].median())
        rt_unit = "sec" if rt_med > 60 else "min"
    qc = qc_samples(X.columns)
    blanks = [c for c in X.columns if BLANK_PATTERN.search(str(c))]
    ready = ["preprocess", "plot"]
    if qc:
        ready.append("qc")
    if X.shape[1] >= 2:
        ready += ["multivariate", "univariate"]
    if "mz" in meta.columns:
        ready.append("annotate")
    return {
        "kind": "feature_table", "format": fmt,
        "orientation": "features-as-rows",
        "n_features": int(X.shape[0]), "n_samples": int(X.shape[1]),
        "samples": [str(c) for c in X.columns],
        "has_mz": "mz" in meta.columns, "has_rt": "rt" in meta.columns,
        "rt_unit": rt_unit,
        "missing_rate": missing, "zero_rate": zeros,
        "qc_samples": [str(c) for c in qc],
        "blank_samples": [str(c) for c in blanks],
        "ready_steps": ready,
    }


def inspect_file(path):
    ext = path.suffix.lower()
    if ext in {".msp", ".mgf"} or ext not in TABLE_EXTS:
        text_kind = classify_text(path)
        if text_kind:
            return text_kind
        if ext in TABLE_EXTS:
            return inspect_table(path)
        return {"kind": "unknown", "note": f"unsupported extension: {ext}"}
    text_kind = classify_text(path)
    if text_kind:
        return text_kind
    return inspect_table(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    args = ap.parse_args()
    root = Path(args.input)
    if not root.exists():
        fail(f"input does not exist: {root}")
    if root.is_dir():
        paths = sorted(p for p in root.iterdir()
                       if p.is_file() and not p.name.startswith(".")
                       and not p.name.endswith(SKIP_SUFFIXES))
    else:
        paths = [root]
    files = []
    for p in paths:
        entry = {"path": str(p)}
        entry.update(inspect_file(p))
        files.append(entry)
    emit({"input": str(root), "files": files})


if __name__ == "__main__":
    main()