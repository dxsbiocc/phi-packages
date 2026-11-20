"""Shared helpers for the lcms-workflow script tools.

Every tool prints exactly one JSON object to stdout; diagnostics go to stderr.
"""
import json
import re
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


MZ_ALIASES = {"mz", "m/z", "mass", "mzmed", "average mz", "monoisotopic mass"}
RT_ALIASES = {
    "rt", "retention time", "retention time (min)", "rtmed", "average rt(min)",
    "average rt", "rt (min)", "retention_time",
}
QC_PATTERN = re.compile(r"(qc|pool)", re.I)
BLANK_PATTERN = re.compile(r"blank", re.I)


def split_matrix(df):
    """Split a feature table into (feature_meta, intensities).

    First column is the feature ID. Columns named mz/rt (vendor aliases
    included) become feature metadata; remaining numeric columns are samples.
    Non-numeric annotation columns are kept in feature_meta, never in stats.
    """
    df = df.copy()
    id_col = df.columns[0]
    df = df.set_index(id_col)
    df.index = df.index.astype(str)
    meta = {}
    samples = []
    for col in df.columns:
        key = str(col).strip().lower()
        if key in MZ_ALIASES:
            meta["mz"] = pd.to_numeric(df[col], errors="coerce")
        elif key in RT_ALIASES:
            meta["rt"] = pd.to_numeric(df[col], errors="coerce")
        elif pd.api.types.is_numeric_dtype(df[col]):
            samples.append(col)
        else:
            meta[str(col)] = df[col]
    intensities = df[samples].apply(pd.to_numeric, errors="coerce")
    meta_df = pd.DataFrame(meta, index=df.index) if meta else pd.DataFrame(index=df.index)
    return meta_df, intensities


def read_matrix(path):
    return split_matrix(read_tsv(path))


def qc_samples(columns):
    return [c for c in columns if QC_PATTERN.search(str(c)) and not BLANK_PATTERN.search(str(c))]


def setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def write_artifact(path, tool, title, *, skill, inputs=(), parameters=None, shape=None):
    """Write <path>.phi-artifact.json beside a produced file (artifact contract 1.0.0).

    Returns str(path) so callers can build the artifacts list inline.
    """
    import hashlib
    import os
    import struct
    from datetime import datetime, timezone
    from pathlib import Path

    p = Path(path)
    ext = p.suffix.lower()
    figures = {".png": "image/png", ".pdf": "application/pdf", ".svg": "image/svg+xml"}
    tables = {".tsv": "text/tab-separated-values", ".csv": "text/csv"}
    if ext in figures:
        kind, media = "figure", figures[ext]
    elif ext in tables:
        kind, media = "table", tables[ext]
    else:
        raise ValueError(f"unsupported artifact extension: {ext}")
    descriptor = {
        "contractVersion": "1.0.0",
        "kind": kind,
        "file": p.name,
        "mediaType": media,
        "title": title[:200],
        "provenance": {
            "createdAt": datetime.now(timezone.utc).isoformat(),
            "tool": tool,
            "skill": skill,
        },
    }
    if kind == "figure":
        block = {"format": ext[1:]}
        if ext == ".png":
            with open(p, "rb") as fh:
                header = fh.read(24)
            if header[:8] == b"\x89PNG\r\n\x1a\n":
                w, h = struct.unpack(">II", header[16:24])
                block["widthPx"], block["heightPx"] = int(w), int(h)
        descriptor["figure"] = block
    elif shape is not None:
        descriptor["table"] = {"rows": int(shape[0]), "columns": int(shape[1])}
    cwd = Path.cwd().resolve()
    prov_inputs = []
    for ip in inputs:
        ip = Path(ip)
        try:
            rel = ip.resolve().relative_to(cwd)
        except ValueError:
            continue
        entry = {"path": str(rel)}
        if ip.is_file():
            entry["sha256"] = hashlib.sha256(ip.read_bytes()).hexdigest()
        prov_inputs.append(entry)
    if prov_inputs:
        descriptor["provenance"]["inputs"] = prov_inputs
    if parameters:
        descriptor["provenance"]["parameters"] = parameters
    target = p.with_name(p.name + ".phi-artifact.json")
    tmp = target.with_name(target.name + ".tmp")
    tmp.write_text(json.dumps(descriptor, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, target)
    return str(p)