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