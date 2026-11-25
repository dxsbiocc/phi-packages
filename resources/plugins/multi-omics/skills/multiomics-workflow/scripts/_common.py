#!/usr/bin/env python3
"""Shared utilities for the multi-omics integration plugin (toolPrefix: mo)."""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def emit(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)


def fail(message: str, code: int = 1) -> None:
    print(json.dumps({"error": message}, indent=2, ensure_ascii=False))
    sys.exit(code)


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    n = int(valid.sum())
    if n == 0:
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / np.arange(1, n + 1, dtype=float)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    res = np.empty(n, dtype=float)
    res[order] = adj
    out[valid] = res
    return out


def parse_labeled_paths(specs: list[str], prefix: str = "modality") -> list[tuple[str, Path]]:
    out: list[tuple[str, Path]] = []
    seen: set[str] = set()
    for idx, raw in enumerate(specs, start=1):
        for token in str(raw).split(","):
            s = token.strip()
            if not s:
                continue
            if "=" in s and not Path(s).is_file():
                lbl, p_str = s.split("=", 1)
                lbl = lbl.strip() or f"{prefix}_{idx}"
            else:
                p_str = s
                lbl = Path(s).stem
            base_lbl = lbl
            k = 2
            while lbl in seen:
                lbl = f"{base_lbl}_{k}"
                k += 1
            seen.add(lbl)
            p = Path(p_str)
            if not p.is_file():
                fail(f"input file not found for '{lbl}': {p_str}")
            out.append((lbl, p))
    return out


def read_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.is_file():
        fail(f"table file not found: {p}")
    sep = "," if p.name.lower().endswith(".csv") else "\t"
    try:
        df = pd.read_csv(p, sep=sep)
    except Exception as exc:
        fail(f"failed to parse {p}: {exc}")
    if df.empty:
        fail(f"table is empty: {p}")
    return df


def read_matrix(path: str | Path) -> pd.DataFrame:
    """Read a feature x sample matrix (first column = feature IDs, remaining columns = samples)."""
    p = Path(path)
    if not p.is_file():
        fail(f"matrix file not found: {p}")
    sep = "," if p.name.lower().endswith(".csv") else "\t"
    try:
        df = pd.read_csv(p, sep=sep, index_col=0)
    except Exception as exc:
        fail(f"failed to parse matrix {p}: {exc}")
    if df.empty or df.shape[1] == 0:
        fail(f"matrix {p} has 0 features or 0 sample columns")
    df.index = df.index.astype(str).str.strip()
    df.columns = [str(c).strip() for c in df.columns]
    num = df.apply(pd.to_numeric, errors="coerce")
    if num.index.duplicated().any():
        num = num.groupby(level=0, sort=False).mean()
    return num


def infer_matrix_scale(df: pd.DataFrame) -> str:
    arr = df.to_numpy(dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return "empty"
    has_neg = bool(np.any(finite < -1e-9))
    is_int = bool(np.all(np.abs(finite - np.round(finite)) < 1e-6)) and not has_neg
    col_sums = np.nansum(arr, axis=0)
    valid_sums = col_sums[col_sums > 0]
    if not has_neg and valid_sums.size > 0:
        med_sum = float(np.median(valid_sums))
        if (0.95 <= med_sum <= 1.05 and float(np.max(finite)) <= 1.01) or (
            95.0 <= med_sum <= 105.0 and float(np.max(finite)) <= 100.5
        ):
            return "relative_compositional"
    if is_int and float(np.max(finite)) >= 20.0:
        return "raw_counts"
    if has_neg or float(np.max(finite)) <= 35.0:
        return "log_continuous"
    return "continuous"


def transform_block(df: pd.DataFrame, method: str = "auto") -> tuple[pd.DataFrame, str]:
    """Return sample-aligned transformed feature x sample DataFrame and applied method name."""
    chosen = method
    if chosen == "auto":
        scale = infer_matrix_scale(df)
        if scale == "relative_compositional":
            chosen = "clr"
        elif scale in ("raw_counts", "continuous"):
            chosen = "log1p"
        else:
            chosen = "none"

    x = df.to_numpy(dtype=float).copy()
    # Impute residual NaNs per feature with row median (or 0 if all NaN)
    row_med = np.nanmedian(np.where(np.isfinite(x), x, np.nan), axis=1)
    row_med = np.where(np.isfinite(row_med), row_med, 0.0)
    inds = np.where(~np.isfinite(x))
    if inds[0].size > 0:
        x[inds] = row_med[inds[0]]

    if chosen == "log1p":
        x = np.log2(np.clip(x, 0.0, None) + 1.0)
    elif chosen == "clr":
        x_pos = np.clip(x, 0.0, None)
        pos_vals = x_pos[x_pos > 0]
        pseudo = min(1e-6, 0.5 * float(np.min(pos_vals))) if pos_vals.size > 0 else 1e-6
        lx = np.log(x_pos + pseudo)
        # CLR across features within each sample (column)
        lx = lx - np.mean(lx, axis=0, keepdims=True)
        x = lx
    elif chosen == "zscore":
        mu = np.mean(x, axis=1, keepdims=True)
        sd = np.std(x, axis=1, ddof=1, keepdims=True)
        sd = np.where(sd > 1e-12, sd, 1.0)
        x = (x - mu) / sd

    return pd.DataFrame(x, index=df.index, columns=df.columns), chosen


def read_metadata_aligned(
    path: str | Path, samples: list[str], sample_col: str | None = None
) -> pd.DataFrame:
    meta = read_table(path)
    scol = sample_col or str(meta.columns[0])
    if scol not in meta.columns:
        fail(f"sample_col '{scol}' not found in metadata columns: {list(meta.columns)}")
    meta = meta.copy()
    meta[scol] = meta[scol].astype(str).str.strip()
    meta = meta.drop_duplicates(subset=[scol], keep="first").set_index(scol)
    missing = [s for s in samples if s not in meta.index]
    if missing:
        fail(f"{len(missing)} shared omics sample(s) missing from metadata '{scol}': {missing[:5]}")
    return meta.loc[samples]


def write_artifact(
    path: str | Path,
    tool: str,
    title: str,
    *,
    inputs: tuple[str | Path | None, ...] | list[str | Path | None] = (),
    parameters: dict[str, Any] | None = None,
    kind: str | None = None,
    shape: tuple[int, int] | None = None,
) -> str:
    p = Path(path)
    ext = p.suffix.lower()
    if kind is None:
        if ext in (".png", ".svg", ".pdf"):
            kind = "figure"
        else:
            kind = "table"

    media_map = {
        ".png": "image/png",
        ".svg": "image/svg+xml",
        ".pdf": "application/pdf",
        ".csv": "text/csv",
        ".tsv": "text/tab-separated-values",
    }
    media = media_map.get(ext, "text/tab-separated-values")
    desc: dict[str, Any] = {
        "contractVersion": "1.0.0",
        "kind": kind,
        "file": p.name,
        "mediaType": media,
        "title": title[:200],
        "provenance": {
            "createdAt": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "tool": tool,
            "skill": "multiomics-workflow",
        },
    }
    if parameters:
        desc["provenance"]["parameters"] = parameters
    if kind == "figure":
        fmt = ext.lstrip(".")
        desc["figure"] = {"format": fmt if fmt in ("png", "pdf", "svg") else "png"}
    elif kind == "table" and shape is not None:
        desc["table"] = {"rows": int(shape[0]), "columns": int(shape[1])}

    cwd = Path.cwd().resolve()
    prov_in = []
    for raw in inputs:
        if not raw:
            continue
        ip = Path(raw)
        try:
            rel = ip.resolve().relative_to(cwd)
        except ValueError:
            continue
        entry: dict[str, Any] = {"path": str(rel)}
        if ip.is_file():
            entry["sha256"] = hashlib.sha256(ip.read_bytes()).hexdigest()
        prov_in.append(entry)
    if prov_in:
        desc["provenance"]["inputs"] = prov_in

    sidecar = p.with_name(p.name + ".phi-artifact.json")
    tmp = sidecar.with_name(sidecar.name + ".tmp")
    tmp.write_text(json.dumps(desc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, sidecar)
    return str(p)
