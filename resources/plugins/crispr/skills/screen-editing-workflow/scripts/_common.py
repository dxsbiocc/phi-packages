"""Shared utilities for CRISPR pooled screen QC, alpha-RRA hit calling, and editing quantification."""
from __future__ import annotations

import hashlib
import json
import os
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def emit(payload: dict) -> None:
    sys.stdout.write(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()


def fail(message: str, code: int = 1) -> "NoReturn":
    sys.stdout.write(json.dumps({"error": message}, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()
    raise SystemExit(code)


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def read_table(path: str | Path, *, header: int | None = 0):
    import io
    import pandas as pd

    p = Path(path)
    if not p.is_file():
        fail(f"Input file not found: {path}")
    name = p.name.lower()
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        sep = ","
    elif name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".txt"):
        sep = "\t"
    else:
        with open(p, "rt", encoding="utf-8", errors="replace") as fh:
            first = ""
            for line in fh:
                if not line.startswith("#") and line.strip():
                    first = line
                    break
        sep = "\t" if "\t" in first else ","
    # Preserve headers containing '#Reads' (CRISPResso2) while skipping leading full-line comments
    if name.endswith(".gz"):
        return pd.read_csv(p, sep=sep, header=header)
    lines: list[str] = []
    with open(p, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not lines and line.lstrip().startswith("#") and "\t" not in line:
                continue
            lines.append(line)
    return pd.read_csv(io.StringIO("".join(lines)), sep=sep, header=header)


def pick_col(
    cols: list[str], requested: str, candidates: list[str], *, required: bool = True
) -> str | None:
    if requested:
        if requested not in cols:
            fail(f"Requested column '{requested}' not found in columns: {cols}")
        return requested
    lower_map = {str(c).lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    if required:
        fail(f"Could not auto-detect column from {candidates} in columns: {cols}")
    return None


def parse_col_list(raw_items: list[str] | None) -> list[str]:
    if not raw_items:
        return []
    out: list[str] = []
    for item in raw_items:
        for tok in str(item).split(","):
            tok = tok.strip()
            if tok and tok not in out:
                out.append(tok)
    return out


def bh_fdr(pvals) -> list[float]:
    import numpy as np

    arr = np.asarray(pvals, dtype=float)
    n = arr.size
    if n == 0:
        return []
    safe = np.where(np.isfinite(arr), np.clip(arr, 0.0, 1.0), 1.0)
    order = np.argsort(safe)
    ranked = safe[order] * n / np.arange(1, n + 1, dtype=float)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out.tolist()


def deseq_size_factors(mat) -> "np.ndarray":
    """Compute DESeq2/MAGeCK median-of-ratios size factors on counts (n_sgrnas x n_samples)."""
    import numpy as np

    arr = np.asarray(mat, dtype=float)
    arr = np.maximum(arr, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_counts = np.where(arr > 0, np.log(arr), np.nan)
    geo_means = np.exp(np.nanmean(log_counts, axis=1))
    valid = np.isfinite(geo_means) & (geo_means > 0) & np.all(arr > 0, axis=1)
    if valid.sum() < 10:
        # Relax all-positive requirement if library has zero dropouts
        valid = np.isfinite(geo_means) & (geo_means > 0)
    if valid.sum() == 0:
        totals = arr.sum(axis=0)
        med = float(np.median(totals[totals > 0])) if np.any(totals > 0) else 1.0
        return np.where(totals > 0, totals / med, 1.0)
    ratios = arr[valid] / geo_means[valid, None]
    ratios = np.where(ratios > 0, ratios, np.nan)
    sf = np.nanmedian(ratios, axis=0)
    sf = np.where(np.isfinite(sf) & (sf > 0), sf, 1.0)
    return sf


def setup_matplotlib():
    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def write_artifact(
    path: str | Path,
    tool: str,
    title: str,
    *,
    skill: str = "screen-editing-workflow",
    inputs=(),
    parameters: dict | None = None,
    shape: tuple[int, int] | dict | None = None,
    scientific_context: dict | None = None,
) -> str:
    p = Path(path)
    ext = p.suffix.lower()
    figures = {".png": "image/png", ".pdf": "application/pdf", ".svg": "image/svg+xml"}
    tables = {".tsv": "text/tab-separated-values", ".csv": "text/csv"}
    if ext in figures:
        kind, media = "figure", figures[ext]
    elif ext in tables:
        kind, media = "table", tables[ext]
    else:
        raise ValueError(f"Unsupported artifact extension: {ext}")

    descriptor: dict = {
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
        block: dict = {"format": ext[1:]}
        if ext == ".png":
            with open(p, "rb") as fh:
                header = fh.read(24)
            if header[:8] == b"\x89PNG\r\n\x1a\n":
                w, h = struct.unpack(">II", header[16:24])
                block["widthPx"], block["heightPx"] = int(w), int(h)
        descriptor["figure"] = block
    elif shape is not None:
        if isinstance(shape, dict):
            rows_n = int(shape.get("rows", 0))
            cols_n = int(shape.get("columns", shape.get("cols", 0)))
        else:
            rows_n, cols_n = int(shape[0]), int(shape[1])
        descriptor["table"] = {"rows": rows_n, "columns": cols_n}

    if scientific_context:
        descriptor["scientificContext"] = scientific_context

    cwd = Path.cwd().resolve()
    prov_inputs = []
    for ip in inputs:
        if not ip:
            continue
        ip_p = Path(ip)
        try:
            rel = ip_p.resolve().relative_to(cwd)
        except ValueError:
            continue
        entry = {"path": str(rel)}
        if ip_p.is_file():
            entry["sha256"] = hashlib.sha256(ip_p.read_bytes()).hexdigest()
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
