"""Shared I/O, statistics, nonlinear-fit, and artifact utilities for wetlab-assay (wa_*)."""
from __future__ import annotations

import hashlib
import json
import math
import os
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd


def emit(payload: dict[str, Any]) -> None:
    """Write JSON result payload to stdout and exit 0."""
    sys.stdout.write(json.dumps(clean_json(payload), ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()
    raise SystemExit(0)


def fail(message: str, code: int = 1) -> None:
    """Write JSON error payload to stdout and exit non-zero."""
    sys.stdout.write(json.dumps({"error": str(message)}, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()
    raise SystemExit(code)


def clean_json(obj: Any) -> Any:
    """Recursively convert numpy types and non-finite floats into strict JSON values."""
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        val = float(obj)
        return val if math.isfinite(val) else None
    if isinstance(obj, dict):
        return {str(k): clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean_json(v) for v in obj]
    return str(obj)


def read_table(path: str | Path) -> pd.DataFrame:
    """Read TSV/CSV input table with delimiter auto-detection."""
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
    try:
        df = pd.read_csv(p, sep=sep, comment="#", low_memory=False)
    except Exception as exc:
        fail(f"Failed to read table '{path}': {exc}")
    if df.empty:
        fail(f"Input table is empty: {path}")
    df.columns = [str(c).strip() for c in df.columns]
    return df


def resolve_column(df: pd.DataFrame, candidates: list[str], label: str, required: bool = True) -> str | None:
    """Resolve column name case-insensitively against a list of candidate names."""
    lower_map = {str(c).strip().lower(): str(c).strip() for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        hit = lower_map.get(cand.lower())
        if hit is not None:
            return hit
    if required:
        fail(f"Missing required {label} column (expected one of {candidates}; found {list(df.columns)})")
    return None


def parse_str_list(raw: list[str] | None) -> list[str]:
    """Flatten repeated or comma-separated CLI string arguments while preserving order."""
    out: list[str] = []
    if not raw:
        return out
    seen: set[str] = set()
    for item in raw:
        for part in str(item).split(","):
            tok = part.strip()
            if tok and tok not in seen:
                seen.add(tok)
                out.append(tok)
    return out


def bh_fdr(pvals: Any) -> list[float]:
    """Benjamini-Hochberg FDR adjustment preserving NaNs as 1.0."""
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
    return [float(x) for x in out]


def compute_aicc(rss: float, n_obs: int, n_params: int) -> float:
    """Compute small-sample corrected AIC (AICc) from residual sum of squares."""
    if n_obs <= 0:
        return float("inf")
    rss_safe = max(float(rss), 1e-30)
    aic = n_obs * math.log(rss_safe / n_obs) + 2.0 * n_params
    denom = n_obs - n_params - 1
    if denom <= 0:
        return float(aic + 2.0 * n_params * (n_params + 1))
    return float(aic + (2.0 * n_params * (n_params + 1)) / denom)


def akaike_weights(aicc_values: list[float]) -> list[float]:
    """Compute normalized Akaike weights from a list of AICc scores."""
    arr = np.asarray(aicc_values, dtype=float)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return [0.0 for _ in aicc_values]
    min_aicc = float(np.min(arr[finite]))
    delta = np.where(finite, arr - min_aicc, np.inf)
    raw = np.where(finite, np.exp(-0.5 * np.clip(delta, 0.0, 700.0)), 0.0)
    total = float(np.sum(raw))
    if total <= 0.0:
        return [0.0 for _ in aicc_values]
    return [float(x / total) for x in raw]


def jacobian_se_ci(res, n_obs: int, alpha: float = 0.05) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Estimate parameter SEs and 95% t-CIs from a scipy.optimize.least_squares result."""
    from scipy import stats

    n_p = len(res.x)
    df = max(1, n_obs - n_p)
    rss = float(np.sum(res.fun ** 2))
    s2 = rss / df
    try:
        jtj = res.jac.T @ res.jac
        cov = np.linalg.pinv(jtj) * s2
        se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    except Exception:
        se = np.full(n_p, np.nan, dtype=float)
    t_crit = float(stats.t.ppf(1.0 - alpha / 2.0, df))
    ci_lo = res.x - t_crit * se
    ci_hi = res.x + t_crit * se
    return se, ci_lo, ci_hi


def setup_matplotlib():
    """Configure headless Agg backend and publication plot defaults."""
    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 180,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.5,
            "legend.fontsize": 8.5,
        }
    )
    return plt


def write_tsv(df: pd.DataFrame, path: str | Path, tool: str, title: str, **kwargs: Any) -> str:
    """Write DataFrame to TSV and emit its .phi-artifact.json descriptor."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, sep="\t", index=False)
    return write_artifact(p, tool=tool, title=title, shape=df.shape, **kwargs)


def write_artifact(
    path: str | Path,
    tool: str,
    title: str,
    *,
    skill: str = "assay-workflow",
    inputs: tuple[str | None, ...] | list[str | None] = (),
    parameters: dict[str, Any] | None = None,
    shape: tuple[int, int] | dict[str, int] | None = None,
    scientific_context: dict[str, Any] | None = None,
) -> str:
    """Write Contract-1.0.0 .phi-artifact.json sidecar for tables and figures."""
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

    descriptor: dict[str, Any] = {
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
        block: dict[str, Any] = {"format": ext[1:]}
        if ext == ".png":
            with open(p, "rb") as fh:
                hdr = fh.read(24)
            if hdr[:8] == b"\x89PNG\r\n\x1a\n":
                w, h = struct.unpack(">II", header_bytes := hdr[16:24])
                del header_bytes
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
        descriptor["scientificContext"] = clean_json(scientific_context)

    cwd = Path.cwd().resolve()
    prov_in = []
    for raw in inputs:
        if not raw:
            continue
        ip = Path(raw)
        try:
            rel = ip.resolve().relative_to(cwd)
        except ValueError:
            rel = Path(ip.name)
        entry: dict[str, Any] = {"path": str(rel)}
        if ip.is_file():
            entry["sha256"] = hashlib.sha256(ip.read_bytes()).hexdigest()
        prov_in.append(entry)
    if prov_in:
        descriptor["provenance"]["inputs"] = prov_in
    if parameters:
        descriptor["provenance"]["parameters"] = clean_json(parameters)

    target = p.with_name(p.name + ".phi-artifact.json")
    tmp = target.with_name(target.name + ".tmp")
    tmp.write_text(json.dumps(descriptor, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, target)
    return str(p)
