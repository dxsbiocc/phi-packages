#!/usr/bin/env python3
"""Shared utilities for the DNA methylation downstream plugin (toolPrefix: meth)."""

from __future__ import annotations

import datetime
import gzip
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma

COORD_RE = re.compile(r"^(chr[0-9A-Za-z_.-]+|[0-9]{1,2}|[XYM]|MT)[:_](\d+)(?:[-_](\d+))?$")
COV_EXTS = (".bismark.cov.gz", ".bismark.cov", ".cov.gz", ".cov", ".bedgraph.gz", ".bedgraph")


def emit(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)


def fail(message: str, code: int = 1) -> None:
    print(json.dumps({"error": message}, indent=2, ensure_ascii=False))
    sys.exit(code)


def open_text(path: str | Path):
    p = str(path)
    if p.lower().endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def clean_sample_stem(path: Path) -> str:
    name = path.name
    for ext in (".bismark.cov.gz", ".bismark.cov", ".cov.gz", ".cov", ".bedGraph.gz", ".bedGraph", ".bedgraph.gz", ".bedgraph", ".tsv", ".csv"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return name


def is_coverage_file(path: Path) -> bool:
    low = path.name.lower()
    return any(low.endswith(ext) for ext in COV_EXTS)


def parse_cpg_coord(cpg_id: str) -> tuple[str, int, int] | None:
    m = COORD_RE.match(str(cpg_id).strip())
    if not m:
        return None
    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3)) if m.group(3) else start
    return chrom, start, end


def parse_coverage_records(path: Path, max_lines: int | None = None) -> tuple[str, list[tuple[str, float, int]]]:
    """Parse Bismark .cov(.gz) or MethylDackel .bedGraph(.gz).
    Returns (detected_format, [(cpg_id 'chr:pos', beta_0_to_1, total_depth), ...]).
    Bismark cov: chrom, start(1-based), end, pct_meth(0..100), count_meth, count_unmeth
    MethylDackel bedGraph: chrom, start(0-based), end(1-based), pct_meth, count_meth, count_unmeth
    """
    records: list[tuple[str, float, int]] = []
    fmt = "bismark_cov"
    low_name = path.name.lower()
    if "bedgraph" in low_name:
        fmt = "methyldackel_bedgraph"

    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.lower().startswith("track"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 6:
                parts = line.rstrip("\r\n").split()
            if len(parts) < 6:
                continue
            chrom = parts[0].strip()
            try:
                s_pos = int(parts[1])
                e_pos = int(parts[2])
                pct = float(parts[3])
                n_m = int(round(float(parts[4])))
                n_u = int(round(float(parts[5])))
            except ValueError:
                continue
            if e_pos == s_pos + 1 and fmt == "methyldackel_bedgraph":
                pos_1based = e_pos
            elif e_pos == s_pos + 1 and s_pos >= 0 and "bismark" not in low_name:
                fmt = "methyldackel_bedgraph"
                pos_1based = e_pos
            else:
                pos_1based = s_pos
            tot = n_m + n_u
            if tot <= 0:
                continue
            beta = n_m / float(tot) if (n_m >= 0 and n_u >= 0) else np.clip(pct / 100.0, 0.0, 1.0)
            records.append((f"{chrom}:{pos_1based}", float(np.clip(beta, 0.0, 1.0)), tot))
            if max_lines is not None and len(records) >= max_lines:
                break
    return fmt, records


def read_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.is_file():
        fail(f"table file not found: {p}")
    sep = "," if p.name.lower().endswith(".csv") else "\t"
    return pd.read_csv(p, sep=sep)


def read_methylation_matrix(path: str | Path) -> tuple[pd.DataFrame, str]:
    """Read CpG x sample matrix and return (numeric_df, 'beta' | 'm_value')."""
    p = Path(path)
    if not p.is_file():
        fail(f"matrix file not found: {p}")
    sep = "," if p.name.lower().endswith(".csv") else "\t"
    df = pd.read_csv(p, sep=sep, index_col=0)
    if df.empty or df.shape[1] == 0:
        fail(f"matrix {p} has 0 CpGs or 0 sample columns")
    df.index = df.index.astype(str).str.strip()
    df.columns = [str(c).strip() for c in df.columns]
    num = df.apply(pd.to_numeric, errors="coerce")
    arr = num.to_numpy(dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        fail(f"matrix {p} contains 0 finite numeric values")
    vmin, vmax = float(np.min(finite)), float(np.max(finite))
    if vmin >= -1e-6 and vmax <= 1.000001:
        scale = "beta"
    elif vmin >= -1e-6 and vmax <= 100.001 and float(np.median(finite)) > 1.5:
        # Percent methylation 0..100 -> rescale to 0..1 beta
        num = num / 100.0
        scale = "beta"
    else:
        scale = "m_value"
    return num, scale


def beta_to_m(beta: np.ndarray, eps: float = 1e-4) -> np.ndarray:
    b = np.clip(beta, eps, 1.0 - eps)
    return np.log2(b / (1.0 - b))


def m_to_beta(m_val: np.ndarray) -> np.ndarray:
    two_m = np.power(2.0, np.clip(m_val, -20.0, 20.0))
    return two_m / (1.0 + two_m)


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


def moderated_ols_contrast(Y_feat_by_sample: np.ndarray, X_design: np.ndarray, col_idx: int = 1) -> dict[str, Any]:
    """Fit row-wise OLS Y = B X^T with Smyth (2004) empirical-Bayes variance moderation
    for column `col_idx` of `X_design` (handling row-wise missing values cleanly).
    """
    n_feat, n_samp = Y_feat_by_sample.shape
    p_cols = X_design.shape[1]
    df_full = n_samp - p_cols
    if df_full < 2:
        fail(f"insufficient residual degrees of freedom (n={n_samp}, p={p_cols})")

    coef = np.full(n_feat, np.nan, dtype=float)
    s2 = np.full(n_feat, np.nan, dtype=float)
    v_unscaled = np.full(n_feat, np.nan, dtype=float)
    df_vec = np.zeros(n_feat, dtype=float)

    no_na_rows = ~np.isnan(Y_feat_by_sample).any(axis=1)
    if np.any(no_na_rows):
        xtx_inv = np.linalg.pinv(X_design.T @ X_design)
        proj = xtx_inv @ X_design.T  # p x n
        Y_sub = Y_feat_by_sample[no_na_rows]  # m x n
        B_sub = Y_sub @ proj.T  # m x p
        resid = Y_sub - B_sub @ X_design.T
        rss = np.sum(resid * resid, axis=1)
        coef[no_na_rows] = B_sub[:, col_idx]
        s2[no_na_rows] = rss / float(df_full)
        v_unscaled[no_na_rows] = float(xtx_inv[col_idx, col_idx])
        df_vec[no_na_rows] = float(df_full)

    na_indices = np.flatnonzero(~no_na_rows)
    for idx in na_indices:
        y_i = Y_feat_by_sample[idx]
        ok = np.isfinite(y_i)
        if int(ok.sum()) <= p_cols + 1:
            continue
        Xi = X_design[ok]
        if np.linalg.matrix_rank(Xi) < p_cols:
            continue
        yi = y_i[ok]
        xtx_i = np.linalg.pinv(Xi.T @ Xi)
        bi = xtx_i @ (Xi.T @ yi)
        ri = yi - Xi @ bi
        df_i = int(ok.sum()) - p_cols
        coef[idx] = float(bi[col_idx])
        s2[idx] = float(np.dot(ri, ri)) / float(df_i)
        v_unscaled[idx] = float(xtx_i[col_idx, col_idx])
        df_vec[idx] = float(df_i)

    valid = np.isfinite(coef) & np.isfinite(s2) & (s2 > 1e-12) & (v_unscaled > 1e-15) & (df_vec >= 2)
    d0 = 0.0
    s0_sq = float(np.nanmedian(s2[valid])) if np.any(valid) else 1.0
    if int(valid.sum()) >= 10:
        z = np.log(s2[valid])
        df_v = df_vec[valid]
        e_z = z - digamma(df_v / 2.0) + np.log(df_v / 2.0)
        mean_e = float(np.mean(e_z))
        var_e = float(np.var(e_z, ddof=1)) - float(np.mean(polygamma(1, df_v / 2.0)))
        if var_e > 1e-6:
            y = max(0.5, 1.0 / var_e)
            for _ in range(30):
                fv = float(polygamma(1, y)) - var_e
                fd = float(polygamma(2, y))
                if abs(fd) < 1e-12:
                    break
                step = fv / fd
                y = max(1e-3, y - step)
                if abs(step) < 1e-7:
                    break
            d0 = min(float(2.0 * y), 1e4)
            s0_sq = float(np.exp(mean_e + digamma(d0 / 2.0) - np.log(d0 / 2.0)))
        else:
            d0 = 50.0
            s0_sq = float(np.exp(mean_e))

    s2_post = np.where(valid, (d0 * s0_sq + df_vec * np.nan_to_num(s2, nan=s0_sq)) / np.maximum(1.0, df_vec + d0), np.nan)
    se = np.sqrt(s2_post * v_unscaled)
    t_mod = np.where(valid & (se > 0), coef / se, np.nan)
    pval = np.where(np.isfinite(t_mod), 2.0 * stats.t.sf(np.abs(t_mod), df=df_vec + d0), np.nan)
    fdr = bh_fdr(pval)

    return {
        "logFC_M": coef,
        "t": t_mod,
        "pvalue": pval,
        "fdr": fdr,
        "d0": round(float(d0), 3),
    }


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
        kind = "figure" if ext in (".png", ".svg", ".pdf") else "table"
    media_map = {
        ".png": "image/png",
        ".svg": "image/svg+xml",
        ".pdf": "application/pdf",
        ".csv": "text/csv",
        ".tsv": "text/tab-separated-values",
    }
    desc: dict[str, Any] = {
        "contractVersion": "1.0.0",
        "kind": kind,
        "file": p.name,
        "mediaType": media_map.get(ext, "text/tab-separated-values"),
        "title": title[:200],
        "provenance": {
            "createdAt": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "tool": tool,
            "skill": "methylation-workflow",
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
