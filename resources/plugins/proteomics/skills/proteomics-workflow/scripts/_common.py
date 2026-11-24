"""Shared utilities for quantitative MS proteomics tools (prot_*)."""

from __future__ import annotations

import json
import math
import os
import re
import sys
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats


CONTAM_PREFIXES = ("CON__", "REV__", "DECOY_", "CONTAM_")
META_FLAG_COLS = {
    "potential contaminant": "+",
    "reverse": "+",
    "only identified by site": "+",
}
KNOWN_ID_COLS = [
    "protein_id",
    "protein",
    "protein.group",
    "majority protein ids",
    "protein ids",
    "genes",
    "gene",
    "site_id",
    "ptm_id",
    "id",
]


def emit(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))


def fail(message: str, code: int = 1) -> None:
    emit({"error": message})
    sys.exit(code)


def info(message: str) -> None:
    print(f"[proteomics] {message}", file=sys.stderr)


def read_table(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        fail(f"file not found: {path}")
    sep = "," if path.lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, low_memory=False)


def detect_format_and_columns(
    df: pd.DataFrame, id_col: str | None = None
) -> dict[str, Any]:
    """Identify export format (MaxQuant / DIA-NN / FragPipe / generic), ID col, and sample cols."""
    cols = list(df.columns)
    lower_map = {c.lower(): c for c in cols}

    # Detect format
    lfq_cols = [c for c in cols if c.lower().startswith("lfq intensity ")]
    maxlfq_cols = [c for c in cols if c.lower().endswith(" maxlfq intensity")]
    if lfq_cols or "majority protein ids" in lower_map:
        fmt = "maxquant"
    elif "protein.group" in lower_map or "protein.ids" in lower_map:
        fmt = "diann"
    elif maxlfq_cols or "combined total spectral count" in lower_map:
        fmt = "fragpipe"
    else:
        fmt = "generic"

    # Choose ID column
    chosen_id = id_col
    if chosen_id is not None:
        if chosen_id not in df.columns:
            fail(f"id_col '{chosen_id}' not found in columns: {', '.join(cols[:15])}")
    else:
        for cand in KNOWN_ID_COLS:
            if cand in lower_map:
                chosen_id = lower_map[cand]
                break
        if chosen_id is None:
            chosen_id = cols[0]

    # Choose quantitative sample columns
    non_sample_lower = {
        chosen_id.lower(),
        "protein",
        "protein.group",
        "protein.ids",
        "protein.names",
        "genes",
        "gene",
        "first.protein.description",
        "majority protein ids",
        "protein ids",
        "gene names",
        "fasta headers",
        "potential contaminant",
        "reverse",
        "only identified by site",
        "peptides",
        "razor + unique peptides",
        "unique peptides",
        "sequence coverage [%]",
        "mol. weight [kda]",
        "q-value",
        "score",
        "protein_id",
        "site_id",
        "ptm_id",
        "residue",
        "position",
        "motif",
    }
    if fmt == "maxquant" and lfq_cols:
        sample_cols = lfq_cols
    elif fmt == "fragpipe" and maxlfq_cols:
        sample_cols = maxlfq_cols
    else:
        sample_cols = []
        for c in cols:
            if c == chosen_id or c.lower() in non_sample_lower:
                continue
            if pd.api.types.is_numeric_dtype(df[c]):
                sample_cols.append(c)
            else:
                coerced = pd.to_numeric(df[c], errors="coerce")
                if coerced.notna().mean() >= 0.5 and coerced.notna().sum() > 0:
                    sample_cols.append(c)

    return {
        "format": fmt,
        "id_col": chosen_id,
        "sample_cols": sample_cols,
    }


def contaminant_mask(df: pd.DataFrame, id_col: str) -> np.ndarray:
    """Boolean mask of contaminant / reverse / decoy / only-identified-by-site rows."""
    mask = np.zeros(len(df), dtype=bool)
    ids = df[id_col].astype(str).str.strip()
    for prefix in CONTAM_PREFIXES:
        mask |= ids.str.upper().str.startswith(prefix).to_numpy()
    for col in df.columns:
        cl = col.lower().strip()
        if cl in META_FLAG_COLS:
            flag_val = META_FLAG_COLS[cl]
            mask |= (df[col].astype(str).str.strip() == flag_val).to_numpy()
    return mask


def infer_scale(mat: pd.DataFrame) -> dict[str, Any]:
    """Infer whether quantitative values are raw linear intensities or already log2-transformed."""
    vals = mat.to_numpy(dtype=float)
    finite = vals[np.isfinite(vals)]
    if len(finite) == 0:
        return {"scale": "empty", "min": None, "median": None, "p95": None, "max": None}
    q50 = float(np.median(finite))
    q95 = float(np.percentile(finite, 95))
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    has_neg = vmin < 0
    if not has_neg and (q95 > 50.0 or vmax > 100.0):
        scale = "raw_intensity"
    else:
        scale = "log2_intensity"
    return {
        "scale": scale,
        "min": round(vmin, 4),
        "median": round(q50, 4),
        "p95": round(q95, 4),
        "max": round(vmax, 4),
    }


def load_metadata(
    path: str, sample_col: str | None = None
) -> tuple[pd.DataFrame, str]:
    meta = read_table(path)
    scol = sample_col or meta.columns[0]
    if scol not in meta.columns:
        fail(f"sample_col '{scol}' not found in metadata columns: {', '.join(meta.columns)}")
    meta[scol] = meta[scol].astype(str).str.strip()
    return meta, scol


def check_collinearity(meta: pd.DataFrame, group_col: str, batch_col: str | None) -> dict[str, Any]:
    if not batch_col or batch_col not in meta.columns or group_col not in meta.columns:
        return {"checked": False, "confounded": False}
    ct = pd.crosstab(meta[group_col].astype(str), meta[batch_col].astype(str))
    # Perfectly confounded if every batch contains at most 1 condition or every condition is in 1 batch.
    batches_pure = bool((ct.gt(0).sum(axis=0) <= 1).all())
    groups_pure = bool((ct.gt(0).sum(axis=1) <= 1).all())
    confounded = batches_pure or groups_pure
    return {
        "checked": True,
        "batch_col": batch_col,
        "group_col": group_col,
        "confounded": confounded,
        "cross_tab": {
            str(g): {str(b): int(ct.loc[g, b]) for b in ct.columns}
            for g in ct.index
        },
    }


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment handling NaNs cleanly."""
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


def moderated_t_test(
    x_case: np.ndarray, x_ctrl: np.ndarray
) -> dict[str, Any]:
    """Empirical-Bayes moderated two-sample t-test (Smyth 2004 limma-style shrinkage).

    Computes per-feature pooled residual variance s2 and df = n1 + n2 - 2,
    estimates prior degrees of freedom d0 and prior variance s0^2 from the
    distribution of log(s2) via method of moments (trigamma / digamma), and
    returns moderated t-statistics, p-values, and BH FDR.
    """
    from scipy.special import digamma, polygamma

    n_feat = x_case.shape[0]
    n1 = np.sum(np.isfinite(x_case), axis=1)
    n2 = np.sum(np.isfinite(x_ctrl), axis=1)
    mean1 = np.nanmean(x_case, axis=1)
    mean2 = np.nanmean(x_ctrl, axis=1)
    log2fc = mean1 - mean2

    ss1 = np.nansum((x_case - mean1[:, None]) ** 2, axis=1)
    ss2 = np.nansum((x_ctrl - mean2[:, None]) ** 2, axis=1)
    df = (n1 + n2 - 2).astype(float)
    valid = (n1 >= 2) & (n2 >= 2) & (df >= 2)
    s2 = np.full(n_feat, np.nan, dtype=float)
    s2[valid] = (ss1[valid] + ss2[valid]) / df[valid]
    pos = valid & np.isfinite(s2) & (s2 > 1e-12)

    d0 = 0.0
    s0_sq = float(np.nanmedian(s2[pos])) if pos.sum() > 0 else 1.0
    shrunk = False

    if int(pos.sum()) >= 10:
        z = np.log(s2[pos])
        df_pos = df[pos]
        e_z = z - digamma(df_pos / 2.0) + np.log(df_pos / 2.0)
        mean_e = float(np.mean(e_z))
        # Empirical variance of e_z minus expected trigamma(df/2)
        var_e = float(np.var(e_z, ddof=1)) - float(np.mean(polygamma(1, df_pos / 2.0)))
        if var_e > 1e-6:
            # Invert trigamma(d0/2) = var_e via Newton's method
            y = max(0.5, 1.0 / var_e)
            for _ in range(30):
                f_val = float(polygamma(1, y)) - var_e
                f_der = float(polygamma(2, y))
                if abs(f_der) < 1e-12:
                    break
                step = f_val / f_der
                y = max(1e-3, y - step)
                if abs(step) < 1e-7:
                    break
            d0 = min(float(2.0 * y), 1e4)
            s0_sq = float(np.exp(mean_e + digamma(d0 / 2.0) - np.log(d0 / 2.0)))
            shrunk = True
        else:
            # Extremely homogeneous variances -> strong moderation toward s0_sq
            d0 = 100.0
            s0_sq = float(np.exp(mean_e))
            shrunk = True

    s2_post = np.full(n_feat, np.nan, dtype=float)
    df_total = df + d0
    if shrunk:
        s2_post[valid] = (d0 * s0_sq + df[valid] * np.nan_to_num(s2[valid], nan=s0_sq)) / df_total[valid]
    else:
        s2_post[valid] = s2[valid]

    se = np.sqrt(s2_post * (1.0 / np.maximum(n1, 1) + 1.0 / np.maximum(n2, 1)))
    t_stat = np.full(n_feat, np.nan, dtype=float)
    ok = valid & np.isfinite(se) & (se > 0)
    t_stat[ok] = log2fc[ok] / se[ok]
    pvals = np.full(n_feat, np.nan, dtype=float)
    pvals[ok] = 2.0 * stats.t.sf(np.abs(t_stat[ok]), df=df_total[ok])
    padj = bh_fdr(pvals)

    return {
        "log2fc": log2fc,
        "mean_case": mean1,
        "mean_control": mean2,
        "ave_expr": 0.5 * (mean1 + mean2),
        "stat": t_stat,
        "pvalue": pvals,
        "padj": padj,
        "d0": round(float(d0), 4),
        "s0": round(float(math.sqrt(max(s0_sq, 0.0))), 4),
        "shrunk": shrunk,
    }


def read_gmt(path: str) -> dict[str, list[str]]:
    if not os.path.exists(path):
        fail(f"GMT file not found: {path}")
    sets: dict[str, list[str]] = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = [p.strip() for p in line.rstrip("\n").split("\t") if p.strip()]
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = list(dict.fromkeys(parts[2:]))
            sets[name] = genes
    return sets


def normalize_gene_id(val: str) -> str:
    """Extract primary gene symbol / accession from semicolon-delimited protein group IDs."""
    text = str(val).strip()
    if ";" in text:
        text = text.split(";")[0].strip()
    return text


def clean_sample_name(col: str) -> str:
    """Strip common MaxQuant/FragPipe intensity prefixes/suffixes for metadata matching."""
    c = str(col).strip()
    c = re.sub(r"^LFQ intensity\s+", "", c, flags=re.IGNORECASE)
    c = re.sub(r"\s+MaxLFQ Intensity$", "", c, flags=re.IGNORECASE)
    c = re.sub(r"\s+Intensity$", "", c, flags=re.IGNORECASE)
    return c.strip()
