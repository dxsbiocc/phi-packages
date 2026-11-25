"""Shared helpers for the bulk-rna-seq plugin (rnaseq-workflow).

Every tool prints exactly one JSON object to stdout; diagnostics go to stderr.
Produced tables and figures carry contract-1.0.0 .phi-artifact.json sidecars.
"""
import hashlib
import json
import os
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path

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


def read_expression_matrix(path):
    df = read_table(path)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    num = df.apply(pd.to_numeric, errors="coerce")
    if num.isna().any().any():
        bad = num.columns[num.isna().any()].tolist()
        fail(f"non-numeric values in expression matrix columns: {bad}")
    return num


def infer_value_kind(df):
    """Classify an expression matrix as raw_counts, log_normalized, or tpm_or_cpm."""
    x = df.to_numpy(dtype=float)
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        return "log_normalized"
    if (finite < 0).any():
        return "log_normalized"
    max_v = float(np.max(finite))
    p99 = float(np.percentile(finite, 99))
    if max_v <= 35.0 and p99 <= 25.0:
        return "log_normalized"
    frac_int = float(np.mean(np.abs(finite - np.round(finite)) < 1e-6))
    if frac_int >= 0.98 and max_v >= 20.0:
        return "raw_counts"
    return "tpm_or_cpm"


def to_log2_cpm(df, value_kind=None):
    """Return a DataFrame on log2(CPM + 1) scale (or unchanged if already log_normalized)."""
    kind = value_kind or infer_value_kind(df)
    if kind == "log_normalized":
        return df.copy()
    x = np.clip(df.to_numpy(dtype=float), 0.0, None)
    col_sums = x.sum(axis=0, keepdims=True)
    safe_sums = np.where(col_sums > 0, col_sums, 1.0)
    cpm = x / safe_sums * 1e6
    return pd.DataFrame(np.log2(cpm + 1.0), index=df.index, columns=df.columns)


def read_metadata(path, samples=None, sample_col=None):
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    if sample_col:
        df = pd.read_csv(path, sep=sep)
        if sample_col not in df.columns:
            fail(f"metadata missing sample_col {sample_col!r}; columns: {df.columns.tolist()}")
        df = df.set_index(sample_col)
    else:
        df = pd.read_csv(path, sep=sep, index_col=0)
    df.index = df.index.astype(str).str.strip()
    if samples is not None:
        missing = [s for s in samples if s not in df.index]
        if missing:
            fail(f"metadata missing samples: {missing}")
        df = df.loc[samples]
    return df


def read_de_table(path):
    """Load a DESeq2, edgeR, or limma differential expression result table and
    normalize key column names (gene, log2fc, pvalue, padj, stat).
    """
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    raw = pd.read_csv(path, sep=sep)
    if raw.empty:
        fail(f"DE table is empty: {path}")

    cols_lower = {c.lower().strip(): c for c in raw.columns}
    gene_col = None
    for cand in ("gene", "gene_id", "geneid", "symbol", "gene_name", "id", "unnamed: 0"):
        if cand in cols_lower:
            gene_col = cols_lower[cand]
            break
    if gene_col is None:
        gene_col = raw.columns[0]

    def _pick(cands):
        for c in cands:
            if c in cols_lower:
                return cols_lower[c]
        return None

    fc_col = _pick(("log2foldchange", "log2fc", "logfc", "lfc"))
    p_col = _pick(("pvalue", "pval", "p.value", "p_val", "p"))
    padj_col = _pick(("padj", "fdr", "adj.p.val", "qvalue", "q_value", "p_adj"))
    stat_col = _pick(("stat", "t", "lr", "f", "wald"))

    if fc_col is None:
        fail(f"{path} missing log2 fold-change column (expected log2FoldChange / logFC / log2FC)")
    if padj_col is None and p_col is None:
        fail(f"{path} missing p-value/padj column")

    out = pd.DataFrame()
    out["gene"] = raw[gene_col].astype(str).str.strip()
    out["log2fc"] = pd.to_numeric(raw[fc_col], errors="coerce")
    out["pvalue"] = pd.to_numeric(raw[p_col], errors="coerce") if p_col else np.nan
    out["padj"] = pd.to_numeric(raw[padj_col], errors="coerce") if padj_col else np.nan
    if stat_col:
        out["stat"] = pd.to_numeric(raw[stat_col], errors="coerce")
    else:
        # Fallback signed ranking statistic: sign(log2fc) * -log10(pvalue)
        p_safe = np.clip(out["pvalue"].fillna(1.0).to_numpy(dtype=float), 1e-300, 1.0)
        out["stat"] = np.sign(out["log2fc"].fillna(0.0).to_numpy(dtype=float)) * (-np.log10(p_safe))

    out = out.loc[out["gene"] != ""].drop_duplicates(subset=["gene"], keep="first").reset_index(drop=True)
    return out


def read_gmt(path):
    """Parse a standard GMT file: <gene_set_name>\\t<description>\\t<gene1>\\t<gene2>..."""
    p = Path(path)
    if not p.is_file():
        fail(f"GMT file not found: {p}")
    gene_sets = {}
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        name = parts[0].strip()
        genes = [g.strip() for g in parts[2:] if g.strip()]
        if name and genes:
            gene_sets[name] = sorted(set(genes))
    if not gene_sets:
        fail(f"no valid gene sets parsed from GMT: {p}")
    return gene_sets


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out


def setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def write_artifact(
    path,
    tool,
    title,
    *,
    skill="rnaseq-workflow",
    inputs=(),
    parameters=None,
    shape=None,
    scientific_context=None,
):
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
