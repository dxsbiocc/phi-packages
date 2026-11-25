#!/usr/bin/env python3
"""Shared utilities for the mendelian-randomization plugin (toolPrefix: mr)."""

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
from scipy import stats

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

SNP_CANDS = ("snp", "rsid", "variant_id", "id", "markername", "hm_rsid")
CHR_CANDS = ("chr", "chrom", "chromosome", "#chrom", "hm_chrom")
POS_CANDS = ("pos", "bp", "position", "base_pair_location", "hm_pos", "start")
EA_CANDS = ("ea", "effect_allele", "a1", "alt", "allele1", "hm_effect_allele")
OA_CANDS = ("oa", "other_allele", "non_effect_allele", "a2", "ref", "allele2", "hm_other_allele")
BETA_CANDS = ("beta", "b", "effect", "log_odds", "hm_beta")
OR_CANDS = ("or", "odds_ratio", "hm_odds_ratio")
SE_CANDS = ("se", "standard_error", "stderr")
P_CANDS = ("p", "pval", "p_value", "pvalue", "p-value")
EAF_CANDS = ("eaf", "effect_allele_frequency", "af", "frq", "freq", "a1freq", "maf", "hm_effect_allele_frequency")
N_CANDS = ("n", "samplesize", "sample_size", "n_total", "neff")


def emit(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)


def fail(message: str, code: int = 1) -> None:
    print(json.dumps({"error": message}, indent=2, ensure_ascii=False))
    sys.exit(code)


def is_palindromic(a1: str, a2: str) -> bool:
    pair = {a1.upper(), a2.upper()}
    return pair in ({"A", "T"}, {"C", "G"})


def strand_flip(allele: str) -> str:
    return "".join(COMPLEMENT.get(b, b) for b in allele.upper())


def read_summary_stats(path: str | Path, default_n: int | None = None) -> tuple[pd.DataFrame, str]:
    """Read and standardize a GWAS/QTL summary-statistics table (.tsv/.csv/.gz)."""
    p = Path(path)
    if not p.is_file():
        fail(f"summary statistics file not found: {p}")
    low_name = p.name.lower()
    sep = "," if (low_name.endswith(".csv") or low_name.endswith(".csv.gz")) else "\t"
    try:
        df = pd.read_csv(p, sep=sep)
        if df.shape[1] == 1:
            df = pd.read_csv(p, sep=r"\s+")
    except Exception as exc:
        fail(f"failed to parse summary statistics {p}: {exc}")
    if df.empty:
        fail(f"summary statistics table is empty: {p}")

    # If already harmonized/standardized with exact standard columns, keep fast path
    lower = {str(c).strip().lower(): c for c in df.columns}
    pick = lambda cands: next((lower[c] for c in cands if c in lower), None)

    snp_col = pick(SNP_CANDS)
    chr_col = pick(CHR_CANDS)
    pos_col = pick(POS_CANDS)
    ea_col = pick(EA_CANDS)
    oa_col = pick(OA_CANDS)
    beta_col = pick(BETA_CANDS)
    or_col = pick(OR_CANDS)
    se_col = pick(SE_CANDS)
    p_col = pick(P_CANDS)
    eaf_col = pick(EAF_CANDS)
    n_col = pick(N_CANDS)

    if snp_col is None and (chr_col is None or pos_col is None):
        fail(f"summary table {p} must have a SNP/rsid column or CHR+POS columns; found {list(df.columns)}")
    if ea_col is None or oa_col is None:
        fail(f"summary table {p} must have effect_allele (EA/A1) and other_allele (OA/A2) columns; found {list(df.columns)}")
    if beta_col is None and or_col is None:
        fail(f"summary table {p} must have BETA or OR column; found {list(df.columns)}")
    if se_col is None:
        fail(f"summary table {p} must have SE (standard_error) column; found {list(df.columns)}")

    out = pd.DataFrame()
    if snp_col is not None:
        out["SNP"] = df[snp_col].astype(str).str.strip()
    else:
        out["SNP"] = df[chr_col].astype(str).str.strip() + ":" + df[pos_col].astype(str).str.strip()

    out["CHR"] = (
        df[chr_col].astype(str).str.strip().str.replace("^chr", "", case=False, regex=True)
        if chr_col is not None
        else pd.Series([None] * len(df), dtype=object)
    )
    out["POS"] = pd.to_numeric(df[pos_col], errors="coerce") if pos_col is not None else np.nan
    out["EA"] = df[ea_col].astype(str).str.strip().str.upper()
    out["OA"] = df[oa_col].astype(str).str.strip().str.upper()

    if beta_col is not None:
        out["BETA"] = pd.to_numeric(df[beta_col], errors="coerce")
        effect_scale = "beta"
    else:
        or_vals = pd.to_numeric(df[or_col], errors="coerce")
        out["BETA"] = np.where(or_vals > 0, np.log(or_vals), np.nan)
        effect_scale = "or_converted_to_log_odds"

    out["SE"] = pd.to_numeric(df[se_col], errors="coerce")
    if p_col is not None:
        p_raw = pd.to_numeric(df[p_col], errors="coerce")
    else:
        p_raw = pd.Series(np.nan, index=out.index)
    # Recompute/fill missing or underflowed P from Z = BETA / SE where needed
    z_val = out["BETA"] / np.where(out["SE"] > 0, out["SE"], np.nan)
    p_from_z = 2.0 * stats.norm.sf(np.abs(z_val))
    out["P"] = np.where(p_raw.notna() & (p_raw > 0) & (p_raw <= 1.0), p_raw, p_from_z)

    out["EAF"] = pd.to_numeric(df[eaf_col], errors="coerce") if eaf_col is not None else np.nan
    if n_col is not None:
        out["N"] = pd.to_numeric(df[n_col], errors="coerce")
    elif default_n is not None:
        out["N"] = float(default_n)
    else:
        out["N"] = np.nan

    valid = (
        out["SNP"].notna()
        & (out["SNP"] != "")
        & out["EA"].isin({"A", "C", "G", "T"})
        & out["OA"].isin({"A", "C", "G", "T"})
        & (out["EA"] != out["OA"])
        & np.isfinite(out["BETA"])
        & np.isfinite(out["SE"])
        & (out["SE"] > 0)
        & np.isfinite(out["P"])
    )
    clean = out.loc[valid].sort_values("P", ascending=True).drop_duplicates(subset=["SNP"], keep="first").reset_index(drop=True)
    return clean, effect_scale


def clump_by_distance(df: pd.DataFrame, p_cutoff: float, clump_kb: int) -> pd.DataFrame:
    """Greedy minimum-P physical window clumping (clump_kb * 1000 bp)."""
    cand = df[df["P"] <= p_cutoff].sort_values("P", ascending=True).reset_index(drop=True)
    if cand.empty or clump_kb <= 0 or cand["CHR"].isna().all() or cand["POS"].isna().all():
        return cand
    win_bp = int(clump_kb) * 1000
    kept_rows = []
    taken_by_chr: dict[str, list[int]] = {}
    for _, row in cand.iterrows():
        ch = str(row["CHR"])
        pos_v = row["POS"]
        if pd.isna(pos_v):
            kept_rows.append(row)
            continue
        pos = int(pos_v)
        existing = taken_by_chr.get(ch, [])
        if any(abs(pos - p0) <= win_bp for p0 in existing):
            continue
        taken_by_chr.setdefault(ch, []).append(pos)
        kept_rows.append(row)
    return pd.DataFrame(kept_rows).reset_index(drop=True)


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
            "skill": "mr-workflow",
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
