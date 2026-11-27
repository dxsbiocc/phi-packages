#!/usr/bin/env python3
"""Shared I/O, allele harmonization, statistics, and Contract-1.0.0 artifact
utilities for the statistical-genetics plugin (toolPrefix: gwas).
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import re
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
from scipy import stats

IUPAC_DNA = {"A", "C", "G", "T"}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

SNP_CANDS = ("snp", "rsid", "variant_id", "id", "markername", "snpid", "markerid", "hm_rsid")
CHR_CANDS = ("chr", "chrom", "chromosome", "#chrom", "#chr", "hm_chrom")
POS_CANDS = ("pos", "bp", "genpos", "position", "base_pair_location", "hm_pos", "start")
EA_CANDS = ("ea", "effect_allele", "hm_effect_allele")
NEA_CANDS = ("nea", "other_allele", "non_effect_allele", "oa", "hm_other_allele")
BETA_CANDS = ("beta", "b", "effect", "log_odds", "hm_beta")
OR_CANDS = ("or", "odds_ratio", "hm_odds_ratio")
SE_CANDS = ("se", "standard_error", "stderr", "log(or)_se", "se_beta")
Z_CANDS = ("z", "zscore", "z_score", "z_stat", "t_stat", "tstat")
P_CANDS = ("p", "pval", "p_value", "pvalue", "p-value", "p.value", "p_bolt_lmm", "p_bolt_lmm_inf")
LOG10P_CANDS = ("log10p", "log10_p", "mlog10p", "neg_log10_p")
EAF_CANDS = ("eaf", "effect_allele_frequency", "a1freq", "a1_freq", "af_allele2", "af", "frq", "freq", "maf", "hm_effect_allele_frequency")
INFO_CANDS = ("info", "imputationinfo", "rsq", "r2", "impute_info", "mach_r2")
N_CANDS = ("n", "obs_ct", "samplesize", "sample_size", "n_total")
NEFF_CANDS = ("n_eff", "neff", "effective_n")


def emit(payload: dict[str, Any]) -> None:
    """Write strict JSON output payload to stdout and exit 0."""
    sys.stdout.write(json.dumps(clean_json(payload), ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()
    raise SystemExit(0)


def fail(message: str, code: int = 1) -> None:
    """Write JSON error payload to stdout and exit non-zero."""
    sys.stdout.write(json.dumps({"error": str(message)}, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()
    raise SystemExit(code)


def clean_json(obj: Any) -> Any:
    """Recursively convert numpy types and non-finite floats into JSON-safe types."""
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, np.bool_):
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


def setup_matplotlib():
    """Configure headless Agg backend and publication-grade matplotlib defaults."""
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


def read_table(path: str | Path) -> pd.DataFrame:
    """Read TSV/CSV/whitespace-delimited text or gzipped table with #-header support."""
    p = Path(path)
    if not p.is_file():
        fail(f"Input file not found: {path}")

    low = p.name.lower()
    if low.endswith(".gz"):
        import gzip

        with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
            raw_lines = fh.readlines()
    else:
        with open(p, "rt", encoding="utf-8", errors="replace") as fh:
            raw_lines = fh.readlines()

    filtered_lines: list[str] = []
    header_seen = False
    for line in raw_lines:
        stripped = line.strip()
        if not stripped:
            continue
        if not header_seen:
            if stripped.startswith("##"):
                continue
            if stripped.startswith("#"):
                after = stripped[1:]
                if "\t" in after or "," in after or re.search(r"\s", after):
                    filtered_lines.append(line.lstrip()[1:])
                    header_seen = True
                continue
            filtered_lines.append(line)
            header_seen = True
        else:
            if stripped.startswith("#"):
                continue
            filtered_lines.append(line)

    if not filtered_lines:
        fail(f"Input table is empty: {path}")

    first = filtered_lines[0]
    if low.endswith(".csv") or low.endswith(".csv.gz"):
        sep = ","
    elif "\t" in first:
        sep = "\t"
    elif "," in first and not re.search(r"\s", first.replace(",", "")):
        sep = ","
    elif "," in first and len(first.split(",")) > 1:
        sep = ","
    else:
        sep = r"\s+"

    try:
        df = pd.read_csv(io.StringIO("".join(filtered_lines)), sep=sep, low_memory=False)
        if df.shape[1] == 1 and sep != r"\s+":
            df = pd.read_csv(io.StringIO("".join(filtered_lines)), sep=r"\s+", low_memory=False)
    except Exception as exc:
        fail(f"Failed to parse table '{path}': {exc}")

    if df.empty:
        fail(f"Input table is empty: {path}")
    df.columns = [str(c).strip() for c in df.columns]
    return df


def resolve_column(df: pd.DataFrame, candidates: tuple[str, ...] | list[str], label: str, required: bool = True) -> str | None:
    """Resolve column name case-insensitively against candidate names."""
    lower_map = {str(c).strip().lower(): str(c).strip() for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        hit = lower_map.get(cand.lower())
        if hit is not None:
            return hit
    if required:
        fail(f"Missing required {label} column (expected one of {list(candidates)}; found {list(df.columns)})")
    return None


def is_valid_iupac_allele(allele: str) -> bool:
    """Return True if allele is a non-empty A/C/G/T string."""
    s = str(allele).strip().upper()
    return bool(s) and all(ch in IUPAC_DNA for ch in s)


def strand_flip(allele: str) -> str:
    """Return base-complemented allele on the opposite strand."""
    return "".join(COMPLEMENT.get(ch, ch) for ch in str(allele).strip().upper())


def is_palindromic(a1: str, a2: str) -> bool:
    """Return True if alleles are self-complementary (e.g. A/T or C/G)."""
    u1 = str(a1).strip().upper()
    u2 = str(a2).strip().upper()
    return bool(u1 and u2 and u1 == strand_flip(u2))


def normalize_chrom(val: Any) -> str:
    """Normalize chromosome string (strip 'chr', uppercase X/Y/MT)."""
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return ""
    s = str(val).strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    if s.endswith(".0") and s[:-2].isdigit():
        s = s[:-2]
    up = s.upper()
    if up in {"23"}:
        return "X"
    if up in {"24"}:
        return "Y"
    if up in {"25", "26", "M", "MT"}:
        return "MT"
    return up if up in {"X", "Y", "MT"} else s


def chrom_sort_key(ch: Any) -> tuple[int, str]:
    """Natural chromosome sort key (1..22, X, Y, MT)."""
    s = normalize_chrom(ch).upper()
    if s.isdigit():
        return (int(s), "")
    order = {"X": 23, "Y": 24, "MT": 25}
    return (order.get(s, 99), s)


def _pick_alleles(df: pd.DataFrame) -> tuple[pd.Series, pd.Series, str | None]:
    """Resolve (EA, NEA, inferred_eaf_col) across PLINK2, REGENIE, SAIGE, BOLT-LMM, GWAS Catalog."""
    lower = {str(c).strip().lower(): str(c).strip() for c in df.columns}

    ea_col = next((lower[c] for c in EA_CANDS if c in lower), None)
    nea_col = next((lower[c] for c in NEA_CANDS if c in lower), None)
    if ea_col is not None and nea_col is not None:
        return (
            df[ea_col].astype(str).str.strip().str.upper(),
            df[nea_col].astype(str).str.strip().str.upper(),
            None,
        )

    # SAIGE convention: Allele1=REF/non-effect, Allele2=ALT/effect, AF_Allele2
    if "allele1" in lower and "allele2" in lower and ("af_allele2" in lower or "p.value" in lower or "tstat" in lower):
        eaf_hint = lower.get("af_allele2")
        return (
            df[lower["allele2"]].astype(str).str.strip().str.upper(),
            df[lower["allele1"]].astype(str).str.strip().str.upper(),
            eaf_hint,
        )

    # REGENIE / BOLT-LMM convention: ALLELE1=effect, ALLELE0=other
    if "allele1" in lower and "allele0" in lower:
        return (
            df[lower["allele1"]].astype(str).str.strip().str.upper(),
            df[lower["allele0"]].astype(str).str.strip().str.upper(),
            None,
        )

    # PLINK / LDSC convention: A1=effect, A2 or AX/OMITTED or non-A1 of REF/ALT
    if "a1" in lower:
        a1 = df[lower["a1"]].astype(str).str.strip().str.upper()
        for alt_nea in ("a2", "ax", "omitted", "allele2", "allele0"):
            if alt_nea in lower:
                return a1, df[lower[alt_nea]].astype(str).str.strip().str.upper(), None
        if "ref" in lower and "alt" in lower:
            ref = df[lower["ref"]].astype(str).str.strip().str.upper()
            alt = df[lower["alt"]].astype(str).str.strip().str.upper()
            nea = pd.Series(np.where(a1 == alt, ref, alt), index=df.index)
            return a1, nea, None

    # Generic ALT/REF (VCF/GWAS): ALT=effect, REF=non-effect
    if "alt" in lower and "ref" in lower:
        return (
            df[lower["alt"]].astype(str).str.strip().str.upper(),
            df[lower["ref"]].astype(str).str.strip().str.upper(),
            None,
        )

    # Fallback: ALLELE1 + ALLELE2 without SAIGE markers (A1=effect, A2=other)
    if "allele1" in lower and "allele2" in lower:
        return (
            df[lower["allele1"]].astype(str).str.strip().str.upper(),
            df[lower["allele2"]].astype(str).str.strip().str.upper(),
            None,
        )

    return pd.Series("", index=df.index), pd.Series("", index=df.index), None


def parse_sumstats(
    path: str | Path,
    *,
    require_alleles: bool = True,
    require_beta_se: bool = True,
    require_n: bool = False,
    default_n: float | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Read and standardize a GWAS summary-statistics file across caller formats."""
    df = read_table(path)
    n_input = int(len(df))
    lower = {str(c).strip().lower(): str(c).strip() for c in df.columns}
    pick = lambda cands: next((lower[c] for c in cands if c in lower), None)

    snp_col = pick(SNP_CANDS)
    chr_col = pick(CHR_CANDS)
    pos_col = pick(POS_CANDS)
    beta_col = pick(BETA_CANDS)
    or_col = pick(OR_CANDS)
    se_col = pick(SE_CANDS)
    z_col = pick(Z_CANDS)
    p_col = pick(P_CANDS)
    log10p_col = pick(LOG10P_CANDS)
    eaf_col = pick(EAF_CANDS)
    info_col = pick(INFO_CANDS)
    n_col = pick(N_CANDS)
    neff_col = pick(NEFF_CANDS)

    if snp_col is None and (chr_col is None or pos_col is None):
        fail(f"Summary table '{path}' requires SNP/rsid or CHR+POS columns; found {list(df.columns)}")

    ea_s, nea_s, eaf_override = _pick_alleles(df)
    if eaf_override is not None:
        eaf_col = eaf_override
    has_alleles = bool((ea_s != "").any() and (nea_s != "").any())
    if require_alleles and not has_alleles:
        fail(f"Summary table '{path}' requires effect/non-effect allele columns (e.g. EA/NEA, A1/A2, ALLELE1/ALLELE0, REF/ALT)")

    out = pd.DataFrame(index=df.index)
    if chr_col is not None:
        out["CHR"] = df[chr_col].map(normalize_chrom)
    else:
        out["CHR"] = ""
    out["POS"] = pd.to_numeric(df[pos_col], errors="coerce") if pos_col is not None else np.nan

    if snp_col is not None:
        snp_vals = df[snp_col].astype(str).str.strip()
        if chr_col is not None and pos_col is not None:
            missing_snp = snp_vals.isin({"", ".", "NA", "nan", "None"})
            fallback_id = out["CHR"].astype(str) + ":" + out["POS"].fillna(-1).astype(int).astype(str)
            snp_vals = snp_vals.where(~missing_snp, fallback_id)
        out["SNP"] = snp_vals
    else:
        out["SNP"] = out["CHR"].astype(str) + ":" + out["POS"].fillna(-1).astype(int).astype(str)

    out["EA"] = ea_s
    out["NEA"] = nea_s

    # Effect size & SE / Z
    effect_scale = "z_only"
    if beta_col is not None:
        out["BETA"] = pd.to_numeric(df[beta_col], errors="coerce")
        effect_scale = "beta"
    elif or_col is not None:
        or_vals = pd.to_numeric(df[or_col], errors="coerce")
        out["BETA"] = np.where(or_vals > 0, np.log(or_vals), np.nan)
        effect_scale = "or_converted_to_log_odds"
    else:
        out["BETA"] = np.nan

    if se_col is not None:
        out["SE"] = pd.to_numeric(df[se_col], errors="coerce")
    elif or_col is not None and "ci_lower" in lower and "ci_upper" in lower:
        lo = pd.to_numeric(df[lower["ci_lower"]], errors="coerce")
        hi = pd.to_numeric(df[lower["ci_upper"]], errors="coerce")
        out["SE"] = np.where((lo > 0) & (hi > lo), (np.log(hi) - np.log(lo)) / (2.0 * 1.959963984540054), np.nan)
    elif beta_col is not None and "ci_lower" in lower and "ci_upper" in lower:
        lo = pd.to_numeric(df[lower["ci_lower"]], errors="coerce")
        hi = pd.to_numeric(df[lower["ci_upper"]], errors="coerce")
        out["SE"] = np.where(hi > lo, (hi - lo) / (2.0 * 1.959963984540054), np.nan)
    else:
        out["SE"] = np.nan

    if z_col is not None:
        out["Z"] = pd.to_numeric(df[z_col], errors="coerce")
    else:
        out["Z"] = np.nan

    # Derive BETA from Z*SE or Z from BETA/SE when one is missing
    if out["BETA"].isna().all() and out["Z"].notna().any() and out["SE"].notna().any():
        out["BETA"] = out["Z"] * out["SE"]
        effect_scale = "z_and_se"
    if out["Z"].isna().all() and out["BETA"].notna().any() and out["SE"].notna().any():
        se_arr = out["SE"].to_numpy(dtype=float, copy=True)
        beta_arr = out["BETA"].to_numpy(dtype=float, copy=True)
        out["Z"] = np.divide(
            beta_arr,
            se_arr,
            out=np.full_like(beta_arr, np.nan, dtype=float),
            where=np.isfinite(se_arr) & (se_arr > 0),
        )

    if require_beta_se and (out["BETA"].isna().all() or out["SE"].isna().all()):
        fail(f"Summary table '{path}' requires BETA (or OR) and SE columns; found {list(df.columns)}")
    if not require_beta_se and out["Z"].isna().all():
        fail(f"Summary table '{path}' requires Z or BETA+SE columns; found {list(df.columns)}")

    # P-value: P or LOG10P or derived from Z
    if p_col is not None:
        p_raw = pd.to_numeric(df[p_col], errors="coerce")
    elif log10p_col is not None:
        lp = pd.to_numeric(df[log10p_col], errors="coerce").to_numpy(dtype=float, copy=True)
        finite_lp = lp[np.isfinite(lp)]
        if finite_lp.size > 0 and np.nanmedian(finite_lp) >= 0:
            p_raw = pd.Series(np.power(10.0, -np.clip(lp, 0.0, 320.0)), index=df.index)
        else:
            p_raw = pd.Series(np.power(10.0, np.clip(lp, -320.0, 0.0)), index=df.index)
        if effect_scale == "beta":
            effect_scale = "beta_with_log10p"
    else:
        p_raw = pd.Series(np.nan, index=df.index)

    z_arr = out["Z"].to_numpy(dtype=float, copy=True)
    p_from_z = np.where(np.isfinite(z_arr), 2.0 * stats.norm.sf(np.abs(z_arr)), np.nan)
    p_raw_arr = p_raw.to_numpy(dtype=float, copy=True)
    p_final = np.where(np.isfinite(p_raw_arr) & (p_raw_arr > 0.0) & (p_raw_arr <= 1.0), p_raw_arr, p_from_z)
    out["P"] = np.clip(p_final, 1e-300, 1.0)

    out["EAF"] = pd.to_numeric(df[eaf_col], errors="coerce") if eaf_col is not None else np.nan
    out["INFO"] = pd.to_numeric(df[info_col], errors="coerce") if info_col is not None else np.nan

    if n_col is not None:
        out["N"] = pd.to_numeric(df[n_col], errors="coerce")
    elif default_n is not None:
        out["N"] = float(default_n)
    else:
        out["N"] = np.nan

    if neff_col is not None:
        out["N_EFF"] = pd.to_numeric(df[neff_col], errors="coerce")
    else:
        out["N_EFF"] = out["N"]

    if require_n and out["N"].isna().all():
        fail(f"Summary table '{path}' requires sample size N column; found {list(df.columns)}")

    meta = {
        "n_input": n_input,
        "effect_scale_input": effect_scale,
        "has_alleles": has_alleles,
        "has_eaf": bool(out["EAF"].notna().any()),
        "has_info": bool(out["INFO"].notna().any()),
        "has_coords": bool((out["CHR"] != "").any() and out["POS"].notna().any()),
    }
    return out.reset_index(drop=True), meta


def clump_by_distance(df: pd.DataFrame, p_cutoff: float, window_kb: int) -> pd.DataFrame:
    """Greedy minimum-P physical-window clumping with minimum separation window_kb."""
    cand = df[df["P"] <= p_cutoff].sort_values("P", ascending=True).copy()
    if cand.empty or window_kb <= 0 or (cand["CHR"] == "").all() or cand["POS"].isna().all():
        return cand.reset_index(drop=True)

    win_bp = int(window_kb) * 1000
    kept_idx: list[Any] = []
    taken_by_chr: dict[str, list[int]] = {}
    for idx, row in cand.iterrows():
        ch = str(row["CHR"])
        pos_v = row["POS"]
        if not ch or pd.isna(pos_v):
            kept_idx.append(idx)
            continue
        pos = int(pos_v)
        existing = taken_by_chr.get(ch, [])
        if any(abs(pos - p0) < win_bp for p0 in existing):
            continue
        taken_by_chr.setdefault(ch, []).append(pos)
        kept_idx.append(idx)
    return cand.loc[kept_idx].reset_index(drop=True)


def write_tsv(df: pd.DataFrame, path: str | Path, tool: str, title: str, **kwargs: Any) -> str:
    """Write DataFrame to TSV and emit its Contract-1.0.0 .phi-artifact.json sidecar."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, sep="\t", index=False)
    return write_artifact(p, tool=tool, title=title, shape=df.shape, **kwargs)


def write_artifact(
    path: str | Path,
    tool: str,
    title: str,
    *,
    skill: str = "gwas-workflow",
    inputs: tuple[str | Path | None, ...] | list[str | Path | None] = (),
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
        if ext == ".png" and p.is_file():
            with open(p, "rb") as fh:
                hdr = fh.read(24)
            if hdr[:8] == b"\x89PNG\r\n\x1a\n":
                w, h = struct.unpack(">II", hdr[16:24])
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
