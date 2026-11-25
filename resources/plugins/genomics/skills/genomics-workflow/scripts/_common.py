"""Shared utilities for somatic cancer genomics tools (gen_*)."""

from __future__ import annotations

import hashlib
import json
import os
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


MANDATORY_MAF_COLS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "Variant_Type",
    "Tumor_Sample_Barcode",
]

NON_SYNONYMOUS_CLASSES = {
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "Splice_Site",
    "Translation_Start_Site",
}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

SIX_CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
TRANSITIONS = {"C>T", "T>C"}

# Curated core member genes of the 10 canonical TCGA oncogenic signaling pathways
# (Sanchez-Vega et al., Cell 2018, 173(2):321-337.e10).
SANCHEZ_VEGA_PATHWAYS: dict[str, list[str]] = {
    "Cell_Cycle": [
        "CDKN2A", "CDKN2B", "CDKN2C", "CDKN1A", "CDKN1B", "RB1", "CCND1", "CCND2",
        "CCND3", "CCNE1", "CDK4", "CDK6", "CDK2", "E2F1", "E2F3",
    ],
    "Hippo": [
        "FAT1", "FAT2", "FAT3", "FAT4", "NF2", "LATS1", "LATS2", "MOB1A", "MOB1B",
        "SAV1", "STK3", "STK4", "YAP1", "WWTR1", "TEAD1", "TEAD2", "CSNK1D", "CSNK1E",
        "DCHS1", "DCHS2", "PTPN14", "TAOK1", "TAOK2", "TAOK3", "WWC1",
    ],
    "Myc": [
        "MYC", "MYCN", "MYCL", "MAX", "MGA", "MXD1", "MXI1", "MXD3", "MXD4",
        "MLX", "MLXIP", "MLXIPL", "MNT",
    ],
    "Notch": [
        "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "FBXW7", "CREBBP", "EP300", "NCOR1",
        "NCOR2", "SPEN", "MAML1", "MAML2", "MAML3", "RBPJ", "JAG1", "JAG2", "DLL1",
        "DLL3", "DLL4", "NUMB", "ARRDC1", "CNTN6", "CUL1", "HES1", "HEY1",
    ],
    "Nrf2": [
        "NFE2L2", "KEAP1", "CUL3",
    ],
    "PI3K": [
        "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2", "PIK3R3", "PTEN", "AKT1",
        "AKT2", "AKT3", "MTOR", "TSC1", "TSC2", "RHEB", "RICTOR", "RPTOR", "STK11",
        "INPP4B", "PDPK1", "PPP2R1A", "DEPTOR",
    ],
    "RTK_RAS": [
        "EGFR", "ERBB2", "ERBB3", "ERBB4", "MET", "FGFR1", "FGFR2", "FGFR3", "FGFR4",
        "PDGFRA", "PDGFRB", "KIT", "FLT3", "ALK", "RET", "ROS1", "NTRK1", "NTRK2",
        "NTRK3", "IGF1R", "KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "ARAF", "NF1",
        "RASA1", "PTPN11", "SOS1", "MAP2K1", "MAP2K2", "MAPK1", "RIT1", "CBL",
        "ERRFI1", "IRS2", "JAK2", "KDR", "ABL1",
    ],
    "TGF_Beta": [
        "TGFBR1", "TGFBR2", "ACVR2A", "ACVR1B", "SMAD2", "SMAD3", "SMAD4",
    ],
    "TP53": [
        "TP53", "MDM2", "MDM4", "ATM", "CHEK2", "RPS6KA3",
    ],
    "Wnt": [
        "APC", "CTNNB1", "AXIN1", "AXIN2", "AMER1", "RNF43", "ZNRF3", "TCF7L2",
        "TCF7", "LEF1", "GSK3B", "LRP5", "LRP6", "DKK1", "DKK2", "DKK3", "DKK4",
        "FBXW11", "CHD8", "SOX9", "LZTR1",
    ],
}


def emit(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))


def fail(message: str, code: int = 1) -> None:
    emit({"error": message})
    sys.exit(code)


def info(message: str) -> None:
    print(f"[genomics] {message}", file=sys.stderr)


def read_maf(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        fail(f"MAF file not found: {path}")
    sep = "," if path.lower().endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, comment="#", low_memory=False)
    missing = [c for c in MANDATORY_MAF_COLS if c not in df.columns]
    if missing:
        fail(f"MAF missing mandatory GDC columns: {', '.join(missing)}")
    df["Hugo_Symbol"] = df["Hugo_Symbol"].astype(str).str.strip()
    df["Tumor_Sample_Barcode"] = df["Tumor_Sample_Barcode"].astype(str).str.strip()
    df["Variant_Classification"] = df["Variant_Classification"].astype(str).str.strip()
    df["Variant_Type"] = df["Variant_Type"].astype(str).str.strip().str.upper()
    df["Chromosome"] = df["Chromosome"].astype(str).str.strip()
    df["Start_Position"] = pd.to_numeric(df["Start_Position"], errors="coerce")
    df["End_Position"] = pd.to_numeric(df["End_Position"], errors="coerce")
    return df


def read_table(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        fail(f"file not found: {path}")
    sep = "," if path.lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, comment="#", low_memory=False)


def read_seg(path: str) -> pd.DataFrame:
    df = read_table(path)
    lower = {c.lower(): c for c in df.columns}
    sample_c = lower.get("sample") or lower.get("id") or lower.get("tumor_sample_barcode") or df.columns[0]
    chr_c = lower.get("chromosome") or lower.get("chrom") or lower.get("chr")
    start_c = lower.get("start") or lower.get("loc.start") or lower.get("start_position")
    end_c = lower.get("end") or lower.get("loc.end") or lower.get("end_position")
    mean_c = lower.get("segment_mean") or lower.get("seg.mean") or lower.get("log2") or lower.get("cn_log2")
    if not all([chr_c, start_c, end_c, mean_c]):
        fail("segment file must contain Sample, Chromosome, Start, End, and Segment_Mean (or log2) columns")
    out = pd.DataFrame(
        {
            "Sample": df[sample_c].astype(str).str.strip(),
            "Chromosome": df[chr_c].astype(str).str.strip().str.replace("^chr", "", regex=True),
            "Start": pd.to_numeric(df[start_c], errors="coerce"),
            "End": pd.to_numeric(df[end_c], errors="coerce"),
            "Segment_Mean": pd.to_numeric(df[mean_c], errors="coerce"),
        }
    ).dropna()
    return out


def extract_vaf(df: pd.DataFrame, vaf_col: str | None = None) -> tuple[pd.Series | None, str | None]:
    """Return normalized [0, 1] VAF series and source description, or (None, None)."""
    lower = {c.lower(): c for c in df.columns}
    chosen = None
    if vaf_col:
        if vaf_col not in df.columns:
            fail(f"vaf_col '{vaf_col}' not found in MAF")
        chosen = vaf_col
    else:
        for cand in ("t_vaf", "vaf", "tumor_vaf", "i_tumorvaf_wu", "af", "tumor_af"):
            if cand in lower:
                chosen = lower[cand]
                break
    if chosen is not None:
        v = pd.to_numeric(df[chosen], errors="coerce")
        if v.notna().any() and float(v.quantile(0.95)) > 1.0:
            v = v / 100.0
        return v.clip(0.0, 1.0), chosen

    alt_c = lower.get("t_alt_count")
    dep_c = lower.get("t_depth")
    ref_c = lower.get("t_ref_count")
    if alt_c and dep_c:
        alt = pd.to_numeric(df[alt_c], errors="coerce")
        dep = pd.to_numeric(df[dep_c], errors="coerce")
        v = alt / dep.replace(0, np.nan)
        return v.clip(0.0, 1.0), f"{alt_c}/{dep_c}"
    if alt_c and ref_c:
        alt = pd.to_numeric(df[alt_c], errors="coerce")
        ref = pd.to_numeric(df[ref_c], errors="coerce")
        v = alt / (alt + ref).replace(0, np.nan)
        return v.clip(0.0, 1.0), f"{alt_c}/({alt_c}+{ref_c})"
    return None, None


def pyrimidine_sub(ref: str, alt: str) -> str | None:
    r = str(ref).strip().upper()
    a = str(alt).strip().upper()
    if len(r) != 1 or len(a) != 1 or r == a or r not in COMPLEMENT or a not in COMPLEMENT:
        return None
    if r in ("C", "T"):
        return f"{r}>{a}"
    return f"{COMPLEMENT[r]}>{COMPLEMENT[a]}"


def canonical_sbs96(ref: str, alt: str, context: str) -> str | None:
    """Fold a 3-mer context (e.g. 'ACG') and SNV (ref->alt) into canonical COSMIC SBS96 key 'A[C>T]G'."""
    r = str(ref).strip().upper()
    a = str(alt).strip().upper()
    ctx = str(context).strip().upper()
    if len(r) != 1 or len(a) != 1 or r == a or r not in COMPLEMENT or a not in COMPLEMENT:
        return None
    if len(ctx) < 3 or any(b not in COMPLEMENT for b in ctx):
        return None
    mid = len(ctx) // 2
    tri = ctx[mid - 1 : mid + 2]
    if tri[1] == r:
        left, right = tri[0], tri[2]
        if r in ("C", "T"):
            return f"{left}[{r}>{a}]{right}"
        return f"{COMPLEMENT[right]}[{COMPLEMENT[r]}>{COMPLEMENT[a]}]{COMPLEMENT[left]}"
    if tri[1] == COMPLEMENT[r]:
        left, right = tri[0], tri[2]
        if tri[1] in ("C", "T"):
            return f"{left}[{tri[1]}>{COMPLEMENT[a]}]{right}"
    return None


def all_sbs96_channels() -> list[str]:
    bases = ["A", "C", "G", "T"]
    channels = []
    for sub in SIX_CLASSES:
        for b5 in bases:
            for b3 in bases:
                channels.append(f"{b5}[{sub}]{b3}")
    return channels


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


def load_clinical(path: str, sample_col: str | None = None) -> tuple[pd.DataFrame, str]:
    clin = read_table(path)
    if sample_col:
        if sample_col not in clin.columns:
            fail(f"sample_col '{sample_col}' not found in clinical table")
        scol = sample_col
    elif "Tumor_Sample_Barcode" in clin.columns:
        scol = "Tumor_Sample_Barcode"
    else:
        scol = clin.columns[0]
    clin[scol] = clin[scol].astype(str).str.strip()
    return clin, scol


def write_artifact(
    path: str | os.PathLike[str],
    tool: str,
    title: str,
    *,
    skill: str = "genomics-workflow",
    inputs: tuple[str | None, ...] | list[str | None] = (),
    parameters: dict[str, Any] | None = None,
    shape: tuple[int, int] | None = None,
) -> str:
    """Write a contract-1.0.0 .phi-artifact.json sidecar beside path and return path."""
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
                header = fh.read(24)
            if header[:8] == b"\x89PNG\r\n\x1a\n":
                w, h = struct.unpack(">II", header[16:24])
                block["widthPx"], block["heightPx"] = int(w), int(h)
        descriptor["figure"] = block
    elif shape is not None:
        descriptor["table"] = {"rows": int(shape[0]), "columns": int(shape[1])}

    cwd = Path.cwd().resolve()
    prov_inputs = []
    for ip_raw in inputs:
        if not ip_raw:
            continue
        ip = Path(ip_raw)
        try:
            rel = ip.resolve().relative_to(cwd)
        except ValueError:
            continue
        entry: dict[str, Any] = {"path": str(rel)}
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
