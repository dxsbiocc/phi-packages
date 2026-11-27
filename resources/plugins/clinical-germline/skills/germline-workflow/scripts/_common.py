"""Shared I/O, VCF/tabular normalization, Mendelian helpers, and artifact utilities for clinical-germline (germ_*)."""
from __future__ import annotations

import gzip
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

LOF_CONSEQUENCES: set[str] = {
    "stop_gained",
    "nonsense",
    "frameshift_variant",
    "frameshift",
    "frameshift_elongation",
    "frameshift_truncation",
    "splice_acceptor_variant",
    "splice_acceptor",
    "splice_donor_variant",
    "splice_donor",
    "start_lost",
    "initiator_codon_variant",
}

INFRAME_INDEL_CONSEQUENCES: set[str] = {
    "inframe_insertion",
    "inframe_deletion",
    "inframe_indel",
    "disruptive_inframe_insertion",
    "disruptive_inframe_deletion",
    "conservative_inframe_insertion",
    "conservative_inframe_deletion",
}

MISSENSE_CONSEQUENCES: set[str] = {
    "missense_variant",
    "missense",
    "nonsynonymous_snv",
    "non_synonymous_codon",
    "stop_lost",
    "protein_altering_variant",
    "protein_altering",
    "coding_sequence_variant",
}

SPLICE_REGION_CONSEQUENCES: set[str] = {
    "splice_region_variant",
    "splice_region",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
}

SYNONYMOUS_CONSEQUENCES: set[str] = {
    "synonymous_variant",
    "synonymous",
    "silent",
    "stop_retained_variant",
}

CODING_OR_SPLICE_CONSEQUENCES: set[str] = (
    LOF_CONSEQUENCES
    | INFRAME_INDEL_CONSEQUENCES
    | MISSENSE_CONSEQUENCES
    | SPLICE_REGION_CONSEQUENCES
)


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


def is_gzipped(path: Path) -> bool:
    """Detect gzipped file via extension or magic bytes."""
    if str(path).lower().endswith(".gz"):
        return True
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def read_text_lines(path: str | Path) -> list[str]:
    """Read plain or gzipped text lines safely."""
    p = Path(path)
    if not p.is_file():
        fail(f"Input file not found: {path}")
    try:
        if is_gzipped(p):
            with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
                return fh.readlines()
        with open(p, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.readlines()
    except Exception as exc:
        fail(f"Failed to read file '{path}': {exc}")
    return []


def is_vcf_file(path: str | Path) -> bool:
    """Return True if path looks like a VCF or VCF.gz file."""
    p = Path(path)
    name = p.name.lower()
    if name.endswith(".vcf") or name.endswith(".vcf.gz"):
        return True
    lines = read_text_lines(p)[:50]
    for line in lines:
        s = line.strip()
        if not s:
            continue
        if s.startswith("##fileformat=VCF") or s.startswith("#CHROM\tPOS\t"):
            return True
        if not s.startswith("#"):
            break
    return False


def read_table(path: str | Path) -> pd.DataFrame:
    """Read TSV/CSV/BED input table with delimiter auto-detection and #-header preservation."""
    p = Path(path)
    raw_lines = read_text_lines(p)

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
                if "\t" in after or "," in after:
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
    name = p.name.lower()
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        sep = ","
    elif "\t" in first:
        sep = "\t"
    else:
        sep = "," if "," in first else "\t"

    try:
        df = pd.read_csv(io.StringIO("".join(filtered_lines)), sep=sep, low_memory=False)
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


def safe_float(val: Any, default: float | None = None) -> float | None:
    """Convert scalar to finite float or return default."""
    if val is None:
        return default
    if isinstance(val, (float, np.floating)):
        v = float(val)
        return v if math.isfinite(v) else default
    if isinstance(val, (int, np.integer)):
        return float(val)
    s = str(val).strip()
    if not s or s.lower() in {".", "na", "nan", "none", "null", ""}:
        return default
    # Handle multi-value fields separated by comma/ampersand by taking first valid number
    for tok in re.split(r"[,;&]+", s):
        tok = tok.strip()
        if not tok or tok == ".":
            continue
        try:
            v = float(tok)
            if math.isfinite(v):
                return v
        except ValueError:
            continue
    return default


def safe_int(val: Any, default: int | None = None) -> int | None:
    """Convert scalar to int or return default."""
    fv = safe_float(val, None)
    if fv is None:
        return default
    return int(round(fv))


def norm_chrom(chrom: Any) -> str:
    """Normalize chromosome name to canonical chr1..chr22, chrX, chrY, chrM."""
    if chrom is None or (isinstance(chrom, float) and math.isnan(chrom)):
        return ""
    s = str(chrom).strip()
    if not s or s.lower() in {".", "na", "nan", "none"}:
        return ""
    u = s.upper()
    if u.startswith("CHR"):
        u = u[3:]
    if u in {"23"}:
        return "chrX"
    if u in {"24"}:
        return "chrY"
    if u in {"MT", "M"}:
        return "chrM"
    if u in {"X", "Y"}:
        return f"chr{u}"
    return f"chr{u}"


def norm_gene(gene: Any) -> str:
    """Normalize gene symbol to uppercase trimmed string or empty string."""
    if gene is None or (isinstance(gene, float) and math.isnan(gene)):
        return ""
    s = str(gene).strip()
    if not s or s.lower() in {".", "na", "nan", "none", "null", "unknown"}:
        return ""
    s = re.split(r"[,;|/\s]+", s)[0].strip()
    return s.upper()


def classify_consequence_group(consequence: Any) -> str:
    """Map VEP/SnpEff consequence string into canonical group."""
    if consequence is None or (isinstance(consequence, float) and math.isnan(consequence)):
        return "other"
    raw = str(consequence).strip().lower()
    if not raw or raw in {".", "na", "nan", "none"}:
        return "other"
    tokens = [t.strip() for t in re.split(r"[&,;|/\s]+", raw) if t.strip()]
    for tok in tokens:
        if tok in LOF_CONSEQUENCES:
            return "lof"
    for tok in tokens:
        if tok in INFRAME_INDEL_CONSEQUENCES:
            return "inframe_indel"
    for tok in tokens:
        if tok in MISSENSE_CONSEQUENCES:
            return "missense"
    for tok in tokens:
        if tok in SPLICE_REGION_CONSEQUENCES:
            return "splice_region"
    for tok in tokens:
        if tok in SYNONYMOUS_CONSEQUENCES:
            return "synonymous"
    return "other"


def is_coding_or_splice_consequence(consequence: Any) -> bool:
    """Return True if consequence is protein-altering or splice-region/canonical splice."""
    if consequence is None or (isinstance(consequence, float) and math.isnan(consequence)):
        return False
    tokens = [t.strip().lower() for t in re.split(r"[&,;|/\s]+", str(consequence)) if t.strip()]
    return any(tok in CODING_OR_SPLICE_CONSEQUENCES for tok in tokens)


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
    skill: str = "germline-workflow",
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
