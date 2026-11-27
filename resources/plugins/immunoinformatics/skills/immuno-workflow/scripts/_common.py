"""Shared I/O, statistics, biochemical indices, and artifact utilities for immunoinformatics (imm_*)."""
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

# Kyte & Doolittle (1982) hydropathy index
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Grantham (1974) amino-acid physicochemical properties: (composition, polarity, volume)
GRANTHAM_PROPS: dict[str, tuple[float, float, float]] = {
    "A": (0.00, 8.1, 31.0),
    "R": (0.65, 10.5, 124.0),
    "N": (1.33, 11.6, 56.0),
    "D": (1.38, 13.0, 54.0),
    "C": (2.75, 5.5, 55.0),
    "Q": (0.89, 10.5, 85.0),
    "E": (0.92, 12.3, 83.0),
    "G": (0.74, 9.0, 3.0),
    "H": (0.58, 10.4, 96.0),
    "I": (0.00, 5.2, 111.0),
    "L": (0.00, 4.9, 111.0),
    "K": (0.33, 11.3, 119.0),
    "M": (0.00, 5.7, 105.0),
    "F": (0.00, 5.2, 132.0),
    "P": (0.39, 8.0, 32.5),
    "S": (1.42, 9.2, 32.0),
    "T": (0.71, 8.6, 61.0),
    "W": (0.13, 5.4, 170.0),
    "Y": (0.20, 6.2, 136.0),
    "V": (0.00, 5.9, 84.0),
}

STANDARD_AAS = set(KYTE_DOOLITTLE.keys())
TCR_LOCI = {"TRA", "TRB", "TRG", "TRD"}
BCR_LOCI = {"IGH", "IGK", "IGL"}


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
    """Read TSV/CSV input table with delimiter auto-detection and #-header preservation."""
    p = Path(path)
    if not p.is_file():
        fail(f"Input file not found: {path}")

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
                # Tools like TRUST4 (#count\tfrequency...) and seq2HLA (#Locus\tAllele 1...)
                # prefix the column header row with '#'. Preserve it if it looks tabular.
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
    elif name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".txt"):
        sep = "\t" if "\t" in first else ","
    else:
        sep = "\t" if "\t" in first else ","

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


def kyte_doolittle_gravy(seq: str) -> float | None:
    """Compute mean Kyte-Doolittle hydropathy across standard amino acids."""
    clean = [ch for ch in str(seq).strip().upper() if ch in KYTE_DOOLITTLE]
    if not clean:
        return None
    return float(sum(KYTE_DOOLITTLE[ch] for ch in clean) / len(clean))


def grantham_distance(aa1: str, aa2: str) -> float:
    """Compute Grantham (1974) physicochemical distance between two standard amino acids."""
    a1 = aa1.upper()
    a2 = aa2.upper()
    if a1 == a2:
        return 0.0
    if a1 not in GRANTHAM_PROPS or a2 not in GRANTHAM_PROPS:
        return 100.0
    c1, p1, v1 = GRANTHAM_PROPS[a1]
    c2, p2, v2 = GRANTHAM_PROPS[a2]
    return math.sqrt(1.833 * (c1 - c2) ** 2 + 0.1018 * (p1 - p2) ** 2 + 0.000399 * (v1 - v2) ** 2)


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
    skill: str = "immuno-workflow",
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


def _clean_gene_hit(raw: Any, strip_allele: bool = False) -> str:
    """Normalize a V/D/J/C gene call, stripping MiXCR scores and multi-hit commas."""
    if raw is None or (isinstance(raw, float) and math.isnan(raw)):
        return ""
    s = str(raw).strip()
    if not s or s.lower() in {"na", "nan", "none", ".", "unresolved"}:
        return ""
    # Take primary hit before comma/whitespace/score
    s = re.split(r"[,;\s]+", s)[0]
    s = re.sub(r"\([^)]*\)$", "", s).strip()
    if strip_allele:
        s = s.split("*")[0]
    return s


def normalize_gene_base(raw: Any) -> str:
    """Return gene symbol without IMGT allele suffix (e.g. TRBV12-3*01 -> TRBV12-3)."""
    return _clean_gene_hit(raw, strip_allele=True)


def infer_locus(locus_val: Any, v_call: str = "", j_call: str = "", c_call: str = "") -> str:
    """Normalize or infer immune receptor locus (TRA/TRB/TRG/TRD/IGH/IGK/IGL)."""
    if locus_val is not None and not (isinstance(locus_val, float) and math.isnan(locus_val)):
        txt = str(locus_val).strip().upper()
        for tok in ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"):
            if txt == tok or txt.startswith(tok):
                return tok
    for gene in (v_call, j_call, c_call):
        g = str(gene).strip().upper()
        for tok in ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"):
            if g.startswith(tok):
                return tok
    return "UNKNOWN"


def locus_to_receptor_type(locus: str) -> str:
    """Map normalized locus to TCR, BCR, or UNKNOWN."""
    if locus in TCR_LOCI:
        return "TCR"
    if locus in BCR_LOCI:
        return "BCR"
    return "UNKNOWN"


def _infer_sample_id_from_path(path: Path) -> str:
    stem = path.name
    for ext in (".tsv.gz", ".csv.gz", ".tsv", ".csv", ".txt"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
            break
    for suf in ("_filtered_contig_annotations", "_all_contig_annotations", "_airr", "_report", "_clones"):
        if stem.lower().endswith(suf):
            stem = stem[: -len(suf)]
            break
    return stem or path.stem


def detect_repertoire_format(df: pd.DataFrame) -> str:
    """Auto-detect repertoire table format: airr, cellranger_vdj, trust4, or mixcr."""
    cols = {str(c).strip().lower() for c in df.columns}
    if {"aaseqcdr3", "clonecount"} & cols or "allvhitswithscore" in cols or "bestvhit" in cols:
        return "mixcr"
    if "cdr3aa" in cols and ("cdr3nt" in cols or "cid" in cols or {"v", "j"}.issubset(cols)):
        return "trust4"
    if "barcode" in cols and ("umis" in cols or "raw_clonotype_id" in cols or "contig_id" in cols or "chain" in cols):
        return "cellranger_vdj"
    if {"junction_aa", "cdr3_aa", "v_call", "j_call", "duplicate_count"} & cols:
        return "airr"
    fail(f"Could not detect repertoire format from columns: {list(df.columns)}")
    return "airr"


def _parse_bool_flag(val: Any) -> bool | None:
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return None
    if isinstance(val, (bool, np.bool_)):
        return bool(val)
    s = str(val).strip().lower()
    if s in {"true", "t", "1", "yes", "y"}:
        return True
    if s in {"false", "f", "0", "no", "n"}:
        return False
    return None


def harmonize_repertoire_file(
    path: str | Path,
    fmt: str = "auto",
    default_sample_id: str | None = None,
) -> tuple[pd.DataFrame, str]:
    """Read and harmonize one AIRR, 10x V(D)J, TRUST4, or MiXCR file to AIRR-compatible schema."""
    p = Path(path)
    df = read_table(p)
    detected = detect_repertoire_format(df) if fmt == "auto" else fmt
    fallback_sid = default_sample_id or _infer_sample_id_from_path(p)

    s_col = resolve_column(df, ["sample_id", "sample", "subject_id", "repertoire_id"], "sample_id", required=False)
    cell_col = resolve_column(df, ["cell_id", "barcode", "cid"], "cell_id", required=False)
    locus_col = resolve_column(df, ["locus", "chain"], "locus", required=False)
    prod_col = resolve_column(df, ["productive"], "productive", required=False)

    if detected == "mixcr":
        aa_col = resolve_column(df, ["aaSeqCDR3", "cdr3_aa", "junction_aa"], "CDR3 AA")
        nt_col = resolve_column(df, ["nSeqCDR3", "targetSequences", "cdr3_nt", "junction"], "CDR3 nt", required=False)
        v_col = resolve_column(df, ["bestVHit", "allVHitsWithScore", "v_call", "V"], "V gene", required=False)
        d_col = resolve_column(df, ["bestDHit", "allDHitsWithScore", "d_call", "D"], "D gene", required=False)
        j_col = resolve_column(df, ["bestJHit", "allJHitsWithScore", "j_call", "J"], "J gene", required=False)
        c_col = resolve_column(df, ["bestCHit", "allCHitsWithScore", "c_call", "C"], "C gene", required=False)
        cnt_col = resolve_column(df, ["cloneCount", "readCount", "duplicate_count", "count"], "count", required=False)
        freq_col = resolve_column(df, ["cloneFraction", "frequency"], "frequency", required=False)
    elif detected == "trust4":
        aa_col = resolve_column(df, ["CDR3aa", "cdr3_aa", "junction_aa"], "CDR3 AA")
        nt_col = resolve_column(df, ["CDR3nt", "cdr3_nt", "junction"], "CDR3 nt", required=False)
        v_col = resolve_column(df, ["V", "v_call"], "V gene", required=False)
        d_col = resolve_column(df, ["D", "d_call"], "D gene", required=False)
        j_col = resolve_column(df, ["J", "j_call"], "J gene", required=False)
        c_col = resolve_column(df, ["C", "c_call"], "C gene", required=False)
        cnt_col = resolve_column(df, ["count", "#count", "duplicate_count", "consensus_count"], "count", required=False)
        freq_col = resolve_column(df, ["frequency"], "frequency", required=False)
    elif detected == "cellranger_vdj":
        aa_col = resolve_column(df, ["cdr3", "cdr3_aa", "junction_aa"], "CDR3 AA")
        nt_col = resolve_column(df, ["cdr3_nt", "junction"], "CDR3 nt", required=False)
        v_col = resolve_column(df, ["v_gene", "v_call"], "V gene", required=False)
        d_col = resolve_column(df, ["d_gene", "d_call"], "D gene", required=False)
        j_col = resolve_column(df, ["j_gene", "j_call"], "J gene", required=False)
        c_col = resolve_column(df, ["c_gene", "c_call"], "C gene", required=False)
        cnt_col = resolve_column(df, ["umis", "reads", "duplicate_count"], "count", required=False)
        freq_col = resolve_column(df, ["frequency"], "frequency", required=False)
    else:  # airr
        aa_col = resolve_column(df, ["cdr3_aa", "junction_aa", "CDR3aa", "aaSeqCDR3"], "CDR3 AA", required=False)
        nt_col = resolve_column(df, ["cdr3_nt", "cdr3", "junction", "CDR3nt"], "CDR3 nt", required=False)
        if aa_col is None and nt_col is None:
            fail(f"AIRR table lacks junction_aa / cdr3_aa column (found {list(df.columns)})")
        v_col = resolve_column(df, ["v_call", "v_gene", "V"], "V gene", required=False)
        d_col = resolve_column(df, ["d_call", "d_gene", "D"], "D gene", required=False)
        j_col = resolve_column(df, ["j_call", "j_gene", "J"], "J gene", required=False)
        c_col = resolve_column(df, ["c_call", "c_gene", "C"], "C gene", required=False)
        cnt_col = resolve_column(
            df,
            ["duplicate_count", "consensus_count", "umi_count", "count", "reads", "umis"],
            "count",
            required=False,
        )
        freq_col = resolve_column(df, ["frequency", "cloneFraction"], "frequency", required=False)

    rows: list[dict[str, Any]] = []
    for idx in range(len(df)):
        row = df.iloc[idx]
        sid_raw = row[s_col] if s_col else None
        sid = str(sid_raw).strip() if sid_raw is not None and str(sid_raw).strip() and str(sid_raw).strip().lower() != "nan" else fallback_sid

        v_call = _clean_gene_hit(row[v_col]) if v_col else ""
        d_call = _clean_gene_hit(row[d_col]) if d_col else ""
        j_call = _clean_gene_hit(row[j_col]) if j_col else ""
        c_call = _clean_gene_hit(row[c_col]) if c_col else ""
        locus = infer_locus(row[locus_col] if locus_col else None, v_call, j_call, c_call)
        rtype = locus_to_receptor_type(locus)

        raw_aa = "" if not aa_col or pd.isna(row[aa_col]) else str(row[aa_col]).strip()
        raw_nt = "" if not nt_col or pd.isna(row[nt_col]) else str(row[nt_col]).strip().upper()
        if raw_aa.lower() in {"nan", "none", "na", "."}:
            raw_aa = ""
        if raw_nt.lower() in {"nan", "none", "na", "."}:
            raw_nt = ""

        aa_upper = raw_aa.upper()
        has_stop = "*" in aa_upper
        is_oof = False
        if aa_upper in {"OUT_OF_FRAME", "FRAMESHIFT", "NONPRODUCTIVE"} or "_" in aa_upper:
            is_oof = True
        elif raw_nt and len(raw_nt) % 3 != 0:
            is_oof = True

        prod_flag = _parse_bool_flag(row[prod_col]) if prod_col else None
        clean_aa = aa_upper if (aa_upper and not is_oof and aa_upper not in {"PARTIAL", "NA", "NONE"}) else ""
        if clean_aa and not all(ch in STANDARD_AAS or ch == "*" for ch in clean_aa):
            # Non-standard ambiguity characters like X/?/~
            is_valid_chars = all(ch in STANDARD_AAS for ch in clean_aa)
        else:
            is_valid_chars = bool(clean_aa) and not has_stop

        if prod_flag is not None:
            is_prod = bool(prod_flag and not has_stop and not is_oof and bool(clean_aa))
            if not is_prod and not has_stop and not is_oof and (not clean_aa or (raw_nt and len(raw_nt) % 3 != 0)):
                is_oof = True
        else:
            is_prod = bool(clean_aa and is_valid_chars and not has_stop and not is_oof)

        cnt_val = 1.0
        if cnt_col:
            cnum = pd.to_numeric(row[cnt_col], errors="coerce")
            if pd.notna(cnum) and float(cnum) > 0:
                cnt_val = float(cnum)
            elif pd.notna(cnum) and float(cnum) == 0:
                cnt_val = 0.0
        freq_val = np.nan
        if freq_col:
            fnum = pd.to_numeric(row[freq_col], errors="coerce")
            if pd.notna(fnum):
                freq_val = float(fnum)

        cell_id = ""
        if cell_col and not pd.isna(row[cell_col]):
            cid_s = str(row[cell_col]).strip()
            if cid_s.lower() not in {"nan", "none", "na"}:
                cell_id = cid_s

        gravy = kyte_doolittle_gravy(clean_aa) if is_prod else None
        aa_len = len(clean_aa) if (clean_aa and not has_stop) else 0

        rows.append(
            {
                "sample_id": sid,
                "locus": locus,
                "receptor_type": rtype,
                "v_call": v_call,
                "d_call": d_call,
                "j_call": j_call,
                "c_call": c_call,
                "cdr3_aa": clean_aa if clean_aa else raw_aa,
                "cdr3_nt": raw_nt,
                "cdr3_length": int(aa_len),
                "kyte_doolittle_gravy": gravy if gravy is not None else np.nan,
                "productive": bool(is_prod),
                "stop_codon": bool(has_stop),
                "out_of_frame": bool(is_oof),
                "duplicate_count": float(cnt_val),
                "frequency": freq_val,
                "cell_id": cell_id,
                "source_format": detected,
            }
        )

    out_df = pd.DataFrame(rows)
    if out_df.empty:
        fail(f"No repertoire records parsed from {path}")

    # Compute per-sample frequency if not provided
    for sid, grp_idx in out_df.groupby("sample_id").groups.items():
        if out_df.loc[grp_idx, "frequency"].isna().all():
            tot = float(out_df.loc[grp_idx, "duplicate_count"].sum())
            if tot > 0:
                out_df.loc[grp_idx, "frequency"] = out_df.loc[grp_idx, "duplicate_count"] / tot
            else:
                out_df.loc[grp_idx, "frequency"] = 0.0
    return out_df, detected


def harmonize_repertoire_inputs(
    paths: list[str | Path],
    fmt: str = "auto",
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Harmonize one or more repertoire input files into one standardized DataFrame."""
    frames: list[pd.DataFrame] = []
    fmt_map: dict[str, str] = {}
    used_sids: set[str] = set()
    for raw in paths:
        p = Path(raw)
        base_sid = _infer_sample_id_from_path(p)
        sid = base_sid
        k = 2
        while sid in used_sids:
            sid = f"{base_sid}_{k}"
            k += 1
        df_i, det_i = harmonize_repertoire_file(p, fmt=fmt, default_sample_id=sid)
        used_sids.update(df_i["sample_id"].astype(str).unique().tolist())
        fmt_map[str(p.name)] = det_i
        frames.append(df_i)
    combined = pd.concat(frames, ignore_index=True)
    return combined, fmt_map
