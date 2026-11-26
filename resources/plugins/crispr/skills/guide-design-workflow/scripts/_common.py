"""Shared utilities for CRISPR guide design, off-target scan, and library QC."""
from __future__ import annotations

import hashlib
import json
import os
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

RESTRICTION_ENZYMES: dict[str, str] = {
    "BsmBI": "CGTCTC",
    "Esp3I": "CGTCTC",
    "BbsI": "GAAGAC",
    "BpiI": "GAAGAC",
    "BsaI": "GGTCTC",
    "SapI": "GCTCTTC",
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
}


def emit(payload: dict) -> None:
    sys.stdout.write(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()


def fail(message: str, code: int = 1) -> "NoReturn":
    sys.stdout.write(json.dumps({"error": message}, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()
    raise SystemExit(code)


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def parse_fasta(path: str | Path) -> list[tuple[str, str]]:
    p = Path(path)
    if not p.is_file():
        fail(f"FASTA file not found: {path}")
    records: list[tuple[str, str]] = []
    cur_id: str | None = None
    chunks: list[str] = []
    with open(p, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    records.append((cur_id, "".join(chunks).upper()))
                cur_id = line[1:].split()[0] or f"seq_{len(records) + 1}"
                chunks = []
            else:
                chunks.append(line)
    if cur_id is not None:
        records.append((cur_id, "".join(chunks).upper()))
    if not records:
        fail(f"No FASTA records found in {path}")
    return records


def read_table(path: str | Path, *, header: int | None = 0):
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
    return pd.read_csv(p, sep=sep, comment="#", header=header)


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


def parse_enzyme_list(raw_items: list[str] | None) -> list[str]:
    if not raw_items:
        return ["BsmBI", "BbsI"]
    out: list[str] = []
    for item in raw_items:
        for tok in str(item).split(","):
            tok = tok.strip()
            if tok and tok not in out:
                out.append(tok)
    return out or ["BsmBI", "BbsI"]


def find_enzyme_collisions(seq: str, enzymes: list[str]) -> list[str]:
    upper = seq.upper()
    hits: list[str] = []
    for enz in enzymes:
        motif = RESTRICTION_ENZYMES.get(enz)
        if not motif:
            for k, v in RESTRICTION_ENZYMES.items():
                if k.lower() == enz.lower():
                    motif = v
                    enz = k
                    break
        if not motif:
            continue
        rc = revcomp(motif)
        if motif in upper or rc in upper:
            hits.append(enz)
    return hits


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
    skill: str = "guide-design-workflow",
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
