#!/usr/bin/env python3
"""Compile Chai-1 prediction outputs (pred.model_idx_*.cif and
scores.model_idx_*.npz or .json) into flat best-model structures and
standardized lineage tables (TSV + JSON).

Pure-Python implementation including a lightweight .npy/.npz scalar reader so
zero third-party packages (not even numpy) are required.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import re
import shutil
import struct
import zipfile
from pathlib import Path
from typing import Any

SCORE_RE = re.compile(r"^scores\.model_idx_(?P<idx>\d+)\.(?:npz|json)$")


def _read_npy_scalar(raw: bytes) -> float | bool | None:
    """Decode a 1-element numeric/bool NumPy .npy payload in pure Python."""
    if not raw.startswith(b"\x93NUMPY"):
        return None
    major, minor = raw[6], raw[7]
    if major == 1:
        hlen = struct.unpack("<H", raw[8:10])[0]
        hdr_start = 10
    elif major >= 2:
        hlen = struct.unpack("<I", raw[8:12])[0]
        hdr_start = 12
    else:
        return None
    header = ast.literal_eval(raw[hdr_start : hdr_start + hlen].decode("latin1").strip())
    descr = str(header.get("descr", ""))
    payload = raw[hdr_start + hlen :]
    endian = ">" if descr.startswith(">") else "<"
    code = descr.lstrip("<>|=")
    if code == "b1" and len(payload) >= 1:
        return bool(payload[0])
    fmt_map = {"f4": "f", "f8": "d", "i4": "i", "i8": "q", "u4": "I", "u8": "Q"}
    if code in fmt_map:
        sz = struct.calcsize(endian + fmt_map[code])
        if len(payload) >= sz:
            return float(struct.unpack(endian + fmt_map[code], payload[:sz])[0])
    return None


def _read_scores_file(path: Path) -> dict[str, Any]:
    if path.suffix.lower() == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    out: dict[str, Any] = {}
    with zipfile.ZipFile(path, "r") as zf:
        for name in zf.namelist():
            if not name.endswith(".npy"):
                continue
            key = name[:-4]
            if key in ("aggregate_score", "ptm", "iptm", "has_inter_chain_clashes"):
                out[key] = _read_npy_scalar(zf.read(name))
    return out


def _maybe_float(val: Any) -> float | None:
    if val is None or isinstance(val, bool):
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile Chai-1 outputs into standardized tables.")
    parser.add_argument("--out-dir", required=True, help="Directory containing per-design Chai-1 output subdirectories.")
    parser.add_argument("--manifest", default="", help="Optional lineage manifest JSON.")
    parser.add_argument("--step-name", default="chai1", help="Pipeline step name.")
    parser.add_argument("--out-structures-dir", required=True, help="Flat output directory for best predicted CIFs.")
    parser.add_argument("--out-tsv", required=True, help="Output standardized TSV table path.")
    parser.add_argument("--out-json", required=True, help="Output standardized JSON list path.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    manifest: dict[str, dict[str, Any]] = {}
    if args.manifest and Path(args.manifest).is_file():
        manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))

    job_dirs = sorted({p.parent for p in out_dir.rglob("pred.model_idx_*.cif") if p.is_file()})
    if not job_dirs:
        raise SystemExit(f"No Chai-1 pred.model_idx_*.cif files found under {out_dir}")

    flat_dir = Path(args.out_structures_dir)
    flat_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for job_dir in job_dirs:
        job_name = job_dir.name
        best_idx = "0"
        best_scores: dict[str, Any] = {}
        best_agg = float("-inf")

        for sf in sorted(job_dir.iterdir()):
            m = SCORE_RE.match(sf.name)
            if not m:
                continue
            data = _read_scores_file(sf)
            agg = _maybe_float(data.get("aggregate_score"))
            if agg is None:
                agg = _maybe_float(data.get("iptm")) or _maybe_float(data.get("ptm")) or float("-inf")
            if agg > best_agg:
                best_agg = agg
                best_idx = m.group("idx")
                best_scores = data

        best_cif = job_dir / f"pred.model_idx_{best_idx}.cif"
        if not best_cif.is_file():
            cifs = sorted(job_dir.glob("pred.model_idx_*.cif"))
            best_cif = cifs[0] if cifs else best_cif

        info = manifest.get(job_name, {})
        design_id = str(info.get("design_id") or job_name)
        parent_id = str(info.get("parent_design_id") or design_id)

        copied_struct = ""
        if best_cif.is_file():
            dst = flat_dir / f"{design_id}.cif"
            shutil.copy2(best_cif, dst)
            copied_struct = str(dst)

        agg_score = _maybe_float(best_scores.get("aggregate_score"))
        ptm = _maybe_float(best_scores.get("ptm"))
        iptm = _maybe_float(best_scores.get("iptm"))
        has_clash = best_scores.get("has_inter_chain_clashes")

        rows.append(
            {
                "design_id": design_id,
                "parent_design_id": parent_id,
                "step": args.step_name,
                "structure_path": copied_struct,
                "parent_structure_path": info.get("parent_structure_path", ""),
                "confidence_score": agg_score if agg_score is not None else iptm,
                "aggregate_score": agg_score,
                "ptm": ptm,
                "iptm": iptm,
                "has_inter_chain_clashes": has_clash,
            }
        )

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})

    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
