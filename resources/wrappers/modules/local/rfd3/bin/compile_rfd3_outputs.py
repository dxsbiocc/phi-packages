#!/usr/bin/env python3
"""Compile RFdiffusion3 (RFD3) output structures and sidecar JSON metadata
into standardized lineage tables (TSV + JSON).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

STRUCTURE_SUFFIXES = (".cif.gz", ".cif", ".pdb")


def _strip_structure_suffix(name: str) -> str:
    lower = name.lower()
    for suf in STRUCTURE_SUFFIXES:
        if lower.endswith(suf):
            return name[: -len(suf)]
    return Path(name).stem


def _find_sidecar_json(struct_path: Path) -> Path | None:
    stem = _strip_structure_suffix(struct_path.name)
    candidate = struct_path.with_name(f"{stem}.json")
    return candidate if candidate.is_file() else None


def _flatten_metrics(obj: Any, prefix: str = "") -> dict[str, Any]:
    flat: dict[str, Any] = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            if isinstance(v, (int, float, bool)) or (isinstance(v, str) and len(v) <= 120):
                flat[key] = v
            elif isinstance(v, dict):
                flat.update(_flatten_metrics(v, key))
    return flat


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile RFD3 output structures and metadata.")
    parser.add_argument("--out-dir", required=True, help="Directory containing RFD3 output files.")
    parser.add_argument("--step-name", default="rfd3", help="Pipeline step name recorded in lineage.")
    parser.add_argument("--out-tsv", required=True, help="Output standardized TSV table path.")
    parser.add_argument("--out-json", required=True, help="Output standardized JSON list path.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    structures = sorted(
        p
        for p in out_dir.rglob("*")
        if p.is_file() and any(p.name.lower().endswith(s) for s in STRUCTURE_SUFFIXES)
    )
    if not structures:
        raise SystemExit(f"No RFD3 output structures (.cif/.cif.gz/.pdb) found under {out_dir}")

    rows: list[dict[str, Any]] = []
    for sp in structures:
        design_id = _strip_structure_suffix(sp.name)
        row: dict[str, Any] = {
            "design_id": design_id,
            "parent_design_id": design_id,
            "step": args.step_name,
            "structure_path": str(sp),
        }
        sidecar = _find_sidecar_json(sp)
        if sidecar is not None:
            try:
                meta = json.loads(sidecar.read_text(encoding="utf-8"))
                row.update(_flatten_metrics(meta))
            except json.JSONDecodeError:
                pass
        rows.append(row)

    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for k in row.keys():
            if k not in seen:
                fieldnames.append(k)
                seen.add(k)

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fieldnames})

    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
