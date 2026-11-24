#!/usr/bin/env python3
"""Compile ESMFold2 predicted PDB/mmCIF structures and optional sidecar JSON
scores into flat structure outputs and standardized lineage tables (TSV + JSON).

Reads per-residue CA pLDDT from the PDB B-factor column when no sidecar JSON
is present, and emits unified `complex_plddt` (0..1) and `confidence_score`
aliases compatible with DESIGN_FILTER.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any


def _mean_ca_bfactor_pdb(pdb_path: Path) -> float | None:
    vals: list[float] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM  ") or len(line) < 66:
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            vals.append(float(line[60:66]))
        except ValueError:
            continue
    return round(sum(vals) / len(vals), 4) if vals else None


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile ESMFold2 outputs into standardized tables.")
    parser.add_argument("--out-dir", required=True, help="Directory containing ESMFold2 predicted PDB/CIF files.")
    parser.add_argument("--manifest", default="", help="Optional lineage manifest JSON.")
    parser.add_argument("--step-name", default="esmfold2", help="Pipeline step name.")
    parser.add_argument("--out-structures-dir", required=True, help="Flat output directory for predicted structures.")
    parser.add_argument("--out-tsv", required=True, help="Output standardized TSV table path.")
    parser.add_argument("--out-json", required=True, help="Output standardized JSON list path.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    manifest: dict[str, dict[str, Any]] = {}
    if args.manifest and Path(args.manifest).is_file():
        manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))

    structs = sorted(
        p for p in out_dir.rglob("*") if p.is_file() and p.suffix.lower() in (".pdb", ".cif", ".mmcif")
    )
    if not structs:
        raise SystemExit(f"No ESMFold2 predicted PDB/CIF files found under {out_dir}")

    flat_dir = Path(args.out_structures_dir)
    flat_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for sp in structs:
        stem = sp.stem
        info = manifest.get(stem, {})
        design_id = str(info.get("design_id") or stem)
        parent_id = str(info.get("parent_design_id") or design_id)

        dst = flat_dir / f"{design_id}{sp.suffix}"
        shutil.copy2(sp, dst)

        plddt: float | None = None
        ptm: float | None = None
        sidecar = sp.with_suffix(".json")
        if sidecar.is_file():
            try:
                meta = json.loads(sidecar.read_text(encoding="utf-8"))
                if meta.get("plddt") is not None:
                    plddt = float(meta["plddt"])
                if meta.get("ptm") is not None:
                    ptm = float(meta["ptm"])
            except (json.JSONDecodeError, ValueError, TypeError):
                pass

        if plddt is None and sp.suffix.lower() == ".pdb":
            plddt = _mean_ca_bfactor_pdb(sp)

        complex_plddt = round(plddt / 100.0, 4) if (plddt is not None and plddt > 1.0) else plddt
        conf_score = ptm if ptm is not None else complex_plddt

        rows.append(
            {
                "design_id": design_id,
                "parent_design_id": parent_id,
                "step": args.step_name,
                "structure_path": str(dst),
                "parent_structure_path": info.get("parent_structure_path", ""),
                "confidence_score": conf_score,
                "plddt": plddt,
                "complex_plddt": complex_plddt,
                "ptm": ptm,
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
