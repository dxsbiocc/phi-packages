#!/usr/bin/env python3
"""Compile Boltz-2 prediction outputs (structures, confidence JSON, affinity JSON)
into flat structure outputs and standardized lineage tables (TSV + JSON).
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any

CONFIDENCE_KEYS = (
    "confidence_score",
    "ptm",
    "iptm",
    "ligand_iptm",
    "protein_iptm",
    "complex_plddt",
    "complex_iplddt",
    "complex_pde",
    "complex_ipde",
)

AFFINITY_KEYS = (
    "affinity_pred_value",
    "affinity_probability_binary",
)


def _find_prediction_dirs(out_dir: Path) -> list[Path]:
    pred_roots = [p for p in out_dir.rglob("predictions") if p.is_dir()]
    job_dirs: list[Path] = []
    for root in sorted(pred_roots):
        for child in sorted(root.iterdir()):
            if child.is_dir():
                job_dirs.append(child)
    if job_dirs:
        return job_dirs
    # Fallback: treat out_dir subdirectories containing confidence_*.json as job dirs
    return sorted({p.parent for p in out_dir.rglob("confidence_*.json")})


def _select_best_model(job_dir: Path, design_id: str) -> tuple[Path | None, Path | None]:
    conf_files = sorted(job_dir.glob("confidence_*_model_*.json"))
    best_conf: Path | None = None
    best_score = float("-inf")
    best_model_idx = "0"

    for cf in conf_files:
        try:
            data = json.loads(cf.read_text(encoding="utf-8"))
            score = float(data.get("confidence_score", float("-inf")))
        except (ValueError, TypeError, json.JSONDecodeError):
            score = float("-inf")
        if best_conf is None or score > best_score:
            best_conf = cf
            best_score = score
            stem_parts = cf.stem.split("_model_")
            best_model_idx = stem_parts[-1] if len(stem_parts) > 1 else "0"

    struct_path: Path | None = None
    for ext in (".cif", ".pdb"):
        cand = job_dir / f"{design_id}_model_{best_model_idx}{ext}"
        if cand.is_file():
            struct_path = cand
            break
    if struct_path is None:
        for ext in (".cif", ".pdb"):
            matches = sorted(job_dir.glob(f"*{ext}"))
            if matches:
                struct_path = matches[0]
                break
    return struct_path, best_conf


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile Boltz-2 outputs into standardized tables.")
    parser.add_argument("--out-dir", required=True, help="Boltz-2 --out_dir directory.")
    parser.add_argument("--manifest", default="", help="Optional lineage manifest JSON.")
    parser.add_argument("--step-name", default="boltz2", help="Pipeline step name.")
    parser.add_argument("--out-structures-dir", required=True, help="Flat output directory for best predicted structures.")
    parser.add_argument("--out-tsv", required=True, help="Output standardized TSV table path.")
    parser.add_argument("--out-json", required=True, help="Output standardized JSON list path.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    manifest: dict[str, dict[str, Any]] = {}
    if args.manifest and Path(args.manifest).is_file():
        manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))

    job_dirs = _find_prediction_dirs(out_dir)
    if not job_dirs:
        raise SystemExit(f"No Boltz-2 prediction directories found under {out_dir}")

    flat_struct_dir = Path(args.out_structures_dir)
    flat_struct_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for job_dir in job_dirs:
        design_id = job_dir.name
        info = manifest.get(design_id, {})
        parent_id = info.get("parent_design_id", design_id)

        best_struct, best_conf = _select_best_model(job_dir, design_id)
        copied_struct = ""
        if best_struct is not None:
            dst = flat_struct_dir / f"{design_id}{best_struct.suffix}"
            shutil.copy2(best_struct, dst)
            copied_struct = str(dst)

        row: dict[str, Any] = {
            "design_id": design_id,
            "parent_design_id": parent_id,
            "step": args.step_name,
            "structure_path": copied_struct,
            "parent_structure_path": info.get("parent_structure_path", ""),
        }

        if best_conf is not None and best_conf.is_file():
            conf_data = json.loads(best_conf.read_text(encoding="utf-8"))
            for key in CONFIDENCE_KEYS:
                row[key] = conf_data.get(key)
        else:
            for key in CONFIDENCE_KEYS:
                row[key] = None

        affinity_files = sorted(job_dir.glob("affinity_*.json"))
        if affinity_files:
            aff_data = json.loads(affinity_files[0].read_text(encoding="utf-8"))
            for key in AFFINITY_KEYS:
                row[key] = aff_data.get(key)
        else:
            for key in AFFINITY_KEYS:
                row[key] = None

        rows.append(row)

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
