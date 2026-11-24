#!/usr/bin/env python3
"""Compile AlphaFold 3 prediction outputs (<name>_model.cif,
<name>_summary_confidences.json, optional <name>_confidences.json) into flat
structures and standardized lineage tables (TSV + JSON).

Emits both AF3-native fields (ranking_score, ptm, iptm, fraction_disordered,
has_clash, min_pae_interaction, plddt) and cross-predictor aliases
(confidence_score, complex_plddt) for seamless use with DESIGN_FILTER.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any


def _compute_min_pae_interaction(matrix: Any) -> float | None:
    if not isinstance(matrix, list) or len(matrix) < 2:
        return None
    off_diag: list[float] = []
    for i, row in enumerate(matrix):
        if not isinstance(row, list):
            continue
        for j, val in enumerate(row):
            if i == j or val is None:
                continue
            try:
                off_diag.append(float(val))
            except (TypeError, ValueError):
                continue
    return round(min(off_diag), 4) if off_diag else None


def _compute_mean_plddt(conf_json: Path) -> float | None:
    if not conf_json.is_file():
        return None
    try:
        data = json.loads(conf_json.read_text(encoding="utf-8"))
        atom_plddts = data.get("atom_plddts")
        if isinstance(atom_plddts, list) and atom_plddts:
            return round(sum(float(v) for v in atom_plddts) / len(atom_plddts), 4)
    except (json.JSONDecodeError, ValueError, TypeError):
        return None
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile AlphaFold 3 outputs into standardized tables.")
    parser.add_argument("--out-dir", required=True, help="AlphaFold 3 --output_dir directory.")
    parser.add_argument("--manifest", default="", help="Optional lineage manifest JSON.")
    parser.add_argument("--step-name", default="alphafold3", help="Pipeline step name.")
    parser.add_argument("--out-structures-dir", required=True, help="Flat output directory for top predicted CIF files.")
    parser.add_argument("--out-tsv", required=True, help="Output standardized TSV table path.")
    parser.add_argument("--out-json", required=True, help="Output standardized JSON list path.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    manifest: dict[str, dict[str, Any]] = {}
    if args.manifest and Path(args.manifest).is_file():
        manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))

    # Locate top-level job directories containing *_summary_confidences.json (excluding seed-* subdirs)
    summary_files = sorted(
        p
        for p in out_dir.rglob("*_summary_confidences.json")
        if p.is_file() and not p.parent.name.startswith("seed-")
    )
    if not summary_files:
        summary_files = sorted(out_dir.rglob("*_summary_confidences.json"))
    if not summary_files:
        raise SystemExit(f"No AlphaFold 3 *_summary_confidences.json files found under {out_dir}")

    flat_struct_dir = Path(args.out_structures_dir)
    flat_struct_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for sf in summary_files:
        job_dir = sf.parent
        af3_name = sf.name[: -len("_summary_confidences.json")]
        info = manifest.get(af3_name, manifest.get(job_dir.name, {}))
        design_id = str(info.get("design_id") or af3_name)
        parent_id = str(info.get("parent_design_id") or design_id)

        model_cif = job_dir / f"{af3_name}_model.cif"
        if not model_cif.is_file():
            cifs = sorted(job_dir.glob("*.cif"))
            model_cif = cifs[0] if cifs else None  # type: ignore[assignment]

        copied_struct = ""
        if model_cif is not None and model_cif.is_file():
            dst = flat_struct_dir / f"{design_id}.cif"
            shutil.copy2(model_cif, dst)
            copied_struct = str(dst)

        summary = json.loads(sf.read_text(encoding="utf-8"))
        ranking_score = summary.get("ranking_score")
        ptm = summary.get("ptm")
        iptm = summary.get("iptm")
        frac_disordered = summary.get("fraction_disordered")
        has_clash = summary.get("has_clash")
        min_pae_inter = _compute_min_pae_interaction(summary.get("chain_pair_pae_min"))

        plddt = _compute_mean_plddt(job_dir / f"{af3_name}_confidences.json")
        complex_plddt = (
            round(plddt / 100.0, 4) if (plddt is not None and plddt > 1.0) else plddt
        )

        rows.append(
            {
                "design_id": design_id,
                "parent_design_id": parent_id,
                "step": args.step_name,
                "structure_path": copied_struct,
                "parent_structure_path": info.get("parent_structure_path", ""),
                "confidence_score": ranking_score,
                "ranking_score": ranking_score,
                "ptm": ptm,
                "iptm": iptm,
                "plddt": plddt,
                "complex_plddt": complex_plddt,
                "min_pae_interaction": min_pae_inter,
                "fraction_disordered": frac_disordered,
                "has_clash": has_clash,
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
