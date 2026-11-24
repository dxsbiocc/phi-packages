#!/usr/bin/env python3
"""Prepare AlphaFold 3 JSON job inputs from either:
  1. Existing AlphaFold 3 JSON file(s) or directory (--inputs), OR
  2. A sequence-native TSV/CSV table (--sequences-tsv, e.g. from ligandmpnn)
     plus optional target sequence and ligand CCD/SMILES.

Pure-Python implementation (no third-party dependencies required).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import string
from pathlib import Path
from typing import Any


def _sanitize_af3_name(name: str) -> str:
    cleaned = re.sub(r"[^a-zA-Z0-9_.-]+", "_", name.strip()).lower()
    return cleaned or "design"


def _parse_seeds(raw: str) -> list[int]:
    seeds: list[int] = []
    for part in raw.split(","):
        p = part.strip()
        if p:
            seeds.append(int(p))
    return seeds or [1]


def _stage_existing_inputs(inputs: list[str], out_dir: Path) -> dict[str, dict[str, Any]]:
    manifest: dict[str, dict[str, Any]] = {}
    for raw in inputs:
        p = Path(raw)
        candidates = (
            [c for c in sorted(p.rglob("*.json")) if c.is_file()]
            if p.is_dir()
            else ([p] if p.is_file() and p.suffix.lower() == ".json" else [])
        )
        for src in candidates:
            dst = out_dir / src.name
            shutil.copy2(src, dst)
            try:
                data = json.loads(src.read_text(encoding="utf-8"))
                job_name = _sanitize_af3_name(str(data.get("name") or src.stem))
            except (json.JSONDecodeError, AttributeError):
                job_name = _sanitize_af3_name(src.stem)
            manifest[job_name] = {
                "design_id": src.stem,
                "parent_design_id": src.stem,
                "source_input": str(src.resolve()),
            }
    return manifest


def _stage_from_sequences_table(
    table_path: Path,
    out_dir: Path,
    target_sequence: str,
    ligand_ccd: str,
    ligand_smiles: str,
    model_seeds: list[int],
    run_data_pipeline: bool,
) -> dict[str, dict[str, Any]]:
    delimiter = "\t" if table_path.suffix.lower() in (".tsv", ".tab") else ","
    manifest: dict[str, dict[str, Any]] = {}
    chain_pool = list(string.ascii_uppercase)

    with table_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for idx, row in enumerate(reader, start=1):
            design_id = (row.get("design_id") or f"design_{idx}").strip()
            parent_id = (row.get("parent_design_id") or design_id).strip()
            raw_seq = (row.get("sequence") or "").strip()
            if not raw_seq:
                continue

            parts = [p.strip() for p in raw_seq.replace("/", ":").split(":") if p.strip()]
            if target_sequence.strip():
                parts.append(target_sequence.strip())

            sequences: list[dict[str, Any]] = []
            used_chains: set[str] = set()
            for i, seq in enumerate(parts):
                cid = chain_pool[i]
                used_chains.add(cid)
                prot_entry: dict[str, Any] = {
                    "id": cid,
                    "sequence": seq,
                }
                if not run_data_pipeline:
                    prot_entry["unpairedMsa"] = ""
                    prot_entry["pairedMsa"] = ""
                    prot_entry["templates"] = []
                sequences.append({"protein": prot_entry})

            if ligand_ccd.strip() or ligand_smiles.strip():
                lig_cid = next((c for c in "LXYZMNOPQRSTUVW" if c not in used_chains), "L")
                if ligand_ccd.strip():
                    ccds = [c.strip() for c in ligand_ccd.split(",") if c.strip()]
                    sequences.append({"ligand": {"id": lig_cid, "ccdCodes": ccds}})
                else:
                    sequences.append({"ligand": {"id": lig_cid, "smiles": ligand_smiles.strip()}})

            af3_name = _sanitize_af3_name(design_id)
            job_json = {
                "name": af3_name,
                "modelSeeds": model_seeds,
                "sequences": sequences,
                "dialect": "alphafold3",
                "version": 1,
            }
            (out_dir / f"{af3_name}.json").write_text(json.dumps(job_json, indent=2) + "\n", encoding="utf-8")
            manifest[af3_name] = {
                "design_id": design_id,
                "parent_design_id": parent_id,
                "parent_structure_path": (row.get("parent_structure_path") or "").strip(),
            }
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare AlphaFold 3 JSON job directory.")
    parser.add_argument("--inputs", nargs="*", default=[], help="Existing AF3 JSON files or directories.")
    parser.add_argument("--sequences-tsv", default="", help="Optional sequence-native TSV/CSV table (e.g. from ligandmpnn).")
    parser.add_argument("--target-sequence", default="", help="Optional target protein sequence appended as second chain.")
    parser.add_argument("--ligand-ccd", default="", help="Optional comma-separated ligand CCD code(s) (e.g. ATP).")
    parser.add_argument("--ligand-smiles", default="", help="Optional ligand SMILES string.")
    parser.add_argument("--model-seeds", default="1", help="Comma-separated integer modelSeeds (default: '1').")
    parser.add_argument("--run-data-pipeline", action="store_true", help="Omit empty MSA/template stubs so AF3 runs MSA search.")
    parser.add_argument("--out-dir", required=True, help="Output directory for staged AF3 JSON files.")
    parser.add_argument("--out-manifest", required=True, help="Output lineage manifest JSON path.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seeds = _parse_seeds(args.model_seeds)

    manifest: dict[str, dict[str, Any]] = {}
    if args.sequences_tsv:
        manifest.update(
            _stage_from_sequences_table(
                table_path=Path(args.sequences_tsv),
                out_dir=out_dir,
                target_sequence=args.target_sequence,
                ligand_ccd=args.ligand_ccd,
                ligand_smiles=args.ligand_smiles,
                model_seeds=seeds,
                run_data_pipeline=args.run_data_pipeline,
            )
        )
    if args.inputs:
        manifest.update(_stage_existing_inputs(args.inputs, out_dir))

    if not manifest:
        raise SystemExit("No AlphaFold 3 JSON inputs or valid sequence rows were staged.")

    out_manifest = Path(args.out_manifest)
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    out_manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
