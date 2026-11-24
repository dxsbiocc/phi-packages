#!/usr/bin/env python3
"""Prepare Protenix JSON input list from either:
  1. Existing Protenix JSON file(s) or directory (--inputs), OR
  2. A sequence-native TSV/CSV table (--sequences-tsv, e.g. from ligandmpnn)
     plus optional target sequence and ligand CCD/SMILES.

Pure-Python implementation (no third-party dependencies required).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any


def _sanitize_name(raw: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", raw.strip())
    return cleaned or "design"


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare Protenix JSON input file.")
    parser.add_argument("--inputs", nargs="*", default=[], help="Existing Protenix JSON file(s) or directory.")
    parser.add_argument("--sequences-tsv", default="", help="Optional sequence-native TSV/CSV table.")
    parser.add_argument("--target-sequence", default="", help="Optional target protein sequence appended as a chain.")
    parser.add_argument("--ligand-ccd", default="", help="Optional ligand CCD code (formatted as CCD_<code>).")
    parser.add_argument("--ligand-smiles", default="", help="Optional ligand SMILES string.")
    parser.add_argument("--out-json", required=True, help="Output staged Protenix JSON list path.")
    parser.add_argument("--out-manifest", required=True, help="Output lineage manifest JSON path.")
    args = parser.parse_args()

    jobs: list[dict[str, Any]] = []
    manifest: dict[str, dict[str, Any]] = {}

    if args.sequences_tsv:
        tpath = Path(args.sequences_tsv)
        delim = "\t" if tpath.suffix.lower() in (".tsv", ".tab") else ","
        with tpath.open("r", encoding="utf-8", newline="") as fh:
            for idx, row in enumerate(csv.DictReader(fh, delimiter=delim), start=1):
                design_id = (row.get("design_id") or f"design_{idx}").strip()
                parent_id = (row.get("parent_design_id") or design_id).strip()
                raw_seq = (row.get("sequence") or "").strip()
                if not raw_seq:
                    continue

                parts = [p.strip() for p in raw_seq.replace("/", ":").split(":") if p.strip()]
                if args.target_sequence.strip():
                    parts.append(args.target_sequence.strip())

                sequences: list[dict[str, Any]] = [
                    {"proteinChain": {"sequence": seq, "count": 1}} for seq in parts
                ]
                if args.ligand_ccd.strip():
                    ccd = args.ligand_ccd.strip()
                    lig_val = ccd if ccd.startswith("CCD_") else f"CCD_{ccd}"
                    sequences.append({"ligand": {"ligand": lig_val, "count": 1}})
                elif args.ligand_smiles.strip():
                    sequences.append({"ligand": {"ligand": args.ligand_smiles.strip(), "count": 1}})

                job_name = _sanitize_name(design_id)
                jobs.append({"name": job_name, "sequences": sequences})
                manifest[job_name] = {
                    "design_id": design_id,
                    "parent_design_id": parent_id,
                    "parent_structure_path": (row.get("parent_structure_path") or "").strip(),
                }

    for raw in args.inputs:
        p = Path(raw)
        candidates = (
            [c for c in sorted(p.rglob("*.json")) if c.is_file()]
            if p.is_dir()
            else ([p] if p.is_file() and p.suffix.lower() == ".json" else [])
        )
        for src in candidates:
            loaded = json.loads(src.read_text(encoding="utf-8"))
            entries = loaded if isinstance(loaded, list) else [loaded]
            for entry in entries:
                if not isinstance(entry, dict):
                    continue
                jname = _sanitize_name(str(entry.get("name") or src.stem))
                entry_copy = dict(entry)
                entry_copy["name"] = jname
                jobs.append(entry_copy)
                manifest[jname] = {
                    "design_id": str(entry.get("name") or src.stem),
                    "parent_design_id": str(entry.get("name") or src.stem),
                    "parent_structure_path": "",
                }

    if not jobs:
        raise SystemExit("No Protenix jobs staged from --inputs or --sequences-tsv.")

    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(jobs, indent=2) + "\n", encoding="utf-8")

    out_man = Path(args.out_manifest)
    out_man.parent.mkdir(parents=True, exist_ok=True)
    out_man.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
