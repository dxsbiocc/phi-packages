#!/usr/bin/env python3
"""Prepare Chai-1 entity-tagged FASTA files (>protein|name=..., >ligand|name=...)
from either:
  1. Existing Chai-1 FASTA file(s) or directory (--inputs), OR
  2. A sequence-native TSV/CSV table (--sequences-tsv, e.g. from ligandmpnn)
     plus optional target sequence and ligand SMILES string.

Pure-Python implementation (no third-party dependencies required).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from pathlib import Path
from typing import Any

FASTA_EXTS = {".fa", ".fasta", ".faa"}


def _sanitize_name(raw: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", raw.strip())
    return cleaned or "design"


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare Chai-1 entity-tagged FASTA directory.")
    parser.add_argument("--inputs", nargs="*", default=[], help="Existing Chai-1 FASTA file(s) or directory.")
    parser.add_argument("--sequences-tsv", default="", help="Optional sequence-native TSV/CSV table.")
    parser.add_argument("--target-sequence", default="", help="Optional target protein sequence appended as target chain.")
    parser.add_argument("--ligand-smiles", default="", help="Optional ligand SMILES string for >ligand|name=ligand.")
    parser.add_argument("--out-dir", required=True, help="Output directory for staged Chai-1 FASTA files.")
    parser.add_argument("--out-manifest", required=True, help="Output lineage manifest JSON path.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
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
                lines: list[str] = []
                for c_idx, seq in enumerate(parts, start=1):
                    lines.append(f">protein|name=design_chain_{c_idx}")
                    lines.append(seq)
                if args.target_sequence.strip():
                    lines.append(">protein|name=target")
                    lines.append(args.target_sequence.strip())
                if args.ligand_smiles.strip():
                    lines.append(">ligand|name=ligand")
                    lines.append(args.ligand_smiles.strip())

                safe_name = _sanitize_name(design_id)
                (out_dir / f"{safe_name}.fasta").write_text("\n".join(lines) + "\n", encoding="utf-8")
                manifest[safe_name] = {
                    "design_id": design_id,
                    "parent_design_id": parent_id,
                    "parent_structure_path": (row.get("parent_structure_path") or "").strip(),
                }

    for raw in args.inputs:
        p = Path(raw)
        candidates = (
            [c for c in sorted(p.rglob("*")) if c.is_file() and c.suffix.lower() in FASTA_EXTS]
            if p.is_dir()
            else ([p] if p.is_file() and p.suffix.lower() in FASTA_EXTS else [])
        )
        for src in candidates:
            safe_name = _sanitize_name(src.stem)
            shutil.copy2(src, out_dir / f"{safe_name}.fasta")
            manifest[safe_name] = {
                "design_id": src.stem,
                "parent_design_id": src.stem,
                "parent_structure_path": "",
            }

    if not manifest:
        raise SystemExit("No Chai-1 FASTA jobs staged from --inputs or --sequences-tsv.")

    out_man = Path(args.out_manifest)
    out_man.parent.mkdir(parents=True, exist_ok=True)
    out_man.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
