#!/usr/bin/env python3
"""Prepare a clean FASTA file and lineage manifest for ESMFold2 from either:
  1. Existing FASTA file(s) (--fasta), OR
  2. A sequence-native TSV/CSV table (--sequences-tsv, e.g. from ligandmpnn).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any


def _sanitize_id(raw: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", raw.strip())
    return cleaned or "design"


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:].strip().split()[0]
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare staged FASTA input for ESMFold2.")
    parser.add_argument("--fasta", nargs="*", default=[], help="Optional input FASTA file(s).")
    parser.add_argument("--sequences-tsv", default="", help="Optional sequence-native TSV/CSV table.")
    parser.add_argument("--out-fasta", required=True, help="Output staged FASTA path.")
    parser.add_argument("--out-manifest", required=True, help="Output lineage manifest JSON path.")
    args = parser.parse_args()

    fasta_lines: list[str] = []
    manifest: dict[str, dict[str, Any]] = {}

    if args.sequences_tsv:
        tpath = Path(args.sequences_tsv)
        delim = "\t" if tpath.suffix.lower() in (".tsv", ".tab") else ","
        with tpath.open("r", encoding="utf-8", newline="") as fh:
            for idx, row in enumerate(csv.DictReader(fh, delimiter=delim), start=1):
                raw_id = (row.get("design_id") or f"design_{idx}").strip()
                parent_id = (row.get("parent_design_id") or raw_id).strip()
                seq = (row.get("sequence") or "").strip().replace("/", ":")
                if not seq:
                    continue
                safe_id = _sanitize_id(raw_id)
                fasta_lines.append(f">{safe_id}")
                fasta_lines.append(seq)
                manifest[safe_id] = {
                    "design_id": raw_id,
                    "parent_design_id": parent_id,
                    "parent_structure_path": (row.get("parent_structure_path") or "").strip(),
                    "sequence": seq,
                }

    for raw_fa in args.fasta:
        fpath = Path(raw_fa)
        if not fpath.is_file():
            continue
        for header, seq in _parse_fasta(fpath):
            if not seq:
                continue
            safe_id = _sanitize_id(header)
            fasta_lines.append(f">{safe_id}")
            fasta_lines.append(seq)
            manifest[safe_id] = {
                "design_id": header,
                "parent_design_id": header,
                "parent_structure_path": "",
                "sequence": seq,
            }

    if not manifest:
        raise SystemExit("No valid sequences staged for ESMFold2.")

    out_fa = Path(args.out_fasta)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    out_fa.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")

    out_man = Path(args.out_manifest)
    out_man.parent.mkdir(parents=True, exist_ok=True)
    out_man.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
