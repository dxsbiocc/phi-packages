#!/usr/bin/env python3
"""Prepare batch inputs (pdb_path_multi.json + lineage manifest) for LigandMPNN."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

STRUCTURE_EXTS = {".pdb", ".cif", ".mmcif"}


def collect_structures(raw_inputs: list[str]) -> list[Path]:
    found: list[Path] = []
    for raw in raw_inputs:
        p = Path(raw)
        if p.is_dir():
            for child in sorted(p.rglob("*")):
                if child.is_file() and child.suffix.lower() in STRUCTURE_EXTS:
                    found.append(child.resolve())
        elif p.is_file() and p.suffix.lower() in STRUCTURE_EXTS:
            found.append(p.resolve())
    return found


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare LigandMPNN multi-PDB JSON input.")
    parser.add_argument("--structures", nargs="+", required=True, help="Input PDB/CIF files or directories.")
    parser.add_argument("--out-multi-json", required=True, help="Output pdb_path_multi.json path.")
    parser.add_argument("--out-manifest", required=True, help="Output parent lineage manifest JSON path.")
    args = parser.parse_args()

    structures = collect_structures(args.structures)
    if not structures:
        raise SystemExit("No PDB or mmCIF structure files found for LigandMPNN.")

    multi_map: dict[str, str] = {}
    manifest: dict[str, str] = {}
    for sp in structures:
        multi_map[str(sp)] = ""
        manifest[sp.stem] = str(sp)

    out_multi = Path(args.out_multi_json)
    out_multi.parent.mkdir(parents=True, exist_ok=True)
    out_multi.write_text(json.dumps(multi_map, indent=2) + "\n", encoding="utf-8")

    out_manifest = Path(args.out_manifest)
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    out_manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
