#!/usr/bin/env python3
"""Prepare and validate an RFdiffusion3 (RFD3) JSON input specification.

Supports two modes:
  1. Existing JSON/YAML config file (--inputs-json): validates keys, resolves
     relative structure paths, and injects a chunk-specific prefix/seed.
  2. Quick CLI mode (--contig): builds a single-entry RFD3 JSON specification
     for binder/scaffold generation without hand-authoring JSON.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare RFD3 input specification JSON.")
    parser.add_argument("--inputs-json", default="", help="Optional existing RFD3 input JSON file.")
    parser.add_argument("--design-name", default="design", help="Design key name when using --contig mode.")
    parser.add_argument("--contig", default="", help="RFD3 contig string (e.g. '50-80' or 'A1-100/0 60-90').")
    parser.add_argument("--input-structure", default="", help="Optional target/motif PDB or CIF path.")
    parser.add_argument("--ligand", default="", help="Optional ligand residue name(s), comma-separated.")
    parser.add_argument("--chunk-id", type=int, default=0, help="Chunk index for deterministic seed/prefix offset.")
    parser.add_argument("--out-json", required=True, help="Output staged RFD3 JSON file.")
    args = parser.parse_args()

    spec: dict[str, Any]
    if args.inputs_json:
        src = Path(args.inputs_json)
        if not src.is_file():
            raise SystemExit(f"RFD3 inputs file not found: {src}")
        loaded = json.loads(src.read_text(encoding="utf-8"))
        if not isinstance(loaded, dict) or not loaded:
            raise SystemExit("RFD3 inputs JSON must be a non-empty JSON object mapping design names to specs.")
        spec = {}
        for key, entry in loaded.items():
            if not isinstance(entry, dict):
                raise SystemExit(f"Invalid RFD3 entry for key {key!r}: expected object.")
            staged_entry = dict(entry)
            if args.input_structure and "input" not in staged_entry:
                staged_entry["input"] = str(Path(args.input_structure).resolve())
            elif "input" in staged_entry and staged_entry["input"]:
                p = Path(str(staged_entry["input"]))
                if not p.is_absolute():
                    candidate = (src.parent / p).resolve()
                    staged_entry["input"] = str(candidate if candidate.exists() else p)
            staged_key = f"{key}_chunk{args.chunk_id}" if args.chunk_id > 0 else key
            spec[staged_key] = staged_entry
    else:
        if not args.contig:
            raise SystemExit("Either --inputs-json or --contig must be provided.")
        entry: dict[str, Any] = {"contig": args.contig}
        if args.input_structure:
            entry["input"] = str(Path(args.input_structure).resolve())
        if args.ligand:
            entry["ligand"] = args.ligand
        key = f"{args.design_name}_chunk{args.chunk_id}" if args.chunk_id > 0 else args.design_name
        spec = {key: entry}

    out_path = Path(args.out_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
