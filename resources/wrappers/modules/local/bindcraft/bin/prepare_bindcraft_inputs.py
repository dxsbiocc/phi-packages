#!/usr/bin/env python3
"""Prepare an absolute-path BindCraft target settings JSON from either:
  1. An existing BindCraft target settings JSON (--settings-json), OR
  2. Direct CLI parameters (--starting-pdb, --binder-name, --chains,
     --target-hotspot-residues, --min-length, --max-length, --num-designs).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare BindCraft settings JSON.")
    parser.add_argument("--settings-json", default="", help="Optional existing BindCraft settings JSON.")
    parser.add_argument("--starting-pdb", default="", help="Target PDB structure path.")
    parser.add_argument("--design-path", required=True, help="Output directory path written into design_path.")
    parser.add_argument("--binder-name", default="binder", help="Binder campaign prefix.")
    parser.add_argument("--chains", default="A", help="Target chain(s) to bind (e.g. 'A' or 'A,B').")
    parser.add_argument("--target-hotspot-residues", default="", help="Optional hotspot residues (e.g. 'A56,A58-62').")
    parser.add_argument("--min-length", type=int, default=65, help="Minimum binder length.")
    parser.add_argument("--max-length", type=int, default=120, help="Maximum binder length.")
    parser.add_argument("--num-designs", type=int, default=10, help="Target number_of_final_designs stopping condition.")
    parser.add_argument("--out-json", required=True, help="Output staged BindCraft settings JSON path.")
    args = parser.parse_args()

    design_path_abs = str(Path(args.design_path).resolve())

    if args.settings_json:
        src = Path(args.settings_json)
        if not src.is_file():
            raise SystemExit(f"BindCraft settings JSON not found: {src}")
        settings: dict[str, Any] = json.loads(src.read_text(encoding="utf-8"))
        settings["design_path"] = design_path_abs
        if args.starting_pdb:
            settings["starting_pdb"] = str(Path(args.starting_pdb).resolve())
        elif settings.get("starting_pdb"):
            p = Path(str(settings["starting_pdb"]))
            if not p.is_absolute():
                cand = (src.parent / p).resolve()
                settings["starting_pdb"] = str(cand if cand.exists() else p)
    else:
        if not args.starting_pdb:
            raise SystemExit("Either --settings-json or --starting-pdb must be provided.")
        if args.max_length < args.min_length:
            raise SystemExit(f"max_length ({args.max_length}) must be >= min_length ({args.min_length}).")
        settings = {
            "design_path": design_path_abs,
            "binder_name": args.binder_name,
            "starting_pdb": str(Path(args.starting_pdb).resolve()),
            "chains": args.chains,
            "target_hotspot_residues": args.target_hotspot_residues,
            "lengths": [args.min_length, args.max_length],
            "number_of_final_designs": args.num_designs,
        }

    out_path = Path(args.out_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
