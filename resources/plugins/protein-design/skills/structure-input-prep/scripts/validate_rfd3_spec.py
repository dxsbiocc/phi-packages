#!/usr/bin/env python3
"""Validate an RFdiffusion3 JSON specification or contig string, optionally
cross-checking chains, residue bounds, and ligands against an input structure.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from inspect_structure import inspect_structure

FIXED_SEG_RE = re.compile(r"^([A-Za-z])(\d+)(?:-(\d+))?$")
DESIGN_SEG_RE = re.compile(r"^(\d+)(?:-(\d+))?$")


def _validate_contig(
    contig: str,
    struct_info: dict[str, Any] | None,
) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    if not contig or not contig.strip():
        errors.append("Contig string is empty.")
        return errors, warnings

    chain_bounds: dict[str, tuple[int, int]] = {}
    if struct_info:
        for ch in struct_info.get("chains", []):
            chain_bounds[ch["chain_id"]] = (ch["min_resseq"], ch["max_resseq"])

    raw_tokens = [t for part in contig.strip().split() for t in part.split("/") if t != ""]
    if not raw_tokens:
        errors.append("Contig contains no segments.")
        return errors, warnings

    has_design_region = False
    for tok in raw_tokens:
        if tok == "0":
            continue
        m_fixed = FIXED_SEG_RE.match(tok)
        if m_fixed:
            cid = m_fixed.group(1)
            start = int(m_fixed.group(2))
            end = int(m_fixed.group(3)) if m_fixed.group(3) else start
            if end < start:
                errors.append(f"Fixed segment {tok!r} has end < start.")
            if chain_bounds:
                if cid not in chain_bounds:
                    errors.append(f"Chain {cid!r} in segment {tok!r} not found in structure (available: {sorted(chain_bounds)}).")
                else:
                    cmin, cmax = chain_bounds[cid]
                    if start < cmin or end > cmax:
                        warnings.append(
                            f"Segment {tok!r} extends beyond observed chain {cid} range ({cmin}-{cmax})."
                        )
            continue
        m_des = DESIGN_SEG_RE.match(tok)
        if m_des:
            has_design_region = True
            lo = int(m_des.group(1))
            hi = int(m_des.group(2)) if m_des.group(2) else lo
            if lo <= 0 or hi < lo:
                errors.append(f"Invalid de novo length segment {tok!r}.")
            continue
        errors.append(f"Unrecognized contig token: {tok!r}.")

    if not has_design_region:
        warnings.append("Contig contains no de novo design length segment (all segments appear fixed).")
    return errors, warnings


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate RFdiffusion3 JSON/contig specification.")
    parser.add_argument("--inputs-json", "--inputs_json", dest="inputs_json", default="", help="Optional RFD3 JSON file.")
    parser.add_argument("--contig", default="", help="Optional contig string to validate directly.")
    parser.add_argument("--input-structure", "--input_structure", dest="input_structure", default="", help="Optional target PDB/mmCIF.")
    parser.add_argument("--ligand", default="", help="Optional comma-separated ligand residue code(s).")
    args = parser.parse_args()

    errors: list[str] = []
    warnings: list[str] = []
    entries_checked: list[str] = []

    global_struct_info: dict[str, Any] | None = None
    if args.input_structure:
        sp = Path(args.input_structure)
        if not sp.is_file():
            errors.append(f"input_structure not found: {sp}")
        else:
            global_struct_info = inspect_structure(sp)

    def _check_ligands(lig_str: str, s_info: dict[str, Any] | None, label: str) -> None:
        if not lig_str or not s_info:
            return
        present = {lig["resname"] for lig in s_info.get("ligands", [])}
        for code in [c.strip() for c in lig_str.split(",") if c.strip()]:
            if code not in present:
                errors.append(f"[{label}] Ligand {code!r} not found in structure HETATMs (present: {sorted(present)}).")

    if args.inputs_json:
        jp = Path(args.inputs_json)
        if not jp.is_file():
            errors.append(f"inputs_json not found: {jp}")
        else:
            try:
                data = json.loads(jp.read_text(encoding="utf-8"))
            except json.JSONDecodeError as exc:
                errors.append(f"Invalid JSON in {jp}: {exc}")
                data = None
            if isinstance(data, dict) and data:
                for key, entry in data.items():
                    entries_checked.append(key)
                    if not isinstance(entry, dict):
                        errors.append(f"[{key}] Entry must be a JSON object.")
                        continue
                    s_info = global_struct_info
                    if entry.get("input"):
                        ep = Path(str(entry["input"]))
                        if not ep.is_absolute():
                            ep = (jp.parent / ep).resolve()
                        if not ep.is_file():
                            errors.append(f"[{key}] Referenced input structure not found: {entry['input']}")
                        else:
                            s_info = inspect_structure(ep)
                    e_err, e_warn = _validate_contig(str(entry.get("contig", "")), s_info)
                    errors.extend(f"[{key}] {m}" for m in e_err)
                    warnings.extend(f"[{key}] {m}" for m in e_warn)
                    _check_ligands(str(entry.get("ligand", "") or args.ligand), s_info, key)
            elif data is not None:
                errors.append("inputs_json must be a non-empty JSON object.")
    elif args.contig:
        entries_checked.append("cli_contig")
        e_err, e_warn = _validate_contig(args.contig, global_struct_info)
        errors.extend(e_err)
        warnings.extend(e_warn)
        _check_ligands(args.ligand, global_struct_info, "cli_contig")
    else:
        errors.append("Provide either --inputs_json or --contig to validate.")

    result = {
        "valid": len(errors) == 0,
        "entries_checked": entries_checked,
        "errors": errors,
        "warnings": warnings,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
