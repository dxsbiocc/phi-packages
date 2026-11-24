#!/usr/bin/env python3
"""Inspect a PDB or mmCIF structure file and emit a structured JSON summary of
protein chains, residue ranges, numbering gaps, and non-water HETATM ligands.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O"}
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def _tokenize_cif_line(line: str) -> list[str]:
    tokens: list[str] = []
    i = 0
    n = len(line)
    while i < n:
        while i < n and line[i].isspace():
            i += 1
        if i >= n or line[i] == "#":
            break
        if line[i] in ("'", '"'):
            quote = line[i]
            i += 1
            start = i
            while i < n and not (line[i] == quote and (i + 1 == n or line[i + 1].isspace())):
                i += 1
            tokens.append(line[start:i])
            i += 1
        else:
            start = i
            while i < n and not line[i].isspace():
                i += 1
            tokens.append(line[start:i])
    return tokens


def _parse_atoms(path: Path) -> list[tuple[str, str, str, str]]:
    """Return list of (record, chain, resseq_str, resname)."""
    suffix = path.suffix.lower()
    records: list[tuple[str, str, str, str]] = []
    if suffix == ".pdb":
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM") or len(line) < 27:
                continue
            alt = line[16].strip()
            if alt not in ("", "A", "1"):
                continue
            resname = line[17:20].strip()
            chain = line[21].strip() or "A"
            resseq = line[22:26].strip()
            records.append((rec, chain, resseq, resname))
        return records
    if suffix in (".cif", ".mmcif"):
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        i = 0
        n = len(lines)
        while i < n:
            if lines[i].strip() == "loop_":
                i += 1
                cols: list[str] = []
                while i < n and lines[i].strip().startswith("_"):
                    cols.append(lines[i].strip())
                    i += 1
                if cols and cols[0].startswith("_atom_site."):
                    idx = {c: k for k, c in enumerate(cols)}
                    g_i = idx.get("_atom_site.group_PDB")
                    c_i = idx.get("_atom_site.auth_asym_id", idx.get("_atom_site.label_asym_id"))
                    s_i = idx.get("_atom_site.auth_seq_id", idx.get("_atom_site.label_seq_id"))
                    r_i = idx.get("_atom_site.label_comp_id", idx.get("_atom_site.auth_comp_id"))
                    a_i = idx.get("_atom_site.label_alt_id")
                    m_i = idx.get("_atom_site.pdbx_PDB_model_num")
                    first_model: str | None = None
                    while i < n:
                        s = lines[i].strip()
                        if not s or s.startswith("#"):
                            i += 1
                            continue
                        if s == "loop_" or s.startswith("_") or s.startswith("data_"):
                            break
                        toks = _tokenize_cif_line(s)
                        i += 1
                        if len(toks) < len(cols):
                            continue
                        if m_i is not None:
                            if first_model is None:
                                first_model = toks[m_i]
                            elif toks[m_i] != first_model:
                                continue
                        if a_i is not None and toks[a_i] not in (".", "?", "A", "1"):
                            continue
                        rec = toks[g_i] if g_i is not None else "ATOM"
                        if rec not in ("ATOM", "HETATM"):
                            continue
                        chain = toks[c_i] if c_i is not None else "A"
                        resseq = toks[s_i] if s_i is not None else "0"
                        resname = toks[r_i] if r_i is not None else "UNK"
                        records.append((rec, chain, resseq, resname))
                    break
            else:
                i += 1
        return records
    raise ValueError(f"Unsupported structure format: {path.suffix}")


def inspect_structure(path: Path) -> dict[str, Any]:
    raw_records = _parse_atoms(path)
    chain_residues: dict[str, dict[int, str]] = {}
    ligands: dict[tuple[str, str, str], int] = {}

    for rec, chain, resseq, resname in raw_records:
        if rec == "ATOM":
            try:
                num = int(resseq)
            except ValueError:
                continue
            chain_residues.setdefault(chain, {}).setdefault(num, resname)
        elif rec == "HETATM":
            if resname in WATER_RESNAMES:
                continue
            key = (chain, resname, resseq)
            ligands[key] = ligands.get(key, 0) + 1

    chains_out: list[dict[str, Any]] = []
    for cid in sorted(chain_residues.keys()):
        res_map = chain_residues[cid]
        nums = sorted(res_map.keys())
        gaps: list[str] = []
        for a, b in zip(nums, nums[1:]):
            if b > a + 1:
                gaps.append(f"{a + 1}-{b - 1}")
        seq = "".join(AA3_TO_1.get(res_map[n], "X") for n in nums)
        chains_out.append(
            {
                "chain_id": cid,
                "num_residues": len(nums),
                "min_resseq": nums[0],
                "max_resseq": nums[-1],
                "gaps": gaps,
                "suggested_contig": f"{cid}{nums[0]}-{nums[-1]}",
                "sequence": seq,
            }
        )

    ligands_out = [
        {"chain_id": c, "resname": rn, "resseq": rs, "atom_count": cnt}
        for (c, rn, rs), cnt in sorted(ligands.items())
    ]

    return {
        "path": str(path.resolve()),
        "format": path.suffix.lower().lstrip("."),
        "num_chains": len(chains_out),
        "chains": chains_out,
        "ligands": ligands_out,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Inspect PDB/mmCIF structure chains, residue ranges, and ligands.")
    parser.add_argument("--input", required=True, help="Input PDB or mmCIF structure file.")
    args = parser.parse_args()

    path = Path(args.input)
    if not path.is_file():
        raise SystemExit(f"Structure file not found: {path}")

    result = inspect_structure(path)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
