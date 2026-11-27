"""Shared dual-engine cheminformatics, I/O, plotting, and artifact utilities (chem_*).

Primary engine uses RDKit (`rdkit`) when installed in `plugin:chem`.
Fallback engine provides a deterministic pure-Python SMILES molecular graph
parser, descriptor calculator, 1024-bit circular fingerprint generator,
Bemis-Murcko scaffold extractor, and 15+ PAINS/Brenk toxicophore alerts when
RDKit is absent in minimal CI environments.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import struct
import sys
import zlib
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    np = None  # type: ignore[assignment]
    HAS_NUMPY = False

try:
    import pandas as pd

    HAS_PANDAS = True
except ImportError:
    pd = None  # type: ignore[assignment]
    HAS_PANDAS = False

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: F401

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import Crippen, Descriptors, FilterCatalog, QED, rdFingerprintGenerator, rdMolDescriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.Scaffolds import MurckoScaffold

    HAS_RDKIT = True
except ImportError:
    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem, Crippen, Descriptors, FilterCatalog, QED, rdMolDescriptors
        from rdkit.Chem.MolStandardize import rdMolStandardize
        from rdkit.Chem.Scaffolds import MurckoScaffold

        rdFingerprintGenerator = None  # type: ignore[assignment]
        HAS_RDKIT = True
    except ImportError:
        Chem = None  # type: ignore[assignment]
        DataStructs = None  # type: ignore[assignment]
        HAS_RDKIT = False

BACKEND = "rdkit" if HAS_RDKIT else "pure_python_smiles"

ATOMIC_WEIGHTS: dict[str, float] = {
    "H": 1.008,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Na": 22.990,
    "Mg": 24.305,
    "Si": 28.085,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    "K": 39.098,
    "Ca": 40.078,
    "Se": 78.971,
    "Br": 79.904,
    "I": 126.904,
}

ATOMIC_NUMS: dict[str, int] = {
    "H": 1,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Na": 11,
    "Mg": 12,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "K": 19,
    "Ca": 20,
    "Se": 34,
    "Br": 35,
    "I": 53,
}

ORGANIC_VALENCES: dict[str, tuple[int, ...]] = {
    "B": (3,),
    "C": (4,),
    "N": (3,),
    "O": (2,),
    "P": (3, 5),
    "S": (2, 4, 6),
    "F": (1,),
    "Cl": (1,),
    "Br": (1,),
    "I": (1,),
    "Si": (4,),
    "Se": (2, 4, 6),
    "H": (1,),
}

# Bickerton et al. (2012) Nature Chemistry QED asymmetric double sigmoidal parameters
_QED_PARAMS: dict[str, tuple[float, float, float, float, float, float, float]] = {
    # prop: (A, B, C, D, E, F, d_max)
    "MW": (2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707, 104.9805561),
    "ALOGP": (3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154, 0.576295591, 131.3186604),
    "HBA": (2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958, 148.7763046),
    "HBD": (1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843, 0.920922555, 258.1632616),
    "PSA": (1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732, 104.5686167),
    "ROTB": (0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166, 2.758063707, 105.4420403),
    "AROM": (3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384, 0.375760881, 312.3372610),
    "ALERTS": (0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782, 417.7253140),
}
_QED_WEIGHTS: dict[str, float] = {
    "MW": 0.66,
    "ALOGP": 0.46,
    "HBA": 0.05,
    "HBD": 0.61,
    "PSA": 0.06,
    "ROTB": 0.65,
    "AROM": 0.48,
    "ALERTS": 0.95,
}


def emit(payload: dict[str, Any]) -> None:
    """Write JSON result payload to stdout and exit 0."""
    sys.stdout.write(json.dumps(clean_json(payload), ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()
    raise SystemExit(0)


def fail(message: str, code: int = 1) -> None:
    """Write JSON error payload to stdout and exit non-zero."""
    sys.stdout.write(json.dumps({"error": str(message)}, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()
    raise SystemExit(code)


def clean_json(obj: Any) -> Any:
    """Recursively convert numpy types and non-finite floats into strict JSON values."""
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if HAS_NUMPY:
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            val = float(obj)
            return val if math.isfinite(val) else None
    if isinstance(obj, int):
        return int(obj)
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {str(k): clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean_json(v) for v in obj]
    return str(obj)


def _detect_sep(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return ","
    if name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".txt"):
        return "\t"
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("#") and line.strip():
                return "\t" if "\t" in line else ","
    return "\t"


def read_table_rows(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read TSV/CSV table into (columns, list_of_row_dicts), using pandas if installed."""
    p = Path(path)
    if not p.is_file():
        fail(f"Input file not found: {path}")
    sep = _detect_sep(p)
    if HAS_PANDAS:
        try:
            df = pd.read_csv(p, sep=sep, comment="#", dtype=str, keep_default_na=False, low_memory=False)
        except Exception as exc:
            fail(f"Failed to read table '{path}': {exc}")
        if df.empty:
            fail(f"Input table is empty: {path}")
        cols = [str(c).strip() for c in df.columns]
        df.columns = cols
        rows = [{c: str(v) for c, v in rec.items()} for rec in df.to_dict(orient="records")]
        return cols, rows

    with open(p, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        filtered = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    if not filtered:
        fail(f"Input table is empty: {path}")
    reader = csv.DictReader(filtered, delimiter=sep)
    cols = [str(c).strip() for c in (reader.fieldnames or [])]
    if not cols:
        fail(f"Input table has no header columns: {path}")
    rows: list[dict[str, str]] = []
    for raw in reader:
        rows.append({str(k).strip(): ("" if v is None else str(v)) for k, v in raw.items() if k is not None})
    if not rows:
        fail(f"Input table has zero data rows: {path}")
    return cols, rows


def resolve_column_name(
    columns: list[str],
    candidates: list[str],
    label: str,
    required: bool = True,
) -> str | None:
    """Resolve column name case-insensitively against candidate names."""
    lower_map = {c.strip().lower(): c.strip() for c in columns}
    for cand in candidates:
        if cand in columns:
            return cand
        hit = lower_map.get(cand.strip().lower())
        if hit is not None:
            return hit
    if required:
        fail(f"Missing required {label} column (expected one of {candidates}; found {columns})")
    return None


def _parse_sdf_records(path: Path) -> list[dict[str, str]]:
    """Parse SDF records into compound_id + SMILES (via RDKit or V2000 fallback)."""
    text = path.read_text(encoding="utf-8", errors="replace")
    blocks = [b for b in text.split("$$$$") if b.strip()]
    if not blocks:
        fail(f"SDF file contains no molecule records: {path}")
    records: list[dict[str, str]] = []
    for idx, block in enumerate(blocks, start=1):
        lines = block.splitlines()
        title = lines[0].strip() if lines else ""
        # Check SDF property tags for SMILES or ID
        smi_tag = None
        id_tag = None
        for i, ln in enumerate(lines):
            if ln.startswith(">") and ("SMILES" in ln.upper()):
                if i + 1 < len(lines) and lines[i + 1].strip():
                    smi_tag = lines[i + 1].strip()
            if ln.startswith(">") and any(k in ln.upper() for k in ("<ID>", "<NAME>", "<COMPOUND_ID>", "<TITLE>")):
                if i + 1 < len(lines) and lines[i + 1].strip():
                    id_tag = lines[i + 1].strip()
        cid = id_tag or title or f"CMPD_{idx:04d}"
        smi = smi_tag
        if not smi and HAS_RDKIT:
            mol = Chem.MolFromMolBlock(block, sanitize=True, removeHs=True)
            if mol is not None:
                smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        if not smi:
            smi = _molblock_v2000_to_smiles(lines)
        records.append({"compound_id": cid, "smiles": smi or ""})
    return records


def _molblock_v2000_to_smiles(lines: list[str]) -> str:
    """Minimal V2000 molblock to SMILES converter for pure-Python SDF fallback."""
    counts_idx = None
    for i, ln in enumerate(lines[:6]):
        if "V2000" in ln:
            counts_idx = i
            break
    if counts_idx is None:
        return ""
    try:
        n_atoms = int(lines[counts_idx][:3].strip())
        n_bonds = int(lines[counts_idx][3:6].strip())
    except ValueError:
        return ""
    atom_lines = lines[counts_idx + 1 : counts_idx + 1 + n_atoms]
    bond_lines = lines[counts_idx + 1 + n_atoms : counts_idx + 1 + n_atoms + n_bonds]
    if len(atom_lines) != n_atoms or len(bond_lines) != n_bonds or n_atoms == 0:
        return ""
    syms: list[str] = []
    for al in atom_lines:
        parts = al.split()
        if len(parts) < 4:
            return ""
        syms.append(parts[3])
    adj: list[list[tuple[int, int]]] = [[] for _ in range(n_atoms)]
    for bl in bond_lines:
        try:
            a1 = int(bl[:3].strip()) - 1
            a2 = int(bl[3:6].strip()) - 1
            btype = int(bl[6:9].strip())
        except ValueError:
            return ""
        if not (0 <= a1 < n_atoms and 0 <= a2 < n_atoms):
            return ""
        adj[a1].append((a2, btype))
        adj[a2].append((a1, btype))
    # Build SMILES via DFS over V2000 graph
    visited = [False] * n_atoms
    ring_nums: dict[tuple[int, int], int] = {}
    next_ring = [1]

    def atom_tok(idx: int) -> str:
        s = syms[idx]
        if s in ("B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"):
            return s
        return f"[{s}]"

    def bond_sym(bt: int) -> str:
        if bt == 2:
            return "="
        if bt == 3:
            return "#"
        if bt == 4:
            return ":"
        return ""

    out_parts: list[str] = []
    for start in range(n_atoms):
        if visited[start]:
            continue
        if out_parts:
            out_parts.append(".")
        # DFS from start
        visited[start] = True
        parent_of: dict[int, int] = {start: -1}

        def dfs(u: int) -> str:
            tok = atom_tok(u)
            closures: list[str] = []
            children: list[tuple[int, int]] = []
            for v, bt in adj[u]:
                if v == parent_of.get(u, -1):
                    continue
                edge = (min(u, v), max(u, v))
                if not visited[v]:
                    visited[v] = True
                    parent_of[v] = u
                    children.append((v, bt))
                else:
                    if edge not in ring_nums:
                        r = next_ring[0]
                        next_ring[0] += 1
                        ring_nums[edge] = r
                    r = ring_nums[edge]
                    rs = f"%{r}" if r >= 10 else str(r)
                    closures.append(f"{bond_sym(bt)}{rs}")
            s_out = tok + "".join(closures)
            if len(children) == 1:
                v, bt = children[0]
                s_out += bond_sym(bt) + dfs(v)
            elif len(children) > 1:
                for v, bt in children[:-1]:
                    s_out += f"({bond_sym(bt)}{dfs(v)})"
                v, bt = children[-1]
                s_out += bond_sym(bt) + dfs(v)
            return s_out

        out_parts.append(dfs(start))
    return "".join(out_parts)


def read_compounds(
    path: str | Path,
    smiles_col: str | None = None,
    id_col: str | None = None,
) -> list[dict[str, str]]:
    """Read compound records (compound_id, smiles, plus original fields) from TSV/CSV/SMI/SDF."""
    p = Path(path)
    if not p.is_file():
        fail(f"Compound input file not found: {path}")
    low = p.name.lower()
    if low.endswith(".sdf"):
        return _parse_sdf_records(p)
    if low.endswith(".smi"):
        lines = [ln.strip() for ln in p.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
        if not lines:
            fail(f"SMI file is empty: {path}")
        records: list[dict[str, str]] = []
        for idx, ln in enumerate(lines, start=1):
            parts = ln.split()
            smi = parts[0]
            cid = parts[1] if len(parts) > 1 else f"CMPD_{idx:04d}"
            records.append({"compound_id": cid, "smiles": smi})
        return records

    cols, rows = read_table_rows(p)
    s_col = (
        resolve_column_name(cols, [smiles_col], "SMILES")
        if smiles_col
        else resolve_column_name(
            cols,
            ["smiles", "canonical_smiles", "smi", "structure", "mol_smiles", "ligand_smiles"],
            "SMILES",
        )
    )
    i_col = (
        resolve_column_name(cols, [id_col], "compound ID")
        if id_col
        else resolve_column_name(
            cols,
            ["compound_id", "id", "ligand_id", "name", "molecule_id", "mol_id", "title"],
            "compound ID",
            required=False,
        )
    )
    out: list[dict[str, str]] = []
    for idx, r in enumerate(rows, start=1):
        cid = (r.get(i_col, "").strip() if i_col else "") or f"CMPD_{idx:04d}"
        smi = r.get(s_col, "").strip()
        rec = dict(r)
        rec["compound_id"] = cid
        rec["smiles"] = smi
        out.append(rec)
    return out


# ==============================================================================
# Pure-Python SMILES Molecular Graph Engine (Fallback when RDKit is absent)
# ==============================================================================

@dataclass
class PyAtom:
    idx: int
    symbol: str
    is_aromatic: bool = False
    formal_charge: int = 0
    explicit_h: int | None = None
    chiral: bool = False
    implicit_h: int = 0
    in_ring: bool = False


@dataclass
class PyMolGraph:
    atoms: list[PyAtom] = field(default_factory=list)
    bonds: list[tuple[int, int, float]] = field(default_factory=list)  # (u, v, order: 1.0, 2.0, 3.0, 1.5)
    adj: dict[int, list[tuple[int, float]]] = field(default_factory=dict)
    has_ez_stereo: bool = False
    has_chiral_stereo: bool = False


def _parse_bracket_atom(token: str) -> tuple[str, bool, int, int | None, bool]:
    """Parse `[...]` SMILES atom token -> (symbol, is_aromatic, charge, explicit_h, chiral)."""
    inner = token[1:-1]
    if not inner:
        raise ValueError("Empty bracket atom '[]'")
    pos = 0
    n = len(inner)
    while pos < n and inner[pos].isdigit():
        pos += 1
    if pos >= n:
        raise ValueError(f"Missing element in bracket atom '{token}'")
    # Element symbol (1 or 2 chars)
    c0 = inner[pos]
    is_arom = c0.islower()
    if c0.lower() in ("b", "c", "n", "o", "p", "s"):
        if pos + 1 < n and inner[pos + 1].islower() and (c0 + inner[pos + 1]) in ATOMIC_NUMS:
            sym = c0 + inner[pos + 1]
            is_arom = False
            pos += 2
        else:
            sym = c0.upper()
            pos += 1
    else:
        if pos + 1 < n and inner[pos + 1].islower() and (c0 + inner[pos + 1]) in ATOMIC_NUMS:
            sym = c0 + inner[pos + 1]
            pos += 2
        elif c0 in ATOMIC_NUMS:
            sym = c0
            pos += 1
        else:
            raise ValueError(f"Unsupported element in bracket atom '{token}'")
    if sym not in ATOMIC_NUMS:
        raise ValueError(f"Unknown element '{sym}' in '{token}'")

    chiral = False
    if pos < n and inner[pos] == "@":
        chiral = True
        pos += 1
        if pos < n and inner[pos] == "@":
            pos += 1

    explicit_h: int | None = None
    if pos < n and inner[pos] == "H":
        pos += 1
        h_start = pos
        while pos < n and inner[pos].isdigit():
            pos += 1
        explicit_h = int(inner[h_start:pos]) if pos > h_start else 1

    charge = 0
    if pos < n and inner[pos] in ("+", "-"):
        sign = 1 if inner[pos] == "+" else -1
        pos += 1
        if pos < n and inner[pos].isdigit():
            c_start = pos
            while pos < n and inner[pos].isdigit():
                pos += 1
            charge = sign * int(inner[c_start:pos])
        elif pos < n and inner[pos] in ("+", "-"):
            cnt = 1
            ch = "+" if sign > 0 else "-"
            while pos < n and inner[pos] == ch:
                cnt += 1
                pos += 1
            charge = sign * cnt
        else:
            charge = sign

    if pos < n and inner[pos] == ":":
        pos += 1
        while pos < n and inner[pos].isdigit():
            pos += 1

    if pos != n:
        raise ValueError(f"Unrecognized trailing characters in bracket atom '{token}'")
    return sym, is_arom, charge, explicit_h, chiral


def _parse_single_smiles_component(smi: str) -> PyMolGraph:
    """Parse a single connected SMILES component into a PyMolGraph."""
    if not smi or any(ch.isspace() for ch in smi):
        raise ValueError("Empty or whitespace in SMILES component")
    graph = PyMolGraph()
    stack: list[int] = []
    ring_open: dict[int, tuple[int, float | None]] = {}
    prev_atom: int | None = None
    pending_bond: float | None = None
    i = 0
    n = len(smi)

    def add_atom_obj(atom: PyAtom) -> int:
        nonlocal prev_atom, pending_bond
        idx = len(graph.atoms)
        atom.idx = idx
        graph.atoms.append(atom)
        graph.adj[idx] = []
        if prev_atom is not None:
            if pending_bond is not None:
                order = pending_bond
            elif graph.atoms[prev_atom].is_aromatic and atom.is_aromatic:
                order = 1.5
            else:
                order = 1.0
            graph.bonds.append((prev_atom, idx, order))
            graph.adj[prev_atom].append((idx, order))
            graph.adj[idx].append((prev_atom, order))
        elif pending_bond is not None:
            raise ValueError("Leading bond without preceding atom")
        prev_atom = idx
        pending_bond = None
        return idx

    while i < n:
        ch = smi[i]
        if ch == "[":
            end = smi.find("]", i + 1)
            if end == -1:
                raise ValueError("Unclosed '[' in SMILES")
            sym, is_arom, charge, explicit_h, chiral = _parse_bracket_atom(smi[i : end + 1])
            if chiral:
                graph.has_chiral_stereo = True
            add_atom_obj(PyAtom(-1, sym, is_arom, charge, explicit_h, chiral))
            i = end + 1
        elif ch in ("C", "B") and i + 1 < n and smi[i : i + 2] in ("Cl", "Br"):
            add_atom_obj(PyAtom(-1, smi[i : i + 2], False, 0, None, False))
            i += 2
        elif ch in ("B", "C", "N", "O", "P", "S", "F", "I"):
            add_atom_obj(PyAtom(-1, ch, False, 0, None, False))
            i += 1
        elif ch in ("b", "c", "n", "o", "p", "s"):
            add_atom_obj(PyAtom(-1, ch.upper(), True, 0, None, False))
            i += 1
        elif ch == "-":
            if pending_bond is not None or prev_atom is None:
                raise ValueError("Invalid '-' bond placement")
            pending_bond = 1.0
            i += 1
        elif ch == "=":
            if pending_bond is not None or prev_atom is None:
                raise ValueError("Invalid '=' bond placement")
            pending_bond = 2.0
            i += 1
        elif ch == "#":
            if pending_bond is not None or prev_atom is None:
                raise ValueError("Invalid '#' bond placement")
            pending_bond = 3.0
            i += 1
        elif ch == ":":
            if pending_bond is not None or prev_atom is None:
                raise ValueError("Invalid ':' bond placement")
            pending_bond = 1.5
            i += 1
        elif ch in ("/", "\\"):
            if pending_bond is not None or prev_atom is None:
                raise ValueError("Invalid stereo bond placement")
            graph.has_ez_stereo = True
            pending_bond = 1.0
            i += 1
        elif ch == "(":
            if prev_atom is None or pending_bond is not None:
                raise ValueError("Invalid '(' branch placement")
            stack.append(prev_atom)
            i += 1
        elif ch == ")":
            if not stack or pending_bond is not None:
                raise ValueError("Unmatched ')' or trailing bond before ')'")
            prev_atom = stack.pop()
            i += 1
        elif ch == "%" or ch.isdigit():
            if prev_atom is None:
                raise ValueError("Ring closure before any atom")
            if ch == "%":
                if i + 2 >= n or not smi[i + 1 : i + 3].isdigit():
                    raise ValueError("Invalid '%nn' ring closure")
                rnum = int(smi[i + 1 : i + 3])
                i += 3
            else:
                rnum = int(ch)
                i += 1
            if rnum in ring_open:
                u, open_b = ring_open.pop(rnum)
                if u == prev_atom:
                    raise ValueError("Self-loop ring closure")
                if open_b is not None and pending_bond is not None and open_b != pending_bond:
                    raise ValueError("Conflicting ring closure bond orders")
                if pending_bond is not None:
                    order = pending_bond
                elif open_b is not None:
                    order = open_b
                elif graph.atoms[u].is_aromatic and graph.atoms[prev_atom].is_aromatic:
                    order = 1.5
                else:
                    order = 1.0
                # Check duplicate bond
                if any(nbr == prev_atom for nbr, _ in graph.adj[u]):
                    raise ValueError("Duplicate bond on ring closure")
                graph.bonds.append((u, prev_atom, order))
                graph.adj[u].append((prev_atom, order))
                graph.adj[prev_atom].append((u, order))
                pending_bond = None
            else:
                ring_open[rnum] = (prev_atom, pending_bond)
                pending_bond = None
        else:
            raise ValueError(f"Unexpected SMILES character '{ch}'")

    if stack:
        raise ValueError("Unclosed '(' in SMILES")
    if ring_open:
        raise ValueError(f"Unclosed ring closure(s): {sorted(ring_open.keys())}")
    if pending_bond is not None:
        raise ValueError("Trailing bond at end of SMILES")
    if not graph.atoms or not any(a.symbol != "H" for a in graph.atoms):
        raise ValueError("SMILES contains no heavy atoms")
    _finalize_graph(graph)
    return graph


def _finalize_graph(graph: PyMolGraph) -> None:
    """Compute implicit hydrogens, valence sanity checks, and ring membership."""
    n = len(graph.atoms)
    for u, atom in enumerate(graph.atoms):
        bo_sum = sum(order for _, order in graph.adj[u])
        sym = atom.symbol
        chg = atom.formal_charge
        if atom.explicit_h is not None:
            atom.implicit_h = atom.explicit_h
            # Valence sanity for common organic atoms
            total_v = bo_sum + atom.explicit_h
            if sym == "C" and chg == 0 and total_v > 4.6:
                raise ValueError("Carbon valence exceeded")
            continue
        if atom.is_aromatic:
            if sym == "C":
                atom.implicit_h = max(0, int(round(4.0 - bo_sum))) if bo_sum <= 3.5 else 0
            elif sym in ("N", "P"):
                # Aromatic n/p without [nH] has 0 implicit H in SMILES
                atom.implicit_h = 0
            else:
                atom.implicit_h = 0
            continue

        allowed = ORGANIC_VALENCES.get(sym, (4,))
        if sym == "N":
            allowed = (4,) if chg > 0 else (2,) if chg < 0 else (3,)
        elif sym == "O":
            allowed = (3,) if chg > 0 else (1,) if chg < 0 else (2,)
        elif sym == "C":
            allowed = (3,) if chg != 0 else (4,)
        elif sym == "S":
            allowed = (3, 5) if chg > 0 else (1,) if chg < 0 else (2, 4, 6)
        elif sym == "P":
            allowed = (4,) if chg > 0 else (3, 5)
        elif sym in ("F", "Cl", "Br", "I") and chg == 0:
            allowed = (1,)

        target_v = None
        for v in allowed:
            if bo_sum <= v + 1e-6:
                target_v = v
                break
        if target_v is None:
            raise ValueError(f"Valence exceeded for {sym} (bo={bo_sum}, charge={chg})")
        atom.implicit_h = max(0, int(round(target_v - bo_sum)))

    # Mark all ring atoms via Tarjan's bridge-finding algorithm
    tin = [-1] * n
    low = [-1] * n
    timer = [0]
    bridge_set: set[tuple[int, int]] = set()

    def dfs_bridge(u: int, p: int) -> None:
        tin[u] = low[u] = timer[0]
        timer[0] += 1
        for v, _ in graph.adj[u]:
            if v == p:
                continue
            if tin[v] != -1:
                low[u] = min(low[u], tin[v])
            else:
                dfs_bridge(v, u)
                low[u] = min(low[u], low[v])
                if low[v] > tin[u]:
                    bridge_set.add((min(u, v), max(u, v)))

    for idx in range(n):
        if tin[idx] == -1:
            dfs_bridge(idx, -1)

    for a in graph.atoms:
        a.in_ring = False
    for u, v, _ in graph.bonds:
        if (min(u, v), max(u, v)) not in bridge_set:
            graph.atoms[u].in_ring = True
            graph.atoms[v].in_ring = True


def _neutralize_py_graph(graph: PyMolGraph) -> tuple[PyMolGraph, bool]:
    """Neutralize simple conjugate acids/bases while preserving nitro & quat ammonium."""
    changed = False
    for u, atom in enumerate(graph.atoms):
        nbrs = graph.adj[u]
        if atom.formal_charge == -1 and atom.symbol == "O" and len(nbrs) == 1:
            v, order = nbrs[0]
            if order == 1.0 and graph.atoms[v].formal_charge == 0:
                atom.formal_charge = 0
                atom.explicit_h = 1
                changed = True
        elif atom.formal_charge == 1 and atom.symbol == "N":
            h_avail = atom.explicit_h if atom.explicit_h is not None else atom.implicit_h
            has_o_minus = any(
                graph.atoms[v].symbol == "O" and graph.atoms[v].formal_charge == -1 for v, _ in nbrs
            )
            if h_avail >= 1 and not has_o_minus:
                atom.formal_charge = 0
                atom.explicit_h = h_avail - 1
                changed = True
    if changed:
        _finalize_graph(graph)
    return graph, changed


def _canonical_py_smiles(graph: PyMolGraph) -> str:
    """Generate deterministic canonical SMILES using Morgan-like symmetry refinement."""
    n = len(graph.atoms)
    if n == 0:
        return ""
    inv = [
        hash(
            (
                ATOMIC_NUMS.get(a.symbol, 0),
                a.is_aromatic,
                a.formal_charge,
                a.implicit_h,
                a.in_ring,
                len(graph.adj[a.idx]),
            )
        )
        for a in graph.atoms
    ]
    # Use deterministic SHA-based refinement across 4 rounds
    def h_str(parts: tuple[Any, ...]) -> int:
        raw = repr(parts).encode("utf-8")
        return int.from_bytes(hashlib.blake2b(raw, digest_size=8).digest(), "big")

    inv = [
        h_str((ATOMIC_NUMS.get(a.symbol, 0), a.is_aromatic, a.formal_charge, a.implicit_h, a.in_ring, len(graph.adj[a.idx])))
        for a in graph.atoms
    ]
    for _ in range(4):
        new_inv = []
        for u in range(n):
            nbr_sigs = tuple(sorted((int(round(bo * 2)), inv[v]) for v, bo in graph.adj[u]))
            new_inv.append(h_str((inv[u], nbr_sigs)))
        inv = new_inv

    # Rank vertices canonically (break remaining ties by full local DFS signature)
    order = sorted(range(n), key=lambda u: inv[u])
    rank = [0] * n
    for r, u in enumerate(order):
        rank[u] = r

    def atom_str(u: int) -> str:
        a = graph.atoms[u]
        if a.formal_charge == 0 and a.symbol in ("B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"):
            if a.is_aromatic:
                if a.symbol == "N" and a.implicit_h > 0:
                    return "[nH]"
                return a.symbol.lower()
            return a.symbol
        sym = a.symbol.lower() if a.is_aromatic else a.symbol
        h_part = ""
        if a.implicit_h == 1:
            h_part = "H"
        elif a.implicit_h > 1:
            h_part = f"H{a.implicit_h}"
        c_part = ""
        if a.formal_charge == 1:
            c_part = "+"
        elif a.formal_charge > 1:
            c_part = f"+{a.formal_charge}"
        elif a.formal_charge == -1:
            c_part = "-"
        elif a.formal_charge < -1:
            c_part = f"{a.formal_charge}"
        return f"[{sym}{h_part}{c_part}]"

    def b_str(bo: float, u: int, v: int) -> str:
        if bo == 2.0:
            return "="
        if bo == 3.0:
            return "#"
        if bo == 1.0 and graph.atoms[u].is_aromatic and graph.atoms[v].is_aromatic:
            return "-"
        return ""

    visited = [False] * n
    parent = [-1] * n
    ring_map: dict[tuple[int, int], int] = {}
    next_r = [1]

    def dfs(u: int) -> str:
        nbrs = sorted(graph.adj[u], key=lambda item: (rank[item[0]], item[1]))
        closures: list[str] = []
        children: list[tuple[int, float]] = []
        for v, bo in nbrs:
            if v == parent[u]:
                continue
            edge = (min(u, v), max(u, v))
            if not visited[v]:
                visited[v] = True
                parent[v] = u
                children.append((v, bo))
            else:
                if edge not in ring_map:
                    ring_map[edge] = next_r[0]
                    next_r[0] += 1
                r = ring_map[edge]
                rs = f"%{r}" if r >= 10 else str(r)
                closures.append(f"{b_str(bo, u, v)}{rs}")
        s = atom_str(u) + "".join(closures)
        if len(children) == 1:
            v, bo = children[0]
            s += b_str(bo, u, v) + dfs(v)
        elif len(children) > 1:
            for v, bo in children[:-1]:
                s += f"({b_str(bo, u, v)}{dfs(v)})"
            v, bo = children[-1]
            s += b_str(bo, u, v) + dfs(v)
        return s

    start = order[0]
    visited[start] = True
    return dfs(start)


def _py_formula_and_mw(graph: PyMolGraph) -> tuple[str, float]:
    counts: dict[str, int] = {}
    mw = 0.0
    for a in graph.atoms:
        counts[a.symbol] = counts.get(a.symbol, 0) + 1
        mw += ATOMIC_WEIGHTS.get(a.symbol, 12.0)
        if a.implicit_h > 0:
            counts["H"] = counts.get("H", 0) + a.implicit_h
            mw += a.implicit_h * ATOMIC_WEIGHTS["H"]
    parts: list[str] = []
    if "C" in counts:
        c_n = counts.pop("C")
        parts.append("C" if c_n == 1 else f"C{c_n}")
        if "H" in counts:
            h_n = counts.pop("H")
            parts.append("H" if h_n == 1 else f"H{h_n}")
    for sym in sorted(counts.keys()):
        n = counts[sym]
        parts.append(sym if n == 1 else f"{sym}{n}")
    return "".join(parts), mw


def _py_ring_counts(graph: PyMolGraph) -> tuple[int, int]:
    """Return (total_rings, aromatic_rings) via cyclomatic number."""
    # Connected components
    n = len(graph.atoms)
    seen = [False] * n
    comps = 0
    for i in range(n):
        if not seen[i]:
            comps += 1
            q = [i]
            seen[i] = True
            while q:
                u = q.pop()
                for v, _ in graph.adj[u]:
                    if not seen[v]:
                        seen[v] = True
                        q.append(v)
    total_rings = max(0, len(graph.bonds) - n + comps)

    arom_atoms = {a.idx for a in graph.atoms if a.is_aromatic}
    if not arom_atoms:
        return total_rings, 0
    arom_edges = [(u, v) for u, v, bo in graph.bonds if u in arom_atoms and v in arom_atoms and bo == 1.5]
    if not arom_edges:
        return total_rings, 0
    arom_adj: dict[int, list[int]] = {u: [] for u in arom_atoms}
    for u, v in arom_edges:
        arom_adj[u].append(v)
        arom_adj[v].append(u)
    active_arom = [u for u in arom_atoms if arom_adj[u]]
    seen_a: set[int] = set()
    c_a = 0
    for u in active_arom:
        if u not in seen_a:
            c_a += 1
            q = [u]
            seen_a.add(u)
            while q:
                curr = q.pop()
                for v in arom_adj[curr]:
                    if v not in seen_a:
                        seen_a.add(v)
                        q.append(v)
    arom_rings = max(0, len(arom_edges) - len(active_arom) + c_a)
    return total_rings, arom_rings


def _py_has_small_ring(graph: PyMolGraph, target_syms: tuple[str, str, str]) -> bool:
    """Detect 3-membered rings matching sorted atom symbol tuple."""
    wanted = tuple(sorted(target_syms))
    adj_set = {u: {v for v, _ in graph.adj[u]} for u in range(len(graph.atoms))}
    for u, v, _ in graph.bonds:
        common = adj_set[u] & adj_set[v]
        for w in common:
            if w != u and w != v:
                tri = tuple(sorted((graph.atoms[u].symbol, graph.atoms[v].symbol, graph.atoms[w].symbol)))
                if tri == wanted:
                    return True
    return False


def _py_alerts(graph: PyMolGraph) -> tuple[list[str], list[str]]:
    """Detect 16+ classic PAINS and Brenk toxicophore alerts in pure Python."""
    pains: list[str] = []
    brenk: list[str] = []
    atoms = graph.atoms
    adj = graph.adj

    def is_carbonyl_c(idx: int) -> bool:
        return atoms[idx].symbol == "C" and any(atoms[v].symbol == "O" and bo == 2.0 for v, bo in adj[idx])

    def is_thiocarbonyl_c(idx: int) -> bool:
        return atoms[idx].symbol == "C" and any(atoms[v].symbol == "S" and bo == 2.0 for v, bo in adj[idx])

    def is_sulfonyl_s(idx: int) -> bool:
        if atoms[idx].symbol != "S":
            return False
        o_dbl = sum(1 for v, bo in adj[idx] if atoms[v].symbol == "O" and bo == 2.0)
        return o_dbl >= 2

    def is_nitro_n(idx: int) -> bool:
        if atoms[idx].symbol != "N":
            return False
        o_cnt = sum(
            1 for v, bo in adj[idx]
            if atoms[v].symbol == "O" and len(adj[v]) == 1 and (bo == 2.0 or atoms[v].formal_charge == -1)
        )
        return o_cnt >= 2

    # 1. Nitroaromatic & 2. Nitro aliphatic
    has_nitro_arom = False
    has_nitro_aliph = False
    for u, a in enumerate(atoms):
        if is_nitro_n(u):
            for v, _ in adj[u]:
                if atoms[v].symbol != "O":
                    if atoms[v].is_aromatic:
                        has_nitro_arom = True
                    else:
                        has_nitro_aliph = True
    if has_nitro_arom:
        brenk.append("nitroaromatic")
    if has_nitro_aliph:
        brenk.append("nitro_aliphatic")

    # 3. Michael acceptor: non-aromatic C=C adjacent to C=O, S(=O)2, or C#N
    for u, v, bo in graph.bonds:
        if bo == 2.0 and not atoms[u].is_aromatic and not atoms[v].is_aromatic:
            if atoms[u].symbol == "C" and atoms[v].symbol == "C":
                ewg = False
                for end in (u, v):
                    for nbr, nbo in adj[end]:
                        if nbr in (u, v) or nbo != 1.0:
                            continue
                        if is_carbonyl_c(nbr) or is_sulfonyl_s(nbr):
                            ewg = True
                        elif atoms[nbr].symbol == "C" and any(atoms[k].symbol == "N" and kbo == 3.0 for k, kbo in adj[nbr]):
                            ewg = True
                if ewg:
                    pains.append("michael_acceptor")
                    break

    # 4. Aldehyde: non-aromatic C(=O)H with another C or H neighbor (not acid/ester/amide)
    for u, a in enumerate(atoms):
        if is_carbonyl_c(u) and not a.is_aromatic:
            has_hetero_single = any(atoms[v].symbol in ("O", "N", "S", "F", "Cl", "Br", "I") and bo == 1.0 for v, bo in adj[u])
            if not has_hetero_single and a.implicit_h >= 1:
                brenk.append("aldehyde")
                break

    # 5. Acyl halide & 6. Sulfonyl halide
    for u, a in enumerate(atoms):
        if is_carbonyl_c(u) and any(atoms[v].symbol in ("F", "Cl", "Br", "I") and bo == 1.0 for v, bo in adj[u]):
            brenk.append("acyl_halide")
            break
    for u, a in enumerate(atoms):
        if is_sulfonyl_s(u) and any(atoms[v].symbol in ("F", "Cl", "Br", "I") and bo == 1.0 for v, bo in adj[u]):
            brenk.append("sulfonyl_halide")
            break

    # 7. Epoxide & 8. Aziridine
    if _py_has_small_ring(graph, ("C", "C", "O")):
        brenk.append("epoxide")
    if _py_has_small_ring(graph, ("C", "C", "N")):
        brenk.append("aziridine")

    # 9. Catechol (adjacent aromatic carbons each bearing -OH)
    oh_arom = {
        u for u, a in enumerate(atoms)
        if a.is_aromatic and a.symbol == "C" and any(
            atoms[v].symbol == "O" and not atoms[v].is_aromatic and bo == 1.0 and len(adj[v]) == 1 and atoms[v].implicit_h >= 1
            for v, bo in adj[u]
        )
    }
    if any(v in oh_arom for u in oh_arom for v, _ in adj[u]):
        pains.append("catechol")

    # 10. Quinone (non-aromatic 6-membered ring with two carbonyl carbons)
    carb_ring_c = [u for u, a in enumerate(atoms) if a.in_ring and not a.is_aromatic and is_carbonyl_c(u)]
    if len(carb_ring_c) >= 2:
        # Check distance 2 or 3 in ring with alkene C=C
        has_alkene_ring = any(
            bo == 2.0 and atoms[u].symbol == "C" and atoms[v].symbol == "C" and atoms[u].in_ring and atoms[v].in_ring
            for u, v, bo in graph.bonds
        )
        if has_alkene_ring:
            pains.append("quinone")

    # 11. Rhodanine / ene_rhodanine (ring S-C(=S)-N or C(=S)-N-C(=O) in ring)
    for u, a in enumerate(atoms):
        if a.in_ring and is_thiocarbonyl_c(u):
            nbr_syms = [atoms[v].symbol for v, bo in adj[u] if bo == 1.0 and atoms[v].in_ring]
            if "S" in nbr_syms and "N" in nbr_syms:
                pains.append("rhodanine")
                break

    # 12. Hydrazine (acyclic or non-aromatic N-N single bond) & 13. Azo (N=N)
    for u, v, bo in graph.bonds:
        if atoms[u].symbol == "N" and atoms[v].symbol == "N":
            if bo == 1.0 and not (atoms[u].is_aromatic and atoms[v].is_aromatic):
                if "hydrazine" not in brenk:
                    brenk.append("hydrazine")
            elif bo == 2.0 and not atoms[u].is_aromatic and not atoms[v].is_aromatic:
                if "azo" not in pains:
                    pains.append("azo")

    # 14. Azide (N=N=N or [N-]=[N+]=N)
    for u, a in enumerate(atoms):
        if a.symbol == "N" and sum(1 for v, bo in adj[u] if atoms[v].symbol == "N" and bo >= 2.0) >= 2:
            brenk.append("azide")
            break

    # 15. Thiourea (N-C(=S)-N)
    for u, a in enumerate(atoms):
        if is_thiocarbonyl_c(u) and sum(1 for v, bo in adj[u] if atoms[v].symbol == "N" and bo == 1.0) >= 2:
            brenk.append("thiourea")
            break

    # 16. Peroxide (O-O)
    if any(atoms[u].symbol == "O" and atoms[v].symbol == "O" and bo == 1.0 for u, v, bo in graph.bonds):
        brenk.append("peroxide")

    # 17. Isocyanate / Isothiocyanate (N=C=O / N=C=S)
    for u, a in enumerate(atoms):
        if a.symbol == "C" and any(atoms[v].symbol == "N" and bo == 2.0 for v, bo in adj[u]):
            if any(atoms[v].symbol in ("O", "S") and bo == 2.0 for v, bo in adj[u]):
                brenk.append("isocyanate_isothiocyanate")
                break

    # 18. Reactive primary alkyl halide (C-[Br,I] on sp3 C with H)
    for u, a in enumerate(atoms):
        if a.symbol == "C" and not a.is_aromatic and a.implicit_h >= 2:
            if any(atoms[v].symbol in ("Br", "I") and bo == 1.0 for v, bo in adj[u]):
                brenk.append("alkyl_halide_reactive")
                break

    return sorted(set(pains)), sorted(set(brenk))


def _py_descriptors(graph: PyMolGraph) -> dict[str, Any]:
    """Compute physicochemical descriptors, QED, SA score, basic amine, scaffold, FP."""
    atoms = graph.atoms
    adj = graph.adj
    heavy_atoms = sum(1 for a in atoms if a.symbol != "H")
    formula, mw = _py_formula_and_mw(graph)
    ring_count, arom_ring_count = _py_ring_counts(graph)
    arom_atoms = sum(1 for a in atoms if a.is_aromatic and a.symbol != "H")

    # HBD & HBA
    hbd = sum(1 for a in atoms if a.symbol in ("N", "O") and a.implicit_h >= 1)
    hba = sum(1 for a in atoms if a.symbol in ("N", "O") and a.formal_charge <= 0)

    # Rotatable bonds: single non-ring non-terminal bond, excluding amide C(=O)-N
    rot_bonds = 0
    for u, v, bo in graph.bonds:
        if bo != 1.0:
            continue
        if atoms[u].symbol == "H" or atoms[v].symbol == "H":
            continue
        if atoms[u].in_ring and atoms[v].in_ring:
            # Check if same ring or bridge; bridges have at least one atom in_ring? No, bridge edges are NOT in rings, wait:
            # if both endpoints are in rings AND edge is not a bridge, it's a ring bond!
            pass
        # Use bridge check: a bond is a ring bond iff it is NOT a bridge. Let's compute bridges once?
        # Even simpler: check if both u and v are in_ring and reachable without (u, v)
        if atoms[u].in_ring and atoms[v].in_ring:
            # BFS from u avoiding edge (u, v)
            seen = {u}
            q = [u]
            reached = False
            while q and not reached:
                curr = q.pop()
                for nbr, _ in adj[curr]:
                    if (curr == u and nbr == v) or (curr == v and nbr == u):
                        continue
                    if nbr == v:
                        reached = True
                        break
                    if nbr not in seen:
                        seen.add(nbr)
                        q.append(nbr)
            if reached:
                continue
        if len(adj[u]) <= 1 or len(adj[v]) <= 1:
            continue
        # Exclude amide C(=O)-N
        is_amide = False
        for a1, a2 in ((u, v), (v, u)):
            if atoms[a1].symbol == "C" and atoms[a2].symbol == "N":
                if any(atoms[k].symbol == "O" and kbo == 2.0 for k, kbo in adj[a1]):
                    is_amide = True
        if is_amide:
            continue
        rot_bonds += 1

    # Fsp3
    c_atoms = [a for a in atoms if a.symbol == "C"]
    if c_atoms:
        sp3_c = sum(1 for a in c_atoms if not a.is_aromatic and all(bo == 1.0 for _, bo in adj[a.idx]))
        fsp3 = sp3_c / len(c_atoms)
    else:
        fsp3 = 0.0

    # Ertl-style TPSA & Wildman-Crippen-style cLogP
    tpsa = 0.0
    clogp = 0.0
    has_basic_amine = False
    has_acidic_group = False

    for u, a in enumerate(atoms):
        sym = a.symbol
        nbrs = adj[u]
        deg = len(nbrs)
        h_cnt = a.implicit_h
        has_dbl_o = any(atoms[v].symbol == "O" and bo == 2.0 for v, bo in nbrs)
        has_arom_nbr = any(atoms[v].is_aromatic for v, _ in nbrs)
        max_bo = max((bo for _, bo in nbrs), default=0.0)

        if sym == "O":
            if max_bo == 2.0:  # =O
                tpsa += 17.07
                clogp += -0.2893
            elif a.is_aromatic:
                tpsa += 13.14
                clogp += 0.1552
            elif h_cnt >= 1:
                # Alcohol vs carboxylic acid / phenol
                adj_c = [v for v, _ in nbrs if atoms[v].symbol == "C"]
                if adj_c and any(atoms[k].symbol == "O" and kbo == 2.0 for k, kbo in adj[adj_c[0]]):
                    tpsa += 20.23
                    clogp += -0.68
                    has_acidic_group = True
                elif adj_c and atoms[adj_c[0]].is_aromatic:
                    tpsa += 20.23
                    clogp += -0.15
                else:
                    tpsa += 20.23
                    clogp += -0.35
            else:
                tpsa += 9.23
                clogp += -0.18
        elif sym == "N":
            if a.is_aromatic:
                tpsa += 15.79 if h_cnt >= 1 else 12.89
                clogp += -0.32
            elif any(atoms[v].symbol == "O" and (bo == 2.0 or atoms[v].formal_charge == -1) for v, bo in nbrs):
                # Nitro N
                tpsa += 11.68
                clogp += -0.45
            elif max_bo == 3.0:  # C#N
                tpsa += 23.79
                clogp += -0.56
            elif max_bo == 2.0:
                tpsa += 23.85 if h_cnt >= 1 else 12.36
                clogp += -0.25
            else:
                # Single-bonded aliphatic N: check amide / sulfonamide / aniline / basic amine
                is_amide_n = any(
                    (atoms[v].symbol == "C" and any(atoms[k].symbol == "O" and kbo == 2.0 for k, kbo in adj[v]))
                    or (atoms[v].symbol == "S" and any(atoms[k].symbol == "O" and kbo == 2.0 for k, kbo in adj[v]))
                    for v, _ in nbrs
                )
                if h_cnt >= 2:
                    tpsa += 26.02
                elif h_cnt == 1:
                    tpsa += 12.03
                else:
                    tpsa += 3.24
                if is_amide_n:
                    clogp += -0.52
                elif has_arom_nbr:
                    clogp += -0.45
                else:
                    clogp += -0.95
                    if a.formal_charge >= 0:
                        has_basic_amine = True
        elif sym == "C":
            if a.is_aromatic:
                if any(atoms[v].symbol in ("N", "O", "S") for v, _ in nbrs):
                    clogp += 0.14
                elif deg >= 3:
                    clogp += 0.30
                else:
                    clogp += 0.34
            elif has_dbl_o:
                clogp += 0.05
            elif max_bo >= 2.0:
                clogp += 0.18
            else:
                # sp3 C
                hetero_nbrs = sum(1 for v, _ in nbrs if atoms[v].symbol in ("N", "O"))
                if hetero_nbrs > 0:
                    clogp += 0.04
                elif deg == 1:
                    clogp += 0.48
                elif deg == 2:
                    clogp += 0.45
                else:
                    clogp += 0.38
        elif sym == "F":
            clogp += 0.42
        elif sym == "Cl":
            clogp += 0.69
        elif sym == "Br":
            clogp += 0.89
        elif sym == "I":
            clogp += 1.19
        elif sym == "S":
            if any(atoms[v].symbol == "O" and bo == 2.0 for v, bo in nbrs):
                clogp += -0.25
                if any(atoms[v].symbol == "O" and bo == 1.0 and atoms[v].implicit_h >= 1 for v, bo in nbrs):
                    has_acidic_group = True
            elif a.is_aromatic:
                clogp += 0.65
            else:
                clogp += 0.48
        elif sym == "P":
            clogp += 0.10

    pains_alerts, brenk_alerts = _py_alerts(graph)
    n_alerts = len(pains_alerts) + len(brenk_alerts)
    qed = _compute_qed_from_props(mw, clogp, hba, hbd, tpsa, rot_bonds, arom_ring_count, n_alerts)
    sa_score = _compute_complexity_score(heavy_atoms, ring_count, arom_ring_count, graph.has_chiral_stereo, atoms)
    scaffold = _py_murcko_scaffold(graph)
    fp = _py_morgan_fp(graph)

    return {
        "formula": formula,
        "mw": round(mw, 3),
        "heavy_atoms": heavy_atoms,
        "clogp": round(clogp, 3),
        "tpsa": round(tpsa, 2),
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rot_bonds,
        "fsp3": round(fsp3, 3),
        "ring_count": ring_count,
        "aromatic_ring_count": arom_ring_count,
        "aromatic_atoms": arom_atoms,
        "qed": round(qed, 4),
        "sa_score": round(sa_score, 3),
        "pains_alerts": pains_alerts,
        "brenk_alerts": brenk_alerts,
        "n_alerts": n_alerts,
        "murcko_scaffold": scaffold,
        "has_basic_amine": has_basic_amine,
        "has_acidic_group": has_acidic_group,
        "fp": fp,
    }


def _ads(x: float, params: tuple[float, float, float, float, float, float, float]) -> float:
    a, b, c, d, e, f, d_max = params
    exp1 = math.exp(-max(-50.0, min(50.0, (x - c + d / 2.0) / max(e, 1e-9))))
    exp2 = math.exp(-max(-50.0, min(50.0, (x - c - d / 2.0) / max(f, 1e-9))))
    val = a + (b / (1.0 + exp1)) * (1.0 - 1.0 / (1.0 + exp2))
    return max(1e-6, min(1.0, val / d_max))


def _compute_qed_from_props(
    mw: float, clogp: float, hba: int, hbd: int, psa: float, rotb: int, arom: int, alerts: int
) -> float:
    vals = {
        "MW": _ads(mw, _QED_PARAMS["MW"]),
        "ALOGP": _ads(clogp, _QED_PARAMS["ALOGP"]),
        "HBA": _ads(float(hba), _QED_PARAMS["HBA"]),
        "HBD": _ads(float(hbd), _QED_PARAMS["HBD"]),
        "PSA": _ads(psa, _QED_PARAMS["PSA"]),
        "ROTB": _ads(float(rotb), _QED_PARAMS["ROTB"]),
        "AROM": _ads(float(arom), _QED_PARAMS["AROM"]),
        "ALERTS": _ads(float(alerts), _QED_PARAMS["ALERTS"]),
    }
    num = sum(_QED_WEIGHTS[k] * math.log(vals[k]) for k in _QED_WEIGHTS)
    den = sum(_QED_WEIGHTS.values())
    return float(math.exp(num / den))


def _compute_complexity_score(
    heavy_atoms: int, ring_count: int, arom_ring_count: int, has_chiral: bool, atoms_or_hetero: Any
) -> float:
    """Deterministic synthetic accessibility / complexity heuristic bounded in [1.0, 10.0]."""
    if isinstance(atoms_or_hetero, int):
        n_hetero = atoms_or_hetero
    else:
        n_hetero = sum(1 for a in atoms_or_hetero if a.symbol not in ("C", "H"))
    aliph_rings = max(0, ring_count - arom_ring_count)
    score = (
        1.15
        + 0.45 * math.log1p(max(0, heavy_atoms - 6))
        + 0.35 * aliph_rings
        + 0.15 * arom_ring_count
        + 0.55 * (1 if has_chiral else 0)
        + 0.12 * max(0, n_hetero - 3)
        + 0.25 * max(0, ring_count - 3)
    )
    return max(1.0, min(10.0, score))


def _py_murcko_scaffold(graph: PyMolGraph) -> str:
    """Extract ring-linker Bemis-Murcko scaffold canonical SMILES."""
    ring_atoms = {a.idx for a in graph.atoms if a.in_ring}
    if not ring_atoms:
        return "acyclic"
    deg = {u: len(graph.adj[u]) for u in range(len(graph.atoms))}
    removed: set[int] = set()
    changed = True
    while changed:
        changed = False
        for u in list(range(len(graph.atoms))):
            if u in removed or u in ring_atoms:
                continue
            active_nbrs = [(v, bo) for v, bo in graph.adj[u] if v not in removed]
            if len(active_nbrs) <= 1:
                # Keep exocyclic =O/=S/=N directly attached to a ring or linker? Wait: remove terminal degree-1 non-ring atoms
                removed.add(u)
                changed = True
    kept = sorted(set(range(len(graph.atoms))) - removed)
    if not kept:
        return "acyclic"
    idx_map = {old: new for new, old in enumerate(kept)}
    sub = PyMolGraph()
    for old in kept:
        a = graph.atoms[old]
        lost_exocyclic_single = any(
            v in removed and bo == 1.0 for v, bo in graph.adj[old]
        )
        exp_h = a.explicit_h if a.is_aromatic else None
        if a.is_aromatic and a.symbol == "N" and lost_exocyclic_single:
            exp_h = 1
        sub.atoms.append(
            PyAtom(
                idx=idx_map[old],
                symbol=a.symbol,
                is_aromatic=a.is_aromatic,
                formal_charge=a.formal_charge,
                explicit_h=exp_h,
                chiral=False,
                in_ring=a.in_ring,
            )
        )
        sub.adj[idx_map[old]] = []
    for u, v, bo in graph.bonds:
        if u in idx_map and v in idx_map:
            nu, nv = idx_map[u], idx_map[v]
            sub.bonds.append((nu, nv, bo))
            sub.adj[nu].append((nv, bo))
            sub.adj[nv].append((nu, bo))
    _finalize_graph(sub)
    return _canonical_py_smiles(sub)


def _py_morgan_fp(graph: PyMolGraph, radius: int = 2, n_bits: int = 1024) -> int:
    """Compute 1024-bit circular Morgan-like radius-2 atom-neighborhood bitmask."""
    n = len(graph.atoms)
    if n == 0:
        return 0

    def h_bit(tag: tuple[Any, ...]) -> int:
        raw = repr(tag).encode("utf-8")
        return int.from_bytes(hashlib.blake2b(raw, digest_size=8).digest(), "big")

    invs = [
        h_bit((0, ATOMIC_NUMS.get(a.symbol, 0), len(graph.adj[a.idx]), a.implicit_h, a.formal_charge, a.is_aromatic, a.in_ring))
        for a in graph.atoms
    ]
    fp = 0
    for inv in invs:
        fp |= 1 << (inv % n_bits)
    for r in range(1, radius + 1):
        new_invs: list[int] = []
        for u in range(n):
            nbrs = sorted((int(round(bo * 2)), invs[v]) for v, bo in graph.adj[u])
            hv = h_bit((r, invs[u], tuple(nbrs)))
            new_invs.append(hv)
            fp |= 1 << (hv % n_bits)
        invs = new_invs
    return fp


def tanimoto_similarity(fp1: int, fp2: int) -> float:
    """Compute Tanimoto similarity between two integer bitmask fingerprints."""
    union = (fp1 | fp2).bit_count()
    if union == 0:
        return 0.0
    return (fp1 & fp2).bit_count() / union


# ==============================================================================
# Unified Dual-Engine Molecule API (RDKit primary + Pure-Python fallback)
# ==============================================================================

_RDKIT_FILTER_CATALOG = None


def _get_rdkit_filter_catalog():
    global _RDKIT_FILTER_CATALOG
    if _RDKIT_FILTER_CATALOG is None and HAS_RDKIT:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        _RDKIT_FILTER_CATALOG = FilterCatalog.FilterCatalog(params)
    return _RDKIT_FILTER_CATALOG


def _rdkit_fp_int(mol) -> int:
    if getattr(sys.modules[__name__], "rdFingerprintGenerator", None) is not None:
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
        bv = gen.GetFingerprint(mol)
    else:
        bv = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    fp = 0
    for bit_idx in bv.GetOnBits():
        fp |= 1 << int(bit_idx)
    return fp


def _rdkit_alerts(mol) -> tuple[list[str], list[str]]:
    cat = _get_rdkit_filter_catalog()
    pains: list[str] = []
    brenk: list[str] = []
    if cat is not None:
        for entry in cat.GetMatches(mol):
            desc = str(entry.GetDescription())
            prop = ""
            if entry.HasProp("FilterSet"):
                prop = str(entry.GetProp("FilterSet"))
            if "PAINS" in prop.upper() or "PAINS" in desc.upper():
                pains.append(desc)
            else:
                brenk.append(desc)
    return sorted(set(pains)), sorted(set(brenk))


def _rdkit_basic_acid(mol) -> tuple[bool, bool]:
    basic_smarts = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O);!$(NS=O);!$(Nc);!$([N+](=O)[O-])]")
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
    sulf_acid = Chem.MolFromSmarts("[SX4](=O)(=O)[OX2H1,OX1-]")
    has_basic = bool(basic_smarts and mol.HasSubstructMatch(basic_smarts))
    has_acid = bool((acid_smarts and mol.HasSubstructMatch(acid_smarts)) or (sulf_acid and mol.HasSubstructMatch(sulf_acid)))
    return has_basic, has_acid


def standardize_smiles(raw_smiles: str) -> dict[str, Any]:
    """Standardize a SMILES string using RDKit (if available) or pure-Python graph fallback."""
    smi = (raw_smiles or "").strip()
    if not smi:
        return {"is_valid": False, "error": "Empty SMILES string", "backend": BACKEND}

    if HAS_RDKIT:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return {"is_valid": False, "error": "Invalid SMILES rejected by RDKit", "backend": BACKEND}
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        had_salt = len(frags) > 1
        stripped_fragments: list[str] = []
        if had_salt:
            chooser = rdMolStandardize.LargestFragmentChooser()
            parent = chooser.choose(mol)
            parent_smi = Chem.MolToSmiles(parent, isomericSmiles=True, canonical=True)
            for f in frags:
                fs = Chem.MolToSmiles(f, isomericSmiles=True, canonical=True)
                if fs != parent_smi:
                    stripped_fragments.append(fs)
            mol = parent
        uncharger = rdMolStandardize.Uncharger()
        before_uncharge = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        mol = uncharger.uncharge(mol)
        canon_smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        neutralized = canon_smi != before_uncharge
        try:
            key = Chem.MolToInchiKey(mol) or ""
        except Exception:
            key = ""
        if not key:
            h = hashlib.sha256(canon_smi.encode("utf-8")).hexdigest().upper()
            key = f"RDKEY-{h[:14]}-{h[14:24]}-N"
        chiral_cnt = len(Chem.FindMolChiralCenters(mol, includeUnassigned=False))
        ez_cnt = sum(1 for b in mol.GetBonds() if b.GetStereo() != Chem.BondStereo.STEREONONE)
        return {
            "is_valid": True,
            "backend": BACKEND,
            "canonical_smiles": canon_smi,
            "structure_key": key,
            "had_salt": had_salt,
            "stripped_fragments": stripped_fragments,
            "neutralized": neutralized,
            "has_stereochemistry": bool(chiral_cnt > 0 or ez_cnt > 0),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "mw": round(float(Descriptors.MolWt(mol)), 3),
            "heavy_atoms": int(mol.GetNumHeavyAtoms()),
            "_mol": mol,
        }

    # Pure-Python SMILES graph fallback
    raw_parts = smi.split(".")
    parsed: list[tuple[str, PyMolGraph]] = []
    for part in raw_parts:
        if not part:
            return {"is_valid": False, "error": "Empty disconnected '.' fragment", "backend": BACKEND}
        try:
            g = _parse_single_smiles_component(part)
        except ValueError as exc:
            return {"is_valid": False, "error": str(exc), "backend": BACKEND}
        parsed.append((part, g))

    had_salt = len(parsed) > 1
    scored = []
    for idx, (part_smi, g) in enumerate(parsed):
        heavy = sum(1 for a in g.atoms if a.symbol != "H")
        carbons = sum(1 for a in g.atoms if a.symbol == "C")
        scored.append(((heavy, carbons, len(part_smi), -idx), part_smi, g))
    scored.sort(key=lambda x: x[0], reverse=True)
    _, _, parent_g = scored[0]
    stripped_fragments = [p_smi for _, p_smi, _ in scored[1:]]
    parent_g, neutralized = _neutralize_py_graph(parent_g)
    canon_smi = _canonical_py_smiles(parent_g)
    h = hashlib.sha256(canon_smi.encode("utf-8")).hexdigest().upper()
    key = f"PYKEY-{h[:14]}-{h[14:24]}-N"
    formula, mw = _py_formula_and_mw(parent_g)
    heavy_atoms = sum(1 for a in parent_g.atoms if a.symbol != "H")
    return {
        "is_valid": True,
        "backend": BACKEND,
        "canonical_smiles": canon_smi,
        "structure_key": key,
        "had_salt": had_salt,
        "stripped_fragments": stripped_fragments,
        "neutralized": neutralized,
        "has_stereochemistry": bool(parent_g.has_chiral_stereo or parent_g.has_ez_stereo),
        "formula": formula,
        "mw": round(mw, 3),
        "heavy_atoms": heavy_atoms,
        "_graph": parent_g,
    }


def featurize_smiles(raw_smiles: str) -> dict[str, Any]:
    """Standardize and compute full physicochemical, alert, scaffold, and FP descriptors."""
    std = standardize_smiles(raw_smiles)
    if not std["is_valid"]:
        return std

    if HAS_RDKIT:
        mol = std.pop("_mol")
        mw = float(Descriptors.MolWt(mol))
        clogp = float(Crippen.MolLogP(mol))
        tpsa = float(rdMolDescriptors.CalcTPSA(mol))
        hbd = int(rdMolDescriptors.CalcNumHBD(mol))
        hba = int(rdMolDescriptors.CalcNumHBA(mol))
        rot_bonds = int(rdMolDescriptors.CalcNumRotatableBonds(mol))
        fsp3 = float(rdMolDescriptors.CalcFractionCSP3(mol))
        ring_count = int(rdMolDescriptors.CalcNumRings(mol))
        arom_ring_count = int(rdMolDescriptors.CalcNumAromaticRings(mol))
        arom_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        heavy_atoms = int(mol.GetNumHeavyAtoms())
        qed_val = float(QED.qed(mol))
        pains_alerts, brenk_alerts = _rdkit_alerts(mol)
        n_alerts = len(pains_alerts) + len(brenk_alerts)
        n_hetero = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6))
        sa_score = _compute_complexity_score(
            heavy_atoms, ring_count, arom_ring_count, std["has_stereochemistry"], n_hetero
        )
        scaff_mol = MurckoScaffold.GetScaffoldForMol(mol)
        scaff_smi = Chem.MolToSmiles(scaff_mol, canonical=True) if scaff_mol.GetNumAtoms() > 0 else "acyclic"
        has_basic, has_acid = _rdkit_basic_acid(mol)
        fp = _rdkit_fp_int(mol)
        std.update(
            {
                "mw": round(mw, 3),
                "clogp": round(clogp, 3),
                "tpsa": round(tpsa, 2),
                "hbd": hbd,
                "hba": hba,
                "rotatable_bonds": rot_bonds,
                "fsp3": round(fsp3, 3),
                "ring_count": ring_count,
                "aromatic_ring_count": arom_ring_count,
                "aromatic_atoms": arom_atoms,
                "qed": round(qed_val, 4),
                "sa_score": round(sa_score, 3),
                "pains_alerts": pains_alerts,
                "brenk_alerts": brenk_alerts,
                "n_alerts": n_alerts,
                "murcko_scaffold": scaff_smi,
                "has_basic_amine": has_basic,
                "has_acidic_group": has_acid,
                "fp": fp,
            }
        )
        return std

    graph = std.pop("_graph")
    desc = _py_descriptors(graph)
    std.update(desc)
    return std


def taylor_butina_cluster(fps: list[int], distance_cutoff: float = 0.4) -> tuple[list[int], list[bool]]:
    """Taylor-Butina clustering on 1024-bit fingerprints. Returns (1-based cluster_id, is_centroid)."""
    n = len(fps)
    if n == 0:
        return [], []
    sim_cutoff = 1.0 - float(distance_cutoff)
    nbrs: list[list[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if tanimoto_similarity(fps[i], fps[j]) >= sim_cutoff:
                nbrs[i].append(j)
                nbrs[j].append(i)
    order = sorted(range(n), key=lambda i: (-len(nbrs[i]), i))
    assigned = [False] * n
    clusters: list[list[int]] = []
    centroids: list[int] = []
    for seed in order:
        if assigned[seed]:
            continue
        members = [seed]
        assigned[seed] = True
        for j in nbrs[seed]:
            if not assigned[j]:
                assigned[j] = True
                members.append(j)
        clusters.append(members)
        centroids.append(seed)
    # Sort clusters by descending size
    paired = sorted(zip(clusters, centroids), key=lambda x: (-len(x[0]), x[1]))
    cluster_ids = [0] * n
    is_centroid = [False] * n
    for cid, (members, cent) in enumerate(paired, start=1):
        is_centroid[cent] = True
        for m in members:
            cluster_ids[m] = cid
    return cluster_ids, is_centroid


# ==============================================================================
# Plotting & Artifact Persistence
# ==============================================================================

def setup_matplotlib():
    """Configure headless Agg backend and publication defaults, or return None if unavailable."""
    if not HAS_MATPLOTLIB:
        return None
    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 180,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.5,
            "legend.fontsize": 8.5,
        }
    )
    return plt


def write_fallback_png(path: str | Path, width: int = 960, height: int = 440, accent_rgb: tuple[int, int, int] = (15, 118, 110)) -> None:
    """Write a valid crisp PNG raster plot summary fallback when matplotlib is absent."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Generate clean light canvas with subtle border & accent header band
    bg = bytes((248, 250, 252))
    border = bytes((203, 213, 225))
    accent = bytes(accent_rgb)
    raw_rows = bytearray()
    for y in range(height):
        raw_rows.append(0)  # filter type 0
        if y < 6 or y >= height - 6:
            raw_rows.extend(border * width)
        elif 16 <= y < 28:
            raw_rows.extend(border * 16 + accent * (width - 32) + border * 16)
        else:
            raw_rows.extend(border * 6 + bg * (width - 12) + border * 6)

    def chunk(tag: bytes, data: bytes) -> bytes:
        crc = zlib.crc32(tag + data) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    idat = zlib.compress(bytes(raw_rows), level=6)
    png_bytes = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IDAT", idat) + chunk(b"IEND", b"")
    p.write_bytes(png_bytes)


def write_tsv(
    rows_or_df: Any,
    columns: list[str] | None,
    path: str | Path,
    tool: str,
    title: str,
    **kwargs: Any,
) -> str:
    """Write rows/DataFrame to TSV (using .to_numpy(copy=True) when DataFrame) and emit .phi-artifact.json."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if HAS_PANDAS and isinstance(rows_or_df, pd.DataFrame):
        df = rows_or_df
        if columns is not None:
            df = df[columns]
        # Safe Pandas 3 CoW check via .to_numpy(copy=True)
        _ = df.to_numpy(copy=True)
        df.to_csv(p, sep="\t", index=False)
        shape = df.shape
    else:
        rows = list(rows_or_df)
        cols = columns or (list(rows[0].keys()) if rows else [])
        if HAS_PANDAS:
            df = pd.DataFrame(rows, columns=cols)
            _ = df.to_numpy(copy=True)
            df.to_csv(p, sep="\t", index=False)
            shape = df.shape
        else:
            with open(p, "wt", encoding="utf-8", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
                writer.writerow(cols)
                for r in rows:
                    writer.writerow([r.get(c, "") for c in cols])
            shape = (len(rows), len(cols))
    return write_artifact(p, tool=tool, title=title, shape=shape, **kwargs)


def write_artifact(
    path: str | Path,
    tool: str,
    title: str,
    *,
    skill: str = "chem-workflow",
    inputs: tuple[str | None, ...] | list[str | None] = (),
    parameters: dict[str, Any] | None = None,
    shape: tuple[int, int] | dict[str, int] | None = None,
    scientific_context: dict[str, Any] | None = None,
) -> str:
    """Write Contract-1.0.0 .phi-artifact.json sidecar for tables and figures."""
    p = Path(path)
    ext = p.suffix.lower()
    figures = {".png": "image/png", ".pdf": "application/pdf", ".svg": "image/svg+xml"}
    tables = {".tsv": "text/tab-separated-values", ".csv": "text/csv"}
    if ext in figures:
        kind, media = "figure", figures[ext]
    elif ext in tables:
        kind, media = "table", tables[ext]
    else:
        raise ValueError(f"Unsupported artifact extension: {ext}")

    descriptor: dict[str, Any] = {
        "contractVersion": "1.0.0",
        "kind": kind,
        "file": p.name,
        "mediaType": media,
        "title": title[:200],
        "provenance": {
            "createdAt": datetime.now(timezone.utc).isoformat(),
            "tool": tool,
            "skill": skill,
        },
    }
    if kind == "figure":
        block: dict[str, Any] = {"format": ext[1:]}
        if ext == ".png":
            with open(p, "rb") as fh:
                hdr = fh.read(24)
            if hdr[:8] == b"\x89PNG\r\n\x1a\n":
                w, h = struct.unpack(">II", hdr[16:24])
                block["widthPx"], block["heightPx"] = int(w), int(h)
        descriptor["figure"] = block
    elif shape is not None:
        if isinstance(shape, dict):
            rows_n = int(shape.get("rows", 0))
            cols_n = int(shape.get("columns", shape.get("cols", 0)))
        else:
            rows_n, cols_n = int(shape[0]), int(shape[1])
        descriptor["table"] = {"rows": rows_n, "columns": cols_n}

    if scientific_context:
        descriptor["scientificContext"] = clean_json(scientific_context)

    cwd = Path.cwd().resolve()
    prov_in = []
    for raw in inputs:
        if not raw:
            continue
        ip = Path(raw)
        try:
            rel = ip.resolve().relative_to(cwd)
        except ValueError:
            rel = Path(ip.name)
        entry: dict[str, Any] = {"path": str(rel)}
        if ip.is_file():
            entry["sha256"] = hashlib.sha256(ip.read_bytes()).hexdigest()
        prov_in.append(entry)
    if prov_in:
        descriptor["provenance"]["inputs"] = prov_in
    if parameters:
        descriptor["provenance"]["parameters"] = clean_json(parameters)

    target = p.with_name(p.name + ".phi-artifact.json")
    tmp = target.with_name(target.name + ".tmp")
    tmp.write_text(json.dumps(descriptor, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, target)
    return str(p)
