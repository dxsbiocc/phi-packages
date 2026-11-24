#!/usr/bin/env python3
"""Compute structural metrics for designed protein structures (PDB / mmCIF).

Pure-Python implementation (no third-party dependencies required) covering
core Chaperone analysis metrics:
  - Kabsch-aligned CA and backbone (N, CA, C, O) RMSD against a reference
  - Inter-chain interface heavy-atom contacts & interface residue counts
  - Cross-interface putative hydrogen-bond counts (N/O donor-acceptor <= cutoff)
  - Protein-ligand heavy-atom contacts (HETATM non-water)
  - Steric clash counts (inter-residue heavy atoms < clash_cutoff)
  - Radius of gyration (CA) and chain/residue composition
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

BACKBONE_ATOMS = {"N", "CA", "C", "O"}
POLAR_ELEMENTS = {"N", "O"}
WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O"}
STRUCTURE_EXTS = {".pdb", ".cif", ".mmcif"}


def _infer_element(atom_name: str, raw_element: str) -> str:
    elem = raw_element.strip().upper()
    if elem:
        return elem
    cleaned = "".join(ch for ch in atom_name.strip() if ch.isalpha())
    return cleaned[0].upper() if cleaned else ""


def parse_pdb(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            if len(line) < 54:
                continue
            altloc = line[16].strip()
            if altloc not in ("", "A", "1"):
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            if rec == "HETATM" and resname in WATER_RESNAMES:
                continue
            chain = line[21].strip() or "A"
            resseq = line[22:26].strip()
            icode = line[26].strip()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            raw_elem = line[76:78] if len(line) >= 78 else ""
            element = _infer_element(atom_name, raw_elem)
            if element == "H":
                continue
            atoms.append(
                {
                    "record": rec,
                    "atom_name": atom_name,
                    "resname": resname,
                    "chain": chain,
                    "res_id": f"{resseq}{icode}",
                    "x": x,
                    "y": y,
                    "z": z,
                    "element": element,
                }
            )
    return atoms


def _tokenize_cif_line(line: str) -> list[str]:
    tokens: list[str] = []
    i = 0
    n = len(line)
    while i < n:
        while i < n and line[i].isspace():
            i += 1
        if i >= n:
            break
        if line[i] == "#":
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


def parse_cif(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
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
                col_idx = {c: idx for idx, c in enumerate(cols)}
                group_i = col_idx.get("_atom_site.group_PDB")
                atom_i = col_idx.get("_atom_site.label_atom_id", col_idx.get("_atom_site.auth_atom_id"))
                comp_i = col_idx.get("_atom_site.label_comp_id", col_idx.get("_atom_site.auth_comp_id"))
                chain_i = col_idx.get("_atom_site.auth_asym_id", col_idx.get("_atom_site.label_asym_id"))
                seq_i = col_idx.get("_atom_site.auth_seq_id", col_idx.get("_atom_site.label_seq_id"))
                alt_i = col_idx.get("_atom_site.label_alt_id")
                x_i = col_idx.get("_atom_site.Cartn_x")
                y_i = col_idx.get("_atom_site.Cartn_y")
                z_i = col_idx.get("_atom_site.Cartn_z")
                elem_i = col_idx.get("_atom_site.type_symbol")
                model_i = col_idx.get("_atom_site.pdbx_PDB_model_num")
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
                    if model_i is not None:
                        m = toks[model_i]
                        if first_model is None:
                            first_model = m
                        elif m != first_model:
                            continue
                    rec = toks[group_i] if group_i is not None else "ATOM"
                    if rec not in ("ATOM", "HETATM"):
                        continue
                    if alt_i is not None and toks[alt_i] not in (".", "?", "A", "1"):
                        continue
                    resname = toks[comp_i] if comp_i is not None else "UNK"
                    if rec == "HETATM" and resname in WATER_RESNAMES:
                        continue
                    atom_name = toks[atom_i] if atom_i is not None else ""
                    chain = toks[chain_i] if chain_i is not None else "A"
                    res_id = toks[seq_i] if seq_i is not None else "0"
                    try:
                        x = float(toks[x_i])  # type: ignore[arg-type]
                        y = float(toks[y_i])  # type: ignore[arg-type]
                        z = float(toks[z_i])  # type: ignore[arg-type]
                    except (TypeError, ValueError):
                        continue
                    raw_elem = toks[elem_i] if elem_i is not None else ""
                    element = _infer_element(atom_name, raw_elem)
                    if element == "H":
                        continue
                    atoms.append(
                        {
                            "record": rec,
                            "atom_name": atom_name,
                            "resname": resname,
                            "chain": chain,
                            "res_id": res_id,
                            "x": x,
                            "y": y,
                            "z": z,
                            "element": element,
                        }
                    )
                break
        else:
            i += 1
    return atoms


def parse_structure(path: Path) -> list[dict[str, Any]]:
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        return parse_pdb(path)
    if suffix in (".cif", ".mmcif"):
        return parse_cif(path)
    raise ValueError(f"Unsupported structure format: {path}")


# --- Pure-Python 3x3 linear algebra & Kabsch RMSD ---

def _mat_mul_3x3(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def _mat_t_3x3(a: list[list[float]]) -> list[list[float]]:
    return [[a[j][i] for j in range(3)] for i in range(3)]


def _det_3x3(m: list[list[float]]) -> float:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def _jacobi_eigen_3x3(sym: list[list[float]], max_iter: int = 50) -> tuple[list[float], list[list[float]]]:
    a = [row[:] for row in sym]
    v = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    for _ in range(max_iter):
        # Find largest off-diagonal element
        p, q = 0, 1
        max_val = abs(a[0][1])
        for i, j in ((0, 2), (1, 2)):
            if abs(a[i][j]) > max_val:
                max_val = abs(a[i][j])
                p, q = i, j
        if max_val < 1e-12:
            break
        diff = a[q][q] - a[p][p]
        if abs(a[p][q]) < 1e-15 * abs(diff):
            t = a[p][q] / diff
        else:
            phi = diff / (2.0 * a[p][q])
            t = 1.0 / (abs(phi) + math.sqrt(phi * phi + 1.0))
            if phi < 0.0:
                t = -t
        c = 1.0 / math.sqrt(1.0 + t * t)
        s = t * c
        tau = s / (1.0 + c)
        apq = a[p][q]
        a[p][q] = a[q][p] = 0.0
        a[p][p] -= t * apq
        a[q][q] += t * apq
        for r in range(3):
            if r != p and r != q:
                arp, arq = a[r][p], a[r][q]
                a[r][p] = a[p][r] = arp - s * (arq + tau * arp)
                a[r][q] = a[q][r] = arq + s * (arp - tau * arq)
        for r in range(3):
            vrp, vrq = v[r][p], v[r][q]
            v[r][p] = vrp - s * (vrq + tau * vrp)
            v[r][q] = vrq + s * (vrp - tau * vrq)
    evals = [a[0][0], a[1][1], a[2][2]]
    return evals, v


def kabsch_rmsd(coords_p: list[tuple[float, float, float]], coords_q: list[tuple[float, float, float]]) -> float | None:
    """Compute optimal Kabsch RMSD between two equal-length 3D point lists."""
    n = len(coords_p)
    if n == 0 or n != len(coords_q):
        return None
    cx_p = sum(c[0] for c in coords_p) / n
    cy_p = sum(c[1] for c in coords_p) / n
    cz_p = sum(c[2] for c in coords_p) / n
    cx_q = sum(c[0] for c in coords_q) / n
    cy_q = sum(c[1] for c in coords_q) / n
    cz_q = sum(c[2] for c in coords_q) / n

    p = [(c[0] - cx_p, c[1] - cy_p, c[2] - cz_p) for c in coords_p]
    q = [(c[0] - cx_q, c[1] - cy_q, c[2] - cz_q) for c in coords_q]

    # Covariance H = P^T Q
    h = [[0.0] * 3 for _ in range(3)]
    for (px, py, pz), (qx, qy, qz) in zip(p, q):
        h[0][0] += px * qx
        h[0][1] += px * qy
        h[0][2] += px * qz
        h[1][0] += py * qx
        h[1][1] += py * qy
        h[1][2] += py * qz
        h[2][0] += pz * qx
        h[2][1] += pz * qy
        h[2][2] += pz * qz

    # Use Kabsch eigenvalue identity: RMSD^2 = (E0 - 2 * sum(s_k * sigma_k)) / N
    # where sigma_k = sqrt(eigenvalues of H^T H), and s_3 = sign(det(H)).
    e0 = sum(px * px + py * py + pz * pz for px, py, pz in p) + sum(
        qx * qx + qy * qy + qz * qz for qx, qy, qz in q
    )
    hth = _mat_mul_3x3(_mat_t_3x3(h), h)
    evals, _ = _jacobi_eigen_3x3(hth)
    evals.sort(reverse=True)
    sigmas = [math.sqrt(max(0.0, ev)) for ev in evals]
    sign = -1.0 if _det_3x3(h) < 0.0 else 1.0
    msd = max(0.0, (e0 - 2.0 * (sigmas[0] + sigmas[1] + sign * sigmas[2])) / n)
    return round(math.sqrt(msd), 4)


def _extract_ordered_coords(
    atoms: list[dict[str, Any]],
    atom_filter: set[str],
    chains: set[str] | None = None,
) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    for a in atoms:
        if a["record"] != "ATOM":
            continue
        if a["atom_name"] not in atom_filter:
            continue
        if chains and a["chain"] not in chains:
            continue
        coords.append((a["x"], a["y"], a["z"]))
    return coords


def _dist_sq(a: dict[str, Any], b: dict[str, Any]) -> float:
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return dx * dx + dy * dy + dz * dz


def compute_metrics_for_structure(
    structure_path: Path,
    reference_atoms: list[dict[str, Any]] | None,
    reference_path: Path | None,
    binder_chains: set[str] | None,
    target_chains: set[str] | None,
    contact_cutoff: float,
    hbond_cutoff: float,
    clash_cutoff: float,
) -> dict[str, Any]:
    atoms = parse_structure(structure_path)
    protein_atoms = [a for a in atoms if a["record"] == "ATOM"]
    ligand_atoms = [a for a in atoms if a["record"] == "HETATM"]

    chain_residues: dict[str, set[str]] = {}
    for a in protein_atoms:
        chain_residues.setdefault(a["chain"], set()).add(a["res_id"])
    all_chains = sorted(chain_residues.keys())

    # Resolve binder / target chains automatically if unspecified
    resolved_binder: set[str]
    resolved_target: set[str]
    if binder_chains:
        resolved_binder = binder_chains & set(all_chains)
    elif all_chains:
        resolved_binder = {all_chains[0]}
    else:
        resolved_binder = set()

    if target_chains:
        resolved_target = target_chains & set(all_chains)
    else:
        resolved_target = set(all_chains) - resolved_binder

    binder_atoms = [a for a in protein_atoms if a["chain"] in resolved_binder]
    target_atoms = [a for a in protein_atoms if a["chain"] in resolved_target]

    # Interface contacts & putative H-bonds
    contact_sq = contact_cutoff * contact_cutoff
    hbond_sq = hbond_cutoff * hbond_cutoff
    interface_contacts = 0
    interface_hbonds = 0
    iface_res_binder: set[tuple[str, str]] = set()
    iface_res_target: set[tuple[str, str]] = set()

    if binder_atoms and target_atoms:
        for ba in binder_atoms:
            b_polar = ba["element"] in POLAR_ELEMENTS
            for ta in target_atoms:
                d2 = _dist_sq(ba, ta)
                if d2 <= contact_sq:
                    interface_contacts += 1
                    iface_res_binder.add((ba["chain"], ba["res_id"]))
                    iface_res_target.add((ta["chain"], ta["res_id"]))
                if b_polar and ta["element"] in POLAR_ELEMENTS and d2 <= hbond_sq:
                    interface_hbonds += 1

    # Protein-ligand contacts
    ligand_contacts = 0
    if ligand_atoms and protein_atoms:
        for pa in protein_atoms:
            for la in ligand_atoms:
                if _dist_sq(pa, la) <= contact_sq:
                    ligand_contacts += 1

    # Inter-residue steric clashes (using simple spatial grid hashing for O(N) speed)
    clash_sq = clash_cutoff * clash_cutoff
    cell_size = max(clash_cutoff, 1.0)
    grid: dict[tuple[int, int, int], list[int]] = {}
    for idx, a in enumerate(protein_atoms):
        key = (
            int(math.floor(a["x"] / cell_size)),
            int(math.floor(a["y"] / cell_size)),
            int(math.floor(a["z"] / cell_size)),
        )
        grid.setdefault(key, []).append(idx)

    clash_count = 0
    for idx, a in enumerate(protein_atoms):
        cx = int(math.floor(a["x"] / cell_size))
        cy = int(math.floor(a["y"] / cell_size))
        cz = int(math.floor(a["z"] / cell_size))
        res_key_a = (a["chain"], a["res_id"])
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((cx + dx, cy + dy, cz + dz), ()):
                        if j <= idx:
                            continue
                        b = protein_atoms[j]
                        if (b["chain"], b["res_id"]) == res_key_a:
                            continue
                        # Skip covalent peptide bond C-N between adjacent residues
                        if a["chain"] == b["chain"] and {a["atom_name"], b["atom_name"]} == {"C", "N"}:
                            continue
                        if _dist_sq(a, b) < clash_sq:
                            clash_count += 1

    # Radius of gyration on binder CA atoms
    binder_ca = _extract_ordered_coords(protein_atoms, {"CA"}, resolved_binder if resolved_binder else None)
    if binder_ca:
        mx = sum(c[0] for c in binder_ca) / len(binder_ca)
        my = sum(c[1] for c in binder_ca) / len(binder_ca)
        mz = sum(c[2] for c in binder_ca) / len(binder_ca)
        rg = round(
            math.sqrt(
                sum((c[0] - mx) ** 2 + (c[1] - my) ** 2 + (c[2] - mz) ** 2 for c in binder_ca)
                / len(binder_ca)
            ),
            4,
        )
    else:
        rg = None

    # Reference Kabsch RMSD (CA & backbone)
    ca_rmsd = None
    backbone_rmsd = None
    if reference_atoms is not None:
        ref_ca = _extract_ordered_coords(reference_atoms, {"CA"}, resolved_binder if binder_chains else None)
        des_ca = _extract_ordered_coords(protein_atoms, {"CA"}, resolved_binder if binder_chains else None)
        if len(ref_ca) == len(des_ca) and len(ref_ca) > 0:
            ca_rmsd = kabsch_rmsd(des_ca, ref_ca)
        ref_bb = _extract_ordered_coords(reference_atoms, BACKBONE_ATOMS, resolved_binder if binder_chains else None)
        des_bb = _extract_ordered_coords(protein_atoms, BACKBONE_ATOMS, resolved_binder if binder_chains else None)
        if len(ref_bb) == len(des_bb) and len(ref_bb) > 0:
            backbone_rmsd = kabsch_rmsd(des_bb, ref_bb)

    return {
        "design_id": structure_path.stem,
        "structure_path": str(structure_path),
        "reference_path": str(reference_path) if reference_path else "",
        "num_chains": len(all_chains),
        "num_residues": sum(len(v) for v in chain_residues.values()),
        "num_atoms": len(protein_atoms),
        "binder_chains": ",".join(sorted(resolved_binder)),
        "target_chains": ",".join(sorted(resolved_target)),
        "binder_length": sum(len(chain_residues[c]) for c in resolved_binder),
        "target_length": sum(len(chain_residues[c]) for c in resolved_target),
        "ca_rmsd": ca_rmsd,
        "backbone_rmsd": backbone_rmsd,
        "interface_contacts": interface_contacts,
        "interface_residues_binder": len(iface_res_binder),
        "interface_residues_target": len(iface_res_target),
        "interface_hbonds": interface_hbonds,
        "ligand_contacts": ligand_contacts,
        "clash_count": clash_count,
        "radius_of_gyration": rg,
    }


def collect_structures(inputs: list[str]) -> list[Path]:
    found: list[Path] = []
    for raw in inputs:
        p = Path(raw)
        if p.is_dir():
            for child in sorted(p.rglob("*")):
                if child.is_file() and child.suffix.lower() in STRUCTURE_EXTS:
                    found.append(child)
        elif p.is_file() and p.suffix.lower() in STRUCTURE_EXTS:
            found.append(p)
    return found


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute protein design structural metrics.")
    parser.add_argument("--structures", nargs="+", required=True, help="Input PDB/CIF files or directories.")
    parser.add_argument("--reference", default=None, help="Optional reference PDB/CIF for Kabsch RMSD.")
    parser.add_argument("--binder-chains", default="", help="Comma-separated binder chain IDs (default: first chain).")
    parser.add_argument("--target-chains", default="", help="Comma-separated target chain IDs (default: remaining chains).")
    parser.add_argument("--contact-cutoff", type=float, default=4.5, help="Heavy-atom interface contact cutoff (Angstrom).")
    parser.add_argument("--hbond-cutoff", type=float, default=3.5, help="Putative polar H-bond distance cutoff (Angstrom).")
    parser.add_argument("--clash-cutoff", type=float, default=2.0, help="Inter-residue heavy-atom clash cutoff (Angstrom).")
    parser.add_argument("--out-tsv", required=True, help="Output TSV metrics table path.")
    parser.add_argument("--out-json", required=True, help="Output JSON metrics path.")
    args = parser.parse_args()

    structures = collect_structures(args.structures)
    if not structures:
        raise SystemExit("No PDB or mmCIF structure files found in --structures.")

    binder_chains = {c.strip() for c in args.binder_chains.split(",") if c.strip()} or None
    target_chains = {c.strip() for c in args.target_chains.split(",") if c.strip()} or None

    ref_path = Path(args.reference) if args.reference else None
    ref_atoms = parse_structure(ref_path) if ref_path else None

    rows = [
        compute_metrics_for_structure(
            structure_path=s,
            reference_atoms=ref_atoms,
            reference_path=ref_path,
            binder_chains=binder_chains,
            target_chains=target_chains,
            contact_cutoff=args.contact_cutoff,
            hbond_cutoff=args.hbond_cutoff,
            clash_cutoff=args.clash_cutoff,
        )
        for s in structures
    ]

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})

    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
