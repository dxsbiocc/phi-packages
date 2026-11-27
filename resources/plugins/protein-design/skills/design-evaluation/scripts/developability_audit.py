#!/usr/bin/env python3
"""Pure-stdlib pre-wetlab biophysical and chemical-liability audit for designed proteins.

Accepts either a protein FASTA (.fa/.fasta) or a DESIGN_FILTER TSV/CSV table
(filtered.tsv / merged.tsv) and writes developability_audit.tsv + its
.phi-artifact.json sidecar into output_dir.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# Average & monoisotopic residue masses (Da), plus H2O = 18.01528 / 18.010564684
AVG_MASS: dict[str, float] = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "Q": 128.1307, "E": 129.1155, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
MONO_MASS: dict[str, float] = {
    "A": 71.037114, "R": 156.101111, "N": 114.042927, "D": 115.026943, "C": 103.009185,
    "Q": 128.058578, "E": 129.042593, "G": 57.021464, "H": 137.058912, "I": 113.084064,
    "L": 113.084064, "K": 128.094963, "M": 131.040485, "F": 147.068414, "P": 97.052764,
    "S": 87.032028, "T": 101.047679, "W": 186.079313, "Y": 163.063329, "V": 99.068414,
}
H2O_AVG = 18.01528
H2O_MONO = 18.010564684

# Kyte-Doolittle (1982) hydropathy scale
KD_HYDRO: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# EMBOSS / Bjellqvist-style pKa values
PKA_NTERM = 8.6
PKA_CTERM = 3.6
PKA_POS = {"K": 10.8, "R": 12.5, "H": 6.5}
PKA_NEG = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}

# Guruprasad, Reddy & Pandit (1990) Dipeptide Instability Weight Values (DIWV)
DIWV: dict[str, dict[str, float]] = {
    "W": {"W": 1.0, "C": 1.0, "M": 24.68, "H": 24.68, "Y": 1.0, "F": 1.0, "Q": 1.0, "N": 13.34, "I": 1.0, "R": 1.0, "D": 1.0, "P": 1.0, "T": -14.03, "K": 1.0, "E": 1.0, "V": -7.49, "S": 1.0, "G": -9.37, "A": -14.03, "L": 13.34},
    "C": {"W": 24.68, "C": 1.0, "M": 33.60, "H": 33.60, "Y": 1.0, "F": 1.0, "Q": -6.54, "N": 1.0, "I": 1.0, "R": 1.0, "D": 20.26, "P": 20.26, "T": 33.60, "K": 1.0, "E": 1.0, "V": -6.54, "S": 1.0, "G": 1.0, "A": 1.0, "L": 20.26},
    "M": {"W": 1.0, "C": 1.0, "M": -1.88, "H": 58.28, "Y": 24.68, "F": 1.0, "Q": -6.54, "N": 1.0, "I": 1.0, "R": -6.54, "D": 1.0, "P": 44.94, "T": -1.88, "K": 1.0, "E": 1.0, "V": 1.0, "S": 44.94, "G": 1.0, "A": 13.34, "L": 1.0},
    "H": {"W": -1.88, "C": 1.0, "M": 1.0, "H": 1.0, "Y": 44.94, "F": -9.37, "Q": 1.0, "N": 24.68, "I": 44.94, "R": 1.0, "D": 1.0, "P": -1.88, "T": -6.54, "K": 24.68, "E": 1.0, "V": 1.0, "S": 1.0, "G": -9.37, "A": 1.0, "L": 1.0},
    "Y": {"W": -9.37, "C": 1.0, "M": 44.94, "H": 13.34, "Y": 13.34, "F": 1.0, "Q": 1.0, "N": 1.0, "I": 1.0, "R": -15.91, "D": 24.68, "P": 13.34, "T": -7.49, "K": 1.0, "E": -6.54, "V": 1.0, "S": 1.0, "G": -7.49, "A": 24.68, "L": 1.0},
    "F": {"W": 1.0, "C": 1.0, "M": 1.0, "H": 1.0, "Y": 33.60, "F": 1.0, "Q": 1.0, "N": 1.0, "I": 1.0, "R": 1.0, "D": 13.34, "P": 20.26, "T": 1.0, "K": -14.03, "E": 1.0, "V": 1.0, "S": 1.0, "G": 1.0, "A": 1.0, "L": 1.0},
    "Q": {"W": 1.0, "C": -6.54, "M": 1.0, "H": 1.0, "Y": -6.54, "F": -6.54, "Q": 20.26, "N": 1.0, "I": 1.0, "R": 1.0, "D": 20.26, "P": 20.26, "T": 1.0, "K": 1.0, "E": 20.26, "V": -6.54, "S": 44.94, "G": 1.0, "A": 1.0, "L": 1.0},
    "N": {"W": -9.37, "C": -1.88, "M": 1.0, "H": 1.0, "Y": 1.0, "F": -14.03, "Q": -6.54, "N": 1.0, "I": 44.94, "R": 1.0, "D": 1.0, "P": -1.88, "T": -7.49, "K": 24.68, "E": 1.0, "V": 1.0, "S": 1.0, "G": -14.03, "A": 1.0, "L": 1.0},
    "I": {"W": 1.0, "C": 1.0, "M": 1.0, "H": 13.34, "Y": 1.0, "F": 1.0, "Q": 1.0, "N": 1.0, "I": 1.0, "R": 1.0, "D": 1.0, "P": -1.88, "T": 1.0, "K": -7.49, "E": 44.94, "V": -7.49, "S": 1.0, "G": 1.0, "A": 1.0, "L": 20.26},
    "R": {"W": 58.28, "C": 1.0, "M": 1.0, "H": 20.26, "Y": -6.54, "F": 1.0, "Q": 20.26, "N": 13.34, "I": 1.0, "R": 58.28, "D": 1.0, "P": 20.26, "T": 1.0, "K": 1.0, "E": 1.0, "V": 1.0, "S": 44.94, "G": -7.49, "A": 1.0, "L": 1.0},
    "D": {"W": 1.0, "C": 1.0, "M": 1.0, "H": 1.0, "Y": 1.0, "F": -6.54, "Q": 1.0, "N": 1.0, "I": 1.0, "R": -6.54, "D": 1.0, "P": 1.0, "T": -14.03, "K": -7.49, "E": 1.0, "V": 1.0, "S": 20.26, "G": 1.0, "A": 1.0, "L": 1.0},
    "P": {"W": -1.88, "C": -6.54, "M": -6.54, "H": 1.0, "Y": 1.0, "F": 20.26, "Q": 20.26, "N": 1.0, "I": 1.0, "R": -6.54, "D": -6.54, "P": 20.26, "T": 1.0, "K": 1.0, "E": 18.38, "V": 20.26, "S": 20.26, "G": 1.0, "A": 20.26, "L": 1.0},
    "T": {"W": -14.03, "C": 1.0, "M": 1.0, "H": 1.0, "Y": 1.0, "F": 13.34, "Q": -6.54, "N": -14.03, "I": 1.0, "R": 1.0, "D": 1.0, "P": 1.0, "T": 1.0, "K": 1.0, "E": 20.26, "V": 1.0, "S": 1.0, "G": -7.49, "A": 1.0, "L": 1.0},
    "K": {"W": 1.0, "C": 1.0, "M": 33.60, "H": 1.0, "Y": 1.0, "F": 1.0, "Q": 24.68, "N": 1.0, "I": -7.49, "R": 33.60, "D": 1.0, "P": -6.54, "T": 1.0, "K": 1.0, "E": 1.0, "V": -7.49, "S": 1.0, "G": -7.49, "A": 1.0, "L": -7.49},
    "E": {"W": -14.03, "C": 44.94, "M": 1.0, "H": -6.54, "Y": 1.0, "F": 1.0, "Q": 20.26, "N": 1.0, "I": 20.26, "R": 1.0, "D": 20.26, "P": 20.26, "T": 1.0, "K": 1.0, "E": 33.60, "V": 1.0, "S": 20.26, "G": 1.0, "A": 1.0, "L": 1.0},
    "V": {"W": 1.0, "C": 1.0, "M": 1.0, "H": 1.0, "Y": -6.54, "F": 1.0, "Q": 1.0, "N": 1.0, "I": 1.0, "R": 1.0, "D": -14.03, "P": 20.26, "T": -7.49, "K": -1.88, "E": 1.0, "V": 1.0, "S": 1.0, "G": -7.49, "A": 1.0, "L": 1.0},
    "S": {"W": 1.0, "C": 33.60, "M": 1.0, "H": 1.0, "Y": 1.0, "F": 1.0, "Q": 20.26, "N": 1.0, "I": 1.0, "R": 20.26, "D": 1.0, "P": 44.94, "T": 1.0, "K": 1.0, "E": 20.26, "V": 1.0, "S": 20.26, "G": 1.0, "A": 1.0, "L": 1.0},
    "G": {"W": 13.34, "C": 1.0, "M": 1.0, "H": 1.0, "Y": -7.49, "F": 1.0, "Q": 1.0, "N": -7.49, "I": -7.49, "R": 1.0, "D": 1.0, "P": 1.0, "T": -7.49, "K": -7.49, "E": -6.54, "V": 1.0, "S": 1.0, "G": 13.34, "A": -7.49, "L": 1.0},
    "A": {"W": 1.0, "C": 44.94, "M": 1.0, "H": -7.49, "Y": 1.0, "F": 1.0, "Q": 1.0, "N": 1.0, "I": 1.0, "R": 1.0, "D": -7.49, "P": 20.26, "T": 1.0, "K": 1.0, "E": 1.0, "V": 1.0, "S": 1.0, "G": 1.0, "A": 1.0, "L": 1.0},
    "L": {"W": 24.68, "C": 1.0, "M": 1.0, "H": 1.0, "Y": 1.0, "F": 1.0, "Q": 33.60, "N": 1.0, "I": 1.0, "R": 20.26, "D": 1.0, "P": 20.26, "T": 1.0, "K": -7.49, "E": 1.0, "V": 1.0, "S": 1.0, "G": 1.0, "A": 1.0, "L": 1.0},
}


def _fail(msg: str) -> None:
    sys.stdout.write(json.dumps({"error": msg}, indent=2) + "\n")
    sys.stdout.flush()
    raise SystemExit(1)


def _net_charge(counts: dict[str, int], ph: float) -> float:
    q = 1.0 / (1.0 + 10.0 ** (ph - PKA_NTERM))
    q -= 1.0 / (1.0 + 10.0 ** (PKA_CTERM - ph))
    for aa, pka in PKA_POS.items():
        n = counts.get(aa, 0)
        if n:
            q += n / (1.0 + 10.0 ** (ph - pka))
    for aa, pka in PKA_NEG.items():
        n = counts.get(aa, 0)
        if n:
            q -= n / (1.0 + 10.0 ** (pka - ph))
    return q


def _isoelectric_point(counts: dict[str, int]) -> float:
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _net_charge(counts, mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _parse_inputs(
    path: Path, seq_col_req: str, id_col_req: str
) -> tuple[str, list[tuple[str, str]]]:
    if not path.is_file():
        _fail(f"sequences_input not found: {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    stripped = text.lstrip()
    ext = path.suffix.lower()
    is_fasta = ext in (".fa", ".fasta", ".faa", ".fna") or stripped.startswith(">")

    entries: list[tuple[str, str]] = []
    if is_fasta:
        cur_id: str | None = None
        chunks: list[str] = []
        for raw in text.splitlines():
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    entries.append((cur_id, "".join(chunks)))
                cur_id = line[1:].split()[0] or f"seq_{len(entries) + 1}"
                chunks = []
            else:
                chunks.append(line)
        if cur_id is not None:
            entries.append((cur_id, "".join(chunks)))
        if not entries:
            _fail(f"No FASTA records found in {path}")
        return "fasta", entries

    sep = "," if ext == ".csv" else "\t"
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        _fail(f"Empty input table: {path}")
    if "\t" not in lines[0] and "," in lines[0]:
        sep = ","
    reader = csv.DictReader(lines, delimiter=sep)
    cols = list(reader.fieldnames or [])
    if not cols:
        _fail(f"No header columns found in {path}")

    lower_map = {c.lower(): c for c in cols}
    seq_col = None
    if seq_col_req:
        if seq_col_req not in cols:
            _fail(f"Requested seq_col '{seq_col_req}' not in columns: {cols}")
        seq_col = seq_col_req
    else:
        for cand in ("sequence", "binder_seq", "seq", "designed_sequence", "aa_seq", "protein_seq"):
            if cand in lower_map:
                seq_col = lower_map[cand]
                break
        if seq_col is None:
            _fail(f"Could not auto-detect sequence column in {cols}; pass --seq_col.")

    id_col = None
    if id_col_req:
        if id_col_req not in cols:
            _fail(f"Requested id_col '{id_col_req}' not in columns: {cols}")
        id_col = id_col_req
    else:
        for cand in ("design_id", "id", "name", "seq_id", "binder_id", "candidate_id"):
            if cand in lower_map:
                id_col = lower_map[cand]
                break

    for idx, row in enumerate(reader, start=1):
        cid = str(row.get(id_col, "")).strip() if id_col else f"seq_{idx}"
        if not cid:
            cid = f"seq_{idx}"
        seq = str(row.get(seq_col, "")).strip()
        if seq:
            entries.append((cid, seq))
    if not entries:
        _fail(f"No sequences extracted from table {path}")
    return "table", entries


def _find_motifs(seq: str) -> dict[str, Any]:
    deam_pos = [m.start() + 1 for m in re.finditer(r"(?=(N[GST]))", seq)]
    isom_pos = [m.start() + 1 for m in re.finditer(r"(?=(D[GS]))", seq)]
    ngly_pos = [m.start() + 1 for m in re.finditer(r"(?=(N[^P][ST]))", seq)]

    nterm_notes: list[str] = []
    if seq.startswith("M") and len(seq) >= 2:
        if seq[1] in "GASPTVC":
            nterm_notes.append(f"M{seq[1]}_met_excision_likely")
        if seq[1] == "P":
            nterm_notes.append("MP_nterm_proline")
    elif seq.startswith("P"):
        nterm_notes.append("P_nterm_proline")
    return {
        "deamidation_positions": deam_pos,
        "isomerization_positions": isom_pos,
        "nglycosylation_positions": ngly_pos,
        "nterm_notes": nterm_notes,
    }


def _audit_sequence(
    cid: str, raw_seq: str, buffer_ph: float, max_hydrophobic_patch: float
) -> dict[str, Any]:
    seq = re.sub(r"\s+", "", raw_seq.upper()).rstrip("*")
    if not seq:
        _fail(f"Empty sequence for candidate '{cid}'")

    non_std = sorted({ch for ch in seq if ch not in AVG_MASS})
    std_seq = "".join(ch for ch in seq if ch in AVG_MASS)
    if not std_seq:
        _fail(f"Candidate '{cid}' contains no standard amino acids: {raw_seq[:30]}")

    L = len(std_seq)
    counts: dict[str, int] = {}
    for ch in std_seq:
        counts[ch] = counts.get(ch, 0) + 1

    mw_avg = sum(AVG_MASS[ch] for ch in std_seq) + H2O_AVG
    mw_mono = sum(MONO_MASS[ch] for ch in std_seq) + H2O_MONO
    pi = _isoelectric_point(counts)
    charge_ph = _net_charge(counts, buffer_ph)
    pi_clash = abs(pi - buffer_ph) < 1.0

    n_w = counts.get("W", 0)
    n_y = counts.get("Y", 0)
    n_c = counts.get("C", 0)
    n_ss = n_c // 2
    eps_red = 5500 * n_w + 1490 * n_y
    eps_ox = eps_red + 125 * n_ss
    abs01_red = eps_red / mw_avg if mw_avg > 0 else 0.0
    abs01_ox = eps_ox / mw_avg if mw_avg > 0 else 0.0
    uv280_blind = (n_w == 0 and n_y == 0)
    unpaired_cys = (n_c % 2 == 1)

    gravy = sum(KD_HYDRO[ch] for ch in std_seq) / L
    if L >= 2:
        ii = (10.0 / L) * sum(DIWV[std_seq[i]][std_seq[i + 1]] for i in range(L - 1))
    else:
        ii = 0.0
    unstable = ii > 40.0

    # Ikai (1980) Aliphatic Index
    x_a = 100.0 * counts.get("A", 0) / L
    x_v = 100.0 * counts.get("V", 0) / L
    x_i = 100.0 * counts.get("I", 0) / L
    x_l = 100.0 * counts.get("L", 0) / L
    aliphatic_idx = x_a + 2.9 * x_v + 3.9 * (x_i + x_l)

    # 9-aa sliding-window max Kyte-Doolittle hydrophobic aggregation patch
    win = min(9, L)
    kd_vals = [KD_HYDRO[ch] for ch in std_seq]
    max_patch = max(sum(kd_vals[i : i + win]) / win for i in range(L - win + 1))
    patch_risk = max_patch > max_hydrophobic_patch

    motifs = _find_motifs(std_seq)
    n_deam = len(motifs["deamidation_positions"])
    n_isom = len(motifs["isomerization_positions"])
    n_ngly = len(motifs["nglycosylation_positions"])

    high_flags: list[str] = []
    caution_flags: list[str] = []
    if uv280_blind:
        high_flags.append("uv280_blind")
    if unpaired_cys:
        high_flags.append("unpaired_cys")
    if pi_clash:
        high_flags.append("pi_buffer_clash")
    if patch_risk:
        high_flags.append("hydrophobic_patch_risk")
    if ii > 45.0:
        high_flags.append("instability_index_high")
    elif unstable:
        caution_flags.append("instability_index_gt40")
    if n_deam > 0:
        caution_flags.append(f"deamidation_NG_NS_NT({n_deam})")
    if n_isom > 0:
        caution_flags.append(f"isomerization_DG_DS({n_isom})")
    if n_ngly > 0:
        caution_flags.append(f"nglycosylation_NXST({n_ngly})")
    if non_std:
        caution_flags.append("nonstandard_aa:" + "".join(non_std))

    if high_flags:
        tier = "high_risk"
    elif caution_flags:
        tier = "caution"
    else:
        tier = "pass"

    all_flags = high_flags + caution_flags
    return {
        "design_id": cid,
        "length_aa": L,
        "mw_avg_da": round(mw_avg, 2),
        "mw_mono_da": round(mw_mono, 4),
        "mw_kda": round(mw_avg / 1000.0, 3),
        "pI": round(pi, 3),
        "net_charge_at_ph": round(charge_ph, 3),
        "pi_buffer_clash": pi_clash,
        "n_trp": n_w,
        "n_tyr": n_y,
        "n_cys": n_c,
        "unpaired_cys": unpaired_cys,
        "eps280_reduced": eps_red,
        "eps280_oxidized": eps_ox,
        "abs280_01pct_reduced": round(abs01_red, 4),
        "abs280_01pct_oxidized": round(abs01_ox, 4),
        "uv280_blind": uv280_blind,
        "gravy": round(gravy, 4),
        "instability_index": round(ii, 2),
        "is_unstable": unstable,
        "aliphatic_index": round(aliphatic_idx, 2),
        "max_9mer_kd": round(max_patch, 3),
        "hydrophobic_patch_risk": patch_risk,
        "n_deamidation_motifs": n_deam,
        "deamidation_positions": ",".join(map(str, motifs["deamidation_positions"])),
        "n_isomerization_motifs": n_isom,
        "isomerization_positions": ",".join(map(str, motifs["isomerization_positions"])),
        "n_nglycosylation_sequons": n_ngly,
        "nglycosylation_positions": ",".join(map(str, motifs["nglycosylation_positions"])),
        "nterm_notes": ";".join(motifs["nterm_notes"]) if motifs["nterm_notes"] else "none",
        "developability_tier": tier,
        "flags": ";".join(all_flags) if all_flags else "none",
    }


def _write_sidecar(
    tsv_path: Path,
    input_path: Path,
    parameters: dict[str, Any],
    rows_n: int,
    cols_n: int,
) -> None:
    cwd = Path.cwd().resolve()
    prov_inputs = []
    try:
        rel = input_path.resolve().relative_to(cwd)
        prov_inputs.append(
            {"path": str(rel), "sha256": hashlib.sha256(input_path.read_bytes()).hexdigest()}
        )
    except ValueError:
        pass

    descriptor: dict[str, Any] = {
        "contractVersion": "1.0.0",
        "kind": "table",
        "file": tsv_path.name,
        "mediaType": "text/tab-separated-values",
        "title": "Pre-Wetlab Protein Developability & Chemical Liability Audit",
        "table": {"rows": rows_n, "columns": cols_n},
        "provenance": {
            "createdAt": datetime.now(timezone.utc).isoformat(),
            "tool": "pd_developability_audit",
            "skill": "design-evaluation",
            "parameters": parameters,
        },
    }
    if prov_inputs:
        descriptor["provenance"]["inputs"] = prov_inputs

    target = tsv_path.with_name(tsv_path.name + ".phi-artifact.json")
    tmp = target.with_name(target.name + ".tmp")
    tmp.write_text(json.dumps(descriptor, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, target)


def main() -> None:
    ap = argparse.ArgumentParser(description="Pre-wetlab biophysical & liability audit for designed proteins")
    ap.add_argument("--sequences_input", "--sequences-input", dest="sequences_input", required=True)
    ap.add_argument("--output_dir", "--output-dir", dest="output_dir", required=True)
    ap.add_argument("--seq_col", "--seq-col", dest="seq_col", default="")
    ap.add_argument("--id_col", "--id-col", dest="id_col", default="")
    ap.add_argument("--buffer_ph", "--buffer-ph", dest="buffer_ph", type=float, default=7.4)
    ap.add_argument(
        "--max_hydrophobic_patch",
        "--max-hydrophobic-patch",
        dest="max_hydrophobic_patch",
        type=float,
        default=2.2,
    )
    args = ap.parse_args()

    if not (0.0 <= args.buffer_ph <= 14.0):
        _fail(f"buffer_ph must be in [0, 14], got {args.buffer_ph}")

    in_path = Path(args.sequences_input)
    input_fmt, entries = _parse_inputs(in_path, args.seq_col, args.id_col)

    rows = [
        _audit_sequence(cid, seq, args.buffer_ph, args.max_hydrophobic_patch)
        for cid, seq in entries
    ]

    tier_counts = {"pass": 0, "caution": 0, "high_risk": 0}
    for r in rows:
        tier_counts[r["developability_tier"]] += 1

    n_uv_blind = sum(1 for r in rows if r["uv280_blind"])
    n_pi_clash = sum(1 for r in rows if r["pi_buffer_clash"])
    n_unpaired_cys = sum(1 for r in rows if r["unpaired_cys"])
    n_instab = sum(1 for r in rows if r["is_unstable"])
    n_patch = sum(1 for r in rows if r["hydrophobic_patch_risk"])

    warnings: list[str] = []
    if n_uv_blind:
        warnings.append(f"{n_uv_blind} candidate(s) are UV280-blind (0 Trp, 0 Tyr); NanoDrop A280 will fail.")
    if n_pi_clash:
        warnings.append(f"{n_pi_clash} candidate(s) have |pI - {args.buffer_ph}| < 1.0 (precipitation risk).")
    if n_unpaired_cys:
        warnings.append(f"{n_unpaired_cys} candidate(s) carry odd/unpaired cysteines (thiol scrambling risk).")
    if n_patch:
        warnings.append(f"{n_patch} candidate(s) exceed 9-mer hydrophobic patch threshold ({args.max_hydrophobic_patch}).")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / "developability_audit.tsv"

    fieldnames = list(rows[0].keys())
    with tsv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    _write_sidecar(
        tsv_path,
        in_path,
        {
            "buffer_ph": args.buffer_ph,
            "max_hydrophobic_patch": args.max_hydrophobic_patch,
        },
        len(rows),
        len(fieldnames),
    )

    tier_order = {"pass": 0, "caution": 1, "high_risk": 2}
    ranked = sorted(rows, key=lambda r: (tier_order[r["developability_tier"]], r["instability_index"]))

    payload = {
        "status": "ok",
        "tool": "pd_developability_audit",
        "input_format": input_fmt,
        "buffer_ph": args.buffer_ph,
        "max_hydrophobic_patch": args.max_hydrophobic_patch,
        "n_sequences": len(rows),
        "tier_counts": tier_counts,
        "n_uv280_blind": n_uv_blind,
        "n_pi_clash": n_pi_clash,
        "n_unpaired_cys": n_unpaired_cys,
        "n_instability_flagged": n_instab,
        "n_hydrophobic_patch_flagged": n_patch,
        "top_candidates": ranked[:10],
        "outputs": {
            "developability_audit_tsv": str(tsv_path),
        },
        "artifacts": [str(tsv_path)],
        "warnings": warnings,
    }
    sys.stdout.write(json.dumps(payload, indent=2) + "\n")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
