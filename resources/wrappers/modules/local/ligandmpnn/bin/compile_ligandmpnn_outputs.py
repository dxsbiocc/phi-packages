#!/usr/bin/env python3
"""Compile LigandMPNN / ProteinMPNN FASTA and optional packed PDB outputs into
standardized sequence-native lineage tables (TSV, FASTA, JSON).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any

KV_PATTERN = re.compile(r"([A-Za-z0-9_]+)\s*=\s*([^\s,]+)")


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def _parse_header_kv(header: str) -> dict[str, str]:
    return {m.group(1): m.group(2) for m in KV_PATTERN.finditer(header)}


def _maybe_float(val: str | None) -> float | None:
    if val is None:
        return None
    try:
        return float(val)
    except ValueError:
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile LigandMPNN outputs into standardized lineage tables.")
    parser.add_argument("--out-dir", required=True, help="LigandMPNN --out_folder directory.")
    parser.add_argument("--manifest", default="", help="Optional parent lineage manifest JSON.")
    parser.add_argument("--model-type", default="ligand_mpnn", help="MPNN model_type used.")
    parser.add_argument("--step-name", default="ligandmpnn", help="Pipeline step name.")
    parser.add_argument("--out-tsv", required=True, help="Output standardized sequences TSV path.")
    parser.add_argument("--out-fasta", required=True, help="Output merged FASTA path.")
    parser.add_argument("--out-json", required=True, help="Output standardized JSON path.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    manifest: dict[str, str] = {}
    if args.manifest and Path(args.manifest).is_file():
        manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))

    seqs_dir = out_dir / "seqs"
    fa_files = sorted(seqs_dir.glob("*.fa")) + sorted(seqs_dir.glob("*.fasta")) if seqs_dir.is_dir() else []
    if not fa_files:
        fa_files = sorted(out_dir.rglob("*.fa")) + sorted(out_dir.rglob("*.fasta"))
    if not fa_files:
        raise SystemExit(f"No LigandMPNN FASTA files found under {out_dir}")

    packed_dir = out_dir / "packed"
    backbones_dir = out_dir / "backbones"

    rows: list[dict[str, Any]] = []
    fasta_lines: list[str] = []

    for fa in fa_files:
        parent_id = fa.stem
        parent_struct = manifest.get(parent_id, "")
        records = _parse_fasta(fa)
        sample_idx = 0
        for rec_i, (header, seq) in enumerate(records):
            kv = _parse_header_kv(header)
            # LigandMPNN / ProteinMPNN writes the native input sequence first (no sample id/T)
            if rec_i == 0 and "id" not in kv and "sample" not in kv and len(records) > 1:
                continue
            sample_idx += 1
            sid = kv.get("id") or kv.get("sample") or str(sample_idx)
            design_id = f"{parent_id}_mpnn_{sid}"

            # Locate optional packed or backbone PDB if output_mode=structures
            struct_path = ""
            for cand_dir in (packed_dir, backbones_dir):
                if cand_dir.is_dir():
                    for pattern in (f"{parent_id}_packed_{sid}_*.pdb", f"{parent_id}_{sid}.pdb"):
                        matches = sorted(cand_dir.glob(pattern))
                        if matches:
                            struct_path = str(matches[0])
                            break
                if struct_path:
                    break

            overall_conf = _maybe_float(kv.get("overall_confidence"))
            ligand_conf = _maybe_float(kv.get("ligand_confidence"))
            score = _maybe_float(kv.get("score"))
            global_score = _maybe_float(kv.get("global_score"))
            seq_rec = _maybe_float(kv.get("seq_rec") or kv.get("seq_recovery"))
            temp = _maybe_float(kv.get("T") or kv.get("temperature"))
            seed = kv.get("seed", "")

            rows.append(
                {
                    "design_id": design_id,
                    "parent_design_id": parent_id,
                    "step": args.step_name,
                    "model_type": args.model_type,
                    "sequence": seq,
                    "sequence_length": len(seq.replace(":", "").replace("/", "")),
                    "structure_path": struct_path,
                    "parent_structure_path": parent_struct,
                    "overall_confidence": overall_conf,
                    "ligand_confidence": ligand_conf,
                    "score": score,
                    "global_score": global_score,
                    "seq_recovery": seq_rec,
                    "temperature": temp,
                    "seed": seed,
                }
            )
            fasta_lines.append(f">{design_id} parent={parent_id} model={args.model_type}")
            fasta_lines.append(seq)

    if not rows:
        raise SystemExit(f"No designed sequences parsed from {out_dir}")

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})

    out_fasta = Path(args.out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    out_fasta.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")

    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
