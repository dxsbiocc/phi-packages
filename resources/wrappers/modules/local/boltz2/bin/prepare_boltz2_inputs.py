#!/usr/bin/env python3
"""Prepare Boltz-2 YAML input files from either:
  1. Existing Boltz YAML/FASTA files or directories (--inputs), OR
  2. A sequence-native TSV/CSV table (--sequences-tsv, e.g. from ligandmpnn)
     plus optional target sequence and ligand CCD/SMILES.

Pure-Python implementation (no third-party YAML library required).
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import string
from pathlib import Path
from typing import Any

BOLTZ_INPUT_EXTS = {".yaml", ".yml", ".fasta", ".fa"}


def _yaml_quote(text: str) -> str:
    escaped = text.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def _write_boltz_yaml(
    out_path: Path,
    protein_chains: list[tuple[str, str]],
    ligand_ccd: str,
    ligand_smiles: str,
    predict_affinity: bool,
    msa_mode: str,
) -> None:
    lines: list[str] = ["version: 1", "sequences:"]
    for chain_id, seq in protein_chains:
        lines.append("  - protein:")
        lines.append(f"      id: {chain_id}")
        lines.append(f"      sequence: {seq}")
        if msa_mode == "empty":
            lines.append("      msa: empty")

    ligand_chain_id: str | None = None
    if ligand_ccd or ligand_smiles:
        used = {cid for cid, _ in protein_chains}
        ligand_chain_id = next((c for c in "LXYZMNOPQRSTUVW" if c not in used), "L")
        lines.append("  - ligand:")
        lines.append(f"      id: {ligand_chain_id}")
        if ligand_ccd:
            lines.append(f"      ccd: {ligand_ccd}")
        else:
            lines.append(f"      smiles: {_yaml_quote(ligand_smiles)}")

    if predict_affinity and ligand_chain_id is not None:
        lines.append("properties:")
        lines.append("  - affinity:")
        lines.append(f"      binder: {ligand_chain_id}")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _stage_existing_inputs(inputs: list[str], out_dir: Path) -> dict[str, dict[str, Any]]:
    manifest: dict[str, dict[str, Any]] = {}
    for raw in inputs:
        p = Path(raw)
        candidates = (
            [c for c in sorted(p.rglob("*")) if c.is_file() and c.suffix.lower() in BOLTZ_INPUT_EXTS]
            if p.is_dir()
            else ([p] if p.is_file() and p.suffix.lower() in BOLTZ_INPUT_EXTS else [])
        )
        for src in candidates:
            dst = out_dir / src.name
            shutil.copy2(src, dst)
            manifest[src.stem] = {"parent_design_id": src.stem, "source_input": str(src.resolve())}
    return manifest


def _stage_from_sequences_table(
    table_path: Path,
    out_dir: Path,
    target_sequence: str,
    ligand_ccd: str,
    ligand_smiles: str,
    predict_affinity: bool,
    msa_mode: str,
) -> dict[str, dict[str, Any]]:
    delimiter = "\t" if table_path.suffix.lower() in (".tsv", ".tab") else ","
    manifest: dict[str, dict[str, Any]] = {}
    chain_pool = list(string.ascii_uppercase)

    with table_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for idx, row in enumerate(reader, start=1):
            design_id = (row.get("design_id") or f"design_{idx}").strip()
            parent_id = (row.get("parent_design_id") or design_id).strip()
            raw_seq = (row.get("sequence") or "").strip()
            if not raw_seq:
                continue

            parts = [p.strip() for p in raw_seq.replace("/", ":").split(":") if p.strip()]
            if target_sequence.strip():
                parts.append(target_sequence.strip())

            protein_chains = [(chain_pool[i], seq) for i, seq in enumerate(parts)]
            yaml_path = out_dir / f"{design_id}.yaml"
            _write_boltz_yaml(
                out_path=yaml_path,
                protein_chains=protein_chains,
                ligand_ccd=ligand_ccd.strip(),
                ligand_smiles=ligand_smiles.strip(),
                predict_affinity=predict_affinity,
                msa_mode=msa_mode,
            )
            manifest[design_id] = {
                "parent_design_id": parent_id,
                "parent_structure_path": (row.get("parent_structure_path") or "").strip(),
            }
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare Boltz-2 prediction input YAML directory.")
    parser.add_argument("--inputs", nargs="*", default=[], help="Existing Boltz YAML/FASTA files or directories.")
    parser.add_argument("--sequences-tsv", default="", help="Optional sequence-native TSV/CSV table (e.g. from ligandmpnn).")
    parser.add_argument("--target-sequence", default="", help="Optional target protein sequence appended as second chain.")
    parser.add_argument("--ligand-ccd", default="", help="Optional ligand CCD code (e.g. ATP, HEM).")
    parser.add_argument("--ligand-smiles", default="", help="Optional ligand SMILES string.")
    parser.add_argument("--predict-affinity", action="store_true", help="Enable Boltz-2 ligand affinity property block.")
    parser.add_argument("--msa-mode", choices=("server", "empty"), default="empty", help="Set 'empty' to write msa: empty in YAML.")
    parser.add_argument("--out-dir", required=True, help="Output directory for staged Boltz-2 YAML files.")
    parser.add_argument("--out-manifest", required=True, help="Output lineage manifest JSON path.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, dict[str, Any]] = {}
    if args.sequences_tsv:
        manifest.update(
            _stage_from_sequences_table(
                table_path=Path(args.sequences_tsv),
                out_dir=out_dir,
                target_sequence=args.target_sequence,
                ligand_ccd=args.ligand_ccd,
                ligand_smiles=args.ligand_smiles,
                predict_affinity=args.predict_affinity,
                msa_mode=args.msa_mode,
            )
        )
    if args.inputs:
        manifest.update(_stage_existing_inputs(args.inputs, out_dir))

    if not manifest:
        raise SystemExit("No Boltz-2 input YAML/FASTA files or valid sequence rows were staged.")

    out_manifest = Path(args.out_manifest)
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    out_manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
