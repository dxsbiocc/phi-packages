#!/usr/bin/env python3
"""rna_id_map: deterministic gene identifier cleaning (Ensembl version stripping),
GTF/GFF3 or custom table ID mapping, biotype filtering, and collision-aware
aggregation for bulk RNA-seq expression matrices and DE result tables.

Pure Python 3 stdlib implementation for fast, dependency-free execution.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any


VERSIONED_ID_RE = re.compile(r"^((?:ENS[A-Z]*[GT]\d+|NM_\d+|NR_\d+|XM_\d+|XR_\d+))\.\d+$", re.IGNORECASE)
DE_FC_CANDS = {"log2foldchange", "log2fc", "logfc", "lfc"}
DE_P_CANDS = {"pvalue", "pval", "p.value", "p_val", "padj", "fdr", "adj.p.val", "qvalue"}
DE_STAT_CANDS = {"stat", "t", "lr", "f", "wald"}


def emit(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)


def fail(message: str, code: int = 1) -> None:
    print(json.dumps({"error": message}, indent=2, ensure_ascii=False))
    sys.exit(code)


def open_text(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def strip_id_version(gene_id: str) -> tuple[str, bool]:
    g = gene_id.strip()
    m = VERSIONED_ID_RE.match(g)
    if m:
        return m.group(1), True
    return g, False


def write_table_artifact(
    path: str,
    title: str,
    *,
    inputs: tuple[str | None, ...],
    parameters: dict[str, Any],
    rows: int,
    columns: int,
) -> str:
    p = Path(path)
    ext = p.suffix.lower()
    media = "text/csv" if ext == ".csv" else "text/tab-separated-values"
    desc: dict[str, Any] = {
        "contractVersion": "1.0.0",
        "kind": "table",
        "file": p.name,
        "mediaType": media,
        "title": title[:200],
        "table": {"rows": int(rows), "columns": int(columns)},
        "provenance": {
            "createdAt": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "tool": "rna_id_map",
            "skill": "rnaseq-workflow",
            "parameters": parameters,
        },
    }
    cwd = Path.cwd().resolve()
    prov_in = []
    for raw in inputs:
        if not raw:
            continue
        ip = Path(raw)
        try:
            rel = ip.resolve().relative_to(cwd)
        except ValueError:
            continue
        entry: dict[str, Any] = {"path": str(rel)}
        if ip.is_file():
            entry["sha256"] = hashlib.sha256(ip.read_bytes()).hexdigest()
        prov_in.append(entry)
    if prov_in:
        desc["provenance"]["inputs"] = prov_in
    sidecar = p.with_name(p.name + ".phi-artifact.json")
    tmp = sidecar.with_name(sidecar.name + ".tmp")
    tmp.write_text(json.dumps(desc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, sidecar)
    return str(p)


def parse_gtf_or_gff(path: str, strip_ver: bool) -> tuple[dict[str, str], dict[str, str]]:
    """Return (id_to_symbol, id_to_biotype) parsed from GTF or GFF3."""
    id_to_sym: dict[str, str] = {}
    id_to_bio: dict[str, str] = {}
    gtf_attr_re = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]+)"')

    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in ("gene", "transcript", "mRNA", "exon"):
                continue
            attr_str = parts[8]
            attrs: dict[str, str] = {}
            if "=" in attr_str and '"' not in attr_str:
                # GFF3 key=value;...
                for item in attr_str.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        attrs[k.strip()] = v.strip()
            else:
                for k, v in gtf_attr_re.findall(attr_str):
                    attrs[k] = v.strip()

            gid = attrs.get("gene_id") or attrs.get("ID", "")
            if gid.startswith("gene:"):
                gid = gid[5:]
            if not gid:
                continue
            gid_clean = strip_id_version(gid)[0] if strip_ver else gid.strip()
            gname = attrs.get("gene_name") or attrs.get("Name") or attrs.get("gene") or ""
            gbio = attrs.get("gene_biotype") or attrs.get("gene_type") or attrs.get("biotype") or ""

            if gname and (gid_clean not in id_to_sym or feature == "gene"):
                id_to_sym[gid_clean] = gname
                id_to_sym[gid.strip()] = gname
            if gbio and (gid_clean not in id_to_bio or feature == "gene"):
                id_to_bio[gid_clean] = gbio
                id_to_bio[gid.strip()] = gbio

    if not id_to_sym:
        fail(f"no gene_id -> gene_name mappings found in GTF/GFF: {path}")
    return id_to_sym, id_to_bio


def parse_mapping_table(
    path: str, from_col: str | None, to_col: str | None, strip_ver: bool
) -> dict[str, str]:
    sep = "," if path.lower().endswith(".csv") else "\t"
    mapping: dict[str, str] = {}
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter=sep)
        header = next(reader, None)
        if not header or len(header) < 2:
            fail(f"mapping_table must have >= 2 columns: {path}")
        lower_map = {h.strip().lower(): idx for idx, h in enumerate(header)}
        f_idx = 0
        t_idx = 1
        if from_col:
            if from_col.lower() not in lower_map:
                fail(f"from_col '{from_col}' not in mapping_table header: {header}")
            f_idx = lower_map[from_col.lower()]
        if to_col:
            if to_col.lower() not in lower_map:
                fail(f"to_col '{to_col}' not in mapping_table header: {header}")
            t_idx = lower_map[to_col.lower()]

        for row in reader:
            if len(row) <= max(f_idx, t_idx):
                continue
            src = row[f_idx].strip()
            dst = row[t_idx].strip()
            if not src or not dst or dst.lower() in ("na", "nan", "null", "none", "."):
                continue
            mapping[src] = dst
            if strip_ver:
                mapping[strip_id_version(src)[0]] = dst
    if not mapping:
        fail(f"mapping_table yielded 0 valid mappings: {path}")
    return mapping


def to_float(val: str) -> float | None:
    s = val.strip()
    if not s or s.lower() in ("na", "nan", "null", "none"):
        return None
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except ValueError:
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Expression matrix or DE table TSV/CSV.")
    parser.add_argument("--output", required=True, help="Output mapped TSV/CSV path.")
    parser.add_argument("--gtf", default=None, help="Optional reference GTF/GFF3 (.gtf/.gtf.gz/.gff3).")
    parser.add_argument("--mapping_table", default=None, help="Optional 2-column mapping TSV/CSV.")
    parser.add_argument("--from_col", default=None, help="Source ID column in --mapping_table.")
    parser.add_argument("--to_col", default=None, help="Target ID column in --mapping_table.")
    parser.add_argument("--id_col", default=None, help="Gene ID column in --input (default first column).")
    parser.add_argument(
        "--keep_version",
        nargs="?",
        const=True,
        type=lambda x: str(x).lower() not in ("false", "0", "no"),
        default=False,
        help="Keep Ensembl/RefSeq '.version' suffixes (default false).",
    )
    parser.add_argument(
        "--strip_version",
        nargs="?",
        const=True,
        type=lambda x: str(x).lower() not in ("false", "0", "no"),
        default=None,
        help="Strip Ensembl/RefSeq '.version' suffixes before mapping (default true unless --keep_version).",
    )
    parser.add_argument("--table_kind", choices=["auto", "matrix", "de"], default="auto")
    parser.add_argument(
        "--agg_method",
        choices=["auto", "sum", "mean", "max_mean", "min_pvalue", "first"],
        default="auto",
    )
    parser.add_argument("--unmapped", choices=["keep", "drop"], default="keep")
    parser.add_argument("--biotype", default=None, help="Optional GTF gene_biotype filter (e.g. protein_coding).")
    args = parser.parse_args()
    if args.strip_version is None:
        args.strip_version = not args.keep_version

    if not os.path.isfile(args.input):
        fail(f"input file not found: {args.input}")
    if args.gtf and not os.path.isfile(args.gtf):
        fail(f"gtf file not found: {args.gtf}")
    if args.mapping_table and not os.path.isfile(args.mapping_table):
        fail(f"mapping_table not found: {args.mapping_table}")

    id_to_sym: dict[str, str] = {}
    id_to_bio: dict[str, str] = {}
    if args.gtf:
        id_to_sym, id_to_bio = parse_gtf_or_gff(args.gtf, args.strip_version)
    if args.mapping_table:
        custom_map = parse_mapping_table(args.mapping_table, args.from_col, args.to_col, args.strip_version)
        id_to_sym.update(custom_map)
    if args.biotype and not id_to_bio:
        fail("--biotype filtering requires --gtf with gene_biotype / gene_type attributes")

    in_sep = "," if args.input.lower().endswith(".csv") else "\t"
    out_sep = "," if args.output.lower().endswith(".csv") else "\t"

    with open_text(args.input) as fh:
        reader = csv.reader(fh, delimiter=in_sep)
        header = next(reader, None)
        if not header:
            fail(f"input table is empty: {args.input}")
        rows_raw = [r for r in reader if any(c.strip() for c in r)]

    if not rows_raw:
        fail(f"input table has 0 data rows: {args.input}")

    # Handle R-style leading unnamed column (len(header) == len(row) - 1)
    if len(header) == len(rows_raw[0]) - 1:
        header = ["gene_id"] + header
    elif not header[0].strip():
        header[0] = "gene_id"

    id_idx = 0
    if args.id_col:
        lower_h = {h.strip().lower(): i for i, h in enumerate(header)}
        if args.id_col.lower() not in lower_h:
            fail(f"id_col '{args.id_col}' not found in header: {header}")
        id_idx = lower_h[args.id_col.lower()]

    other_indices = [i for i in range(len(header)) if i != id_idx]
    other_names = [header[i].strip() for i in other_indices]
    lower_other = {c.lower(): idx for idx, c in enumerate(other_names)}

    # Detect table_kind
    has_fc = any(c in lower_other for c in DE_FC_CANDS)
    has_p = any(c in lower_other for c in DE_P_CANDS)
    detected_kind = args.table_kind
    if detected_kind == "auto":
        detected_kind = "de" if (has_fc and has_p) else "matrix"

    # Inspect numeric values for scale / agg_method selection
    is_integer_counts = True
    has_negative = False
    max_val = 0.0
    for r in rows_raw[:200]:
        for idx in other_indices:
            if idx >= len(r):
                continue
            v = to_float(r[idx])
            if v is None:
                continue
            if v < 0:
                has_negative = True
                is_integer_counts = False
            elif abs(v - round(v)) > 1e-6:
                is_integer_counts = False
            if v > max_val:
                max_val = v

    agg = args.agg_method
    if agg == "auto":
        if detected_kind == "de":
            agg = "min_pvalue"
        else:
            agg = "sum" if is_integer_counts else "max_mean"

    # Map and filter rows
    n_input_rows = len(rows_raw)
    n_versions_stripped = 0
    n_mapped = 0
    n_unmapped = 0
    n_biotype_filtered = 0
    grouped: dict[str, list[tuple[str, list[str]]]] = {}

    for r in rows_raw:
        if id_idx >= len(r):
            continue
        raw_id = r[id_idx].strip()
        if not raw_id:
            continue
        clean_id, stripped = strip_id_version(raw_id) if args.strip_version else (raw_id, False)
        if stripped:
            n_versions_stripped += 1

        if args.biotype:
            bio = id_to_bio.get(clean_id) or id_to_bio.get(raw_id)
            if bio != args.biotype:
                n_biotype_filtered += 1
                continue

        mapped_target = id_to_sym.get(clean_id) or id_to_sym.get(raw_id)
        if mapped_target:
            target_id = mapped_target
            n_mapped += 1
        else:
            n_unmapped += 1
            if args.unmapped == "drop" and (args.gtf or args.mapping_table):
                continue
            target_id = clean_id

        vals = [r[i] if i < len(r) else "" for i in other_indices]
        grouped.setdefault(target_id, []).append((raw_id, vals))

    if not grouped:
        fail("0 genes remained after ID mapping / biotype filtering")

    # Resolve P-value / stat indices for DE deduplication
    p_col_pos = next((lower_other[c] for c in ("pvalue", "pval", "p.value", "padj", "fdr") if c in lower_other), None)
    stat_col_pos = next((lower_other[c] for c in DE_STAT_CANDS if c in lower_other), None)

    out_rows: list[list[str]] = []
    collision_rows: list[list[str]] = []
    n_collisions = 0

    for target_id, items in grouped.items():
        if len(items) == 1:
            out_rows.append([target_id] + items[0][1])
            continue

        n_collisions += 1
        src_ids = [it[0] for it in items]
        collision_rows.append([target_id, str(len(items)), ";".join(src_ids)])

        if agg == "first":
            chosen = items[0][1]
        elif agg == "min_pvalue":
            def de_key(item: tuple[str, list[str]]) -> tuple[float, float]:
                vals = item[1]
                pv = to_float(vals[p_col_pos]) if p_col_pos is not None and p_col_pos < len(vals) else None
                st = to_float(vals[stat_col_pos]) if stat_col_pos is not None and stat_col_pos < len(vals) else None
                return (pv if pv is not None else 2.0, -(abs(st) if st is not None else 0.0))

            chosen = min(items, key=de_key)[1]
        elif agg == "max_mean":
            def row_mean(item: tuple[str, list[str]]) -> float:
                nums = [v for c in item[1] if (v := to_float(c)) is not None]
                return sum(nums) / len(nums) if nums else -1e30

            chosen = max(items, key=row_mean)[1]
        elif agg in ("sum", "mean"):
            n_cols = len(other_indices)
            agg_vals: list[str] = []
            for c_i in range(n_cols):
                nums = [v for _, row_v in items if c_i < len(row_v) and (v := to_float(row_v[c_i])) is not None]
                if not nums:
                    agg_vals.append(items[0][1][c_i] if c_i < len(items[0][1]) else "")
                elif agg == "sum":
                    s = sum(nums)
                    agg_vals.append(str(int(round(s))) if is_integer_counts else f"{s:.6g}")
                else:
                    m = sum(nums) / len(nums)
                    agg_vals.append(f"{m:.6g}")
            chosen = agg_vals
        else:
            chosen = items[0][1]

        out_rows.append([target_id] + chosen)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_header = [header[id_idx].strip() or "gene_id"] + other_names
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter=out_sep)
        writer.writerow(out_header)
        writer.writerows(out_rows)

    params_meta = {
        "strip_version": args.strip_version,
        "table_kind": detected_kind,
        "agg_method": agg,
        "unmapped": args.unmapped,
        "biotype": args.biotype,
    }
    in_files = (args.input, args.gtf, args.mapping_table)
    artifacts = [
        write_table_artifact(
            str(out_path),
            f"Mapped {detected_kind} table ({len(out_rows)} unique genes, agg={agg})",
            inputs=in_files,
            parameters=params_meta,
            rows=len(out_rows),
            columns=len(out_header),
        )
    ]

    collision_path_str: str | None = None
    if collision_rows:
        col_path = out_path.with_name(out_path.stem + "-id-collisions.tsv")
        with open(col_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["target_gene", "n_source_ids", "source_ids"])
            writer.writerows(collision_rows)
        collision_path_str = str(col_path)
        artifacts.append(
            write_table_artifact(
                collision_path_str,
                f"Many-to-one gene ID collision report ({n_collisions} merged symbols)",
                inputs=in_files,
                parameters=params_meta,
                rows=len(collision_rows),
                columns=3,
            )
        )

    if detected_kind == "de":
        sem = "de_stats"
        sci_ctx = {str(out_path): {"valueSemantics": sem}}
    else:
        if is_integer_counts:
            sem = "raw_counts"
        elif has_negative or max_val < 30.0:
            sem = "log_transformed"
        else:
            sem = "normalized"
        sci_ctx = {str(out_path): {"valueSemantics": sem, "sampleIds": other_names}}

    warnings: list[str] = []
    if (args.gtf or args.mapping_table) and n_mapped == 0:
        warnings.append(
            "0 gene IDs matched the supplied mapping source; check identifier namespace (e.g. ENSG vs ENSMUSG vs Entrez)"
        )
    elif (args.gtf or args.mapping_table) and n_unmapped / max(1, n_mapped + n_unmapped) > 0.3:
        warnings.append(
            f"{n_unmapped}/{n_mapped + n_unmapped} rows did not match the mapping source (unmapped={args.unmapped})"
        )

    emit(
        {
            "tool": "rna_id_map",
            "input": args.input,
            "output": str(out_path),
            "collision_report": collision_path_str,
            "table_kind": detected_kind,
            "agg_method": agg,
            "strip_version": args.strip_version,
            "biotype_filter": args.biotype,
            "n_input_rows": n_input_rows,
            "n_versions_stripped": n_versions_stripped,
            "n_biotype_filtered": n_biotype_filtered,
            "n_mapped": n_mapped,
            "n_unmapped": n_unmapped,
            "n_collisions_merged": n_collisions,
            "n_output_genes": len(out_rows),
            "scientificContext": sci_ctx,
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
