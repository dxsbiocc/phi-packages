#!/usr/bin/env python3
"""rna_geo_import: parse an NCBI GEO Series Matrix file (GSE*_series_matrix.txt/.gz)
into a clean expression matrix TSV and structured samplesheet TSV, with optional
GPL platform probe-to-gene mapping and value-scale inference.

Pure Python 3 stdlib implementation.
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


GENE_SYMBOL_COL_CANDS = (
    "gene symbol",
    "gene_symbol",
    "symbol",
    "ilmn_gene",
    "orf",
    "gene_assignment",
    "gene",
)


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


def slug_col(name: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9]+", "_", name.strip().lower()).strip("_")
    return s or "characteristic"


def unquote_tsv_line(line: str) -> list[str]:
    reader = csv.reader([line.rstrip("\r\n")], delimiter="\t", quotechar='"')
    return [c.strip() for c in next(reader, [])]


def to_float(val: str) -> float | None:
    s = val.strip()
    if not s or s.lower() in ("na", "nan", "null", "none"):
        return None
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except ValueError:
        return None


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
    desc: dict[str, Any] = {
        "contractVersion": "1.0.0",
        "kind": "table",
        "file": p.name,
        "mediaType": "text/tab-separated-values",
        "title": title[:200],
        "table": {"rows": int(rows), "columns": int(columns)},
        "provenance": {
            "createdAt": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "tool": "rna_geo_import",
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


def parse_platform_table(
    path: str, probe_col: str | None, gene_col: str | None, multi_mode: str
) -> dict[str, str]:
    sep = "," if path.lower().endswith(".csv") else "\t"
    probe_map: dict[str, str] = {}
    with open_text(path) as fh:
        header: list[str] | None = None
        for line in fh:
            if not line.strip() or line.startswith(("#", "!", "^")):
                continue
            header = next(csv.reader([line.rstrip("\r\n")], delimiter=sep))
            break
        if not header or len(header) < 2:
            fail(f"platform_table has no valid header: {path}")

        lower_map = {h.strip().lower(): i for i, h in enumerate(header)}
        p_idx = 0
        if probe_col:
            if probe_col.lower() not in lower_map:
                fail(f"probe_col '{probe_col}' not found in platform_table: {header}")
            p_idx = lower_map[probe_col.lower()]
        elif "id" in lower_map:
            p_idx = lower_map["id"]

        g_idx = None
        if gene_col:
            if gene_col.lower() not in lower_map:
                fail(f"gene_col '{gene_col}' not found in platform_table: {header}")
            g_idx = lower_map[gene_col.lower()]
        else:
            for cand in GENE_SYMBOL_COL_CANDS:
                if cand in lower_map and lower_map[cand] != p_idx:
                    g_idx = lower_map[cand]
                    break
            if g_idx is None:
                g_idx = 1 if p_idx == 0 else 0

        reader = csv.reader(fh, delimiter=sep)
        for row in reader:
            if not row or row[0].startswith(("!", "#")):
                continue
            if len(row) <= max(p_idx, g_idx):
                continue
            pid = row[p_idx].strip()
            raw_g = row[g_idx].strip()
            if not pid or not raw_g or raw_g in ("---", "NA", "nan", ""):
                continue
            # Handle Affymetrix / GEO 'GENE1 /// GENE2' or GenBank 'NM_001 // TP53 // ...'
            if "///" in raw_g:
                syms = [s.strip() for s in raw_g.split("///") if s.strip() and s.strip() != "---"]
                if not syms:
                    continue
                uniq_syms = list(dict.fromkeys(syms))
                if len(uniq_syms) > 1 and multi_mode == "drop":
                    continue
                gene = uniq_syms[0]
            elif "//" in raw_g and "gene_assignment" == header[g_idx].strip().lower():
                parts = [p.strip() for p in raw_g.split("//")]
                gene = parts[1] if len(parts) > 1 else parts[0]
            else:
                gene = raw_g
            if gene and gene != "---":
                probe_map[pid] = gene

    if not probe_map:
        fail(f"platform_table yielded 0 valid probe -> gene mappings: {path}")
    return probe_map


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--series_matrix", required=True, help="GEO GSE*_series_matrix.txt or .txt.gz.")
    parser.add_argument("--output_dir", required=True, help="Output directory for samplesheet & matrix.")
    parser.add_argument("--platform_table", default=None, help="Optional GPL probe annotation TSV/CSV.")
    parser.add_argument("--probe_col", default=None, help="Probe ID column in --platform_table.")
    parser.add_argument("--gene_col", default=None, help="Gene symbol column in --platform_table.")
    parser.add_argument(
        "--multi_gene_probes",
        choices=["primary", "drop"],
        default="primary",
        help="How to treat multi-gene 'A /// B' probes (default primary).",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.series_matrix):
        fail(f"series_matrix file not found: {args.series_matrix}")
    if args.platform_table and not os.path.isfile(args.platform_table):
        fail(f"platform_table not found: {args.platform_table}")

    series_meta: dict[str, str] = {}
    sample_headers: dict[str, list[str]] = {}
    characteristics_rows: list[list[str]] = []
    matrix_header: list[str] = []
    matrix_rows: list[tuple[str, list[str]]] = []
    in_table = False

    with open_text(args.series_matrix) as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                break
            if in_table:
                if not line.strip():
                    continue
                cols = unquote_tsv_line(line)
                if not cols:
                    continue
                if not matrix_header:
                    matrix_header = cols
                else:
                    matrix_rows.append((cols[0], cols[1:]))
                continue

            if line.startswith("!Series_"):
                cols = unquote_tsv_line(line)
                if len(cols) >= 2:
                    key = cols[0][len("!Series_") :]
                    series_meta.setdefault(key, cols[1])
            elif line.startswith("!Sample_characteristics_ch1"):
                cols = unquote_tsv_line(line)
                if len(cols) >= 2:
                    characteristics_rows.append(cols[1:])
            elif line.startswith("!Sample_"):
                cols = unquote_tsv_line(line)
                if len(cols) >= 2:
                    key = cols[0][len("!Sample_") :]
                    if key not in sample_headers:
                        sample_headers[key] = cols[1:]

    gsms = sample_headers.get("geo_accession") or (matrix_header[1:] if len(matrix_header) > 1 else [])
    if not gsms:
        fail("no !Sample_geo_accession or matrix sample columns found in series_matrix")
    n_samples = len(gsms)

    # Parse !Sample_characteristics_ch1 into structured columns (key: value)
    char_cols: dict[str, list[str]] = {}
    unnamed_idx = 0
    for ch_idx, vals in enumerate(characteristics_rows):
        padded = [vals[i] if i < len(vals) else "" for i in range(n_samples)]
        # Infer dominant 'key: value' prefix across samples in this row
        keys_found = []
        for v in padded:
            if ":" in v:
                keys_found.append(slug_col(v.split(":", 1)[0]))
        if keys_found:
            # Handle mixed characteristics rows where different samples put keys on different rows
            distinct_keys = list(dict.fromkeys(keys_found))
            if len(distinct_keys) == 1:
                col_name = distinct_keys[0]
                if col_name in char_cols:
                    col_name = f"{col_name}_{ch_idx + 1}"
                char_cols[col_name] = [
                    v.split(":", 1)[1].strip() if ":" in v else v.strip() for v in padded
                ]
            else:
                for i, v in enumerate(padded):
                    if ":" in v:
                        k_raw, v_raw = v.split(":", 1)
                        cname = slug_col(k_raw)
                        char_cols.setdefault(cname, [""] * n_samples)
                        if not char_cols[cname][i]:
                            char_cols[cname][i] = v_raw.strip()
        else:
            unnamed_idx += 1
            char_cols[f"characteristic_{unnamed_idx}"] = padded

    # Build samplesheet
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    meta_cols = ["sample", "title", "source_name", "organism"] + list(char_cols.keys())
    samplesheet_path = outdir / "samplesheet.tsv"
    titles = sample_headers.get("title", [""] * n_samples)
    sources = sample_headers.get("source_name_ch1", [""] * n_samples)
    organisms = sample_headers.get("organism_ch1", [""] * n_samples)

    with open(samplesheet_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(meta_cols)
        for i in range(n_samples):
            row = [
                gsms[i],
                titles[i] if i < len(titles) else "",
                sources[i] if i < len(sources) else "",
                organisms[i] if i < len(organisms) else "",
            ] + [char_cols[c][i] for c in char_cols]
            writer.writerow(row)

    warnings: list[str] = []
    artifacts: list[str] = []
    in_files = (args.series_matrix, args.platform_table)
    params_meta = {"multi_gene_probes": args.multi_gene_probes, "mapped_platform": bool(args.platform_table)}

    artifacts.append(
        write_table_artifact(
            str(samplesheet_path),
            f"GEO samplesheet ({series_meta.get('geo_accession', 'GSE')}, {n_samples} samples)",
            inputs=in_files,
            parameters=params_meta,
            rows=n_samples,
            columns=len(meta_cols),
        )
    )

    # Process expression matrix if non-empty (note: RNA-seq GEO Series Matrix files often omit the matrix and attach raw counts in supplementary files!)
    matrix_path_str: str | None = None
    inferred_scale = "absent_supplementary_only"
    n_features_raw = len(matrix_rows)
    n_features_out = 0
    sci_ctx: dict[str, Any] = {}

    if matrix_rows and len(matrix_header) > 1:
        mat_samples = matrix_header[1:]
        probe_map = (
            parse_platform_table(args.platform_table, args.probe_col, args.gene_col, args.multi_gene_probes)
            if args.platform_table
            else None
        )

        grouped: dict[str, list[list[str]]] = {}
        is_int = True
        has_neg = False
        max_val = 0.0

        for feat_id, vals in matrix_rows:
            fid = feat_id.strip()
            if not fid:
                continue
            if probe_map is not None:
                gene = probe_map.get(fid)
                if not gene:
                    continue
            else:
                gene = fid

            padded_vals = [vals[j] if j < len(vals) else "" for j in range(len(mat_samples))]
            for cell in padded_vals[:20]:
                v = to_float(cell)
                if v is None:
                    continue
                if v < 0:
                    has_neg = True
                    is_int = False
                elif abs(v - round(v)) > 1e-6:
                    is_int = False
                if v > max_val:
                    max_val = v
            grouped.setdefault(gene, []).append(padded_vals)

        out_mat_rows: list[list[str]] = []
        for gene, item_list in grouped.items():
            if len(item_list) == 1:
                chosen = item_list[0]
            else:
                # Pick highest-mean probe (standard microarray / probe-level collapse)
                def probe_mean(vals: list[str]) -> float:
                    nums = [v for c in vals if (v := to_float(c)) is not None]
                    return sum(nums) / len(nums) if nums else -1e30

                chosen = max(item_list, key=probe_mean)
            out_mat_rows.append([gene] + chosen)

        n_features_out = len(out_mat_rows)
        mat_path = outdir / "expression_matrix.tsv"
        with open(mat_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["gene_id"] + mat_samples)
            writer.writerows(out_mat_rows)
        matrix_path_str = str(mat_path)

        if is_int and max_val >= 20:
            inferred_scale = "raw_counts"
        elif has_neg or max_val <= 30.0:
            inferred_scale = "log_transformed"
            warnings.append(
                "matrix values appear log-transformed (typical microarray / normalized RNA-seq); "
                "use limma or rna_pathway_score rather than DESeq2/edgeR"
            )
        else:
            inferred_scale = "normalized"
            warnings.append(
                "matrix values are non-integer linear intensities/TPM/FPKM; do NOT feed into DESeq2/edgeR"
            )

        sci_ctx[matrix_path_str] = {
            "valueSemantics": inferred_scale,
            "sampleIds": mat_samples,
        }
        artifacts.insert(
            0,
            write_table_artifact(
                matrix_path_str,
                f"GEO expression matrix ({n_features_out} features x {len(mat_samples)} samples, {inferred_scale})",
                inputs=in_files,
                parameters=params_meta,
                rows=n_features_out,
                columns=len(mat_samples) + 1,
            ),
        )
    else:
        warnings.append(
            "Series Matrix contains 0 inline expression rows (common for high-throughput RNA-seq GSE entries); "
            "samplesheet.tsv was extracted — pair it with the GSE supplementary raw count matrix"
        )

    emit(
        {
            "tool": "rna_geo_import",
            "series_accession": series_meta.get("geo_accession"),
            "series_title": series_meta.get("title"),
            "platform_id": series_meta.get("platform_id"),
            "n_samples": n_samples,
            "metadata_columns": meta_cols,
            "n_features_raw": n_features_raw,
            "n_features_output": n_features_out,
            "inferred_scale": inferred_scale,
            "outputs": {
                "samplesheet": str(samplesheet_path),
                "expression_matrix": matrix_path_str,
            },
            **({"scientificContext": sci_ctx} if sci_ctx else {}),
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
