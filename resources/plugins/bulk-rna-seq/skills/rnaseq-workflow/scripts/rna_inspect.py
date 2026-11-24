#!/usr/bin/env python3
"""rna_inspect: inspect bulk RNA-seq expression matrices, DE tables, sample
metadata, and GMT gene-set files. Detects value kind (raw_counts vs
normalized_tpm_or_fpkm vs log_transformed), duplicate gene IDs, sample
alignment, low-expression fraction, and batch-condition collinearity.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import emit, fail, read_de_table, read_expression_matrix, read_gmt, read_metadata


def infer_matrix_kind(df):
    x = df.to_numpy(dtype=float)
    min_v = float(np.nanmin(x))
    max_v = float(np.nanmax(x))
    if min_v < 0:
        return "log_or_vst_centered", "Contains negative values (likely VST/rlog/z-score); do NOT pass to DESeq2/edgeR."
    is_int = bool(np.all(np.abs(x - np.round(x)) < 1e-6))
    col_sums = x.sum(axis=0)
    median_sum = float(np.median(col_sums))
    if is_int and max_v >= 20:
        return "raw_counts", "Non-negative integers compatible with DESeq2 / edgeR raw count wrappers."
    if 9e5 <= median_sum <= 1.1e6 and not is_int:
        return "tpm_or_cpm", "Column sums ~1e6 with fractional values (TPM/CPM); do NOT pass to DESeq2/edgeR raw count models."
    if max_v <= 30.0 and not is_int:
        return "log_normalized", "Low-range fractional values (likely log2-CPM/log2-TPM/VST); use limma or rna_pathway_score."
    return "fractional_normalized", "Non-integer non-negative values (FPKM/TPM/estimated counts); round only if RSEM/Salmon expected counts, never if TPM/FPKM."


def check_collinearity(meta, group_col, batch_col):
    if group_col not in meta.columns or batch_col not in meta.columns:
        return None
    ct = pd.crosstab(meta[group_col].astype(str), meta[batch_col].astype(str))
    # Each batch only appears in one group OR each group only appears in one batch
    batch_confounded = bool((ct > 0).sum(axis=0).max() == 1)
    group_confounded = bool((ct > 0).sum(axis=1).max() == 1)
    return {
        "group_col": group_col,
        "batch_col": batch_col,
        "confounded": bool(batch_confounded or group_confounded),
        "crosstab": {str(g): {str(b): int(ct.loc[g, b]) for b in ct.columns} for g in ct.index},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", help="Gene x sample expression/count matrix TSV/CSV.")
    ap.add_argument("--metadata", help="Sample metadata TSV/CSV.")
    ap.add_argument("--sample_col", default=None)
    ap.add_argument("--group_col", default="condition")
    ap.add_argument("--batch_col", default=None)
    ap.add_argument("--de_table", help="Optional existing DESeq2/edgeR/limma DE result table.")
    ap.add_argument("--gmt", help="Optional GMT gene-set file.")
    args = ap.parse_args()

    if not (args.matrix or args.de_table or args.metadata or args.gmt):
        fail("provide at least one of --matrix, --de_table, --metadata, or --gmt")

    report = {"tool": "rna_inspect", "warnings": [], "ready_steps": []}

    matrix_samples = None
    if args.matrix:
        p = Path(args.matrix)
        if not p.is_file():
            fail(f"matrix file not found: {p}")
        df = read_expression_matrix(p)
        matrix_samples = df.columns.tolist()
        dup_genes = int(df.index.duplicated().sum())
        kind, kind_note = infer_matrix_kind(df)
        zero_all = int((df.sum(axis=1) == 0).sum())
        col_sums = df.sum(axis=0)
        report["matrix"] = {
            "path": str(p),
            "n_genes": int(df.shape[0]),
            "n_samples": int(df.shape[1]),
            "samples": matrix_samples,
            "duplicate_genes": dup_genes,
            "all_zero_genes": zero_all,
            "inferred_kind": kind,
            "kind_note": kind_note,
            "min_col_sum": float(col_sums.min()),
            "median_col_sum": float(col_sums.median()),
            "max_col_sum": float(col_sums.max()),
        }
        if dup_genes > 0:
            report["warnings"].append(f"Matrix has {dup_genes} duplicate gene IDs; deduplicate before DE/pathway analysis.")
        if kind == "raw_counts":
            report["ready_steps"].extend(["rna_qc", "local/modules/differential-expression-deseq2", "local/modules/differential-expression-edger"])
        else:
            report["warnings"].append(kind_note)
            report["ready_steps"].extend(["rna_qc", "rna_pathway_score", "local/modules/differential-expression-limma"])

    if args.metadata:
        mp = Path(args.metadata)
        if not mp.is_file():
            fail(f"metadata file not found: {mp}")
        meta = read_metadata(mp, sample_col=args.sample_col)
        meta_info = {
            "path": str(mp),
            "n_rows": int(len(meta)),
            "columns": meta.columns.tolist(),
        }
        if matrix_samples is not None:
            missing_in_meta = [s for s in matrix_samples if s not in meta.index]
            extra_in_meta = [s for s in meta.index if s not in matrix_samples]
            meta_info["missing_matrix_samples"] = missing_in_meta
            meta_info["extra_metadata_samples"] = extra_in_meta
            if missing_in_meta:
                report["warnings"].append(f"Metadata missing matrix samples: {missing_in_meta}")
        if args.group_col in meta.columns:
            vc = meta[args.group_col].astype(str).value_counts().to_dict()
            meta_info["group_counts"] = {str(k): int(v) for k, v in vc.items()}
            min_rep = min(vc.values()) if vc else 0
            if min_rep < 3:
                report["warnings"].append(
                    f"Group column {args.group_col!r} has a group with only {min_rep} replicate(s); "
                    "DESeq2/edgeR statistical power and dispersion estimates are weak below n=3."
                )
        if args.batch_col:
            collin = check_collinearity(meta, args.group_col, args.batch_col)
            if collin is not None:
                meta_info["batch_collinearity"] = collin
                if collin["confounded"]:
                    report["warnings"].append(
                        f"CONFOUNDED DESIGN: {args.batch_col!r} is completely collinear with {args.group_col!r}; "
                        "batch effect cannot be separated from biological condition in ~ batch + condition."
                    )
        report["metadata"] = meta_info

    if args.de_table:
        dp = Path(args.de_table)
        if not dp.is_file():
            fail(f"DE table not found: {dp}")
        de = read_de_table(dp)
        sig_05 = int(((de["padj"] < 0.05) & (de["log2fc"].abs() >= 1.0)).sum()) if de["padj"].notna().any() else 0
        report["de_table"] = {
            "path": str(dp),
            "n_genes": int(len(de)),
            "n_sig_padj005_lfc1": sig_05,
            "n_up_padj005_lfc1": int(((de["padj"] < 0.05) & (de["log2fc"] >= 1.0)).sum()) if de["padj"].notna().any() else 0,
            "n_down_padj005_lfc1": int(((de["padj"] < 0.05) & (de["log2fc"] <= -1.0)).sum()) if de["padj"].notna().any() else 0,
        }
        report["ready_steps"].extend(["rna_compare", "rna_ora", "rna_gsea"])

    if args.gmt:
        gs = read_gmt(args.gmt)
        sizes = [len(v) for v in gs.values()]
        report["gmt"] = {
            "path": str(args.gmt),
            "n_gene_sets": int(len(gs)),
            "min_set_size": int(min(sizes)),
            "median_set_size": float(np.median(sizes)),
            "max_set_size": int(max(sizes)),
        }

    report["ready_steps"] = sorted(set(report["ready_steps"]))
    emit(report)


if __name__ == "__main__":
    main()
