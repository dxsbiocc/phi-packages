#!/usr/bin/env python3
"""meth_annotate: annotate DMRs or DMCs against a reference GTF/GFF3
(Promoter / GeneBody / Downstream / Distal_Intergenic + nearest TSS) and
optional CpG Island BED (Island / Shore / Shelf / OpenSea), plus optional GMT
hypergeometric ORA against the assayed CpG gene universe.
"""

from __future__ import annotations

import argparse
import bisect
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import bh_fdr, emit, fail, open_text, parse_cpg_coord, read_table, write_artifact

GTF_ATTR_RE = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]+)"')


def parse_gtf_genes(path: str) -> dict[str, list[tuple[int, int, int, str, str]]]:
    """Return chrom -> sorted list of (tss, start, end, strand, gene_symbol)."""
    by_chrom: dict[str, list[tuple[int, int, int, str, str]]] = {}
    seen: set[tuple[str, str]] = set()
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 9:
                continue
            feat = parts[2]
            if feat not in ("gene", "transcript", "mRNA"):
                continue
            chrom, start_s, end_s, strand, attr = parts[0], parts[3], parts[4], parts[6], parts[8]
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                continue
            kv = dict(GTF_ATTR_RE.findall(attr))
            if not kv and "=" in attr:
                for item in attr.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        kv[k.strip()] = v.strip().strip('"')
            sym = kv.get("gene_name") or kv.get("Name") or kv.get("gene_id") or kv.get("ID")
            if not sym:
                continue
            if feat != "gene" and (chrom, sym) in seen:
                continue
            seen.add((chrom, sym))
            tss = start if strand != "-" else end
            by_chrom.setdefault(chrom, []).append((tss, start, end, strand, sym))

    if not by_chrom:
        fail(f"0 gene/transcript records parsed from GTF/GFF3: {path}")
    for ch in by_chrom:
        by_chrom[ch].sort(key=lambda x: x[0])
    return by_chrom


def parse_bed_intervals(path: str) -> dict[str, tuple[list[int], list[tuple[int, int]]]]:
    raw: dict[str, list[tuple[int, int]]] = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.lower().startswith("track"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                s, e = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            raw.setdefault(parts[0].strip(), []).append((s, e))
    out = {}
    for ch, ivs in raw.items():
        ivs.sort()
        out[ch] = ([x[0] for x in ivs], ivs)
    return out


def classify_cgi(chrom: str, mid: int, cgi_index: dict[str, tuple[list[int], list[tuple[int, int]]]]) -> str:
    if chrom not in cgi_index:
        return "OpenSea"
    starts, ivs = cgi_index[chrom]
    idx = bisect.bisect_right(starts, mid)
    min_dist = 10**9
    for k in range(max(0, idx - 3), min(len(ivs), idx + 3)):
        s, e = ivs[k]
        if s <= mid <= e:
            return "Island"
        dist = (s - mid) if mid < s else (mid - e)
        if dist < min_dist:
            min_dist = dist
    if min_dist <= 2000:
        return "Shore"
    if min_dist <= 4000:
        return "Shelf"
    return "OpenSea"


def annotate_locus(
    chrom: str,
    start: int,
    end: int,
    genes_by_chrom: dict[str, list[tuple[int, int, int, str, str]]],
    prom_up: int,
    prom_down: int,
) -> tuple[str, str, int | None]:
    genes = genes_by_chrom.get(chrom)
    if not genes and not chrom.startswith("chr"):
        genes = genes_by_chrom.get("chr" + chrom)
    elif not genes and chrom.startswith("chr"):
        genes = genes_by_chrom.get(chrom[3:])
    if not genes:
        return "Distal_Intergenic", "Unmapped", None

    mid = (start + end) // 2
    tss_list = [g[0] for g in genes]
    pos_i = bisect.bisect_left(tss_list, mid)

    # Search nearby genes for Promoter -> GeneBody -> Downstream -> Nearest TSS
    cand_slice = genes[max(0, pos_i - 25) : min(len(genes), pos_i + 25)]
    best_gene = None
    best_abs_dist = 10**12
    best_signed_dist = 0
    in_body_gene = None
    in_down_gene = None

    for tss, gs, ge, strand, sym in cand_slice:
        signed_d = (mid - tss) if strand != "-" else (tss - mid)
        p_lo = (tss - prom_up) if strand != "-" else (tss - prom_down)
        p_hi = (tss + prom_down) if strand != "-" else (tss + prom_up)
        if end >= p_lo and start <= p_hi:
            return "Promoter", sym, int(signed_d)
        if end >= gs and start <= ge and in_body_gene is None:
            in_body_gene = (sym, int(signed_d))
        if 0 < signed_d <= (ge - tss + 3000 if strand != "-" else tss - gs + 3000) and in_down_gene is None:
            if (strand != "-" and ge < mid <= ge + 3000) or (strand == "-" and gs - 3000 <= mid < gs):
                in_down_gene = (sym, int(signed_d))
        if abs(mid - tss) < best_abs_dist:
            best_abs_dist = abs(mid - tss)
            best_gene = sym
            best_signed_dist = int(signed_d)

    if in_body_gene is not None:
        return "GeneBody", in_body_gene[0], in_body_gene[1]
    if in_down_gene is not None:
        return "Downstream", in_down_gene[0], in_down_gene[1]
    return "Distal_Intergenic", best_gene or "Unmapped", best_signed_dist


def read_gmt(path: str) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with open_text(path) as fh:
        for line in fh:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) >= 3:
                genes = {g.strip() for g in parts[2:] if g.strip()}
                if genes:
                    sets[parts[0].strip()] = genes
    return sets


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--cpg_islands", default=None)
    parser.add_argument("--gmt", default=None)
    parser.add_argument("--all_tested_table", default=None)
    parser.add_argument("--promoter_upstream", type=int, default=2000)
    parser.add_argument("--promoter_downstream", type=int, default=500)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    df = read_table(args.input)
    genes_by_chrom = parse_gtf_genes(args.gtf)
    cgi_idx = parse_bed_intervals(args.cpg_islands) if args.cpg_islands else None

    # Resolve coordinates per row
    has_dmr_cols = {"chrom", "start", "end"}.issubset(set(df.columns))
    has_dmc_cols = {"chrom", "pos"}.issubset(set(df.columns))

    feat_classes = []
    nearest_genes = []
    tss_dists = []
    cgi_classes = []

    for _, row in df.iterrows():
        ch, st, en = None, None, None
        if has_dmr_cols and pd.notna(row["chrom"]) and pd.notna(row["start"]):
            ch, st, en = str(row["chrom"]), int(row["start"]), int(row["end"])
        elif has_dmc_cols and pd.notna(row["chrom"]) and pd.notna(row["pos"]):
            ch, st = str(row["chrom"]), int(row["pos"])
            en = st
        elif "cpg_id" in df.columns:
            parsed = parse_cpg_coord(str(row["cpg_id"]))
            if parsed:
                ch, st, en = parsed
        if ch is None or st is None or en is None:
            feat_classes.append("Unmapped")
            nearest_genes.append("Unmapped")
            tss_dists.append(None)
            cgi_classes.append("Unknown")
            continue

        fc, gsym, d_tss = annotate_locus(
            ch, st, en, genes_by_chrom, args.promoter_upstream, args.promoter_downstream
        )
        feat_classes.append(fc)
        nearest_genes.append(gsym)
        tss_dists.append(d_tss)
        if cgi_idx is not None:
            cgi_classes.append(classify_cgi(ch, (st + en) // 2, cgi_idx))

    out_df = df.copy()
    out_df["genomic_feature"] = feat_classes
    out_df["nearest_gene"] = nearest_genes
    out_df["distance_to_tss"] = tss_dists
    if cgi_idx is not None:
        out_df["cpg_island_context"] = cgi_classes

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    ann_path = outdir / "annotated_methylation.tsv"
    out_df.to_csv(ann_path, sep="\t", index=False)

    feat_counts = {str(k): int(v) for k, v in out_df["genomic_feature"].value_counts().items()}
    cgi_counts = (
        {str(k): int(v) for k, v in out_df["cpg_island_context"].value_counts().items()}
        if cgi_idx is not None
        else None
    )

    # Bar plot of genomic context (and CpG island context if available)
    n_panels = 2 if cgi_counts else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(5.8 * n_panels, 4.2))
    ax0 = axes[0] if n_panels == 2 else axes
    f_keys = list(feat_counts.keys())
    f_vals = [feat_counts[k] for k in f_keys]
    ax0.barh(f_keys[::-1], f_vals[::-1], color="#2563EB", alpha=0.85)
    ax0.set_xlabel("Count")
    ax0.set_title("Genomic Feature Distribution")
    ax0.spines["top"].set_visible(False)
    ax0.spines["right"].set_visible(False)

    if cgi_counts:
        ax1 = axes[1]
        c_keys = list(cgi_counts.keys())
        c_vals = [cgi_counts[k] for k in c_keys]
        ax1.barh(c_keys[::-1], c_vals[::-1], color="#10B981", alpha=0.85)
        ax1.set_xlabel("Count")
        ax1.set_title("CpG Island Context Distribution")
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)

    fig.tight_layout()
    fig_path = outdir / "meth-genomic-context.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    in_files = [args.input, args.gtf, args.cpg_islands, args.gmt, args.all_tested_table]
    outputs: dict[str, str | None] = {"annotated_table": str(ann_path), "ora_table": None}
    artifacts = [
        write_artifact(fig_path, "meth_annotate", "Methylation genomic feature distribution", inputs=in_files),
        write_artifact(ann_path, "meth_annotate", "Genomically annotated DMR/DMC table", inputs=in_files, shape=out_df.shape),
    ]

    ora_summary = None
    warnings: list[str] = []
    if args.gmt:
        gmt_sets = read_gmt(args.gmt)
        # Foreground genes: significant rows if 'significant' or 'status' column exists, else all rows
        if "significant" in out_df.columns:
            fg_sub = out_df[out_df["significant"].astype(bool)]
        elif "status" in out_df.columns:
            fg_sub = out_df[out_df["status"].isin(["hyper", "hypo"])]
        else:
            fg_sub = out_df
        fg_genes = {g for g in fg_sub["nearest_gene"].astype(str) if g and g != "Unmapped"}

        # Universe: all genes in GTF or all genes mapped from --all_tested_table
        if args.all_tested_table:
            bg_df = read_table(args.all_tested_table)
            bg_genes = set()
            for _, r in bg_df.iterrows():
                if "chrom" in bg_df.columns and "pos" in bg_df.columns and pd.notna(r["chrom"]) and pd.notna(r["pos"]):
                    _, gs, _ = annotate_locus(
                        str(r["chrom"]), int(r["pos"]), int(r["pos"]), genes_by_chrom, args.promoter_upstream, args.promoter_downstream
                    )
                    if gs != "Unmapped":
                        bg_genes.add(gs)
        else:
            bg_genes = {g[4] for ch_list in genes_by_chrom.values() for g in ch_list}
            warnings.append(
                "ORA universe defaulted to all GTF genes; pass --all_tested_table (dmc_results.tsv) to restrict universe to assayed CpG genes"
            )

        fg_in_bg = fg_genes & bg_genes
        M = len(bg_genes)
        N = len(fg_in_bg)
        ora_rows = []
        if N >= 2 and M >= 10:
            for pw, pgenes in gmt_sets.items():
                pset = pgenes & bg_genes
                K = len(pset)
                if K < 3 or K > 500:
                    continue
                overlap = sorted(fg_in_bg & pset)
                k = len(overlap)
                if k == 0:
                    continue
                pval = float(stats.hypergeom.sf(k - 1, M, K, N))
                ora_rows.append(
                    {
                        "pathway": pw,
                        "overlap": k,
                        "set_size": K,
                        "fold_enrichment": round((k / N) / (K / M), 4),
                        "pvalue": pval,
                        "genes": ";".join(overlap[:25]),
                    }
                )
        if ora_rows:
            ora_df = pd.DataFrame(ora_rows)
            ora_df["fdr"] = bh_fdr(ora_df["pvalue"].to_numpy(dtype=float))
            ora_df = ora_df.sort_values("pvalue", ascending=True)
            ora_path = outdir / "meth_ora.tsv"
            ora_df.to_csv(ora_path, sep="\t", index=False)
            outputs["ora_table"] = str(ora_path)
            artifacts.append(write_artifact(ora_path, "meth_annotate", "Methylation target-gene ORA table", inputs=in_files, shape=ora_df.shape))
            ora_summary = {
                "n_foreground_genes": N,
                "n_universe_genes": M,
                "n_significant_pathways": int((ora_df["fdr"] < 0.05).sum()),
            }

    emit(
        {
            "tool": "meth_annotate",
            "n_regions_annotated": int(len(out_df)),
            "genomic_feature_counts": feat_counts,
            "cpg_island_context_counts": cgi_counts,
            "ora_summary": ora_summary,
            "outputs": outputs,
            "figures": [str(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
