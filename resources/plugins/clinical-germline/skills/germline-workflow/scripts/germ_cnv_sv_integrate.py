#!/usr/bin/env python3
"""germ_cnv_sv_integrate: Score germline CNV/SV calls (DEL, DUP, INV, BND)
under ClinGen/ACMG 2020 (Riggs et al.) dosage-sensitivity criteria and detect
trans-acting SNV + exon-DEL compound heterozygotes in recessive disease genes.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    is_vcf_file,
    norm_chrom,
    norm_gene,
    read_table,
    read_text_lines,
    resolve_column,
    safe_float,
    safe_int,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)

# Built-in ClinGen dosage-sensitive & rare-disease genes:
# (chrom, start_hg38, end_hg38, HI_score_0_3, TS_score_0_3, inheritance_mode)
CLINGEN_GENES: dict[str, tuple[str, int, int, int, int, str]] = {
    "SCN1A": ("chr2", 165984641, 166149161, 3, 0, "AD"),
    "NRXN1": ("chr2", 49918503, 51225575, 3, 0, "AD"),
    "MSH2": ("chr2", 47403067, 47663146, 3, 0, "AD"),
    "MLH1": ("chr3", 36993330, 37050918, 3, 0, "AD"),
    "VHL": ("chr3", 10141635, 10153670, 3, 0, "AD"),
    "USH2A": ("chr1", 215622894, 216423448, 0, 0, "AR"),
    "ABCA4": ("chr1", 93992834, 94121149, 0, 0, "AR"),
    "GBA1": ("chr1", 155234452, 155244699, 0, 0, "AR"),
    "RPE65": ("chr1", 68428822, 68450013, 0, 0, "AR"),
    "SMN1": ("chr5", 70924941, 70953015, 0, 0, "AR"),
    "EYS": ("chr6", 63719980, 65707225, 0, 0, "AR"),
    "PMS2": ("chr7", 5970925, 6009106, 3, 0, "AD"),
    "CFTR": ("chr7", 117465784, 117715971, 0, 0, "AR"),
    "SLC26A4": ("chr7", 107660635, 107717809, 0, 0, "AR"),
    "TSC1": ("chr9", 132891348, 132946874, 3, 0, "AD"),
    "PTEN": ("chr10", 87863113, 87971930, 3, 0, "AD"),
    "KCNQ1": ("chr11", 2444684, 2849110, 3, 0, "AD"),
    "MYBPC3": ("chr11", 47331406, 47352702, 3, 0, "AD"),
    "MYO7A": ("chr11", 77128255, 77215240, 0, 0, "AR"),
    "PAH": ("chr12", 102836885, 102958441, 0, 0, "AR"),
    "ATP7B": ("chr13", 51930436, 52012130, 0, 0, "AR"),
    "RB1": ("chr13", 48303744, 48481890, 3, 0, "AD"),
    "GJB2": ("chr13", 20187463, 20192938, 0, 0, "AR"),
    "FBN1": ("chr15", 48408313, 48645709, 3, 0, "AD"),
    "UBE3A": ("chr15", 25333728, 25439056, 3, 0, "AD"),
    "PKD1": ("chr16", 2088708, 2135898, 3, 0, "AD"),
    "TSC2": ("chr16", 2047985, 2089491, 3, 0, "AD"),
    "TP53": ("chr17", 7661779, 7687550, 3, 0, "AD"),
    "NF1": ("chr17", 31094927, 31377677, 3, 0, "AD"),
    "BRCA1": ("chr17", 43044295, 43170245, 3, 0, "AD"),
    "PMP22": ("chr17", 15229785, 15265326, 3, 3, "AD"),
    "GAA": ("chr17", 80101535, 80119881, 0, 0, "AR"),
    "BRCA2": ("chr13", 32315086, 32400268, 3, 0, "AD"),
    "LDLR": ("chr19", 11089362, 11133820, 3, 0, "AD"),
    "STK11": ("chr19", 1205798, 1228434, 3, 0, "AD"),
    "SHANK3": ("chr22", 50674415, 50733298, 3, 2, "AD"),
    "DMD": ("chrX", 31097677, 33339609, 3, 0, "XL"),
    "MECP2": ("chrX", 154021573, 154097731, 3, 3, "XL"),
    "PLP1": ("chrX", 103776506, 103792619, 2, 3, "XL"),
    "COL4A5": ("chrX", 108439800, 108697545, 3, 0, "XL"),
}

SV_TIER_ORDER = ["Pathogenic", "Likely_pathogenic", "VUS", "Benign"]


def _load_gene_intervals(gene_bed: str | None) -> list[tuple[str, int, int, str]]:
    intervals: list[tuple[str, int, int, str]] = [
        (chrom, s, e, g) for g, (chrom, s, e, _, _, _) in CLINGEN_GENES.items()
    ]
    if not gene_bed:
        return intervals
    p = Path(gene_bed)
    if not p.is_file():
        fail(f"gene_bed file not found: {gene_bed}")
    lines = read_text_lines(p)
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        parts = re.split(r"[\t,]+", s)
        if len(parts) < 4:
            continue
        if parts[0].lower() in {"chrom", "chr", "#chrom"}:
            continue
        c = norm_chrom(parts[0])
        st = safe_int(parts[1])
        en = safe_int(parts[2])
        g = norm_gene(parts[3])
        if c and st is not None and en is not None and g:
            intervals.append((c, st, en, g))
    return intervals


def _parse_sv_vcf(path: Path) -> pd.DataFrame:
    lines = read_text_lines(path)
    samples: list[str] = []
    rows: list[dict[str, Any]] = []
    for line in lines:
        s = line.rstrip("\r\n")
        if not s or s.startswith("##"):
            continue
        if s.startswith("#CHROM"):
            cols = s.split("\t")
            samples = [c.strip() for c in cols[9:]] if len(cols) > 9 else ["SAMPLE_1"]
            continue
        cols = s.split("\t")
        if len(cols) < 8:
            continue
        chrom = norm_chrom(cols[0])
        start = safe_int(cols[1], 0) or 0
        info_map: dict[str, str] = {}
        for item in cols[7].split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                info_map[k.strip()] = v.strip()
        sv_type = str(info_map.get("SVTYPE", "DEL")).upper()
        end = safe_int(info_map.get("END"))
        svlen = safe_int(info_map.get("SVLEN"))
        if end is None:
            end = start + abs(svlen) if svlen is not None else start + 1
        genes = info_map.get("GENE") or info_map.get("GENES") or info_map.get("SYMBOL") or ""
        cn = safe_float(info_map.get("CN"))
        sid = samples[0] if samples else "SAMPLE_1"
        if len(cols) > 9 and "CN" in cols[8].split(":"):
            fmt = cols[8].split(":")
            vals = cols[9].split(":")
            cn_idx = fmt.index("CN")
            if cn_idx < len(vals):
                cn_v = safe_float(vals[cn_idx])
                if cn_v is not None:
                    cn = cn_v
        rows.append(
            {
                "sample": sid,
                "chrom": chrom,
                "start": start,
                "end": end,
                "sv_type": sv_type,
                "copy_number": cn if cn is not None else np.nan,
                "genes": genes,
            }
        )
    if not rows:
        fail(f"No SV/CNV records found in VCF: {path}")
    return pd.DataFrame(rows)


def _parse_sv_table(path: Path) -> pd.DataFrame:
    lines = [ln.strip() for ln in read_text_lines(path) if ln.strip() and not ln.strip().startswith("##")]
    if not lines:
        fail(f"SV/CNV file is empty: {path}")
    first = lines[0].lstrip("#")
    parts = re.split(r"[\t,]+", first)
    # Detect headerless BED (col1/col2 are ints)
    if len(parts) >= 3 and safe_int(parts[1]) is not None and safe_int(parts[2]) is not None:
        rows = []
        for ln in lines:
            if ln.startswith("#"):
                continue
            p = re.split(r"[\t,]+", ln)
            if len(p) < 3:
                continue
            rows.append(
                {
                    "sample": p[5] if len(p) > 5 else "SAMPLE_1",
                    "chrom": norm_chrom(p[0]),
                    "start": safe_int(p[1], 0) or 0,
                    "end": safe_int(p[2], 0) or 0,
                    "sv_type": p[3].upper() if len(p) > 3 else "DEL",
                    "copy_number": safe_float(p[4]) if len(p) > 4 else np.nan,
                    "genes": p[6] if len(p) > 6 else "",
                }
            )
        return pd.DataFrame(rows)

    df = read_table(path)
    chrom_col = resolve_column(df, ["chrom", "chr", "#chrom", "chromosome"], "chrom")
    start_col = resolve_column(df, ["start", "pos", "begin"], "start")
    end_col = resolve_column(df, ["end", "stop"], "end")
    type_col = resolve_column(df, ["sv_type", "svtype", "type", "cnv_type", "call"], "sv_type")
    samp_col = resolve_column(df, ["sample", "sample_id", "patient_id"], "sample", required=False)
    cn_col = resolve_column(df, ["copy_number", "cn", "copies"], "copy_number", required=False)
    gene_col = resolve_column(df, ["genes", "gene", "symbol", "gene_symbol"], "genes", required=False)

    out = pd.DataFrame(
        {
            "sample": df[samp_col].fillna("SAMPLE_1").astype(str) if samp_col else "SAMPLE_1",
            "chrom": df[chrom_col].map(norm_chrom),
            "start": df[start_col].map(lambda x: safe_int(x, 0) or 0),
            "end": df[end_col].map(lambda x: safe_int(x, 0) or 0),
            "sv_type": df[type_col].fillna("DEL").astype(str).str.upper().str.strip(),
            "copy_number": df[cn_col].map(lambda x: safe_float(x)) if cn_col else np.nan,
            "genes": df[gene_col].fillna("").astype(str) if gene_col else "",
        }
    )
    return out


def _normalize_sv_type(sv_type: str, cn: float | None) -> str:
    s = sv_type.upper()
    if s in {"DEL", "DELETION", "LOSS"}:
        return "DEL"
    if s in {"DUP", "DUPLICATION", "GAIN"}:
        return "DUP"
    if s in {"INV", "INVERSION"}:
        return "INV"
    if s in {"BND", "TRA", "TRANSLOCATION", "BREAKEND"}:
        return "BND"
    if cn is not None:
        if cn < 2.0:
            return "DEL"
        if cn > 2.0:
            return "DUP"
    return s


def _score_sv_riggs(sv_type: str, length_bp: int, cn: float | None, genes: list[str]) -> dict[str, Any]:
    """Score CNV/SV under ClinGen/ACMG 2020 Riggs dosage-sensitivity criteria."""
    hi3_genes = [g for g in genes if CLINGEN_GENES.get(g, ("", 0, 0, 0, 0, ""))[3] == 3]
    hi2_genes = [g for g in genes if CLINGEN_GENES.get(g, ("", 0, 0, 0, 0, ""))[3] == 2]
    ts3_genes = [g for g in genes if CLINGEN_GENES.get(g, ("", 0, 0, 0, 0, ""))[4] == 3]
    ar_genes = [g for g in genes if CLINGEN_GENES.get(g, ("", 0, 0, 0, 0, ""))[5] == "AR"]

    score = 0.0
    rationale: list[str] = []

    if sv_type == "DEL":
        if hi3_genes:
            score = 1.00
            rationale.append(f"ClinGen HI=3 haploinsufficient gene(s): {','.join(hi3_genes)} (+1.00)")
        elif ar_genes and cn is not None and cn == 0.0:
            score = 0.95
            rationale.append(f"Homozygous deletion (CN=0) of AR gene(s): {','.join(ar_genes)} (+0.95)")
        elif hi2_genes or (genes and length_bp >= 1_000_000):
            score = 0.90
            rationale.append("Emerging HI=2 or >=1 Mb gene-disrupting deletion (+0.90)")
        elif ar_genes:
            score = 0.45
            rationale.append(f"Heterozygous deletion of AR gene(s) {','.join(ar_genes)}; check trans SNV (+0.45)")
        elif genes:
            score = 0.30
            rationale.append("Coding gene deletion without established HI=3 score (+0.30)")
        else:
            score = -0.95 if length_bp < 50_000 else 0.0
            rationale.append("Intergenic/non-coding deletion")
    elif sv_type == "DUP":
        if ts3_genes:
            score = 1.00
            rationale.append(f"ClinGen TS=3 triplosensitive gene(s): {','.join(ts3_genes)} (+1.00)")
        elif hi3_genes and length_bp < 2_000_000:
            score = 0.45
            rationale.append("Duplication overlapping HI gene; uncertain reading-frame disruption (+0.45)")
        elif genes and length_bp >= 2_000_000:
            score = 0.90
            rationale.append("Large (>=2 Mb) coding duplication (+0.90)")
        elif genes:
            score = 0.15
            rationale.append("Coding duplication without established TS=3 score (+0.15)")
        else:
            score = -0.95 if length_bp < 50_000 else 0.0
            rationale.append("Intergenic/non-coding duplication")
    else:  # INV / BND
        if hi3_genes:
            score = 0.90
            rationale.append(f"{sv_type} disrupting HI=3 gene(s): {','.join(hi3_genes)} (+0.90)")
        elif genes:
            score = 0.30
            rationale.append(f"{sv_type} overlapping coding gene(s) (+0.30)")
        else:
            score = 0.0
            rationale.append(f"Intergenic {sv_type}")

    if score >= 0.99:
        acmg_class = "Pathogenic"
    elif score >= 0.90:
        acmg_class = "Likely_pathogenic"
    elif score <= -0.90:
        acmg_class = "Benign"
    else:
        acmg_class = "VUS"

    return {
        "riggs_score": round(score, 2),
        "acmg_cnv_class": acmg_class,
        "clingen_hi_genes": ",".join(hi3_genes) if hi3_genes else "",
        "clingen_ts_genes": ",".join(ts3_genes) if ts3_genes else "",
        "ar_genes_overlapped": ",".join(ar_genes) if ar_genes else "",
        "dosage_rationale": "; ".join(rationale),
    }


def _plot_cnv_summary(sv_df: pd.DataFrame, tier_counts: dict[str, int], n_comp: int, out_png: Path) -> None:
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.5))

    tier_cols = {"Pathogenic": "#dc2626", "Likely_pathogenic": "#ea580c", "VUS": "#d97706", "Benign": "#0f766e"}
    counts = np.array([tier_counts.get(t, 0) for t in SV_TIER_ORDER], dtype=float)
    x = np.arange(len(SV_TIER_ORDER))
    axes[0].bar(x, counts, color=[tier_cols[t] for t in SV_TIER_ORDER], edgecolor="#1e293b", linewidth=0.6)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(SV_TIER_ORDER, rotation=20, ha="right")
    axes[0].set_ylabel("CNV/SV count")
    axes[0].set_title(f"ClinGen Riggs CNV/SV Tiers (SNV+DEL hits: {n_comp})")
    max_c = max(float(np.max(counts)) if counts.size else 1.0, 1.0)
    axes[0].set_ylim(0, max_c * 1.2)
    for xi, c in zip(x, counts):
        axes[0].text(xi, c + max_c * 0.03, f"{int(c)}", ha="center", fontsize=8.5)

    if not sv_df.empty:
        by_type = sv_df["sv_type"].value_counts()
        tx = np.arange(len(by_type))
        tvals = by_type.to_numpy(dtype=float, copy=True)
        axes[1].bar(tx, tvals, color="#0284c7", edgecolor="#1e293b", linewidth=0.6)
        axes[1].set_xticks(tx)
        axes[1].set_xticklabels(by_type.index.tolist())
        axes[1].set_ylabel("Variant count")
        axes[1].set_title("Germline CNV/SV Type Breakdown")
        max_t = max(float(np.max(tvals)) if tvals.size else 1.0, 1.0)
        axes[1].set_ylim(0, max_t * 1.2)
        for xi, c in zip(tx, tvals):
            axes[1].text(xi, c + max_t * 0.03, f"{int(c)}", ha="center", fontsize=8.5)
    else:
        axes[1].text(0.5, 0.5, "No SV calls", ha="center", va="center")

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sv_calls", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--snv_candidates", default=None)
    ap.add_argument("--gene_bed", default=None)
    args = ap.parse_args()

    sv_path = Path(args.sv_calls)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sv_raw = _parse_sv_vcf(sv_path) if is_vcf_file(sv_path) else _parse_sv_table(sv_path)
    if sv_raw.empty:
        fail(f"No SV/CNV calls parsed from {sv_path}")

    intervals = _load_gene_intervals(args.gene_bed)
    warnings: list[str] = []

    sv_rows: list[dict[str, Any]] = []
    for idx in range(len(sv_raw)):
        r = sv_raw.iloc[idx]
        chrom = str(r["chrom"])
        start = int(r["start"])
        end = max(start, int(r["end"]))
        cn = safe_float(r["copy_number"])
        sv_type = _normalize_sv_type(str(r["sv_type"]), cn)
        length_bp = max(1, end - start)

        explicit_genes = [norm_gene(g) for g in re.split(r"[,;|\s]+", str(r["genes"])) if norm_gene(g)]
        ov_genes = [
            g for ic, ist, ien, g in intervals if ic == chrom and not (end < ist or start > ien)
        ]
        all_genes = sorted(set(explicit_genes) | set(ov_genes))

        scored = _score_sv_riggs(sv_type, length_bp, cn, all_genes)
        sv_rows.append(
            {
                "sv_id": f"{chrom}:{start}-{end}:{sv_type}",
                "sample": str(r["sample"]),
                "chrom": chrom,
                "start": start,
                "end": end,
                "sv_length_bp": length_bp,
                "sv_type": sv_type,
                "copy_number": cn if cn is not None else np.nan,
                "genes": ",".join(all_genes),
                **scored,
            }
        )

    sv_df = pd.DataFrame(sv_rows)

    # Detect SNV + CNV Deletion compound heterozygotes
    comp_rows: list[dict[str, Any]] = []
    if args.snv_candidates:
        snv_df = read_table(args.snv_candidates)
        g_col = resolve_column(snv_df, ["gene", "symbol", "gene_name", "hgnc_symbol"], "gene")
        vid_col = resolve_column(snv_df, ["variant_id", "id", "rsid"], "variant_id", required=False)
        csq_col = resolve_column(snv_df, ["consequence", "effect", "annotation"], "consequence", required=False)
        gt_col = resolve_column(snv_df, ["proband_gt", "gt", "genotype"], "proband_gt", required=False)
        acmg_col = resolve_column(snv_df, ["acmg_class"], "acmg_class", required=False)
        samp_col = resolve_column(snv_df, ["sample", "sample_id"], "sample", required=False)
        af_col = resolve_column(snv_df, ["gnomad_af", "af"], "gnomad_af", required=False)

        snv_samples = set(snv_df[samp_col].dropna().astype(str)) if samp_col else set()
        sv_samples = set(sv_df["sample"].dropna().astype(str))

        for _, sv_r in sv_df.iterrows():
            if sv_r["sv_type"] not in {"DEL", "INV", "BND"}:
                continue
            sv_genes = {norm_gene(g) for g in str(sv_r["genes"]).split(",") if norm_gene(g)}
            if not sv_genes:
                continue
            for s_idx in range(len(snv_df)):
                snv_r = snv_df.iloc[s_idx]
                sg = norm_gene(snv_r[g_col])
                if sg not in sv_genes:
                    continue
                if samp_col and len(sv_samples) > 1 and len(snv_samples) > 1:
                    if str(snv_r[samp_col]).strip() != str(sv_r["sample"]).strip():
                        continue
                af_v = safe_float(snv_r[af_col]) if af_col else None
                if af_v is not None and af_v > 0.01:
                    continue
                snv_gt = str(snv_r[gt_col]).strip() if gt_col and pd.notna(snv_r[gt_col]) else "0/1"
                # Apparent 1/1 SNV in a hemizygous deletion region can reflect allelic dropout / LOH over DEL
                mechanism = (
                    "apparent_hom_snv_over_del" if snv_gt in {"1/1", "1|1"} else "trans_snv_plus_del_compound_het"
                )
                vid = (
                    str(snv_r[vid_col]).strip()
                    if vid_col and pd.notna(snv_r[vid_col])
                    else f"{sg}_snv_{s_idx+1}"
                )
                comp_rows.append(
                    {
                        "sample": str(sv_r["sample"]),
                        "gene": sg,
                        "compound_mechanism": mechanism,
                        "snv_variant_id": vid,
                        "snv_consequence": str(snv_r[csq_col]) if csq_col and pd.notna(snv_r[csq_col]) else "",
                        "snv_gt": snv_gt,
                        "snv_acmg_class": str(snv_r[acmg_col]) if acmg_col and pd.notna(snv_r[acmg_col]) else "",
                        "sv_id": str(sv_r["sv_id"]),
                        "sv_type": str(sv_r["sv_type"]),
                        "sv_copy_number": sv_r["copy_number"],
                        "sv_acmg_class": str(sv_r["acmg_cnv_class"]),
                        "clinical_note": (
                            f"Biallelic compound hit in {sg}: {vid} ({snv_gt}) + {sv_r['sv_id']} ({sv_r['sv_type']}); confirm trans phase if parental samples available."
                        ),
                    }
                )

    comp_cols = [
        "sample", "gene", "compound_mechanism", "snv_variant_id", "snv_consequence",
        "snv_gt", "snv_acmg_class", "sv_id", "sv_type", "sv_copy_number", "sv_acmg_class", "clinical_note",
    ]
    comp_df = pd.DataFrame(comp_rows, columns=comp_cols) if comp_rows else pd.DataFrame(columns=comp_cols)

    # Elevate SV rows involved in compound hits if currently VUS
    if not comp_df.empty:
        comp_sv_ids = set(comp_df["sv_id"].astype(str))
        for i in range(len(sv_df)):
            if str(sv_df.loc[i, "sv_id"]) in comp_sv_ids and str(sv_df.loc[i, "acmg_cnv_class"]) == "VUS":
                sv_df.loc[i, "acmg_cnv_class"] = "Likely_pathogenic"
                sv_df.loc[i, "riggs_score"] = max(float(sv_df.loc[i, "riggs_score"]), 0.90)
                sv_df.loc[i, "dosage_rationale"] = (
                    str(sv_df.loc[i, "dosage_rationale"]) + "; Elevated via SNV+DEL compound heterozygote match"
                )
        warnings.append(
            f"Detected {len(comp_df)} SNV+CNV compound heterozygote hit(s) in: {', '.join(sorted(comp_df['gene'].unique()))}."
        )

    sv_df = sv_df.sort_values(by=["riggs_score", "sv_length_bp"], ascending=[False, False]).reset_index(drop=True)
    tier_counts = {t: int((sv_df["acmg_cnv_class"] == t).sum()) for t in SV_TIER_ORDER}
    n_plp = int(tier_counts["Pathogenic"] + tier_counts["Likely_pathogenic"])

    sv_tsv = out_dir / "classified_germline_cnv_sv.tsv"
    comp_tsv = out_dir / "snv_cnv_compound_hits.tsv"
    plot_png = out_dir / "germline_cnv_dosage_summary.png"

    _plot_cnv_summary(sv_df, tier_counts, int(len(comp_df)), plot_png)

    sci_ctx = {"tier_counts": tier_counts, "n_snv_cnv_compound_hits": int(len(comp_df))}
    write_tsv(
        sv_df,
        sv_tsv,
        tool="germ_cnv_sv_integrate",
        title="ClinGen Riggs Classified Germline CNV/SV Calls",
        inputs=[str(sv_path), args.snv_candidates, args.gene_bed],
        scientific_context=sci_ctx,
    )
    write_tsv(
        comp_df,
        comp_tsv,
        tool="germ_cnv_sv_integrate",
        title="SNV + CNV/DEL Compound Heterozygote Hits",
        inputs=[str(sv_path), args.snv_candidates, args.gene_bed],
        scientific_context=sci_ctx,
    )
    write_artifact(
        plot_png,
        tool="germ_cnv_sv_integrate",
        title="Germline CNV/SV Dosage Classification & Type Summary",
        inputs=[str(sv_path), args.snv_candidates, args.gene_bed],
        scientific_context=sci_ctx,
    )

    emit(
        {
            "status": "ok",
            "tool": "germ_cnv_sv_integrate",
            "n_sv_calls": int(len(sv_df)),
            "tier_counts": tier_counts,
            "n_pathogenic_or_likely": n_plp,
            "n_snv_cnv_compound_hits": int(len(comp_df)),
            "compound_hits": comp_df.to_dict(orient="records"),
            "outputs": {
                "classified_germline_cnv_sv": str(sv_tsv),
                "snv_cnv_compound_hits": str(comp_tsv),
                "germline_cnv_dosage_summary": str(plot_png),
            },
            "artifacts": [str(sv_tsv), str(comp_tsv), str(plot_png)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
