#!/usr/bin/env python3
"""germ_vcf_filter: Filter single-sample or trio annotated VCF/VCF.gz (INFO/CSQ or
INFO/ANN) or variant TSV/CSV by DP, GQ, coding/splice consequence, and gnomAD AF,
and evaluate Mendelian inheritance fit (de_novo, autosomal_recessive_hom,
phased in-trans compound_heterozygous vs singleton unphased_candidate_comp_het,
autosomal_dominant, and x_linked).
"""
from __future__ import annotations

import argparse
import math
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
    classify_consequence_group,
    emit,
    fail,
    is_coding_or_splice_consequence,
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

DEFAULT_CSQ_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
    "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
    "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE",
    "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "gnomAD_AF", "AF", "CLIN_SIG",
    "REVEL", "CADD_PHRED", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "pLI", "LOEUF",
]

DEFAULT_ANN_FIELDS = [
    "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID",
    "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p",
    "cDNA.pos", "CDS.pos", "AA.pos", "Distance", "ERRORS",
]

IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "": 0}


def _parse_csq_format(header_line: str) -> list[str] | None:
    m = re.search(r"Format:\s*([^\"'>]+)", header_line)
    if not m:
        return None
    fields = [f.strip() for f in m.group(1).split("|") if f.strip()]
    return fields or None


def _parse_spliceai(val: Any) -> float | None:
    if val is None:
        return None
    s = str(val).strip()
    if not s or s in {".", "NA", "nan", "None"}:
        return None
    # Handle pipe-delimited VEP/Illumina SpliceAI string: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|...
    if "|" in s:
        parts = s.split("|")
        nums = [safe_float(x) for x in parts[2:6]]
        valid = [x for x in nums if x is not None]
        if valid:
            return float(max(valid))
    nums = [safe_float(x) for x in re.split(r"[,;|\s]+", s)]
    valid = [x for x in nums if x is not None]
    return float(max(valid)) if valid else None


def _norm_gt(raw_gt: Any) -> str:
    """Normalize GT string to '0/0', '0/1', '1/1', '0', '1', or './.'."""
    if raw_gt is None or (isinstance(raw_gt, float) and math.isnan(raw_gt)):
        return "./."
    s = str(raw_gt).strip()
    if not s or s in {".", "./.", ".|.", "NA", "nan", "None", "missing"}:
        return "./."
    low = s.lower()
    if low in {"hom_ref", "ref"}:
        return "0/0"
    if low in {"het", "heterozygous"}:
        return "0/1"
    if low in {"hom_alt", "hom", "homozygous"}:
        return "1/1"
    if low in {"hemi", "hemizygous"}:
        return "1"
    sep = "/" if "/" in s else ("|" if "|" in s else None)
    if sep is None:
        if s == "0":
            return "0"
        if s.isdigit() and int(s) > 0:
            return "1"
        return "./."
    alleles = s.split(sep)
    if len(alleles) != 2 or "." in alleles:
        return "./."
    try:
        a1 = int(alleles[0])
        a2 = int(alleles[1])
    except ValueError:
        return "./."
    if a1 == 0 and a2 == 0:
        return "0/0"
    if (a1 == 0 and a2 > 0) or (a2 == 0 and a1 > 0):
        return "0/1"
    return "1/1"


def _gt_carries_alt(gt: str) -> bool:
    return gt in {"0/1", "1/1", "1"}


def _gt_is_hom_ref(gt: str) -> bool:
    return gt in {"0/0", "0"}


def _gt_is_called(gt: str) -> bool:
    return gt != "./."


def _parse_sample_format(fmt_keys: list[str], sample_raw: str) -> dict[str, Any]:
    vals = sample_raw.split(":")
    fmap = {k: vals[i] if i < len(vals) else "." for i, k in enumerate(fmt_keys)}
    gt = _norm_gt(fmap.get("GT", "./."))
    dp = safe_int(fmap.get("DP"))
    gq = safe_int(fmap.get("GQ"))
    vaf = None
    ad_raw = fmap.get("AD", "")
    if ad_raw and ad_raw != ".":
        ad_parts = [safe_int(x) for x in ad_raw.split(",")]
        valid_ad = [x for x in ad_parts if x is not None and x >= 0]
        if valid_ad:
            tot_ad = sum(valid_ad)
            if dp is None and tot_ad > 0:
                dp = int(tot_ad)
            if len(valid_ad) >= 2 and tot_ad > 0:
                vaf = float(sum(valid_ad[1:]) / tot_ad)
    if vaf is None and "AF" in fmap:
        vaf = safe_float(fmap.get("AF"))
    return {"gt": gt, "dp": dp, "gq": gq, "vaf": vaf}


def _select_best_transcript(entries: list[dict[str, str]], alt_allele: str) -> dict[str, str]:
    if not entries:
        return {}
    best = entries[0]
    best_score = (-1, -1, -1)
    for ent in entries:
        impact = str(ent.get("IMPACT") or ent.get("Annotation_Impact") or "").upper()
        gene = norm_gene(ent.get("SYMBOL") or ent.get("Gene_Name") or ent.get("Gene"))
        csq = ent.get("Consequence") or ent.get("Annotation") or ""
        allele_match = 1 if (not ent.get("Allele") or ent.get("Allele") == alt_allele) else 0
        coding_bonus = 1 if is_coding_or_splice_consequence(csq) else 0
        score = (
            allele_match,
            IMPACT_RANK.get(impact, 0) * 2 + coding_bonus,
            1 if gene else 0,
        )
        if score > best_score:
            best_score = score
            best = ent
    return best


def _parse_vcf(
    path: Path,
    proband_id: str | None,
    mother_id: str | None,
    father_id: str | None,
) -> tuple[pd.DataFrame, str, str | None, str | None, list[str]]:
    lines = read_text_lines(path)
    csq_fields: list[str] = list(DEFAULT_CSQ_FIELDS)
    ann_fields: list[str] = list(DEFAULT_ANN_FIELDS)
    samples: list[str] = []
    records: list[dict[str, Any]] = []
    warnings: list[str] = []

    header_found = False
    for line in lines:
        s = line.rstrip("\r\n")
        if not s:
            continue
        if s.startswith("##"):
            if "ID=CSQ" in s:
                parsed = _parse_csq_format(s)
                if parsed:
                    csq_fields = parsed
            elif "ID=ANN" in s:
                parsed = _parse_csq_format(s)
                if parsed:
                    ann_fields = parsed
            continue
        if s.startswith("#CHROM"):
            cols = s.split("\t")
            header_found = True
            if len(cols) > 9:
                samples = [c.strip() for c in cols[9:]]
            break

    if not header_found:
        fail(f"Invalid VCF file (missing #CHROM header): {path}")

    if not samples:
        resolved_proband = proband_id or "PROBAND"
        warnings.append("VCF header contains no sample columns; treating all records as proband carriers.")
    else:
        if proband_id:
            if proband_id not in samples:
                fail(f"proband_id '{proband_id}' not found in VCF samples {samples}")
            resolved_proband = proband_id
        else:
            resolved_proband = samples[0]
            if len(samples) > 1:
                warnings.append(
                    f"Multi-sample VCF without --proband_id; automatically selected first sample '{resolved_proband}'."
                )
    if mother_id and mother_id not in samples:
        fail(f"mother_id '{mother_id}' not found in VCF samples {samples}")
    if father_id and father_id not in samples:
        fail(f"father_id '{father_id}' not found in VCF samples {samples}")

    sample_idx = {name: 9 + idx for idx, name in enumerate(samples)}

    for line in lines:
        s = line.rstrip("\r\n")
        if not s or s.startswith("#"):
            continue
        cols = s.split("\t")
        if len(cols) < 8:
            continue
        chrom = norm_chrom(cols[0])
        pos = safe_int(cols[1], 0) or 0
        vid = cols[2].strip() if cols[2].strip() != "." else ""
        ref = cols[3].strip()
        alt_raw = cols[4].strip()
        alt = alt_raw.split(",")[0].strip()
        flt = cols[6].strip()
        info_raw = cols[7].strip()

        info_map: dict[str, str] = {}
        if info_raw and info_raw != ".":
            for item in info_raw.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_map[k.strip()] = v.strip()
                else:
                    info_map[item.strip()] = "True"

        ann_best: dict[str, str] = {}
        if "CSQ" in info_map:
            entries = []
            for rec in info_map["CSQ"].split(","):
                parts = rec.split("|")
                entries.append({csq_fields[i]: parts[i] if i < len(parts) else "" for i in range(len(csq_fields))})
            ann_best = _select_best_transcript(entries, alt)
        elif "ANN" in info_map:
            entries = []
            for rec in info_map["ANN"].split(","):
                parts = rec.split("|")
                entries.append({ann_fields[i]: parts[i] if i < len(parts) else "" for i in range(len(ann_fields))})
            ann_best = _select_best_transcript(entries, alt)

        gene = norm_gene(
            ann_best.get("SYMBOL")
            or ann_best.get("Gene_Name")
            or info_map.get("SYMBOL")
            or info_map.get("GENE")
            or info_map.get("Gene")
        )
        consequence = (
            ann_best.get("Consequence")
            or ann_best.get("Annotation")
            or info_map.get("Consequence")
            or info_map.get("FUNC")
            or ""
        )
        impact = (
            ann_best.get("IMPACT")
            or ann_best.get("Annotation_Impact")
            or info_map.get("IMPACT")
            or ""
        )
        hgvsc = ann_best.get("HGVSc") or ann_best.get("HGVS.c") or info_map.get("HGVSc") or ""
        hgvsp = ann_best.get("HGVSp") or ann_best.get("HGVS.p") or info_map.get("HGVSp") or ""

        # Population AF preference: gnomAD / MAX_AF over generic cohort AF
        gnomad_af = None
        for af_key in ("gnomAD_AF", "gnomADe_AF", "gnomADg_AF", "MAX_AF", "gnomad_af", "AF"):
            cand_v = ann_best.get(af_key) if af_key in ann_best else info_map.get(af_key)
            val = safe_float(cand_v)
            if val is not None:
                gnomad_af = val
                break

        clinvar_sig = (
            ann_best.get("CLIN_SIG")
            or info_map.get("CLIN_SIG")
            or info_map.get("CLNSIG")
            or info_map.get("clinvar_sig")
            or ""
        )
        revel = safe_float(ann_best.get("REVEL") or info_map.get("REVEL"))
        cadd = safe_float(ann_best.get("CADD_PHRED") or info_map.get("CADD_PHRED") or info_map.get("CADD"))
        splice_vals = [
            _parse_spliceai(ann_best.get(k) or info_map.get(k))
            for k in (
                "SpliceAI_max",
                "SPLICEAI",
                "SpliceAI",
                "SpliceAI_pred_DS_AG",
                "SpliceAI_pred_DS_AL",
                "SpliceAI_pred_DS_DG",
                "SpliceAI_pred_DS_DL",
            )
        ]
        valid_sp = [x for x in splice_vals if x is not None]
        spliceai_max = float(max(valid_sp)) if valid_sp else None
        pli = safe_float(ann_best.get("pLI") or info_map.get("pLI") or info_map.get("PLI"))
        loeuf = safe_float(ann_best.get("LOEUF") or info_map.get("LOEUF") or info_map.get("loeuf"))

        fmt_keys = [k.strip() for k in cols[8].split(":")] if len(cols) > 8 else []
        if samples and resolved_proband in sample_idx and sample_idx[resolved_proband] < len(cols):
            p_info = _parse_sample_format(fmt_keys, cols[sample_idx[resolved_proband]])
        else:
            p_info = {"gt": "0/1", "dp": safe_int(info_map.get("DP")), "gq": None, "vaf": None}

        m_info = (
            _parse_sample_format(fmt_keys, cols[sample_idx[mother_id]])
            if mother_id and mother_id in sample_idx and sample_idx[mother_id] < len(cols)
            else {"gt": "./.", "dp": None, "gq": None, "vaf": None}
        )
        f_info = (
            _parse_sample_format(fmt_keys, cols[sample_idx[father_id]])
            if father_id and father_id in sample_idx and sample_idx[father_id] < len(cols)
            else {"gt": "./.", "dp": None, "gq": None, "vaf": None}
        )

        var_id = vid if vid else f"{chrom}:{pos}:{ref}>{alt}"
        records.append(
            {
                "variant_id": var_id,
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "filter": flt,
                "gene": gene,
                "consequence": consequence,
                "impact": impact,
                "hgvsc": hgvsc,
                "hgvsp": hgvsp,
                "gnomad_af": gnomad_af if gnomad_af is not None else np.nan,
                "proband_gt": p_info["gt"],
                "proband_dp": p_info["dp"] if p_info["dp"] is not None else np.nan,
                "proband_gq": p_info["gq"] if p_info["gq"] is not None else np.nan,
                "proband_vaf": p_info["vaf"] if p_info["vaf"] is not None else np.nan,
                "mother_gt": m_info["gt"],
                "mother_dp": m_info["dp"] if m_info["dp"] is not None else np.nan,
                "mother_gq": m_info["gq"] if m_info["gq"] is not None else np.nan,
                "father_gt": f_info["gt"],
                "father_dp": f_info["dp"] if f_info["dp"] is not None else np.nan,
                "father_gq": f_info["gq"] if f_info["gq"] is not None else np.nan,
                "clinvar_sig": clinvar_sig,
                "revel": revel if revel is not None else np.nan,
                "cadd_phred": cadd if cadd is not None else np.nan,
                "spliceai_max": spliceai_max if spliceai_max is not None else np.nan,
                "pli": pli if pli is not None else np.nan,
                "loeuf": loeuf if loeuf is not None else np.nan,
            }
        )

    if not records:
        fail(f"VCF file contains 0 variant records: {path}")
    return pd.DataFrame(records), resolved_proband, mother_id, father_id, warnings


def _parse_tabular(
    path: Path,
    proband_id: str | None,
    mother_id: str | None,
    father_id: str | None,
) -> tuple[pd.DataFrame, str, str | None, str | None, list[str]]:
    df = read_table(path)
    warnings: list[str] = []

    gene_col = resolve_column(df, ["gene", "symbol", "gene_name", "hgnc_symbol"], "gene")
    csq_col = resolve_column(df, ["consequence", "effect", "annotation", "variant_class"], "consequence")
    af_col = resolve_column(df, ["gnomad_af", "af", "popmax_af", "max_af"], "gnomad_af", required=False)

    chrom_col = resolve_column(df, ["chrom", "chr", "#chrom", "chromosome"], "chrom", required=False)
    pos_col = resolve_column(df, ["pos", "position", "start"], "pos", required=False)
    ref_col = resolve_column(df, ["ref", "reference"], "ref", required=False)
    alt_col = resolve_column(df, ["alt", "alternate"], "alt", required=False)
    vid_col = resolve_column(df, ["variant_id", "id", "rsid"], "variant_id", required=False)
    flt_col = resolve_column(df, ["filter"], "filter", required=False)
    imp_col = resolve_column(df, ["impact", "annotation_impact"], "impact", required=False)
    hgvsc_col = resolve_column(df, ["hgvsc", "hgvs_c", "cdna_change"], "hgvsc", required=False)
    hgvsp_col = resolve_column(df, ["hgvsp", "hgvs_p", "aa_change"], "hgvsp", required=False)
    clin_col = resolve_column(df, ["clinvar_sig", "clnsig", "clin_sig", "clinvar"], "clinvar_sig", required=False)
    rev_col = resolve_column(df, ["revel", "revel_score"], "revel", required=False)
    cadd_col = resolve_column(df, ["cadd_phred", "cadd", "cadd_score"], "cadd_phred", required=False)
    sp_col = resolve_column(df, ["spliceai_max", "spliceai", "spliceai_score"], "spliceai_max", required=False)
    pli_col = resolve_column(df, ["pli"], "pli", required=False)
    loeuf_col = resolve_column(df, ["loeuf"], "loeuf", required=False)

    sample_col = resolve_column(df, ["sample_id", "sample"], "sample_id", required=False)
    if proband_id:
        resolved_proband = proband_id
    elif sample_col and df[sample_col].dropna().nunique() == 1:
        resolved_proband = str(df[sample_col].dropna().iloc[0]).strip()
    else:
        resolved_proband = "PROBAND"

    pgt_cands = [f"{resolved_proband}_gt", "proband_gt", "gt", "genotype"]
    pdp_cands = [f"{resolved_proband}_dp", "proband_dp", "dp", "depth"]
    pgq_cands = [f"{resolved_proband}_gq", "proband_gq", "gq"]
    pgt_col = resolve_column(df, pgt_cands, "proband_gt", required=False)
    pdp_col = resolve_column(df, pdp_cands, "proband_dp", required=False)
    pgq_col = resolve_column(df, pgq_cands, "proband_gq", required=False)
    pvaf_col = resolve_column(df, [f"{resolved_proband}_vaf", "proband_vaf", "vaf", "af_sample"], "proband_vaf", required=False)

    mgt_cands = ([f"{mother_id}_gt"] if mother_id else []) + ["mother_gt", "mat_gt"]
    mdp_cands = ([f"{mother_id}_dp"] if mother_id else []) + ["mother_dp", "mat_dp"]
    mgq_cands = ([f"{mother_id}_gq"] if mother_id else []) + ["mother_gq", "mat_gq"]
    mgt_col = resolve_column(df, mgt_cands, "mother_gt", required=False)
    mdp_col = resolve_column(df, mdp_cands, "mother_dp", required=False)
    mgq_col = resolve_column(df, mgq_cands, "mother_gq", required=False)

    fgt_cands = ([f"{father_id}_gt"] if father_id else []) + ["father_gt", "pat_gt"]
    fdp_cands = ([f"{father_id}_dp"] if father_id else []) + ["father_dp", "pat_dp"]
    fgq_cands = ([f"{father_id}_gq"] if father_id else []) + ["father_gq", "pat_gq"]
    fgt_col = resolve_column(df, fgt_cands, "father_gt", required=False)
    fdp_col = resolve_column(df, fdp_cands, "father_dp", required=False)
    fgq_col = resolve_column(df, fgq_cands, "father_gq", required=False)

    if pgt_col is None:
        warnings.append("Tabular input lacks genotype (GT) column; assuming proband heterozygous carrier (0/1).")

    records: list[dict[str, Any]] = []
    for idx in range(len(df)):
        row = df.iloc[idx]
        chrom = norm_chrom(row[chrom_col]) if chrom_col else ""
        pos = safe_int(row[pos_col], idx + 1) if pos_col else idx + 1
        ref = str(row[ref_col]).strip() if ref_col and pd.notna(row[ref_col]) else "N"
        alt = str(row[alt_col]).strip() if alt_col and pd.notna(row[alt_col]) else "N"
        gene = norm_gene(row[gene_col])
        csq = "" if pd.isna(row[csq_col]) else str(row[csq_col]).strip()
        vid = str(row[vid_col]).strip() if vid_col and pd.notna(row[vid_col]) else f"{chrom or 'chr?'}:{pos}:{ref}>{alt}"
        pgt = _norm_gt(row[pgt_col]) if pgt_col else "0/1"
        mgt = _norm_gt(row[mgt_col]) if mgt_col else "./."
        fgt = _norm_gt(row[fgt_col]) if fgt_col else "./."

        records.append(
            {
                "variant_id": vid,
                "chrom": chrom,
                "pos": int(pos or 0),
                "ref": ref,
                "alt": alt,
                "filter": str(row[flt_col]).strip() if flt_col and pd.notna(row[flt_col]) else "PASS",
                "gene": gene,
                "consequence": csq,
                "impact": str(row[imp_col]).strip() if imp_col and pd.notna(row[imp_col]) else "",
                "hgvsc": str(row[hgvsc_col]).strip() if hgvsc_col and pd.notna(row[hgvsc_col]) else "",
                "hgvsp": str(row[hgvsp_col]).strip() if hgvsp_col and pd.notna(row[hgvsp_col]) else "",
                "gnomad_af": safe_float(row[af_col]) if af_col else np.nan,
                "proband_gt": pgt,
                "proband_dp": safe_int(row[pdp_col]) if pdp_col else np.nan,
                "proband_gq": safe_int(row[pgq_col]) if pgq_col else np.nan,
                "proband_vaf": safe_float(row[pvaf_col]) if pvaf_col else np.nan,
                "mother_gt": mgt,
                "mother_dp": safe_int(row[mdp_col]) if mdp_col else np.nan,
                "mother_gq": safe_int(row[mgq_col]) if mgq_col else np.nan,
                "father_gt": fgt,
                "father_dp": safe_int(row[fdp_col]) if fdp_col else np.nan,
                "father_gq": safe_int(row[fgq_col]) if fgq_col else np.nan,
                "clinvar_sig": str(row[clin_col]).strip() if clin_col and pd.notna(row[clin_col]) else "",
                "revel": safe_float(row[rev_col]) if rev_col else np.nan,
                "cadd_phred": safe_float(row[cadd_col]) if cadd_col else np.nan,
                "spliceai_max": _parse_spliceai(row[sp_col]) if sp_col else np.nan,
                "pli": safe_float(row[pli_col]) if pli_col else np.nan,
                "loeuf": safe_float(row[loeuf_col]) if loeuf_col else np.nan,
            }
        )
    return pd.DataFrame(records), resolved_proband, mother_id or ("MOTHER" if mgt_col else None), father_id or ("FATHER" if fgt_col else None), warnings


def _passes_qc(row: pd.Series, min_dp: int, min_gq: int) -> bool:
    flt = str(row["filter"]).strip().upper()
    if flt not in {"PASS", ".", "", "NONE", "NAN"}:
        return False
    dp = safe_float(row["proband_dp"])
    gq = safe_float(row["proband_gq"])
    if dp is not None and dp < min_dp:
        return False
    if gq is not None and gq < min_gq:
        return False
    return True


def _parent_qc_ok(row: pd.Series, min_dp: int, min_gq: int) -> bool:
    for pfx in ("mother", "father"):
        dp = safe_float(row[f"{pfx}_dp"])
        gq = safe_float(row[f"{pfx}_gq"])
        if dp is not None and dp < min_dp:
            return False
        if gq is not None and gq < min_gq:
            return False
    return True


def _plot_funnel(funnel_counts: list[tuple[str, int]], inh_counts: dict[str, int], out_png: Path) -> None:
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

    labels = [x[0] for x in funnel_counts][::-1]
    vals = np.array([x[1] for x in funnel_counts][::-1], dtype=float)
    y_pos = np.arange(len(labels))
    colors = ["#0f766e", "#0d9488", "#0284c7", "#3b82f6", "#64748b"][: len(labels)]
    axes[0].barh(y_pos, vals, color=colors, edgecolor="#1e293b", linewidth=0.6)
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(labels)
    axes[0].set_xlabel("Variant count")
    axes[0].set_title("Germline Variant Filtering Funnel")
    max_v = max(float(np.max(vals)) if vals.size else 1.0, 1.0)
    axes[0].set_xlim(0, max_v * 1.18)
    for y, v in zip(y_pos, vals):
        axes[0].text(v + max_v * 0.02, y, f"{int(v)}", va="center", fontsize=8.5)

    if inh_counts:
        modes = list(inh_counts.keys())
        counts = np.array([inh_counts[m] for m in modes], dtype=float)
        x_pos = np.arange(len(modes))
        palette = {
            "de_novo": "#dc2626",
            "autosomal_recessive_hom": "#7c3aed",
            "compound_heterozygous": "#0f766e",
            "unphased_candidate_comp_het": "#d97706",
            "autosomal_dominant": "#2563eb",
            "x_linked": "#db2777",
        }
        bar_cols = [palette.get(m, "#475569") for m in modes]
        axes[1].bar(x_pos, counts, color=bar_cols, edgecolor="#1e293b", linewidth=0.6)
        axes[1].set_xticks(x_pos)
        axes[1].set_xticklabels(modes, rotation=28, ha="right")
        axes[1].set_ylabel("Candidate variants")
        axes[1].set_title("Mendelian Inheritance Fit")
        max_c = max(float(np.max(counts)) if counts.size else 1.0, 1.0)
        axes[1].set_ylim(0, max_c * 1.2)
        for x, c in zip(x_pos, counts):
            axes[1].text(x, c + max_c * 0.03, f"{int(c)}", ha="center", fontsize=8.5)
    else:
        axes[1].text(0.5, 0.5, "No Mendelian candidates retained", ha="center", va="center")
        axes[1].set_title("Mendelian Inheritance Fit")

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--proband_id", default=None)
    ap.add_argument("--mother_id", default=None)
    ap.add_argument("--father_id", default=None)
    ap.add_argument("--proband_sex", choices=["unknown", "male", "female"], default="unknown")
    ap.add_argument("--max_af_dominant", type=float, default=0.0001)
    ap.add_argument("--max_af_recessive", type=float, default=0.01)
    ap.add_argument("--min_dp", type=int, default=10)
    ap.add_argument("--min_gq", type=int, default=20)
    args = ap.parse_args()

    vpath = Path(args.variants)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if is_vcf_file(vpath):
        input_format = "vcf"
        df, proband_id, mother_id, father_id, warnings = _parse_vcf(
            vpath, args.proband_id, args.mother_id, args.father_id
        )
    else:
        input_format = "tabular"
        df, proband_id, mother_id, father_id, warnings = _parse_tabular(
            vpath, args.proband_id, args.mother_id, args.father_id
        )

    total_variants = int(len(df))
    carrier_mask = df["proband_gt"].map(_gt_carries_alt)
    df_car = df[carrier_mask].copy().reset_index(drop=True)
    proband_carrier_n = int(len(df_car))

    qc_mask = df_car.apply(lambda r: _passes_qc(r, args.min_dp, args.min_gq), axis=1)
    df_qc = df_car[qc_mask].copy().reset_index(drop=True)
    qc_pass_n = int(len(df_qc))

    def _coding_splice_row(r: pd.Series) -> bool:
        if is_coding_or_splice_consequence(r["consequence"]):
            return True
        sp = safe_float(r["spliceai_max"])
        if sp is not None and sp >= 0.20:
            return True
        csig = str(r["clinvar_sig"]).lower()
        if "pathogenic" in csig and "conflicting" not in csig and "benign" not in csig:
            return True
        return False

    cs_mask = df_qc.apply(_coding_splice_row, axis=1) if not df_qc.empty else pd.Series(dtype=bool)
    df_cs = df_qc[cs_mask].copy().reset_index(drop=True) if not df_qc.empty else df_qc.copy()
    coding_splice_n = int(len(df_cs))

    # Rare gate (<= max_af_recessive; missing AF treated as 0.0)
    af_eff = df_cs["gnomad_af"].apply(lambda x: 0.0 if safe_float(x) is None else float(safe_float(x)))
    df_cs["gnomad_af_effective"] = af_eff
    df_rare = df_cs[df_cs["gnomad_af_effective"] <= args.max_af_recessive].copy().reset_index(drop=True)

    # Detect trio availability
    has_mother = bool((df_rare["mother_gt"].map(_gt_is_called)).any()) if not df_rare.empty else bool(mother_id)
    has_father = bool((df_rare["father_gt"].map(_gt_is_called)).any()) if not df_rare.empty else bool(father_id)
    trio_mode = bool(has_mother and has_father)

    # First pass: determine parental origin & per-variant individual modes
    rows_meta: list[dict[str, Any]] = []
    for idx in range(len(df_rare)):
        r = df_rare.iloc[idx]
        chrom = str(r["chrom"])
        pgt = str(r["proband_gt"])
        mgt = str(r["mother_gt"])
        fgt = str(r["father_gt"])
        af = float(r["gnomad_af_effective"])
        both_parents_called = _gt_is_called(mgt) and _gt_is_called(fgt)
        parent_qc = _parent_qc_ok(r, args.min_dp, args.min_gq)

        origin = "unknown"
        if both_parents_called:
            m_alt = _gt_carries_alt(mgt)
            f_alt = _gt_carries_alt(fgt)
            if _gt_is_hom_ref(mgt) and _gt_is_hom_ref(fgt):
                origin = "de_novo" if parent_qc else "low_qc_de_novo"
            elif m_alt and not f_alt:
                origin = "maternal_only"
            elif f_alt and not m_alt:
                origin = "paternal_only"
            elif m_alt and f_alt:
                origin = "biparental"

        compat: list[str] = []
        is_x = chrom == "chrX"
        if origin == "de_novo" and pgt == "0/1" and af <= args.max_af_dominant:
            compat.append("de_novo")
        if is_x:
            if args.proband_sex == "male" and pgt in {"1", "1/1", "0/1"} and af <= args.max_af_recessive:
                compat.append("x_linked")
            elif args.proband_sex == "female" and af <= args.max_af_dominant:
                compat.append("x_linked")
            elif args.proband_sex == "unknown" and af <= args.max_af_recessive:
                compat.append("x_linked")
        if not is_x and pgt == "1/1" and af <= args.max_af_recessive:
            if not (both_parents_called and _gt_is_hom_ref(mgt) and _gt_is_hom_ref(fgt)):
                compat.append("autosomal_recessive_hom")
        if not is_x and af <= args.max_af_dominant:
            compat.append("autosomal_dominant")

        rows_meta.append({"origin": origin, "compat": compat, "parent_qc": parent_qc, "both_parents_called": both_parents_called})

    # Second pass: gene-level compound heterozygote phasing
    gene_groups: dict[str, list[int]] = {}
    for idx in range(len(df_rare)):
        g = str(df_rare.iloc[idx]["gene"]).strip()
        if g and str(df_rare.iloc[idx]["chrom"]) != "chrX":
            gene_groups.setdefault(g, []).append(idx)

    comp_het_status: dict[int, tuple[str, str]] = {}  # idx -> (mode, phasing_status)
    unphased_warned = False
    cis_warned = False
    for gene, idxs in gene_groups.items():
        if len(idxs) < 2:
            continue
        if trio_mode:
            mat_idxs = [i for i in idxs if rows_meta[i]["origin"] == "maternal_only"]
            pat_idxs = [i for i in idxs if rows_meta[i]["origin"] == "paternal_only"]
            dn_idxs = [i for i in idxs if rows_meta[i]["origin"] == "de_novo"]
            # Phased in trans if at least one maternal-only AND one paternal-only (or de novo paired with parental)
            if (mat_idxs and pat_idxs) or (mat_idxs and dn_idxs) or (pat_idxs and dn_idxs):
                trans_set = set(mat_idxs + pat_idxs + dn_idxs)
                for i in trans_set:
                    comp_het_status[i] = ("compound_heterozygous", "phased_in_trans")
                    if "compound_heterozygous" not in rows_meta[i]["compat"]:
                        rows_meta[i]["compat"].append("compound_heterozygous")
            else:
                all_mat = all(rows_meta[i]["origin"] == "maternal_only" for i in idxs)
                all_pat = all(rows_meta[i]["origin"] == "paternal_only" for i in idxs)
                if all_mat or all_pat:
                    cis_warned = True
                    for i in idxs:
                        comp_het_status[i] = ("in_cis_excluded", "in_cis_excluded")
                else:
                    unphased_warned = True
                    for i in idxs:
                        comp_het_status[i] = ("unphased_candidate_comp_het", "unphased_trio_ambiguous")
                        if "unphased_candidate_comp_het" not in rows_meta[i]["compat"]:
                            rows_meta[i]["compat"].append("unphased_candidate_comp_het")
        else:
            unphased_warned = True
            for i in idxs:
                comp_het_status[i] = ("unphased_candidate_comp_het", "unphased_singleton")
                if "unphased_candidate_comp_het" not in rows_meta[i]["compat"]:
                    rows_meta[i]["compat"].append("unphased_candidate_comp_het")

    if unphased_warned:
        warnings.append(
            "One or more multi-hit genes lack full parental phasing and are flagged as 'unphased_candidate_comp_het'; parental or long-read phasing is required to confirm trans configuration."
        )
    if cis_warned:
        warnings.append("Excluded same-parent in-cis variant pair(s) from compound_heterozygous calls.")

    # Assign primary Mendelian inheritance_mode & phasing_status and keep candidates
    priority_order = [
        "de_novo",
        "x_linked",
        "autosomal_recessive_hom",
        "compound_heterozygous",
        "unphased_candidate_comp_het",
        "autosomal_dominant",
    ]
    kept_rows: list[dict[str, Any]] = []
    for idx in range(len(df_rare)):
        r = df_rare.iloc[idx].to_dict()
        meta = rows_meta[idx]
        compat = meta["compat"]
        if not compat:
            continue
        primary = next((m for m in priority_order if m in compat), compat[0])
        if primary == "de_novo":
            phasing = "trio_de_novo"
        elif primary == "compound_heterozygous":
            phasing = "phased_in_trans"
        elif primary == "unphased_candidate_comp_het":
            phasing = comp_het_status[idx][1]
        elif primary == "autosomal_recessive_hom":
            phasing = "trio_biparental" if meta["origin"] == "biparental" else ("unphased_singleton" if not trio_mode else meta["origin"])
        elif idx in comp_het_status and comp_het_status[idx][0] == "in_cis_excluded":
            phasing = "in_cis_excluded"
        else:
            phasing = meta["origin"] if trio_mode else "unphased_singleton"

        r["consequence_group"] = classify_consequence_group(r["consequence"])
        r["inheritance_mode"] = primary
        r["phasing_status"] = phasing
        r["all_compatible_modes"] = ";".join(compat)
        kept_rows.append(r)

    cand_cols = [
        "variant_id", "chrom", "pos", "ref", "alt", "gene", "consequence", "consequence_group",
        "impact", "hgvsc", "hgvsp", "gnomad_af", "proband_gt", "proband_dp", "proband_gq",
        "proband_vaf", "mother_gt", "father_gt", "inheritance_mode", "phasing_status",
        "all_compatible_modes", "clinvar_sig", "revel", "cadd_phred", "spliceai_max", "pli", "loeuf",
    ]
    cand_df = pd.DataFrame(kept_rows, columns=cand_cols) if kept_rows else pd.DataFrame(columns=cand_cols)

    # Build gene-level inheritance summary
    gene_rows: list[dict[str, Any]] = []
    if not cand_df.empty:
        for gene, gdf in cand_df.groupby("gene", sort=True):
            if not str(gene).strip():
                continue
            modes = sorted(gdf["inheritance_mode"].astype(str).unique().tolist())
            phasings = sorted(gdf["phasing_status"].astype(str).unique().tolist())
            afs = gdf["gnomad_af"].dropna().to_numpy(dtype=float, copy=True)
            gene_rows.append(
                {
                    "gene": gene,
                    "n_rare_variants": int(len(gdf)),
                    "inheritance_modes": ";".join(modes),
                    "phasing_status": ";".join(phasings),
                    "best_consequence": str(gdf.iloc[0]["consequence"]),
                    "min_gnomad_af": float(np.min(afs)) if afs.size else 0.0,
                    "has_de_novo": bool("de_novo" in modes),
                    "has_ar_hom": bool("autosomal_recessive_hom" in modes),
                    "has_phased_comp_het": bool("compound_heterozygous" in modes),
                    "has_unphased_comp_het": bool("unphased_candidate_comp_het" in modes),
                }
            )
    gene_cols = [
        "gene", "n_rare_variants", "inheritance_modes", "phasing_status",
        "best_consequence", "min_gnomad_af", "has_de_novo", "has_ar_hom",
        "has_phased_comp_het", "has_unphased_comp_het",
    ]
    gene_df = pd.DataFrame(gene_rows, columns=gene_cols) if gene_rows else pd.DataFrame(columns=gene_cols)

    inh_counts = (
        {str(k): int(v) for k, v in cand_df["inheritance_mode"].value_counts().items()}
        if not cand_df.empty
        else {}
    )

    cand_tsv = out_dir / "filtered_germline_candidates.tsv"
    gene_tsv = out_dir / "inheritance_gene_summary.tsv"
    funnel_png = out_dir / "germline_inheritance_funnel.png"

    funnel_steps = [
        ("Total input", total_variants),
        ("Proband carrier", proband_carrier_n),
        ("QC pass (DP/GQ)", qc_pass_n),
        ("Coding / splice", coding_splice_n),
        ("Mendelian rare", int(len(cand_df))),
    ]
    _plot_funnel(funnel_steps, inh_counts, funnel_png)

    params = {
        "proband_id": proband_id,
        "mother_id": mother_id,
        "father_id": father_id,
        "proband_sex": args.proband_sex,
        "max_af_dominant": args.max_af_dominant,
        "max_af_recessive": args.max_af_recessive,
        "min_dp": args.min_dp,
        "min_gq": args.min_gq,
    }
    sci_ctx = {"trio_mode": trio_mode, "inheritance_counts": inh_counts}
    write_tsv(
        cand_df,
        cand_tsv,
        tool="germ_vcf_filter",
        title="Filtered Mendelian Germline Candidate Variants",
        inputs=[str(vpath)],
        parameters=params,
        scientific_context=sci_ctx,
    )
    write_tsv(
        gene_df,
        gene_tsv,
        tool="germ_vcf_filter",
        title="Gene-Level Mendelian Inheritance Summary",
        inputs=[str(vpath)],
        parameters=params,
        scientific_context=sci_ctx,
    )
    write_artifact(
        funnel_png,
        tool="germ_vcf_filter",
        title="Germline Variant Filtering Funnel & Inheritance Breakdown",
        inputs=[str(vpath)],
        parameters=params,
        scientific_context=sci_ctx,
    )

    emit(
        {
            "status": "ok",
            "tool": "germ_vcf_filter",
            "input_format": input_format,
            "trio_mode": trio_mode,
            "proband_id": proband_id,
            "mother_id": mother_id,
            "father_id": father_id,
            "proband_sex": args.proband_sex,
            "total_variants_input": total_variants,
            "proband_carrier_variants": proband_carrier_n,
            "qc_pass_variants": qc_pass_n,
            "coding_or_splice_variants": coding_splice_n,
            "n_candidates": int(len(cand_df)),
            "n_candidate_genes": int(len(gene_df)),
            "inheritance_counts": inh_counts,
            "outputs": {
                "filtered_germline_candidates": str(cand_tsv),
                "inheritance_gene_summary": str(gene_tsv),
                "germline_inheritance_funnel": str(funnel_png),
            },
            "artifacts": [str(cand_tsv), str(gene_tsv), str(funnel_png)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
