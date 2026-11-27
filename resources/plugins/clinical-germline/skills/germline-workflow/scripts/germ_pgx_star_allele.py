#!/usr/bin/env python3
"""germ_pgx_star_allele: Call/resolve CPIC star-allele diplotypes, CPIC Activity
Scores, metabolizer/transporter phenotypes, and CPIC Level A clinical dosing
alerts across CYP2C19, CYP2D6, CYP2C9, VKORC1, SLCO1B1, DPYD, TPMT, NUDT15,
and UGT1A1 from diplotype tables or rsID/VCF genotypes.
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
    norm_gene,
    read_table,
    read_text_lines,
    resolve_column,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)

# Star-allele function & CPIC Activity Score tables (published CPIC specifications)
ALLELE_ACTIVITY: dict[str, dict[str, float]] = {
    "CYP2C19": {"*1": 1.0, "*2": 0.0, "*3": 0.0, "*4": 0.0, "*17": 1.5},
    "CYP2D6": {"*1": 1.0, "*2": 1.0, "*4": 0.0, "*5": 0.0, "*6": 0.0, "*10": 0.25, "*41": 0.5},
    "CYP2C9": {"*1": 1.0, "*2": 0.5, "*3": 0.0, "*5": 0.0, "*6": 0.0, "*8": 0.5, "*11": 0.5},
    "VKORC1": {"G": 1.0, "A": 0.0, "*1": 1.0, "*2": 0.0},
    "SLCO1B1": {"*1": 1.0, "*1A": 1.0, "*1B": 1.0, "*5": 0.0, "*15": 0.0, "*17": 0.0},
    "DPYD": {"*1": 1.0, "*2A": 0.0, "*13": 0.0, "c.2846A>T": 0.5, "HapB3": 0.5},
    "TPMT": {"*1": 1.0, "*2": 0.0, "*3A": 0.0, "*3B": 0.0, "*3C": 0.0},
    "NUDT15": {"*1": 1.0, "*2": 0.0, "*3": 0.0, "*9": 0.0},
    "UGT1A1": {"*1": 1.0, "*28": 0.5, "*80": 0.5, "*6": 0.5, "*37": 0.5},
}

# Defining rsIDs -> (gene, allele, ref_bases, alt_bases, priority)
RSID_DEFINITIONS: dict[str, dict[str, Any]] = {
    "rs4244285": {"gene": "CYP2C19", "allele": "*2", "ref": {"G"}, "alt": {"A"}, "prio": 10},
    "rs4986893": {"gene": "CYP2C19", "allele": "*3", "ref": {"G"}, "alt": {"A"}, "prio": 10},
    "rs12248560": {"gene": "CYP2C19", "allele": "*17", "ref": {"C"}, "alt": {"T"}, "prio": 5},
    "rs3892097": {"gene": "CYP2D6", "allele": "*4", "ref": {"C", "G"}, "alt": {"T", "A"}, "prio": 10},
    "rs5030655": {"gene": "CYP2D6", "allele": "*6", "ref": {"A", "T"}, "alt": {"-", "DEL"}, "prio": 10},
    "rs1065852": {"gene": "CYP2D6", "allele": "*10", "ref": {"G", "C"}, "alt": {"A", "T"}, "prio": 6},
    "rs28371725": {"gene": "CYP2D6", "allele": "*41", "ref": {"C", "G"}, "alt": {"T", "A"}, "prio": 7},
    "rs1799853": {"gene": "CYP2C9", "allele": "*2", "ref": {"C"}, "alt": {"T"}, "prio": 8},
    "rs1057910": {"gene": "CYP2C9", "allele": "*3", "ref": {"A"}, "alt": {"C"}, "prio": 10},
    "rs9923231": {"gene": "VKORC1", "allele": "A", "ref": {"G", "C"}, "alt": {"A", "T"}, "prio": 10},
    "rs4149056": {"gene": "SLCO1B1", "allele": "*5", "ref": {"T"}, "alt": {"C"}, "prio": 10},
    "rs3918290": {"gene": "DPYD", "allele": "*2A", "ref": {"C", "G"}, "alt": {"T", "A"}, "prio": 10},
    "rs55886062": {"gene": "DPYD", "allele": "*13", "ref": {"A", "T"}, "alt": {"C", "G"}, "prio": 10},
    "rs67376798": {"gene": "DPYD", "allele": "c.2846A>T", "ref": {"T", "A"}, "alt": {"A", "T"}, "prio": 8},
    "rs75017182": {"gene": "DPYD", "allele": "HapB3", "ref": {"G", "C"}, "alt": {"C", "G"}, "prio": 8},
    "rs1800462": {"gene": "TPMT", "allele": "*2", "ref": {"C", "G"}, "alt": {"G", "C"}, "prio": 9},
    "rs1800460": {"gene": "TPMT", "allele": "*3B", "ref": {"C", "G"}, "alt": {"T", "A"}, "prio": 9},
    "rs1142345": {"gene": "TPMT", "allele": "*3C", "ref": {"T", "A"}, "alt": {"C", "G"}, "prio": 9},
    "rs116855232": {"gene": "NUDT15", "allele": "*3", "ref": {"C"}, "alt": {"T"}, "prio": 10},
    "rs887829": {"gene": "UGT1A1", "allele": "*80", "ref": {"C"}, "alt": {"T"}, "prio": 9},
    "rs4148323": {"gene": "UGT1A1", "allele": "*6", "ref": {"G"}, "alt": {"A"}, "prio": 9},
}

DEFAULT_GENE_ORDER = [
    "CYP2C19", "CYP2D6", "CYP2C9", "VKORC1", "SLCO1B1", "DPYD", "TPMT", "NUDT15", "UGT1A1"
]


def _parse_diplotype_alleles(gene: str, diplotype: str) -> tuple[str, str]:
    s = str(diplotype).strip()
    if "/" in s:
        a1, a2 = [x.strip() for x in s.split("/", 1)]
    elif "|" in s:
        a1, a2 = [x.strip() for x in s.split("|", 1)]
    else:
        a1, a2 = s, "*1"
    if gene == "VKORC1":
        a1 = "A" if a1.upper() in {"A", "T", "*2"} else "G"
        a2 = "A" if a2.upper() in {"A", "T", "*2"} else "G"
    return a1, a2


def _diplotype_to_phenotype(gene: str, a1: str, a2: str) -> tuple[float, str, str]:
    tbl = ALLELE_ACTIVITY.get(gene, {})
    s1 = tbl.get(a1, tbl.get(a1.upper(), 1.0))
    s2 = tbl.get(a2, tbl.get(a2.upper(), 1.0))
    as_val = round(float(s1 + s2), 2)

    if gene == "CYP2C19":
        if s1 == 0.0 and s2 == 0.0:
            pheno = "Poor_Metabolizer"
        elif (s1 == 0.0 and s2 >= 1.0) or (s2 == 0.0 and s1 >= 1.0):
            pheno = "Intermediate_Metabolizer"
        elif s1 > 1.0 and s2 > 1.0:
            pheno = "Ultrarapid_Metabolizer"
        elif (s1 > 1.0 and s2 == 1.0) or (s2 > 1.0 and s1 == 1.0):
            pheno = "Rapid_Metabolizer"
        else:
            pheno = "Normal_Metabolizer"
    elif gene == "CYP2D6":
        if as_val == 0.0:
            pheno = "Poor_Metabolizer"
        elif as_val <= 1.0:
            pheno = "Intermediate_Metabolizer"
        elif as_val <= 2.25:
            pheno = "Normal_Metabolizer"
        else:
            pheno = "Ultrarapid_Metabolizer"
    elif gene in {"DPYD", "CYP2C9"}:
        if as_val <= 0.5:
            pheno = "Poor_Metabolizer"
        elif as_val <= 1.5:
            pheno = "Intermediate_Metabolizer"
        else:
            pheno = "Normal_Metabolizer"
    elif gene in {"TPMT", "NUDT15", "UGT1A1"}:
        if as_val == 0.0 or (gene == "UGT1A1" and as_val <= 1.0):
            pheno = "Poor_Metabolizer"
        elif as_val < 2.0:
            pheno = "Intermediate_Metabolizer"
        else:
            pheno = "Normal_Metabolizer"
    elif gene in {"SLCO1B1", "VKORC1"}:
        if as_val < 2.0:
            pheno = "Decreased_Function"
        else:
            pheno = "Normal_Metabolizer"
    else:
        pheno = "Normal_Metabolizer" if as_val >= 1.5 else "Intermediate_Metabolizer"

    dip_str = f"{a1}/{a2}"
    return as_val, pheno, dip_str


def _cpic_recommendations(gene: str, pheno: str, act_score: float, diplotype: str) -> list[dict[str, str]]:
    recs: list[dict[str, str]] = []

    def add(drug: str, tier: str, guidance: str) -> None:
        recs.append(
            {
                "gene": gene,
                "diplotype": diplotype,
                "phenotype": pheno,
                "drug": drug,
                "cpic_level": "A",
                "action_tier": tier,
                "recommendation": guidance,
            }
        )

    if gene == "CYP2C19":
        if pheno in {"Poor_Metabolizer", "Intermediate_Metabolizer"}:
            add("Clopidogrel", "alternative_therapy", "Avoid clopidogrel due to reduced bioactivation and MACE risk; use prasugrel or ticagrelor.")
        else:
            add("Clopidogrel", "standard_dosing", "Use standard clopidogrel dosing per clinical indication.")
        if pheno == "Poor_Metabolizer":
            add("Escitalopram", "dose_reduction_required", "Consider 50% starting dose reduction or alternative SSRI not metabolized by CYP2C19.")
            add("Voriconazole", "alternative_therapy", "Choose an alternative antifungal not dependent on CYP2C19 metabolism (e.g., isavuconazole).")
        elif pheno in {"Ultrarapid_Metabolizer", "Rapid_Metabolizer"}:
            add("Voriconazole", "alternative_therapy", "Avoid voriconazole due to subtherapeutic trough concentrations; select alternative agent.")
        else:
            add("Escitalopram", "standard_dosing", "Initiate therapy with standard recommended starting dose.")
    elif gene == "CYP2D6":
        if pheno in {"Poor_Metabolizer", "Ultrarapid_Metabolizer"}:
            msg = (
                "Avoid codeine/tramadol: risk of life-threatening respiratory depression from ultra-rapid morphine formation."
                if pheno == "Ultrarapid_Metabolizer"
                else "Avoid codeine/tramadol: lack of bioactivation to active morphine/O-desmethyltramadol."
            )
            add("Codeine", "contraindicated_or_avoid", msg)
        else:
            add("Codeine", "standard_dosing", "Use standard age/weight-based dosing; monitor analgesic response.")
        if pheno in {"Poor_Metabolizer", "Intermediate_Metabolizer"}:
            add("Tamoxifen", "alternative_therapy", "Consider alternative hormonal therapy (e.g., aromatase inhibitor) due to reduced endoxifen formation.")
        else:
            add("Tamoxifen", "standard_dosing", "Use standard tamoxifen 20 mg/day dosing.")
    elif gene == "CYP2C9":
        if pheno in {"Poor_Metabolizer", "Intermediate_Metabolizer"}:
            add("Warfarin", "dose_reduction_required", "Reduce initial warfarin dose using a CPIC/IWPC pharmacogenetic dosing algorithm.")
            if pheno == "Poor_Metabolizer":
                add("Phenytoin", "dose_reduction_required", "Reduce phenytoin maintenance dose by 50% and monitor serum levels.")
        else:
            add("Warfarin", "standard_dosing", "Standard CYP2C9 component for warfarin dosing (check VKORC1 status).")
    elif gene == "VKORC1":
        if pheno == "Decreased_Function":
            add("Warfarin", "dose_reduction_required", "Increased warfarin sensitivity (VKORC1 -1639A); reduce starting warfarin dose per CPIC algorithm.")
        else:
            add("Warfarin", "standard_dosing", "Normal VKORC1 warfarin sensitivity.")
    elif gene == "SLCO1B1":
        if pheno == "Decreased_Function":
            add("Simvastatin", "alternative_therapy", "Increased statin myopathy risk (rs4149056); prescribe lower dose (<=20 mg) or switch to rosuvastatin/pravastatin.")
        else:
            add("Simvastatin", "standard_dosing", "Prescribe standard simvastatin starting dose.")
    elif gene == "DPYD":
        if act_score <= 0.5:
            add("Fluorouracil (5-FU) / Capecitabine", "contraindicated_or_avoid", "Avoid fluoropyrimidines completely due to high risk of severe/fatal toxicity.")
        elif act_score <= 1.5:
            add("Fluorouracil (5-FU) / Capecitabine", "dose_reduction_required", "Reduce starting fluoropyrimidine dose by 50% and titrate based on toxicity/TDM.")
        else:
            add("Fluorouracil (5-FU) / Capecitabine", "standard_dosing", "Use standard fluoropyrimidine dosing.")
    elif gene in {"TPMT", "NUDT15"}:
        if pheno == "Poor_Metabolizer":
            add("Azathioprine / Mercaptopurine", "contraindicated_or_avoid", "Avoid thiopurines in non-malignant conditions or reduce dose by >=90% with intensive monitoring.")
        elif pheno == "Intermediate_Metabolizer":
            add("Azathioprine / Mercaptopurine", "dose_reduction_required", "Start thiopurine at 30-80% of standard dose to prevent severe myelosuppression.")
        else:
            add("Azathioprine / Mercaptopurine", "standard_dosing", "Use standard thiopurine starting dose.")
    elif gene == "UGT1A1":
        if pheno == "Poor_Metabolizer":
            add("Irinotecan", "dose_reduction_required", "Reduce initial irinotecan dose by >=30% due to heightened risk of severe neutropenia and diarrhea.")
        else:
            add("Irinotecan", "standard_dosing", "Use standard irinotecan dosing.")
    return recs


def _alt_allele_copies(gt_raw: Any, rsid: str) -> int | None:
    """Return 0, 1, or 2 alt copies for an rsID genotype string."""
    if gt_raw is None:
        return None
    s = str(gt_raw).strip()
    if not s or s in {".", "./.", ".|.", "NA", "nan", "None"}:
        return None
    if s in {"0/0", "0|0", "0"}:
        return 0
    if s in {"0/1", "1/0", "0|1", "1|0"}:
        return 1
    if s in {"1/1", "1|1", "1"}:
        return 2
    low = s.lower()
    if low in {"ref", "hom_ref", "wt"}:
        return 0
    if low in {"het", "heterozygous"}:
        return 1
    if low in {"hom_alt", "hom", "homozygous"}:
        return 2

    defn = RSID_DEFINITIONS.get(rsid, {})
    alt_set = {a.upper() for a in defn.get("alt", set())}
    ref_set = {r.upper() for r in defn.get("ref", set())}
    # Parse base pair like G/A, GA, A|A, G A
    parts = [p for p in re.split(r"[/|\s,]+", s.upper()) if p]
    if len(parts) == 1 and len(parts[0]) == 2 and parts[0][0] in "ACGT" and parts[0][1] in "ACGT":
        parts = [parts[0][0], parts[0][1]]
    if len(parts) == 2:
        alts = sum(1 for p in parts if p in alt_set)
        if alts > 0:
            return min(2, alts)
        if all(p in ref_set for p in parts):
            return 0
    return None


def _resolve_rsid_hits_to_diplotype(gene: str, hits: list[tuple[str, int, int]]) -> tuple[str, str]:
    """hits: list of (allele_name, copies_1_or_2, priority). Returns (allele1, allele2)."""
    default_ref = "G" if gene == "VKORC1" else "*1"
    if not hits:
        return default_ref, default_ref

    # Special TPMT *3A cis haplotype check: rs1800460 (*3B) + rs1142345 (*3C)
    if gene == "TPMT":
        alleles_map = {a: c for a, c, _ in hits}
        c3b = alleles_map.get("*3B", 0)
        c3c = alleles_map.get("*3C", 0)
        if c3b > 0 and c3c > 0:
            if c3b >= 2 and c3c >= 2:
                return "*3A", "*3A"
            return "*1", "*3A"

    # Sort descending by priority and copy count
    sorted_hits = sorted(hits, key=lambda x: (x[1], x[2]), reverse=True)
    called: list[str] = []
    for allele, copies, _ in sorted_hits:
        for _ in range(copies):
            if len(called) < 2:
                called.append(allele)
    while len(called) < 2:
        called.append(default_ref)
    # Canonical ordering: *1 first when het
    if called[1] == default_ref and called[0] != default_ref:
        called = [called[1], called[0]]
    return called[0], called[1]


def _parse_pgx_vcf(path: Path) -> tuple[list[dict[str, Any]], list[str]]:
    lines = read_text_lines(path)
    samples: list[str] = []
    rows: list[dict[str, Any]] = []
    warnings: list[str] = []

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
        rsid = cols[2].strip().lower()
        if rsid not in RSID_DEFINITIONS:
            continue
        fmt_keys = cols[8].split(":") if len(cols) > 8 else ["GT"]
        gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else 0
        for i, sid in enumerate(samples):
            raw_sample = cols[9 + i] if 9 + i < len(cols) else "0/0"
            gt = raw_sample.split(":")[gt_idx] if gt_idx < len(raw_sample.split(":")) else "./."
            copies = _alt_allele_copies(gt, rsid)
            if copies is not None:
                rows.append({"sample_id": sid, "rsid": rsid, "copies": copies})
    if not rows:
        warnings.append("No CPIC defining rsIDs found in VCF ID column; defaulting tested genes to *1/*1.")
        for sid in samples or ["SAMPLE_1"]:
            for g in DEFAULT_GENE_ORDER:
                rows.append({"sample_id": sid, "rsid": None, "gene": g, "copies": 0})
    return rows, warnings


def _plot_pgx_overview(rec_df: pd.DataFrame, pheno_df: pd.DataFrame, out_png: Path) -> None:
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5))

    tier_order = ["contraindicated_or_avoid", "dose_reduction_required", "alternative_therapy", "standard_dosing"]
    tier_colors = {
        "contraindicated_or_avoid": "#dc2626",
        "dose_reduction_required": "#ea580c",
        "alternative_therapy": "#d97706",
        "standard_dosing": "#0f766e",
    }
    if not rec_df.empty:
        counts = np.array([int((rec_df["action_tier"] == t).sum()) for t in tier_order], dtype=float)
        x = np.arange(len(tier_order))
        axes[0].bar(x, counts, color=[tier_colors[t] for t in tier_order], edgecolor="#1e293b", linewidth=0.6)
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(["Avoid /\nContraind.", "Dose\nReduction", "Alternative\nTherapy", "Standard\nDosing"])
        axes[0].set_ylabel("CPIC Level A Recommendations")
        axes[0].set_title("CPIC Clinical Actionability Tiers")
        max_c = max(float(np.max(counts)) if counts.size else 1.0, 1.0)
        axes[0].set_ylim(0, max_c * 1.2)
        for xi, c in zip(x, counts):
            axes[0].text(xi, c + max_c * 0.03, f"{int(c)}", ha="center", fontsize=8.5)
    else:
        axes[0].text(0.5, 0.5, "No drug recommendations", ha="center", va="center")

    if not pheno_df.empty:
        sub = pheno_df.head(12).iloc[::-1]
        y = np.arange(len(sub))
        vals = sub["activity_score"].to_numpy(dtype=float, copy=True)
        labels = [f"{r['sample_id']} | {r['gene']} ({r['diplotype']})" for _, r in sub.iterrows()]
        pcols = [
            "#dc2626" if "Poor" in str(p) or "Ultrarapid" in str(p)
            else ("#d97706" if "Intermediate" in str(p) or "Decreased" in str(p) else "#0f766e")
            for p in sub["phenotype"]
        ]
        axes[1].barh(y, vals, color=pcols, edgecolor="#1e293b", linewidth=0.6)
        axes[1].set_yticks(y)
        axes[1].set_yticklabels(labels, fontsize=8.5)
        axes[1].set_xlabel("CPIC Activity Score")
        axes[1].set_title("Pharmacogene Diplotypes & Activity Scores")
        axes[1].set_xlim(0, max(float(np.max(vals)) if vals.size else 2.0, 2.0) * 1.15)
        for yp, v in zip(y, vals):
            axes[1].text(v + 0.05, yp, f"{v:.2g}", va="center", fontsize=8.5)
    else:
        axes[1].text(0.5, 0.5, "No gene calls", ha="center", va="center")

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    gpath = Path(args.genotypes)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []
    gene_calls: list[dict[str, Any]] = []

    if is_vcf_file(gpath):
        input_mode = "vcf"
        raw_hits, vcf_warns = _parse_pgx_vcf(gpath)
        warnings.extend(vcf_warns)
        by_sample_gene: dict[tuple[str, str], list[tuple[str, int, int]]] = {}
        samples_seen: set[str] = set()
        genes_seen: set[str] = set()
        for item in raw_hits:
            sid = item["sample_id"]
            samples_seen.add(sid)
            if item.get("rsid"):
                defn = RSID_DEFINITIONS[item["rsid"]]
                g = defn["gene"]
                genes_seen.add(g)
                by_sample_gene.setdefault((sid, g), [])
                if item["copies"] > 0:
                    by_sample_gene[(sid, g)].append((defn["allele"], int(item["copies"]), int(defn["prio"])))
            elif item.get("gene"):
                g = item["gene"]
                genes_seen.add(g)
                by_sample_gene.setdefault((sid, g), [])
        for sid in sorted(samples_seen):
            for g in [x for x in DEFAULT_GENE_ORDER if x in genes_seen] + sorted(genes_seen - set(DEFAULT_GENE_ORDER)):
                a1, a2 = _resolve_rsid_hits_to_diplotype(g, by_sample_gene.get((sid, g), []))
                act_score, pheno, dip_str = _diplotype_to_phenotype(g, a1, a2)
                gene_calls.append(
                    {"sample_id": sid, "gene": g, "diplotype": dip_str, "activity_score": act_score, "phenotype": pheno}
                )
    else:
        df = read_table(gpath)
        sample_col = resolve_column(df, ["sample_id", "sample", "patient_id"], "sample_id", required=False)
        dip_col = resolve_column(df, ["diplotype", "star_alleles", "alleles"], "diplotype", required=False)
        rsid_col = resolve_column(df, ["rsid", "snp", "variant_id", "id"], "rsid", required=False)
        gt_col = resolve_column(df, ["genotype", "gt", "call"], "genotype", required=False)
        gene_col = resolve_column(df, ["gene", "pgx_gene", "symbol"], "gene", required=False)

        if dip_col and gene_col:
            input_mode = "diplotype_table"
            for idx in range(len(df)):
                row = df.iloc[idx]
                sid = str(row[sample_col]).strip() if sample_col and pd.notna(row[sample_col]) else "SAMPLE_1"
                g = norm_gene(row[gene_col])
                if not g or pd.isna(row[dip_col]):
                    continue
                a1, a2 = _parse_diplotype_alleles(g, str(row[dip_col]))
                act_score, pheno, dip_str = _diplotype_to_phenotype(g, a1, a2)
                gene_calls.append(
                    {"sample_id": sid, "gene": g, "diplotype": dip_str, "activity_score": act_score, "phenotype": pheno}
                )
        elif rsid_col and gt_col:
            input_mode = "rsid_table"
            by_sg: dict[tuple[str, str], list[tuple[str, int, int]]] = {}
            samples_seen = set()
            genes_seen = set()
            for idx in range(len(df)):
                row = df.iloc[idx]
                sid = str(row[sample_col]).strip() if sample_col and pd.notna(row[sample_col]) else "SAMPLE_1"
                samples_seen.add(sid)
                rsid = str(row[rsid_col]).strip().lower()
                if rsid not in RSID_DEFINITIONS:
                    continue
                defn = RSID_DEFINITIONS[rsid]
                g = defn["gene"]
                genes_seen.add(g)
                by_sg.setdefault((sid, g), [])
                copies = _alt_allele_copies(row[gt_col], rsid)
                if copies is not None and copies > 0:
                    by_sg[(sid, g)].append((defn["allele"], int(copies), int(defn["prio"])))
            if not genes_seen:
                fail("No recognized CPIC rsIDs found in genotype table.")
            for sid in sorted(samples_seen):
                for g in [x for x in DEFAULT_GENE_ORDER if x in genes_seen] + sorted(genes_seen - set(DEFAULT_GENE_ORDER)):
                    a1, a2 = _resolve_rsid_hits_to_diplotype(g, by_sg.get((sid, g), []))
                    act_score, pheno, dip_str = _diplotype_to_phenotype(g, a1, a2)
                    gene_calls.append(
                        {"sample_id": sid, "gene": g, "diplotype": dip_str, "activity_score": act_score, "phenotype": pheno}
                    )
        else:
            fail(f"Unrecognized PGx input columns {list(df.columns)}; provide (gene, diplotype) or (rsid, genotype).")

    if not gene_calls:
        fail("No valid pharmacogene calls produced from input.")

    pheno_df = pd.DataFrame(gene_calls)
    rec_rows: list[dict[str, Any]] = []
    for gc in gene_calls:
        recs = _cpic_recommendations(gc["gene"], gc["phenotype"], float(gc["activity_score"]), gc["diplotype"])
        for r in recs:
            rec_rows.append({"sample_id": gc["sample_id"], "activity_score": gc["activity_score"], **r})
    rec_df = pd.DataFrame(rec_rows)

    tier_counts = (
        {str(k): int(v) for k, v in rec_df["action_tier"].value_counts().items()}
        if not rec_df.empty
        else {}
    )
    n_alerts = int(
        rec_df["action_tier"].isin({"contraindicated_or_avoid", "dose_reduction_required", "alternative_therapy"}).sum()
    ) if not rec_df.empty else 0

    if tier_counts.get("contraindicated_or_avoid", 0) > 0:
        crit = rec_df[rec_df["action_tier"] == "contraindicated_or_avoid"]
        drugs = sorted((crit["gene"] + " -> " + crit["drug"]).unique().tolist())
        warnings.append(f"HIGH-PRIORITY CPIC contraindicated_or_avoid alert(s): {', '.join(drugs)}.")

    pheno_tsv = out_dir / "cpic_diplotype_phenotypes.tsv"
    rec_tsv = out_dir / "cpic_drug_recommendations.tsv"
    plot_png = out_dir / "cpic_actionability_overview.png"

    _plot_pgx_overview(rec_df, pheno_df, plot_png)

    sci_ctx = {"input_mode": input_mode, "action_tier_counts": tier_counts, "n_actionable_alerts": n_alerts}
    write_tsv(
        pheno_df,
        pheno_tsv,
        tool="germ_pgx_star_allele",
        title="CPIC Star-Allele Diplotypes, Activity Scores & Phenotypes",
        inputs=[str(gpath)],
        scientific_context=sci_ctx,
    )
    write_tsv(
        rec_df,
        rec_tsv,
        tool="germ_pgx_star_allele",
        title="CPIC Level A Drug Dosing Recommendations",
        inputs=[str(gpath)],
        scientific_context=sci_ctx,
    )
    write_artifact(
        plot_png,
        tool="germ_pgx_star_allele",
        title="CPIC Pharmacogenomic Actionability Overview",
        inputs=[str(gpath)],
        scientific_context=sci_ctx,
    )

    emit(
        {
            "status": "ok",
            "tool": "germ_pgx_star_allele",
            "input_mode": input_mode,
            "n_samples": int(pheno_df["sample_id"].nunique()),
            "n_gene_calls": int(len(pheno_df)),
            "n_recommendations": int(len(rec_df)),
            "n_actionable_alerts": n_alerts,
            "action_tier_counts": tier_counts,
            "gene_calls": gene_calls,
            "outputs": {
                "cpic_diplotype_phenotypes": str(pheno_tsv),
                "cpic_drug_recommendations": str(rec_tsv),
                "cpic_actionability_overview": str(plot_png),
            },
            "artifacts": [str(pheno_tsv), str(rec_tsv), str(plot_png)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
