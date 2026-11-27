#!/usr/bin/env python3
"""germ_acmg_prioritize: Classify germline variants under ACMG/AMP 2015 +
ClinGen Pejaver (2022) calibrated computational cutoffs + Tavtigian (2018/2020)
Bayesian point scale (P0=0.10, O_PVS=350) and flag ACMG SF v3.2 actionable
secondary findings. 100% original MIT code.
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
    classify_consequence_group,
    emit,
    fail,
    norm_gene,
    read_table,
    read_text_lines,
    resolve_column,
    safe_float,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)

# ACMG SF v3.2 medically actionable genes (Miller et al. 2023) mapped to clinical domain
ACMG_SF_V32_GENES: dict[str, str] = {
    # Hereditary Cancer
    "BRCA1": "Hereditary Breast/Ovarian Cancer",
    "BRCA2": "Hereditary Breast/Ovarian Cancer",
    "PALB2": "Hereditary Breast/Pancreatic Cancer",
    "TP53": "Li-Fraumeni Syndrome",
    "MLH1": "Lynch Syndrome",
    "MSH2": "Lynch Syndrome",
    "MSH6": "Lynch Syndrome",
    "PMS2": "Lynch Syndrome",
    "APC": "Familial Adenomatous Polyposis",
    "MUTYH": "MUTYH-Associated Polyposis",
    "BMPR1A": "Juvenile Polyposis Syndrome",
    "SMAD4": "Juvenile Polyposis / HHT",
    "STK11": "Peutz-Jeghers Syndrome",
    "PTEN": "PTEN Hamartoma Tumor Syndrome",
    "VHL": "Von Hippel-Lindau Syndrome",
    "RET": "MEN2 / Familial Medullary Thyroid Cancer",
    "MEN1": "Multiple Endocrine Neoplasia Type 1",
    "SDHD": "Hereditary Paraganglioma/Pheochromocytoma",
    "SDHAF2": "Hereditary Paraganglioma/Pheochromocytoma",
    "SDHC": "Hereditary Paraganglioma/Pheochromocytoma",
    "SDHB": "Hereditary Paraganglioma/Pheochromocytoma",
    "MAX": "Hereditary Paraganglioma/Pheochromocytoma",
    "TMEM127": "Hereditary Paraganglioma/Pheochromocytoma",
    "NF2": "Neurofibromatosis Type 2",
    "RB1": "Retinoblastoma",
    "WT1": "WT1-Related Wilms Tumor",
    "TSC1": "Tuberous Sclerosis Complex",
    "TSC2": "Tuberous Sclerosis Complex",
    # Aortopathy / Connective Tissue
    "FBN1": "Marfan Syndrome",
    "TGFBR1": "Loeys-Dietz Syndrome",
    "TGFBR2": "Loeys-Dietz Syndrome",
    "SMAD3": "Loeys-Dietz Syndrome",
    "ACTA2": "Familial Thoracic Aortic Aneurysm",
    "MYH11": "Familial Thoracic Aortic Aneurysm",
    "COL3A1": "Vascular Ehlers-Danlos Syndrome",
    # Cardiomyopathy & Arrhythmia
    "MYBPC3": "Hypertrophic Cardiomyopathy",
    "MYH7": "Hypertrophic/Dilated Cardiomyopathy",
    "TNNT2": "Hypertrophic/Dilated Cardiomyopathy",
    "TNNI3": "Hypertrophic Cardiomyopathy",
    "TPM1": "Hypertrophic Cardiomyopathy",
    "MYL3": "Hypertrophic Cardiomyopathy",
    "ACTC1": "Hypertrophic Cardiomyopathy",
    "PRKAG2": "PRKAG2 Cardiac Syndrome",
    "GLA": "Fabry Disease",
    "LAMP2": "Danon Disease",
    "LMNA": "Dilated Cardiomyopathy / Laminopathy",
    "FLNC": "Cardiomyopathy",
    "TTN": "Dilated Cardiomyopathy",
    "RBM20": "Dilated Cardiomyopathy",
    "BAG3": "Dilated Cardiomyopathy",
    "DES": "Desminopathy / Cardiomyopathy",
    "PKP2": "Arrhythmogenic Right Ventricular Cardiomyopathy",
    "DSP": "Arrhythmogenic Cardiomyopathy",
    "DSC2": "Arrhythmogenic Right Ventricular Cardiomyopathy",
    "DSG2": "Arrhythmogenic Right Ventricular Cardiomyopathy",
    "TMEM43": "Arrhythmogenic Right Ventricular Cardiomyopathy",
    "KCNQ1": "Long QT Syndrome Type 1",
    "KCNH2": "Long QT Syndrome Type 2",
    "SCN5A": "Long QT Type 3 / Brugada Syndrome",
    "CALM1": "Calmodulinopathy / Long QT",
    "CALM2": "Calmodulinopathy / Long QT",
    "CALM3": "Calmodulinopathy / Long QT",
    "TRDN": "CPVT",
    "RYR2": "Catecholaminergic Polymorphic VT",
    "CASQ2": "CPVT",
    "TTR": "Hereditary Transthyretin Amyloidosis",
    # Familial Hypercholesterolemia & Miscellaneous Actionable
    "LDLR": "Familial Hypercholesterolemia",
    "APOB": "Familial Hypercholesterolemia",
    "PCSK9": "Familial Hypercholesterolemia",
    "RYR1": "Malignant Hyperthermia Susceptibility",
    "CACNA1S": "Malignant Hyperthermia Susceptibility",
    "ATP7B": "Wilson Disease",
    "HFE": "Hereditary Hemochromatosis (p.Cys282Tyr)",
    "OTC": "Ornithine Transcarbamylase Deficiency",
    "GAA": "Pompe Disease",
    "BTD": "Biotinidase Deficiency",
    "HNF1A": "MODY3 Maturity-Onset Diabetes",
    "RPE65": "RPE65-Related Retinal Dystrophy",
    "ENG": "Hereditary Hemorrhagic Telangiectasia",
    "ACVRL1": "Hereditary Hemorrhagic Telangiectasia",
    "FLCN": "Birt-Hogg-Dube Syndrome",
}

# Curated haploinsufficient / LoF-disease mechanism rare-disease genes (90+ genes)
CURATED_LOF_HI_GENES: set[str] = set(ACMG_SF_V32_GENES.keys()) | {
    "SCN1A", "MECP2", "STXBP1", "CDKL5", "ARID1B", "CHD8", "SYNGAP1", "SHANK3",
    "FOXG1", "KCNQ2", "GRIN2A", "GRIN2B", "ANKRD11", "KMT2D", "NSD1", "CHD7",
    "CREBBP", "EP300", "COL1A1", "COL2A1", "PKD1", "PKD2", "PAX6", "NF1",
    "DMD", "UBE3A", "NRXN1", "PAH", "GBA1", "HEXA", "CFTR", "BTK", "IL2RG",
    "STAT3", "RAG1", "RAG2", "COL4A5", "BBS1", "SMN1", "DYRK1A", "FOXP1",
    "MED13L", "CTNNB1", "ASXL1", "KDM5C", "ATRX", "SLC2A1", "ALDH7A1",
}

TIER_ORDER = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]


def _load_custom_lof_genes(path: str | None) -> set[str]:
    if not path:
        return set()
    p = Path(path)
    if not p.is_file():
        fail(f"lof_genes file not found: {path}")
    lines = read_text_lines(p)
    non_empty = [ln.strip() for ln in lines if ln.strip() and not ln.strip().startswith("#")]
    if not non_empty:
        return set()
    if "\t" in non_empty[0] or "," in non_empty[0]:
        df = read_table(p)
        gcol = resolve_column(df, ["gene", "symbol", "gene_symbol", "hgnc_symbol"], "gene")
        return {norm_gene(g) for g in df[gcol].tolist() if norm_gene(g)}
    genes: set[str] = set()
    for ln in non_empty:
        tok = norm_gene(ln)
        if tok and tok not in {"GENE", "SYMBOL", "GENE_SYMBOL"}:
            genes.add(tok)
    return genes


def _tavtigian_posterior(points: int) -> float:
    """Tavtigian et al. (2018/2020) Bayesian posterior probability (P0=0.10, O_PVS=350)."""
    p0 = 0.10
    prior_odds = p0 / (1.0 - p0)
    o_pvs = 350.0
    combined_odds = prior_odds * (o_pvs ** (points / 8.0))
    return float(combined_odds / (1.0 + combined_odds))


def _evaluate_variant(row: pd.Series, lof_gene_set: set[str]) -> dict[str, Any]:
    gene = str(row["gene"])
    csq_group = str(row["consequence_group"])
    af = safe_float(row["gnomad_af"])
    revel = safe_float(row["revel"])
    cadd = safe_float(row["cadd_phred"])
    spliceai = safe_float(row["spliceai_max"])
    pli = safe_float(row["pli"])
    loeuf = safe_float(row["loeuf"])
    clinvar = str(row["clinvar_sig"]).strip().lower()
    inh = str(row["inheritance_mode"]).strip().lower()
    phasing = str(row["phasing_status"]).strip().lower()

    evidence: list[tuple[str, int]] = []
    is_lof_intolerant = (
        gene in lof_gene_set
        or (pli is not None and pli >= 0.90)
        or (loeuf is not None and loeuf <= 0.35)
    )

    # 1. PVS1 (+8): Null variant in LoF-intolerant gene
    if csq_group == "lof" and is_lof_intolerant:
        evidence.append(("PVS1", 8))

    # 2. PS1/PS5_ClinVar (+4): Established ClinVar Pathogenic without conflicts
    clin_tokens = [t.strip() for t in re.split(r"[,;/|\s]+", clinvar) if t.strip()]
    if (
        "pathogenic" in clin_tokens
        and "conflicting" not in clinvar
        and "benign" not in clinvar
        and "uncertain" not in clinvar
    ):
        evidence.append(("PS5_ClinVar", 4))

    # 3. PS2 (+4): Trio-confirmed de novo
    if inh == "de_novo" and phasing != "unphased_singleton":
        evidence.append(("PS2", 4))

    # 4. BA1 (-8) / BS1 (-4) / PM2_Supporting (+1) population AF rules
    ba1_hit = False
    if af is not None and af > 0.05:
        evidence.append(("BA1", -8))
        ba1_hit = True
    elif af is not None and af > 0.01:
        evidence.append(("BS1", -4))
    else:
        ar_mode = inh in {"autosomal_recessive_hom", "compound_heterozygous", "unphased_candidate_comp_het"}
        pm2_thresh = 1e-4 if ar_mode else 1e-5
        if af is None or af < pm2_thresh:
            evidence.append(("PM2_Supporting", 1))

    # 5. PM3 (+2): Phased in-trans compound heterozygous
    if inh == "compound_heterozygous" and phasing == "phased_in_trans":
        evidence.append(("PM3", 2))

    # 6. PM4 (+2): In-frame indel in non-repeat coding region
    if csq_group == "inframe_indel":
        evidence.append(("PM4", 2))

    # 7. Computational deleteriousness (Pejaver 2022 calibrated PP3 vs BP4; max non-redundant tier)
    pp3_pts = 0
    pp3_label = ""
    if (revel is not None and revel >= 0.932) or (spliceai is not None and spliceai >= 0.80) or (revel is None and cadd is not None and cadd >= 29.8):
        pp3_pts, pp3_label = 4, "PP3_Strong"
    elif (revel is not None and revel >= 0.773) or (spliceai is not None and spliceai >= 0.50) or (revel is None and cadd is not None and cadd >= 28.1):
        pp3_pts, pp3_label = 2, "PP3_Moderate"
    elif (revel is not None and revel >= 0.644) or (spliceai is not None and spliceai >= 0.20) or (cadd is not None and cadd >= 25.3):
        pp3_pts, pp3_label = 1, "PP3"

    if pp3_pts > 0 and csq_group != "lof":
        evidence.append((pp3_label, pp3_pts))
    elif pp3_pts > 0 and csq_group == "lof" and not is_lof_intolerant:
        evidence.append((pp3_label, pp3_pts))
    elif csq_group in {"missense", "synonymous", "other"} and (spliceai is None or spliceai < 0.10):
        if revel is not None and revel <= 0.183:
            evidence.append(("BP4_Moderate", -2))
        elif (revel is not None and revel <= 0.290) or (revel is None and cadd is not None and cadd <= 15.0):
            evidence.append(("BP4", -1))

    # 8. BP7 (-1): Synonymous with low splice score
    if csq_group == "synonymous" and (spliceai is None or spliceai < 0.10):
        evidence.append(("BP7", -1))

    points = int(sum(pts for _, pts in evidence))
    posterior = _tavtigian_posterior(points)

    if ba1_hit or points <= -7:
        acmg_class = "Benign"
    elif points >= 10:
        acmg_class = "Pathogenic"
    elif 6 <= points <= 9:
        acmg_class = "Likely_pathogenic"
    elif 0 <= points <= 5:
        acmg_class = "VUS"
    else:
        acmg_class = "Likely_benign"

    sf_condition = ACMG_SF_V32_GENES.get(gene, "")
    is_sf = bool(sf_condition)
    sf_reportable = bool(is_sf and acmg_class in {"Pathogenic", "Likely_pathogenic"})
    ev_codes = ";".join(f"{code}({pts:+d})" for code, pts in evidence)

    return {
        "acmg_points": points,
        "tavtigian_posterior": round(posterior, 5),
        "acmg_class": acmg_class,
        "acmg_evidence": ev_codes if ev_codes else "None",
        "is_lof_intolerant_gene": bool(is_lof_intolerant),
        "is_acmg_sf_v32": is_sf,
        "acmg_sf_condition": sf_condition,
        "acmg_sf_reportable": sf_reportable,
    }


def _plot_acmg_summary(df: pd.DataFrame, tier_counts: dict[str, int], out_png: Path) -> None:
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))

    tier_colors = {
        "Pathogenic": "#dc2626",
        "Likely_pathogenic": "#ea580c",
        "VUS": "#d97706",
        "Likely_benign": "#0284c7",
        "Benign": "#0f766e",
    }
    counts = np.array([tier_counts.get(t, 0) for t in TIER_ORDER], dtype=float)
    x = np.arange(len(TIER_ORDER))
    axes[0].bar(x, counts, color=[tier_colors[t] for t in TIER_ORDER], edgecolor="#1e293b", linewidth=0.6)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(TIER_ORDER, rotation=22, ha="right")
    axes[0].set_ylabel("Variant count")
    axes[0].set_title("ACMG/AMP 5-Tier Variant Classification")
    max_c = max(float(np.max(counts)) if counts.size else 1.0, 1.0)
    axes[0].set_ylim(0, max_c * 1.2)
    for xi, c in zip(x, counts):
        axes[0].text(xi, c + max_c * 0.03, f"{int(c)}", ha="center", fontsize=8.5)

    top = df.head(12).iloc[::-1]
    if not top.empty:
        labels = [
            f"{r['gene']} ({r['acmg_class'][:3]})" + (" [SF]" if r["acmg_sf_reportable"] else "")
            for _, r in top.iterrows()
        ]
        pts = top["acmg_points"].to_numpy(dtype=float, copy=True)
        cols = [tier_colors.get(str(c), "#64748b") for c in top["acmg_class"]]
        y = np.arange(len(labels))
        axes[1].barh(y, pts, color=cols, edgecolor="#1e293b", linewidth=0.6)
        axes[1].axvline(10, color="#dc2626", linestyle="--", linewidth=1.0, label="Pathogenic (>=10)")
        axes[1].axvline(6, color="#ea580c", linestyle=":", linewidth=1.0, label="Likely Pathogenic (>=6)")
        axes[1].set_yticks(y)
        axes[1].set_yticklabels(labels)
        axes[1].set_xlabel("Tavtigian Bayesian Points")
        axes[1].set_title("Top Prioritized Variants by ACMG Points")
        axes[1].legend(loc="lower right", frameon=False)
    else:
        axes[1].text(0.5, 0.5, "No variants to display", ha="center", va="center")
        axes[1].set_title("Top Prioritized Variants")

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--lof_genes", default=None)
    args = ap.parse_args()

    vpath = Path(args.variants)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = read_table(vpath)
    gene_col = resolve_column(df, ["gene", "symbol", "gene_name", "hgnc_symbol"], "gene")
    csq_col = resolve_column(df, ["consequence", "effect", "annotation", "variant_class"], "consequence")
    af_col = resolve_column(df, ["gnomad_af", "af", "popmax_af", "max_af"], "gnomad_af", required=False)
    vid_col = resolve_column(df, ["variant_id", "id", "rsid"], "variant_id", required=False)
    clin_col = resolve_column(df, ["clinvar_sig", "clnsig", "clin_sig", "clinvar"], "clinvar_sig", required=False)
    rev_col = resolve_column(df, ["revel", "revel_score"], "revel", required=False)
    cadd_col = resolve_column(df, ["cadd_phred", "cadd", "cadd_score"], "cadd_phred", required=False)
    sp_col = resolve_column(df, ["spliceai_max", "spliceai", "spliceai_score"], "spliceai_max", required=False)
    pli_col = resolve_column(df, ["pli"], "pli", required=False)
    loeuf_col = resolve_column(df, ["loeuf"], "loeuf", required=False)
    inh_col = resolve_column(df, ["inheritance_mode", "inheritance"], "inheritance_mode", required=False)
    ph_col = resolve_column(df, ["phasing_status", "phasing"], "phasing_status", required=False)
    chrom_col = resolve_column(df, ["chrom", "chr", "#chrom"], "chrom", required=False)
    pos_col = resolve_column(df, ["pos", "position", "start"], "pos", required=False)

    lof_set = set(CURATED_LOF_HI_GENES) | _load_custom_lof_genes(args.lof_genes)

    work = df.copy()
    work["gene"] = work[gene_col].map(norm_gene)
    work["consequence"] = work[csq_col].fillna("").astype(str)
    work["consequence_group"] = work["consequence"].map(classify_consequence_group)
    work["gnomad_af"] = work[af_col].map(lambda x: safe_float(x)) if af_col else np.nan
    work["clinvar_sig"] = work[clin_col].fillna("").astype(str) if clin_col else ""
    work["revel"] = work[rev_col].map(lambda x: safe_float(x)) if rev_col else np.nan
    work["cadd_phred"] = work[cadd_col].map(lambda x: safe_float(x)) if cadd_col else np.nan
    work["spliceai_max"] = work[sp_col].map(lambda x: safe_float(x)) if sp_col else np.nan
    work["pli"] = work[pli_col].map(lambda x: safe_float(x)) if pli_col else np.nan
    work["loeuf"] = work[loeuf_col].map(lambda x: safe_float(x)) if loeuf_col else np.nan
    work["inheritance_mode"] = work[inh_col].fillna("").astype(str) if inh_col else ""
    work["phasing_status"] = work[ph_col].fillna("").astype(str) if ph_col else ""
    if vid_col:
        work["variant_id"] = work[vid_col].fillna("").astype(str)
    elif chrom_col and pos_col:
        work["variant_id"] = work[chrom_col].astype(str) + ":" + work[pos_col].astype(str) + ":" + work["gene"]
    else:
        work["variant_id"] = [f"var_{i+1}:{g}" for i, g in enumerate(work["gene"])]

    eval_rows = [ _evaluate_variant(work.iloc[i], lof_set) for i in range(len(work)) ]
    eval_df = pd.DataFrame(eval_rows)
    for col in eval_df.columns:
        work[col] = eval_df[col]

    work = work.sort_values(by=["acmg_points", "tavtigian_posterior"], ascending=[False, False]).reset_index(drop=True)

    tier_counts = {t: int((work["acmg_class"] == t).sum()) for t in TIER_ORDER}
    n_plp = int(tier_counts["Pathogenic"] + tier_counts["Likely_pathogenic"])
    n_sf_hit = int(work["is_acmg_sf_v32"].sum())
    n_sf_rep = int(work["acmg_sf_reportable"].sum())

    warnings: list[str] = []
    if n_sf_rep > 0:
        sf_genes = sorted(work.loc[work["acmg_sf_reportable"], "gene"].unique().tolist())
        warnings.append(
            f"ACMG SF v3.2 actionable secondary finding(s) detected in: {', '.join(sf_genes)}."
        )

    out_tsv = out_dir / "acmg_classified_variants.tsv"
    out_png = out_dir / "acmg_evidence_summary.png"
    _plot_acmg_summary(work, tier_counts, out_png)

    params = {"lof_genes_file": args.lof_genes}
    sci_ctx = {
        "tier_counts": tier_counts,
        "n_acmg_sf_reportable": n_sf_rep,
        "bayesian_prior_p0": 0.10,
        "bayesian_o_pvs": 350.0,
    }
    write_tsv(
        work,
        out_tsv,
        tool="germ_acmg_prioritize",
        title="ACMG/AMP 2015 + Tavtigian Bayesian Classified Variants",
        inputs=[str(vpath), args.lof_genes],
        parameters=params,
        scientific_context=sci_ctx,
    )
    write_artifact(
        out_png,
        tool="germ_acmg_prioritize",
        title="ACMG 5-Tier Classification & Evidence Point Summary",
        inputs=[str(vpath), args.lof_genes],
        parameters=params,
        scientific_context=sci_ctx,
    )

    top_variants = []
    for _, r in work.head(10).iterrows():
        top_variants.append(
            {
                "variant_id": str(r["variant_id"]),
                "gene": str(r["gene"]),
                "consequence": str(r["consequence"]),
                "acmg_class": str(r["acmg_class"]),
                "acmg_points": int(r["acmg_points"]),
                "tavtigian_posterior": float(r["tavtigian_posterior"]),
                "acmg_evidence": str(r["acmg_evidence"]),
                "is_acmg_sf_v32": bool(r["is_acmg_sf_v32"]),
                "acmg_sf_reportable": bool(r["acmg_sf_reportable"]),
            }
        )

    emit(
        {
            "status": "ok",
            "tool": "germ_acmg_prioritize",
            "n_variants": int(len(work)),
            "tier_counts": tier_counts,
            "n_pathogenic_or_likely": n_plp,
            "n_acmg_sf_v32_genes_hit": n_sf_hit,
            "n_acmg_sf_reportable": n_sf_rep,
            "top_variants": top_variants,
            "outputs": {
                "acmg_classified_variants": str(out_tsv),
                "acmg_evidence_summary": str(out_png),
            },
            "artifacts": [str(out_tsv), str(out_png)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
