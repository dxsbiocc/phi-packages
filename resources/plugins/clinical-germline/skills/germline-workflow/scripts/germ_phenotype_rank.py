#!/usr/bin/env python3
"""germ_phenotype_rank: Compute Information-Content-weighted HPO phenotypic
similarity between patient HPO terms/keywords and rare-disease genes, integrating
ACMG points and Mendelian inheritance fit into an Exomiser-inspired
combined_diagnostic_score (0..100).
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
    norm_gene,
    parse_str_list,
    read_table,
    resolve_column,
    safe_float,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)

# Canonical HPO ID dictionary with English/Chinese synonyms & clinical domain
HPO_DICT: dict[str, dict[str, Any]] = {
    "HP:0001250": {"label": "Seizure", "domain": "epilepsy", "syns": ["seizure", "seizures", "epilepsy", "epileptic", "convulsion", "癫痫", "抽搐", "惊厥"]},
    "HP:0002373": {"label": "Febrile seizure", "domain": "epilepsy", "syns": ["febrile seizure", "febrile seizures", "热性惊厥"]},
    "HP:0002133": {"label": "Status epilepticus", "domain": "epilepsy", "syns": ["status epilepticus", "癫痫持续状态"]},
    "HP:0011097": {"label": "Epileptic spasm", "domain": "epilepsy", "syns": ["epileptic spasm", "infantile spasms", "west syndrome", "婴儿痉挛"]},
    "HP:0001263": {"label": "Global developmental delay", "domain": "neurodev", "syns": ["developmental delay", "global developmental delay", "gdd", "发育迟缓", "全面发育迟缓"]},
    "HP:0001249": {"label": "Intellectual disability", "domain": "neurodev", "syns": ["intellectual disability", "mental retardation", "cognitive impairment", "智力障碍"]},
    "HP:0000729": {"label": "Autistic behavior", "domain": "neurodev", "syns": ["autism", "autistic", "asd", "自闭症", "孤独症"]},
    "HP:0001252": {"label": "Hypotonia", "domain": "neurodev", "syns": ["hypotonia", "floppy", "muscle weakness", "肌张力低下"]},
    "HP:0000252": {"label": "Microcephaly", "domain": "neurodev", "syns": ["microcephaly", "小头畸形"]},
    "HP:0000256": {"label": "Macrocephaly", "domain": "neurodev", "syns": ["macrocephaly", "大头畸形"]},
    "HP:0001251": {"label": "Ataxia", "domain": "neurodev", "syns": ["ataxia", "gait ataxia", "共济失调"]},
    "HP:0002376": {"label": "Developmental regression", "domain": "neurodev", "syns": ["regression", "developmental regression", "发育倒退"]},
    "HP:0001638": {"label": "Cardiomyopathy", "domain": "cardiac", "syns": ["cardiomyopathy", "心肌病"]},
    "HP:0001639": {"label": "Hypertrophic cardiomyopathy", "domain": "cardiac", "syns": ["hypertrophic cardiomyopathy", "hcm", "肥厚型心肌病"]},
    "HP:0001644": {"label": "Dilated cardiomyopathy", "domain": "cardiac", "syns": ["dilated cardiomyopathy", "dcm", "扩张型心肌病"]},
    "HP:0011675": {"label": "Arrhythmia", "domain": "cardiac", "syns": ["arrhythmia", "ventricular tachycardia", "brugada", "心律失常"]},
    "HP:0001657": {"label": "Prolonged QT interval", "domain": "cardiac", "syns": ["long qt", "prolonged qt", "lqts", "qt延长"]},
    "HP:0001645": {"label": "Sudden cardiac death", "domain": "cardiac", "syns": ["sudden cardiac death", "syncope", "心源性猝死", "晕厥"]},
    "HP:0002616": {"label": "Aortic root aneurysm", "domain": "connective", "syns": ["aortic aneurysm", "aortic root dilation", "aortic dissection", "主动脉瘤", "主动脉夹层"]},
    "HP:0001083": {"label": "Ectopia lentis", "domain": "connective", "syns": ["ectopia lentis", "lens dislocation", "晶状体脱位"]},
    "HP:0001166": {"label": "Arachnodactyly", "domain": "connective", "syns": ["arachnodactyly", "marfan", "marfanoid", "蜘蛛指", "马凡"]},
    "HP:0002659": {"label": "Increased susceptibility to fractures", "domain": "connective", "syns": ["fracture", "osteogenesis imperfecta", "brittle bone", "病理性骨折", "成骨不全"]},
    "HP:0000977": {"label": "Soft skin", "domain": "connective", "syns": ["hyperextensible skin", "ehlers-danlos", "joint hypermobility", "关节松弛"]},
    "HP:0003496": {"label": "Short stature", "domain": "skeletal", "syns": ["short stature", "dwarfism", "achondroplasia", "身材矮小", "软骨发育不全"]},
    "HP:0002240": {"label": "Hepatomegaly", "domain": "metabolic", "syns": ["hepatomegaly", "hepatosplenomegaly", "肝肿大", "肝脾肿大"]},
    "HP:0001987": {"label": "Hyperammonemia", "domain": "metabolic", "syns": ["hyperammonemia", "高氨血症"]},
    "HP:0004923": {"label": "Hyperphenylalaninemia", "domain": "metabolic", "syns": ["hyperphenylalaninemia", "phenylketonuria", "pku", "苯丙酮尿症"]},
    "HP:0001395": {"label": "Hepatic fibrosis", "domain": "metabolic", "syns": ["liver failure", "cholestasis", "cirrhosis", "kayser-fleischer", "肝衰竭", "胆汁淤积"]},
    "HP:0003198": {"label": "Myopathy", "domain": "metabolic", "syns": ["myopathy", "pompe", "elevated ck", "肌病"]},
    "HP:0002719": {"label": "Recurrent infections", "domain": "immune", "syns": ["recurrent infections", "immunodeficiency", "scid", "agammaglobulinemia", "反复感染", "免疫缺陷"]},
    "HP:0004432": {"label": "Agammaglobulinemia", "domain": "immune", "syns": ["agammaglobulinemia", "hypogammaglobulinemia", "低丙种球蛋白血症"]},
    "HP:0000964": {"label": "Eczema", "domain": "immune", "syns": ["eczema", "hyper-ige", "atopic dermatitis", "湿疹"]},
    "HP:0000107": {"label": "Renal cyst", "domain": "renal", "syns": ["renal cyst", "polycystic kidney", "肾囊肿", "多囊肾"]},
    "HP:0000093": {"label": "Proteinuria", "domain": "renal", "syns": ["proteinuria", "nephrotic", "hematuria", "蛋白尿", "肾病综合征", "血尿"]},
    "HP:0000365": {"label": "Hearing impairment", "domain": "renal", "syns": ["hearing loss", "sensorineural deafness", "alport", "听力损失", "耳聋"]},
    "HP:0000510": {"label": "Rod-cone dystrophy", "domain": "ciliopathy", "syns": ["retinitis pigmentosa", "retinal dystrophy", "polydactyly", "视网膜色素变性", "多指"]},
}

# Built-in 61-gene multi-system rare-disease HPO mapping
GENE_HPO_BASE: dict[str, list[str]] = {
    # Epilepsy / Neurodevelopmental (15)
    "SCN1A": ["HP:0001250", "HP:0002373", "HP:0002133", "HP:0001263", "HP:0001251"],
    "MECP2": ["HP:0001250", "HP:0001263", "HP:0002376", "HP:0000729", "HP:0000252", "HP:0001252"],
    "STXBP1": ["HP:0001250", "HP:0011097", "HP:0001263", "HP:0001249", "HP:0001252"],
    "CDKL5": ["HP:0001250", "HP:0011097", "HP:0001263", "HP:0001252", "HP:0000252"],
    "ARID1B": ["HP:0001263", "HP:0001249", "HP:0001252", "HP:0001250"],
    "CHD8": ["HP:0000729", "HP:0000256", "HP:0001263", "HP:0001249"],
    "KCNQ2": ["HP:0001250", "HP:0002133", "HP:0001263", "HP:0001252"],
    "SCN2A": ["HP:0001250", "HP:0000729", "HP:0001263", "HP:0001251"],
    "SCN8A": ["HP:0001250", "HP:0002133", "HP:0001263", "HP:0001252"],
    "SYNGAP1": ["HP:0001250", "HP:0001263", "HP:0001249", "HP:0000729"],
    "SHANK3": ["HP:0000729", "HP:0001263", "HP:0001252", "HP:0001249"],
    "TSC1": ["HP:0001250", "HP:0011097", "HP:0001249", "HP:0000729"],
    "TSC2": ["HP:0001250", "HP:0011097", "HP:0001249", "HP:0000107"],
    "GRIN2A": ["HP:0001250", "HP:0001263", "HP:0001249"],
    "FOXG1": ["HP:0001250", "HP:0000252", "HP:0001263", "HP:0001252"],
    # Cardiomyopathy / Arrhythmia (11)
    "MYH7": ["HP:0001638", "HP:0001639", "HP:0001644", "HP:0001645"],
    "MYBPC3": ["HP:0001638", "HP:0001639", "HP:0001645"],
    "TTN": ["HP:0001638", "HP:0001644", "HP:0011675", "HP:0003198"],
    "SCN5A": ["HP:0011675", "HP:0001657", "HP:0001645", "HP:0001644"],
    "KCNQ1": ["HP:0001657", "HP:0011675", "HP:0001645", "HP:0000365"],
    "LMNA": ["HP:0001638", "HP:0001644", "HP:0011675", "HP:0003198"],
    "KCNH2": ["HP:0001657", "HP:0011675", "HP:0001645"],
    "RYR2": ["HP:0011675", "HP:0001645"],
    "PKP2": ["HP:0001638", "HP:0011675", "HP:0001645"],
    "DSP": ["HP:0001638", "HP:0001644", "HP:0011675"],
    "TNNT2": ["HP:0001638", "HP:0001639", "HP:0001644", "HP:0001645"],
    # Connective / Skeletal (9)
    "FBN1": ["HP:0002616", "HP:0001083", "HP:0001166", "HP:0000977"],
    "COL1A1": ["HP:0002659", "HP:0000365", "HP:0003496"],
    "FGFR3": ["HP:0003496", "HP:0000256"],
    "COL1A2": ["HP:0002659", "HP:0000365", "HP:0000977"],
    "COL2A1": ["HP:0003496", "HP:0000365"],
    "COL3A1": ["HP:0002616", "HP:0000977"],
    "TGFBR1": ["HP:0002616", "HP:0001166"],
    "TGFBR2": ["HP:0002616", "HP:0001166"],
    "RUNX2": ["HP:0003496", "HP:0002659"],
    # Metabolic / Lysosomal (10)
    "PAH": ["HP:0004923", "HP:0001263", "HP:0001249", "HP:0001250", "HP:0000252"],
    "GAA": ["HP:0003198", "HP:0001252", "HP:0001638", "HP:0001639", "HP:0002240"],
    "GBA1": ["HP:0002240", "HP:0001250", "HP:0001251"],
    "HEXA": ["HP:0002376", "HP:0001250", "HP:0001252"],
    "ATP7B": ["HP:0001395", "HP:0002240", "HP:0001251"],
    "OTC": ["HP:0001987", "HP:0001250", "HP:0001263"],
    "BTD": ["HP:0001250", "HP:0001252", "HP:0000365", "HP:0000964"],
    "GLA": ["HP:0001638", "HP:0001639", "HP:0000093"],
    "IDUA": ["HP:0002240", "HP:0003496", "HP:0001249"],
    "SLC25A13": ["HP:0001987", "HP:0001395", "HP:0002240"],
    # Immunodeficiency (8)
    "BTK": ["HP:0002719", "HP:0004432"],
    "IL2RG": ["HP:0002719", "HP:0004432", "HP:0003496"],
    "STAT3": ["HP:0002719", "HP:0000964", "HP:0002659"],
    "RAG1": ["HP:0002719", "HP:0004432"],
    "RAG2": ["HP:0002719", "HP:0004432"],
    "ADA": ["HP:0002719", "HP:0004432"],
    "CYBB": ["HP:0002719"],
    "FOXP3": ["HP:0002719", "HP:0000964"],
    # Renal / Ciliopathy (8)
    "PKD1": ["HP:0000107", "HP:0002616"],
    "COL4A5": ["HP:0000093", "HP:0000365"],
    "BBS1": ["HP:0000510", "HP:0000107", "HP:0001249"],
    "PKD2": ["HP:0000107"],
    "PKHD1": ["HP:0000107", "HP:0001395"],
    "NPHS1": ["HP:0000093"],
    "NPHS2": ["HP:0000093"],
    "PAX2": ["HP:0000093", "HP:0000107"],
}


def _load_custom_hpo(path: str | None) -> tuple[dict[str, set[str]], dict[str, str]]:
    gene_map: dict[str, set[str]] = {g: set(terms) for g, terms in GENE_HPO_BASE.items()}
    extra_labels: dict[str, str] = {}
    if not path:
        return gene_map, extra_labels
    df = read_table(path)
    gcol = resolve_column(df, ["gene_symbol", "gene", "symbol", "ncbi_gene_symbol"], "gene_symbol")
    hcol = resolve_column(df, ["hpo_id", "hpo", "term_id"], "hpo_id")
    ncol = resolve_column(df, ["hpo_name", "term_name", "phenotype", "name"], "hpo_name", required=False)
    for idx in range(len(df)):
        row = df.iloc[idx]
        g = norm_gene(row[gcol])
        hid = str(row[hcol]).strip().upper()
        if not g or not hid.startswith("HP:"):
            continue
        gene_map.setdefault(g, set()).add(hid)
        if ncol and pd.notna(row[ncol]):
            extra_labels[hid] = str(row[ncol]).strip()
    return gene_map, extra_labels


def _resolve_query_terms(
    raw_terms: list[str],
    gene_map: dict[str, set[str]],
    extra_labels: dict[str, str],
) -> tuple[list[str], list[str], list[str]]:
    """Resolve HPO IDs and keyword terms into HPO IDs + fallback keyword tokens."""
    all_corpus_hpos = set().union(*gene_map.values()) if gene_map else set()
    resolved_ids: list[str] = []
    keyword_tokens: list[str] = []
    warnings: list[str] = []

    for term in raw_terms:
        t = term.strip()
        if not t:
            continue
        m = re.match(r"^(HP:\d{7})$", t.upper())
        if m:
            hid = m.group(1)
            if hid in all_corpus_hpos or hid in HPO_DICT:
                if hid not in resolved_ids:
                    resolved_ids.append(hid)
            else:
                warnings.append(f"HPO term '{hid}' not present in built-in or custom annotation map.")
                if hid not in resolved_ids:
                    resolved_ids.append(hid)
            continue

        # Match keyword against HPO_DICT synonyms & extra_labels
        tl = t.lower()
        matched_any = False
        for hid, meta in HPO_DICT.items():
            label_l = meta["label"].lower()
            if tl in label_l or label_l in tl or any(tl in s.lower() or s.lower() in tl for s in meta["syns"]):
                if hid not in resolved_ids:
                    resolved_ids.append(hid)
                matched_any = True
        for hid, lbl in extra_labels.items():
            if tl in lbl.lower():
                if hid not in resolved_ids:
                    resolved_ids.append(hid)
                matched_any = True
        if not matched_any:
            keyword_tokens.append(tl)
            warnings.append(f"Phenotype keyword '{t}' did not match a built-in HPO entry directly.")

    return resolved_ids, keyword_tokens, warnings


def _compute_ic_map(gene_map: dict[str, set[str]]) -> dict[str, float]:
    n_genes = max(len(gene_map), 1)
    freq: dict[str, int] = {}
    for terms in gene_map.values():
        for t in terms:
            freq[t] = freq.get(t, 0) + 1
    ic: dict[str, float] = {}
    for t, cnt in freq.items():
        ic[t] = -math.log((cnt + 0.5) / (n_genes + 1.0))
    return ic


def _gene_phenotype_sim(
    gene: str,
    query_ids: list[str],
    gene_map: dict[str, set[str]],
    ic_map: dict[str, float],
    extra_labels: dict[str, str],
) -> tuple[float, list[str]]:
    g_terms = gene_map.get(gene, set())
    if not g_terms or not query_ids:
        return 0.0, []

    default_ic = math.log(len(gene_map) + 1.0)
    q_domains = {HPO_DICT[q]["domain"] for q in query_ids if q in HPO_DICT}
    g_domains = {HPO_DICT[t]["domain"] for t in g_terms if t in HPO_DICT}

    matched_ic = 0.0
    query_ic = 0.0
    matched_labels: list[str] = []

    for q in query_ids:
        q_ic = ic_map.get(q, default_ic)
        query_ic += q_ic
        if q in g_terms:
            matched_ic += q_ic
            lbl = HPO_DICT.get(q, {}).get("label") or extra_labels.get(q) or q
            matched_labels.append(f"{q} ({lbl})")
        elif q in HPO_DICT and HPO_DICT[q]["domain"] in g_domains:
            # Partial semantic-domain credit for closely related phenotypes in the same organ system
            matched_ic += 0.45 * q_ic

    recall_ic = matched_ic / max(query_ic, 1e-9)
    union_terms = g_terms | set(query_ids)
    union_ic = sum(ic_map.get(t, default_ic) for t in union_terms)
    jaccard_ic = matched_ic / max(union_ic, 1e-9)

    sim = min(1.0, 0.72 * recall_ic + 0.28 * jaccard_ic)
    if q_domains and (q_domains & g_domains) and sim == 0.0:
        sim = 0.20
    return round(float(sim), 4), matched_labels


def _variant_pathogenicity_component(row: pd.Series) -> float:
    """Return 0..50 score from acmg_points if present, or fallback variant severity."""
    pts = safe_float(row.get("acmg_points"))
    if pts is not None:
        return float(np.clip((pts + 4.0) / 16.0, 0.0, 1.0) * 50.0)
    csq_grp = classify_consequence_group(row.get("consequence"))
    base = {"lof": 34.0, "inframe_indel": 24.0, "missense": 20.0, "splice_region": 18.0, "synonymous": 4.0}.get(csq_grp, 8.0)
    rev = safe_float(row.get("revel"))
    if rev is not None and rev >= 0.75:
        base += 8.0
    csig = str(row.get("clinvar_sig", "")).lower()
    if "pathogenic" in csig and "conflicting" not in csig and "benign" not in csig:
        base += 12.0
    return float(min(50.0, base))


def _inheritance_bonus(mode: Any) -> float:
    m = str(mode or "").strip().lower()
    if m in {"de_novo", "compound_heterozygous", "autosomal_recessive_hom"}:
        return 10.0
    if m in {"x_linked"}:
        return 8.0
    if m in {"unphased_candidate_comp_het", "autosomal_dominant"}:
        return 6.0
    return 2.0


def _plot_prioritization(gene_df: pd.DataFrame, out_png: Path) -> None:
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.6))

    if not gene_df.empty:
        x = gene_df["phenotype_sim"].to_numpy(dtype=float, copy=True)
        y = gene_df["max_acmg_points"].to_numpy(dtype=float, copy=True)
        scores = gene_df["combined_diagnostic_score"].to_numpy(dtype=float, copy=True)
        sizes = 40.0 + 3.0 * scores
        sc = axes[0].scatter(x, y, s=sizes, c=scores, cmap="viridis", alpha=0.85, edgecolors="#1e293b", linewidth=0.7)
        fig.colorbar(sc, ax=axes[0], label="Diagnostic score (0-100)")
        for _, r in gene_df.head(6).iterrows():
            axes[0].annotate(
                str(r["gene"]),
                (float(r["phenotype_sim"]), float(r["max_acmg_points"])),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8.5,
                fontweight="bold",
            )
        axes[0].set_xlim(-0.05, 1.08)
        axes[0].set_xlabel("HPO Phenotype Similarity (0..1)")
        axes[0].set_ylabel("Max ACMG Points / Variant Score")
        axes[0].set_title("Genotype vs. Phenotype Diagnostic Space")

        top = gene_df.head(10).iloc[::-1]
        y_pos = np.arange(len(top))
        vals = top["combined_diagnostic_score"].to_numpy(dtype=float, copy=True)
        axes[1].barh(y_pos, vals, color="#0f766e", edgecolor="#1e293b", linewidth=0.6)
        axes[1].set_yticks(y_pos)
        axes[1].set_yticklabels(top["gene"].tolist())
        axes[1].set_xlim(0, 110)
        axes[1].set_xlabel("Combined Diagnostic Score (0-100)")
        axes[1].set_title("Top Ranked Candidate Causal Genes")
        for yp, v in zip(y_pos, vals):
            axes[1].text(v + 1.5, yp, f"{v:.1f}", va="center", fontsize=8.5)
    else:
        axes[0].text(0.5, 0.5, "No genes to plot", ha="center", va="center")
        axes[1].text(0.5, 0.5, "No genes to plot", ha="center", va="center")

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--hpo_terms", nargs="+", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--hpo_annotations", default=None)
    args = ap.parse_args()

    vpath = Path(args.variants)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_terms = parse_str_list(args.hpo_terms)
    if not raw_terms:
        fail("At least one HPO term or phenotype keyword is required.")

    df = read_table(vpath)
    gene_col = resolve_column(df, ["gene", "symbol", "gene_name", "hgnc_symbol"], "gene")
    csq_col = resolve_column(df, ["consequence", "effect", "annotation"], "consequence", required=False)
    inh_col = resolve_column(df, ["inheritance_mode", "inheritance"], "inheritance_mode", required=False)
    acmg_pts_col = resolve_column(df, ["acmg_points"], "acmg_points", required=False)
    acmg_cls_col = resolve_column(df, ["acmg_class"], "acmg_class", required=False)
    vid_col = resolve_column(df, ["variant_id", "id", "rsid"], "variant_id", required=False)

    gene_map, extra_labels = _load_custom_hpo(args.hpo_annotations)
    resolved_ids, _, warnings = _resolve_query_terms(raw_terms, gene_map, extra_labels)
    if not resolved_ids:
        warnings.append("No HPO terms matched the HPO dictionary; phenotypic similarity will be 0.0.")
    ic_map = _compute_ic_map(gene_map)

    work = df.copy()
    work["gene"] = work[gene_col].map(norm_gene)
    if csq_col and "consequence" not in work.columns:
        work["consequence"] = work[csq_col].fillna("").astype(str)
    elif "consequence" not in work.columns:
        work["consequence"] = ""
    if inh_col and "inheritance_mode" not in work.columns:
        work["inheritance_mode"] = work[inh_col].fillna("").astype(str)
    elif "inheritance_mode" not in work.columns:
        work["inheritance_mode"] = ""
    if acmg_pts_col and "acmg_points" not in work.columns:
        work["acmg_points"] = work[acmg_pts_col]
    if acmg_cls_col and "acmg_class" not in work.columns:
        work["acmg_class"] = work[acmg_cls_col]
    if vid_col and "variant_id" not in work.columns:
        work["variant_id"] = work[vid_col].fillna("").astype(str)
    elif "variant_id" not in work.columns:
        work["variant_id"] = [f"var_{i+1}:{g}" for i, g in enumerate(work["gene"])]

    pheno_sims: list[float] = []
    matched_terms_col: list[str] = []
    diag_scores: list[float] = []
    var_pts_col: list[float] = []

    for idx in range(len(work)):
        r = work.iloc[idx]
        g = str(r["gene"])
        sim, m_labels = _gene_phenotype_sim(g, resolved_ids, gene_map, ic_map, extra_labels)
        v_score = _variant_pathogenicity_component(r)
        i_bonus = _inheritance_bonus(r.get("inheritance_mode"))
        comb = round(float(np.clip(v_score + 40.0 * sim + i_bonus, 0.0, 100.0)), 2)

        pheno_sims.append(sim)
        matched_terms_col.append("; ".join(m_labels) if m_labels else "None")
        diag_scores.append(comb)
        apts = safe_float(r.get("acmg_points"))
        var_pts_col.append(float(apts) if apts is not None else round(v_score / 4.0, 2))

    work["phenotype_sim"] = pheno_sims
    work["matched_hpo_terms"] = matched_terms_col
    work["combined_diagnostic_score"] = diag_scores
    work["_plot_var_points"] = var_pts_col

    work = work.sort_values(
        by=["combined_diagnostic_score", "phenotype_sim"], ascending=[False, False]
    ).reset_index(drop=True)

    # Aggregate gene ranking
    gene_rows: list[dict[str, Any]] = []
    for gene, gdf in work.groupby("gene", sort=False):
        if not str(gene).strip():
            continue
        best = gdf.iloc[0]
        # Bonus +3 for phased compound het / multiple rare hits in same gene (capped at 100)
        n_vars = int(len(gdf))
        multi_bonus = 3.0 if (n_vars >= 2 and str(best.get("inheritance_mode", "")) in {"compound_heterozygous", "unphased_candidate_comp_het"}) else 0.0
        g_score = round(min(100.0, float(best["combined_diagnostic_score"]) + multi_bonus), 2)
        gene_rows.append(
            {
                "gene": gene,
                "combined_diagnostic_score": g_score,
                "phenotype_sim": float(best["phenotype_sim"]),
                "max_acmg_points": float(gdf["_plot_var_points"].max()),
                "best_acmg_class": str(best.get("acmg_class", "")),
                "inheritance_modes": ";".join(sorted(set(gdf["inheritance_mode"].astype(str)) - {""})),
                "n_variants": n_vars,
                "top_variant_id": str(best["variant_id"]),
                "matched_hpo_terms": str(best["matched_hpo_terms"]),
            }
        )
    gene_df = pd.DataFrame(gene_rows)
    if not gene_df.empty:
        gene_df = gene_df.sort_values(
            by=["combined_diagnostic_score", "phenotype_sim", "max_acmg_points"],
            ascending=[False, False, False],
        ).reset_index(drop=True)
        gene_df.insert(0, "rank", np.arange(1, len(gene_df) + 1))

    work_out = work.drop(columns=["_plot_var_points"])

    gene_tsv = out_dir / "diagnostic_gene_ranking.tsv"
    var_tsv = out_dir / "phenotype_matched_variants.tsv"
    plot_png = out_dir / "diagnostic_prioritization_plot.png"

    _plot_prioritization(gene_df, plot_png)

    params = {"hpo_terms": raw_terms, "resolved_hpo_ids": resolved_ids, "hpo_annotations": args.hpo_annotations}
    sci_ctx = {"n_genes": int(len(gene_df)), "top_gene": str(gene_df.iloc[0]["gene"]) if not gene_df.empty else None}
    write_tsv(
        gene_df,
        gene_tsv,
        tool="germ_phenotype_rank",
        title="HPO Phenotype-Driven Diagnostic Gene Ranking",
        inputs=[str(vpath), args.hpo_annotations],
        parameters=params,
        scientific_context=sci_ctx,
    )
    write_tsv(
        work_out,
        var_tsv,
        tool="germ_phenotype_rank",
        title="Phenotype-Matched Germline Candidate Variants",
        inputs=[str(vpath), args.hpo_annotations],
        parameters=params,
        scientific_context=sci_ctx,
    )
    write_artifact(
        plot_png,
        tool="germ_phenotype_rank",
        title="Exomiser-Inspired Genotype-Phenotype Diagnostic Prioritization",
        inputs=[str(vpath), args.hpo_annotations],
        parameters=params,
        scientific_context=sci_ctx,
    )

    top_genes = gene_df.head(10).to_dict(orient="records") if not gene_df.empty else []

    emit(
        {
            "status": "ok",
            "tool": "germ_phenotype_rank",
            "query_hpo_terms": raw_terms,
            "resolved_hpo_ids": resolved_ids,
            "n_variants": int(len(work_out)),
            "n_genes": int(len(gene_df)),
            "top_genes": top_genes,
            "outputs": {
                "diagnostic_gene_ranking": str(gene_tsv),
                "phenotype_matched_variants": str(var_tsv),
                "diagnostic_prioritization_plot": str(plot_png),
            },
            "artifacts": [str(gene_tsv), str(var_tsv), str(plot_png)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
