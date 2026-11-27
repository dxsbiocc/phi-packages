#!/usr/bin/env python3
"""imm_neoantigen_prioritize: Prioritize candidate neoantigens from upstream
MHC binding predictions or somatic peptide tables (with built-in anchor-motif
PSSM IC50 estimation), computing Differential Agretopicity Index (DAI), anchor
vs TCR-contact mutation status, foreignness/hydropathy, TPM/VAF/CCF gates, and
a composite immunogenicity score (0..100).
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
    STANDARD_AAS,
    emit,
    fail,
    grantham_distance,
    kyte_doolittle_gravy,
    parse_str_list,
    read_table,
    resolve_column,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)

# Curated MHC-I allele primary anchor preferences (P2 and C-terminus) + secondary
# groove weights for built-in PSSM IC50 approximation when NetMHCpan/MHCflurry
# scores are not pre-computed.
HLA_MOTIFS: dict[str, dict[str, dict[str, float]]] = {
    "HLA-A*01:01": {
        "p2": {"T": 2.4, "S": 2.1, "L": 0.6, "I": 0.5, "V": 0.4, "D": -0.8, "E": -0.8, "R": -1.5, "K": -1.5},
        "pc": {"Y": 2.8, "F": 1.4, "W": 0.9, "L": 0.2, "R": -1.6, "K": -1.6, "D": -1.2, "E": -1.2},
        "p3": {"D": 1.3, "E": 1.1, "L": 0.5},
    },
    "HLA-A*02:01": {
        "p2": {"L": 2.6, "M": 2.4, "I": 1.6, "V": 1.2, "Q": 1.1, "A": 0.8, "T": 0.6, "D": -1.8, "E": -1.8, "R": -1.6, "K": -1.6, "P": -1.4},
        "pc": {"V": 2.6, "L": 2.4, "I": 1.9, "A": 1.2, "M": 1.1, "T": 0.4, "R": -1.8, "K": -1.8, "D": -1.7, "E": -1.7, "P": -1.6},
        "p3": {"Y": 0.6, "F": 0.6, "W": 0.5, "A": 0.3},
    },
    "HLA-A*03:01": {
        "p2": {"L": 2.0, "V": 1.8, "I": 1.7, "M": 1.6, "T": 1.2, "A": 0.8, "D": -1.5, "E": -1.5},
        "pc": {"K": 2.7, "Y": 2.0, "R": 1.8, "F": 1.1, "D": -1.8, "E": -1.8},
        "p3": {"F": 0.6, "Y": 0.5, "L": 0.4},
    },
    "HLA-A*11:01": {
        "p2": {"V": 2.1, "T": 2.0, "I": 1.7, "L": 1.4, "S": 1.1, "M": 1.0, "D": -1.5, "E": -1.5},
        "pc": {"K": 2.8, "R": 2.1, "Y": 0.9, "D": -1.8, "E": -1.8},
        "p3": {"F": 0.6, "Y": 0.6},
    },
    "HLA-A*24:02": {
        "p2": {"Y": 2.7, "F": 2.2, "W": 1.5, "M": 0.8, "D": -1.6, "E": -1.6, "K": -1.4},
        "pc": {"F": 2.4, "L": 2.2, "I": 1.8, "W": 1.5, "R": -1.6, "K": -1.6},
        "p3": {"I": 0.5, "L": 0.5, "V": 0.4},
    },
    "HLA-B*07:02": {
        "p2": {"P": 2.8, "A": 1.0, "S": 0.4, "L": -0.4, "D": -1.6, "E": -1.6, "R": -1.4},
        "pc": {"L": 2.5, "F": 1.8, "M": 1.6, "I": 1.3, "V": 1.1, "A": 0.7, "D": -1.7, "E": -1.7},
        "p3": {"R": 0.9, "K": 0.7, "L": 0.4},
    },
    "HLA-B*08:01": {
        "p2": {"L": 1.6, "I": 1.4, "V": 1.2, "A": 0.8, "P": 0.6},
        "pc": {"L": 2.4, "I": 1.8, "V": 1.5, "F": 1.2, "D": -1.5, "E": -1.5},
        "p3": {"K": 1.8, "R": 1.8},
    },
    "HLA-B*15:01": {
        "p2": {"Q": 2.3, "L": 1.9, "M": 1.6, "I": 1.1, "V": 1.0, "D": -1.5},
        "pc": {"Y": 2.5, "F": 2.4, "W": 1.6, "L": 1.1},
        "p3": {"K": 0.5, "V": 0.4},
    },
    "HLA-B*27:05": {
        "p2": {"R": 3.0, "K": 1.5, "Q": 0.8, "D": -2.0, "E": -2.0},
        "pc": {"L": 2.1, "F": 1.8, "K": 1.5, "R": 1.5, "Y": 1.4, "I": 1.2},
        "p3": {"L": 0.6, "F": 0.5},
    },
    "HLA-B*35:01": {
        "p2": {"P": 2.6, "A": 1.0, "S": 0.6, "D": -1.6},
        "pc": {"Y": 2.6, "F": 2.0, "L": 1.5, "M": 1.4, "W": 1.2},
        "p3": {"A": 0.4, "L": 0.4},
    },
    "HLA-B*40:01": {
        "p2": {"E": 2.8, "D": 1.8, "Q": 0.5, "R": -1.8, "K": -1.8},
        "pc": {"L": 2.5, "I": 1.6, "V": 1.4, "A": 1.0},
        "p3": {"F": 0.5, "I": 0.5},
    },
    "HLA-B*44:02": {
        "p2": {"E": 2.9, "D": 1.7, "R": -1.8, "K": -1.8, "P": -1.5},
        "pc": {"F": 2.4, "Y": 2.3, "W": 2.0, "L": 1.4},
        "p3": {"I": 0.5, "L": 0.4},
    },
    "HLA-B*57:01": {
        "p2": {"A": 2.1, "T": 2.0, "S": 1.9, "V": 1.0},
        "pc": {"W": 2.9, "F": 2.2, "Y": 1.6, "L": 0.8},
        "p3": {"K": 0.5, "R": 0.4},
    },
    "HLA-C*04:01": {
        "p2": {"Y": 2.4, "F": 2.2, "W": 1.4, "P": 0.8},
        "pc": {"F": 2.3, "L": 2.0, "I": 1.5, "M": 1.3},
        "p3": {"D": 0.8, "E": 0.6},
    },
    "HLA-C*07:01": {
        "p2": {"Y": 2.0, "R": 1.7, "F": 1.5, "L": 0.8},
        "pc": {"L": 2.2, "F": 2.0, "Y": 1.8, "M": 1.2},
        "p3": {"S": 0.5, "N": 0.4},
    },
    "HLA-C*07:02": {
        "p2": {"Y": 2.2, "R": 1.6, "F": 1.5, "K": 0.9},
        "pc": {"L": 2.3, "F": 2.1, "Y": 1.7, "I": 1.2},
        "p3": {"S": 0.5, "N": 0.4},
    },
}


def _normalize_allele_name(raw: str) -> str:
    s = str(raw).strip().upper().replace("'", "")
    m = re.match(r"^(?:HLA-)?([A-Z0-9]+)\*?(\d{2,3}):(\d{2,3})", s)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"
    if not s.startswith("HLA-"):
        return f"HLA-{s}"
    return s


def _resolve_motif(allele: str) -> dict[str, dict[str, float]]:
    norm = _normalize_allele_name(allele)
    if norm in HLA_MOTIFS:
        return HLA_MOTIFS[norm]
    prefix = norm.split(":")[0]
    for k, m in HLA_MOTIFS.items():
        if k.startswith(prefix + ":"):
            return m
    loc = prefix.split("*")[0]
    for k, m in HLA_MOTIFS.items():
        if k.startswith(loc + "*"):
            return m
    return HLA_MOTIFS["HLA-A*02:01"]


def estimate_pssm_ic50(peptide: str, hla_allele: str) -> float:
    """Estimate MHC-I binding IC50 (nM) from curated anchor-motif PSSM."""
    pep = str(peptide).strip().upper()
    if len(pep) < 8 or not all(c in STANDARD_AAS for c in pep):
        return 25000.0
    motif = _resolve_motif(hla_allele)
    p2_sc = motif["p2"].get(pep[1], -0.4)
    pc_sc = motif["pc"].get(pep[-1], -0.4)
    p3_sc = motif.get("p3", {}).get(pep[2], 0.0)
    # Length preference centered on 9-mer (8/10/11-mer carry mild structural penalty)
    len_pen = {8: -0.4, 9: 0.25, 10: -0.15, 11: -0.5}.get(len(pep), -1.0)
    # Mild penalty for multiple prolines or extreme negative charge at P1/P_C-1
    p1_bonus = 0.25 if pep[0] in {"Y", "F", "K", "R", "A", "M", "L", "S"} else (-0.35 if pep[0] in {"D", "E", "P"} else 0.0)

    score = p2_sc + pc_sc + 0.6 * p3_sc + len_pen + p1_bonus
    # Calibrated logistic mapping:
    # score >= 4.5 -> ~12-40 nM (strong binder)
    # score ~ 2.5..4.0 -> ~80-450 nM (weak binder)
    # score <= 1.0 -> >2000 nM (non-binder)
    prob = 1.0 / (1.0 + math.exp(-0.95 * (score - 1.9)))
    ic50 = 5.0 * (10000.0 ** (1.0 - prob))
    return float(round(np.clip(ic50, 5.0, 50000.0), 2))


def _select_best_9mer(
    mt_pep: str,
    wt_pep: str,
    alleles: list[str],
) -> tuple[str, str, str, float, float | None]:
    """Select best-binding 9-mer (or full minimal peptide) across alleles."""
    mt = str(mt_pep).strip().upper()
    wt = str(wt_pep).strip().upper() if wt_pep else ""

    if len(mt) <= 11:
        best_res: tuple[str, str, str, float, float | None] | None = None
        for al in alleles:
            ic_mt = estimate_pssm_ic50(mt, al)
            ic_wt = estimate_pssm_ic50(wt, al) if wt and len(wt) >= 8 else None
            if best_res is None or ic_mt < best_res[3]:
                best_res = (al, mt, wt, ic_mt, ic_wt)
        assert best_res is not None
        return best_res

    # Long flanking sequence (>11 aa): enumerate 9-mer windows covering differences
    k = 9
    cand_indices: list[int] = []
    if wt and len(wt) == len(mt):
        mut_pos = [i for i, (a, b) in enumerate(zip(mt, wt)) if a != b]
        for start in range(len(mt) - k + 1):
            if any(start <= mp < start + k for mp in mut_pos):
                cand_indices.append(start)
    if not cand_indices:
        cand_indices = list(range(len(mt) - k + 1))

    best_win: tuple[str, str, str, float, float | None] | None = None
    for start in cand_indices:
        sub_mt = mt[start : start + k]
        sub_wt = wt[start : start + k] if (wt and len(wt) == len(mt)) else ""
        if sub_wt == sub_mt:
            continue
        for al in alleles:
            ic_mt = estimate_pssm_ic50(sub_mt, al)
            ic_wt = estimate_pssm_ic50(sub_wt, al) if sub_wt else None
            if best_win is None or ic_mt < best_win[3]:
                best_win = (al, sub_mt, sub_wt, ic_mt, ic_wt)

    if best_win is None:
        sub_mt = mt[:k]
        al = alleles[0]
        return (al, sub_mt, wt[:k] if wt else "", estimate_pssm_ic50(sub_mt, al), None)
    return best_win


def _analyze_mutation_position(mt_ep: str, wt_ep: str) -> tuple[bool, str, str, float]:
    """Classify anchor vs TCR-contact mutation and compute foreignness score (0..1)."""
    if not wt_ep:
        return False, "unknown", "", 0.65
    if len(mt_ep) != len(wt_ep):
        return False, "frameshift_or_indel", "indel", 1.0

    diffs = [i for i, (a, b) in enumerate(zip(mt_ep, wt_ep)) if a != b]
    if not diffs:
        return False, "none", "", 0.0

    pos_str = ",".join(str(i + 1) for i in diffs)
    anchor_idx = {1, len(mt_ep) - 1}
    hit_anchor = any(i in anchor_idx for i in diffs)
    hit_tcr = any(2 <= i <= len(mt_ep) - 2 for i in diffs)

    if hit_anchor and hit_tcr:
        role = "both"
    elif hit_anchor:
        role = "anchor"
    elif hit_tcr:
        role = "tcr_contact"
    else:
        role = "n_term"

    # Physicochemical Grantham shift + TCR-contact aromatic/hydrophobic foreignness
    g_vals = [grantham_distance(mt_ep[i], wt_ep[i]) / 180.0 for i in diffs]
    mean_g = float(np.clip(np.mean(g_vals), 0.0, 1.0))
    tcr_hydro_bonus = 0.0
    for i in diffs:
        if 2 <= i <= len(mt_ep) - 2 and mt_ep[i] in {"F", "W", "Y", "L", "I", "V"}:
            tcr_hydro_bonus += 0.15
    foreignness = float(np.clip(0.35 + 0.50 * mean_g + tcr_hydro_bonus, 0.0, 1.0))
    return bool(hit_anchor), role, pos_str, round(foreignness, 4)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--hla_alleles", nargs="+", default=None)
    ap.add_argument("--ic50_cutoff_nm", type=float, default=500.0)
    ap.add_argument("--min_tpm", type=float, default=1.0)
    args = ap.parse_args()

    df = read_table(args.candidates)
    mt_col = resolve_column(df, ["mt_epitope", "mt_peptide", "mutant_peptide", "epitope", "peptide"], "mutant peptide/epitope")
    wt_col = resolve_column(df, ["wt_epitope", "wt_peptide", "wildtype_peptide"], "wildtype peptide/epitope", required=False)
    gene_col = resolve_column(df, ["gene", "hugo_symbol", "gene_symbol", "symbol"], "gene", required=False)
    hla_col = resolve_column(df, ["hla_allele", "allele", "mhc_allele", "hla"], "HLA allele", required=False)
    ic50_mt_col = resolve_column(
        df,
        ["ic50_mt", "ic50_mt_nm", "mt_ic50", "best_mt_ic50", "ic50", "aff_mt", "affinity_nm"],
        "mutant IC50",
        required=False,
    )
    ic50_wt_col = resolve_column(
        df,
        ["ic50_wt", "ic50_wt_nm", "wt_ic50", "best_wt_ic50", "aff_wt"],
        "wildtype IC50",
        required=False,
    )
    tpm_col = resolve_column(df, ["tpm", "rna_tpm", "expression_tpm", "gene_expression", "expression"], "TPM", required=False)
    vaf_col = resolve_column(df, ["vaf", "dna_vaf", "tumor_vaf"], "VAF", required=False)
    ccf_col = resolve_column(df, ["ccf", "cancer_cell_fraction"], "CCF", required=False)
    mut_col = resolve_column(df, ["mutation_type", "variant_class", "consequence"], "mutation_type", required=False)

    cli_alleles = [_normalize_allele_name(a) for a in parse_str_list(args.hla_alleles)]
    warnings: list[str] = []

    if ic50_mt_col is None and not cli_alleles and hla_col is None:
        cli_alleles = ["HLA-A*02:01"]
        warnings.append(
            "Neither pre-computed ic50_mt nor hla_alleles/hla_allele was supplied; defaulted to HLA-A*02:01 for built-in PSSM estimation."
        )

    used_provided = 0
    used_pssm = 0
    out_rows: list[dict[str, Any]] = []

    for idx in range(len(df)):
        row = df.iloc[idx]
        raw_mt = str(row[mt_col]).strip().upper() if pd.notna(row[mt_col]) else ""
        if not raw_mt or raw_mt in {"NAN", "NONE", "NA"}:
            continue
        raw_wt = str(row[wt_col]).strip().upper() if wt_col and pd.notna(row[wt_col]) else ""
        if raw_wt in {"NAN", "NONE", "NA"}:
            raw_wt = ""

        gene = str(row[gene_col]).strip() if gene_col and pd.notna(row[gene_col]) else f"CAND_{idx + 1}"
        mut_type = str(row[mut_col]).strip() if mut_col and pd.notna(row[mut_col]) else "missense"

        row_alleles: list[str] = []
        if hla_col and pd.notna(row[hla_col]) and str(row[hla_col]).strip():
            row_alleles = [_normalize_allele_name(x) for x in re.split(r"[,;/\s]+", str(row[hla_col]).strip()) if x]
        if not row_alleles:
            row_alleles = cli_alleles if cli_alleles else ["HLA-A*02:01"]

        ic_mt_val = pd.to_numeric(row[ic50_mt_col], errors="coerce") if ic50_mt_col else np.nan
        ic_wt_val = pd.to_numeric(row[ic50_wt_col], errors="coerce") if ic50_wt_col else np.nan

        if pd.notna(ic_mt_val) and float(ic_mt_val) > 0:
            ic_mt = float(ic_mt_val)
            chosen_allele = row_alleles[0]
            mt_ep = raw_mt
            wt_ep = raw_wt
            if pd.notna(ic_wt_val) and float(ic_wt_val) > 0:
                ic_wt = float(ic_wt_val)
            elif wt_ep and len(wt_ep) >= 8:
                ic_wt = estimate_pssm_ic50(wt_ep, chosen_allele)
            else:
                ic_wt = None
            src_mode = "provided"
            used_provided += 1
        else:
            chosen_allele, mt_ep, wt_ep, ic_mt, ic_wt = _select_best_9mer(raw_mt, raw_wt, row_alleles)
            src_mode = "pssm_estimate"
            used_pssm += 1

        if ic_mt <= 50.0:
            tier = "strong_binder"
        elif ic_mt <= args.ic50_cutoff_nm:
            tier = "weak_binder"
        else:
            tier = "non_binder"

        dai = float(round(ic_wt / max(ic_mt, 1e-3), 4)) if (ic_wt is not None and ic_wt > 0) else None
        anchor_mut, mut_role, mut_pos, foreignness = _analyze_mutation_position(mt_ep, wt_ep)
        gravy = kyte_doolittle_gravy(mt_ep)

        tpm = float(pd.to_numeric(row[tpm_col], errors="coerce")) if tpm_col else np.nan
        vaf = float(pd.to_numeric(row[vaf_col], errors="coerce")) if vaf_col else np.nan
        ccf = float(pd.to_numeric(row[ccf_col], errors="coerce")) if ccf_col else np.nan

        pass_bind = bool(ic_mt <= args.ic50_cutoff_nm)
        pass_expr = bool(tpm >= args.min_tpm) if np.isfinite(tpm) else True
        if np.isfinite(ccf):
            pass_clonal = bool(ccf >= 0.5)
        elif np.isfinite(vaf):
            pass_clonal = bool(vaf >= 0.1)
        else:
            pass_clonal = True
        pass_all = bool(pass_bind and pass_expr and pass_clonal)

        # Composite Immunogenicity Score (0..100)
        # 1) Binding affinity (0..40)
        bind_pts = 40.0 * float(np.clip((math.log10(5000.0) - math.log10(max(ic_mt, 5.0))) / (math.log10(5000.0) - math.log10(5.0)), 0.0, 1.0))
        # 2) Agretopicity DAI & anchor/TCR balance (0..25)
        if dai is not None:
            dai_pts = 18.0 * float(np.clip(math.log10(max(dai, 0.5) * 2.0) / math.log10(40.0), 0.0, 1.0))
        else:
            dai_pts = 8.0
        role_bonus = 7.0 if mut_role in {"both", "tcr_contact", "frameshift_or_indel"} else (5.0 if mut_role == "anchor" else 2.0)
        dai_total = min(25.0, dai_pts + role_bonus)
        # 3) Foreignness & peptide hydropathy (0..15)
        hydro_bonus = 3.0 if (gravy is not None and -0.5 <= gravy <= 2.0) else 1.0
        foreign_pts = min(15.0, 12.0 * foreignness + hydro_bonus)
        # 4) Expression & Clonality (0..20)
        expr_pts = 10.0 * float(np.clip(math.log2(max(tpm, 0.0) + 1.0) / math.log2(51.0), 0.0, 1.0)) if np.isfinite(tpm) else 5.0
        if np.isfinite(ccf):
            clonal_pts = 10.0 * float(np.clip(ccf, 0.0, 1.0))
        elif np.isfinite(vaf):
            clonal_pts = 10.0 * float(np.clip(2.0 * vaf, 0.0, 1.0))
        else:
            clonal_pts = 5.0
        comp_score = round(float(np.clip(bind_pts + dai_total + foreign_pts + expr_pts + clonal_pts, 0.0, 100.0)), 2)

        out_rows.append(
            {
                "gene": gene,
                "mutation_type": mut_type,
                "hla_allele": chosen_allele,
                "mt_epitope": mt_ep,
                "wt_epitope": wt_ep,
                "ic50_mt_nm": round(ic_mt, 2),
                "ic50_wt_nm": round(ic_wt, 2) if ic_wt is not None else np.nan,
                "ic50_source": src_mode,
                "binding_tier": tier,
                "agretopicity_dai": dai if dai is not None else np.nan,
                "anchor_mutation": anchor_mut,
                "mutation_role": mut_role,
                "mutated_positions": mut_pos,
                "mt_kd_hydropathy": round(gravy, 3) if gravy is not None else np.nan,
                "foreignness_score": foreignness,
                "tpm": round(tpm, 2) if np.isfinite(tpm) else np.nan,
                "vaf": round(vaf, 4) if np.isfinite(vaf) else np.nan,
                "ccf": round(ccf, 4) if np.isfinite(ccf) else np.nan,
                "passes_binding_gate": pass_bind,
                "passes_expression_gate": pass_expr,
                "passes_clonality_gate": pass_clonal,
                "passes_all_gates": pass_all,
                "immunogenicity_score": comp_score,
            }
        )

    if not out_rows:
        fail("No valid neoantigen candidates parsed from input table.")

    if used_pssm > 0:
        warnings.append(
            f"Built-in anchor-motif PSSM IC50 estimation was used for {used_pssm} candidate(s); validate top hits with NetMHCpan or MHCflurry before clinical selection."
        )

    ic50_source_summary = (
        "provided" if used_pssm == 0 else ("pssm_estimate" if used_provided == 0 else "mixed")
    )

    out_df = pd.DataFrame(out_rows)
    out_df = out_df.sort_values(["passes_all_gates", "immunogenicity_score"], ascending=[False, False]).reset_index(drop=True)
    out_df.insert(0, "rank", np.arange(1, len(out_df) + 1))

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    params = {
        "hla_alleles": cli_alleles,
        "ic50_cutoff_nm": float(args.ic50_cutoff_nm),
        "min_tpm": float(args.min_tpm),
        "ic50_source": ic50_source_summary,
    }
    sci_ctx = {
        "n_candidates": int(len(out_df)),
        "n_passed_all_gates": int(out_df["passes_all_gates"].sum()),
    }

    tsv_path = write_tsv(
        out_df,
        out_dir / "prioritized_neoantigens.tsv",
        tool="imm_neoantigen_prioritize",
        title="Prioritized Neoantigen Candidates (IC50, DAI, TESLA Gates & Immunogenicity Score)",
        inputs=[args.candidates],
        parameters=params,
        scientific_context=sci_ctx,
    )

    # 2-panel Neoantigen Landscape figure
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.9))

    log_ic50 = np.log10(np.clip(out_df["ic50_mt_nm"].to_numpy(dtype=float, copy=True), 1.0, 50000.0))
    dai_vals = out_df["agretopicity_dai"].fillna(1.0).to_numpy(dtype=float, copy=True)
    log_dai = np.log10(np.clip(dai_vals, 0.05, 1000.0))
    scores = out_df["immunogenicity_score"].to_numpy(dtype=float, copy=True)
    passed = out_df["passes_all_gates"].to_numpy(dtype=bool, copy=True)

    sc = ax1.scatter(
        log_ic50[passed],
        log_dai[passed],
        c=scores[passed],
        cmap="viridis",
        vmin=0,
        vmax=100,
        s=70,
        edgecolors="#0f172a",
        linewidth=0.7,
        label=f"Passed gates ({int(passed.sum())})",
    )
    if np.any(~passed):
        ax1.scatter(
            log_ic50[~passed],
            log_dai[~passed],
            c="#94a3b8",
            s=36,
            alpha=0.65,
            marker="x",
            label=f"Filtered ({int((~passed).sum())})",
        )
    ax1.axvline(math.log10(50.0), color="#dc2626", linestyle="--", linewidth=1.0, label="50 nM (Strong)")
    ax1.axvline(math.log10(max(args.ic50_cutoff_nm, 1.0)), color="#f59e0b", linestyle=":", linewidth=1.2, label=f"{args.ic50_cutoff_nm:g} nM cutoff")
    ax1.axhline(math.log10(2.0), color="#0284c7", linestyle="--", linewidth=0.9, alpha=0.7, label="DAI = 2.0")

    for _, r in out_df.head(6).iterrows():
        lx = math.log10(max(float(r["ic50_mt_nm"]), 1.0))
        ly = math.log10(np.clip(float(r["agretopicity_dai"]) if pd.notna(r["agretopicity_dai"]) else 1.0, 0.05, 1000.0))
        ax1.annotate(f"{r['gene']} ({r['mt_epitope'][:6]})", (lx, ly), textcoords="offset points", xytext=(4, 4), fontsize=7.5)

    ax1.set_xlabel("log10 Mutant IC50 (nM)")
    ax1.set_ylabel("log10 Agretopicity (DAI = IC50_WT / IC50_MT)")
    ax1.set_title("Neoantigen Binding vs. Agretopicity Landscape")
    ax1.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax1.grid(alpha=0.2, linestyle="--")
    fig.colorbar(sc, ax=ax1, label="Immunogenicity Score", fraction=0.046, pad=0.04)

    top_plot = out_df.head(15).iloc[::-1]
    labels = [f"{r['gene']} | {r['mt_epitope']}" for _, r in top_plot.iterrows()]
    y_pos = np.arange(len(labels))
    tier_colors = {
        "strong_binder": "#0f766e",
        "weak_binder": "#0284c7",
        "non_binder": "#94a3b8",
    }
    bar_cols = [tier_colors.get(t, "#94a3b8") for t in top_plot["binding_tier"].tolist()]
    ax2.barh(y_pos, top_plot["immunogenicity_score"].to_numpy(dtype=float, copy=True), color=bar_cols)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(labels, fontsize=8)
    ax2.set_xlim(0, 105)
    ax2.set_xlabel("Composite Immunogenicity Score (0–100)")
    ax2.set_title("Top Prioritized Neoantigen Candidates")
    ax2.grid(axis="x", alpha=0.2, linestyle="--")

    fig.tight_layout()
    fig_file = out_dir / "neoantigen_landscape.png"
    fig.savefig(fig_file, bbox_inches="tight")
    plt.close(fig)
    fig_path = write_artifact(
        fig_file,
        tool="imm_neoantigen_prioritize",
        title="Neoantigen Binding, Agretopicity & Immunogenicity Prioritization",
        inputs=[args.candidates],
        parameters=params,
        scientific_context=sci_ctx,
    )

    n_strong = int((out_df["binding_tier"] == "strong_binder").sum())
    n_weak = int((out_df["binding_tier"] == "weak_binder").sum())
    n_passed = int(out_df["passes_all_gates"].sum())
    n_high_dai = int((out_df["agretopicity_dai"].fillna(0.0) >= 2.0).sum())
    n_anchor = int(out_df["anchor_mutation"].sum())

    emit(
        {
            "status": "ok",
            "tool": "imm_neoantigen_prioritize",
            "ic50_source": ic50_source_summary,
            "ic50_cutoff_nm": float(args.ic50_cutoff_nm),
            "min_tpm": float(args.min_tpm),
            "n_candidates": int(len(out_df)),
            "n_strong_binders": n_strong,
            "n_weak_binders": n_weak,
            "n_passed_all_gates": n_passed,
            "n_high_dai": n_high_dai,
            "n_anchor_mutations": n_anchor,
            "top_candidates": out_df.head(15).to_dict(orient="records"),
            "outputs": {
                "prioritized_neoantigens": tsv_path,
                "neoantigen_landscape": fig_path,
            },
            "artifacts": [tsv_path, fig_path],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
