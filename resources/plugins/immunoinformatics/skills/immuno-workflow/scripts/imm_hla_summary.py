#!/usr/bin/env python3
"""imm_hla_summary: Standardize OptiType, arcasHLA, seq2HLA, HLA-LA, and
clinical HLA tables to 2-field (4-digit) alleles across Class I (HLA-A/B/C) and
Class II (HLA-DRB1/DQB1/DPB1), flag Class I homozygosity, and compute
approximate pocket-Grantham HLA Evolutionary Divergence (HED).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    grantham_distance,
    parse_str_list,
    read_table,
    resolve_column,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)

CLASS1_LOCI = ["HLA-A", "HLA-B", "HLA-C"]
CLASS2_LOCI = ["HLA-DRB1", "HLA-DQB1", "HLA-DPB1"]
ALL_CANONICAL_LOCI = CLASS1_LOCI + CLASS2_LOCI

# Curated 12-residue peptide-binding groove pocket pseudo-sequences (B/C/D/E/F
# pocket polymorphic positions) per major HLA allele family, derived from
# canonical MHC supertype / pocket physicochemical classifications.
POCKET_PSEUDOSEQS: dict[str, str] = {
    # HLA-A families
    "A*01": "YFDAMNVKTQKY",
    "A*02": "YFYAVHVKTHTY",
    "A*03": "YFDAVNVDTQQD",
    "A*11": "YYDAVNVDAQRY",
    "A*23": "SFDAVNVEERLY",
    "A*24": "SFDAVNVEEHLY",
    "A*25": "YFDAMNVEKRSY",
    "A*26": "YFDAMNVDKRNY",
    "A*29": "TLDAMNVKTQRY",
    "A*30": "SFDAINVETQQY",
    "A*31": "TFDAINVDTQRY",
    "A*32": "FFDAMNVEEQLY",
    "A*33": "TFDAINVDTQRY",
    "A*68": "YFDAVNVDTHRY",
    "A*69": "YFYAVHVDTHRY",
    # HLA-B families
    "B*07": "YSEYRNQETQIS",
    "B*08": "YSDYRNQETQNS",
    "B*13": "YYTYRETEIQLS",
    "B*14": "YSDYRNVETQTS",
    "B*15": "YYAYRETETQLY",
    "B*18": "YHAYREKETQTS",
    "B*27": "YHEYCKTETQLS",
    "B*35": "YYAYRNQETQLY",
    "B*37": "YHDYREKETQLS",
    "B*38": "YSEYICVETQTS",
    "B*39": "YSEYICVETQTS",
    "B*40": "YHTYKEKETQTS",
    "B*41": "YHTYKEKETQTS",
    "B*44": "YYTYKEKETQLY",
    "B*45": "YHTYKEKETQTS",
    "B*49": "YHTYKEKETQLW",
    "B*50": "YHTYKEKETQTS",
    "B*51": "YYAYRNTEIQLY",
    "B*52": "YYAYRNTEIQLY",
    "B*53": "YYAYRNQETQLY",
    "B*55": "YSEYRNQETQIS",
    "B*56": "YSEYRNQETQIS",
    "B*57": "YYAYRNMKIQLY",
    "B*58": "YYAYRNMKIQLY",
    # HLA-C families
    "C*01": "YCSYKYRKTQLS",
    "C*02": "YYSYKYRKTQLW",
    "C*03": "YYSYKYRDTQLY",
    "C*04": "YFSYKYRKTQDF",
    "C*05": "YYSYKYRKTQKF",
    "C*06": "YDSYKYRKTQLW",
    "C*07": "YDSYKYREKQLS",
    "C*08": "YCSYKYRKTQKF",
    "C*12": "YYSYKYREKQLW",
    "C*14": "YFSYKYRETQLF",
    "C*15": "YYSYKYRETQLW",
    "C*16": "YYSYKYRKTQLS",
    "C*17": "YYSYKYRKTQLS",
    # Class II DRB1 / DQB1 / DPB1 pocket pseudo-sequences
    "DRB1*01": "FLWQFEYDRVYG",
    "DRB1*03": "YLYSFEYNRVYD",
    "DRB1*04": "FLHQYEYDRVYQ",
    "DRB1*07": "WLWQFKYDKVYG",
    "DRB1*08": "YLYSFEYDRVYN",
    "DRB1*09": "WLKQFEYDKVYD",
    "DRB1*10": "VLWQFEYDRRYG",
    "DRB1*11": "YLYSFEYDRVYD",
    "DRB1*12": "YLYSFEYDRVYN",
    "DRB1*13": "YLYSFEYDRVYD",
    "DRB1*14": "YLYSFEYDRVYR",
    "DRB1*15": "PLWQFEYARVYG",
    "DRB1*16": "PLWQFEYDRVYG",
    "DQB1*02": "YGYQFAYDRVYD",
    "DQB1*03": "YAYQFAYDRVYQ",
    "DQB1*04": "YLYQFAYDRVYD",
    "DQB1*05": "YGVQFAYQRVYG",
    "DQB1*06": "YGVQFAYDRVYG",
    "DPB1*01": "YGYQFAEDRVYD",
    "DPB1*02": "YGFQFAEDRVYE",
    "DPB1*03": "YLYQFAEDRVYD",
    "DPB1*04": "YGFQFAEDRVYG",
    "DPB1*05": "YLYQFAEDRVYE",
    "DPB1*06": "YLYQFAEDRVYD",
}

TOOL_PRIORITY = {
    "optitype": 4,
    "hlala": 4,
    "arcashla": 3,
    "seq2hla": 2,
    "clinical_table": 1,
}


def normalize_hla_allele(raw: Any, hint_locus: str | None = None) -> tuple[str | None, str | None, bool]:
    """Normalize an HLA allele string to canonical 'HLA-<LOCUS>*<F1>:<F2>' (2-field).

    Returns (canonical_locus, normalized_allele, is_two_field).
    """
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return None, None, False
    s = str(raw).strip().replace("'", "").replace('"', "")
    if not s or s.lower() in {"na", "nan", "none", "no", "null", ".", "-", "unresolved"}:
        return None, None, False

    s_up = s.upper()
    # Handle format like HLA-A*02:01:01G or A*02:01 or A02:01 or A*02
    m = re.match(r"^(?:HLA-)?([A-Z0-9]+)\*?(\d{2,3}):(\d{2,3})(?::\d+)*[A-Z]*$", s_up)
    if m:
        loc, f1, f2 = m.group(1), m.group(2), m.group(3)
        canon_loc = f"HLA-{loc}"
        return canon_loc, f"{canon_loc}*{f1}:{f2}", True

    m1 = re.match(r"^(?:HLA-)?([A-Z0-9]+)\*?(\d{2,3})[A-Z]*$", s_up)
    if m1:
        loc, f1 = m1.group(1), m1.group(2)
        canon_loc = f"HLA-{loc}"
        return canon_loc, f"{canon_loc}*{f1}", False

    if hint_locus:
        loc = hint_locus.upper().replace("HLA-", "")
        m_num = re.match(r"^(\d{2,3}):(\d{2,3})(?::\d+)*[A-Z]*$", s_up)
        if m_num:
            canon_loc = f"HLA-{loc}"
            return canon_loc, f"{canon_loc}*{m_num.group(1)}:{m_num.group(2)}", True
    return None, None, False


def _allele_pocket_seq(allele: str) -> tuple[str, int]:
    """Return 12-aa curated pocket pseudo-sequence and 2nd-field integer for micro-divergence."""
    body = allele.replace("HLA-", "")
    parts = body.split("*")
    loc = parts[0]
    fields = parts[1].split(":") if len(parts) > 1 else ["01"]
    f1 = fields[0]
    f2 = int(re.sub(r"\D", "", fields[1])) if len(fields) > 1 and re.sub(r"\D", "", fields[1]) else 1
    key = f"{loc}*{f1}"
    if key in POCKET_PSEUDOSEQS:
        return POCKET_PSEUDOSEQS[key], f2
    # Fallback to first available family at same locus
    for k, seq in POCKET_PSEUDOSEQS.items():
        if k.startswith(f"{loc}*"):
            return seq, f2
    return "YFYAVHVKTHTY", f2


def approximate_hed(allele_1: str | None, allele_2: str | None) -> float | None:
    """Compute approximate pocket-Grantham HLA Evolutionary Divergence (HED)."""
    if not allele_1 or not allele_2:
        return None
    if allele_1 == allele_2:
        return 0.0
    seq1, f2_1 = _allele_pocket_seq(allele_1)
    seq2, f2_2 = _allele_pocket_seq(allele_2)
    # Mean Grantham distance across the 12 curated pocket contact residues (scaled by 0.15
    # to match published ~0-16 HLA Evolutionary Divergence units across full groove)
    dists = [grantham_distance(a, b) for a, b in zip(seq1, seq2)]
    mean_g = float(sum(dists) / len(dists)) * 0.15
    if seq1 == seq2 and f2_1 != f2_2:
        # Distinct 2nd-field alleles within the same allele group share primary pocket motif
        # but carry minor alpha-helix / floor polymorphisms (~0.8 to 2.5 HED units)
        delta_f2 = abs(f2_1 - f2_2)
        mean_g = float(min(3.0, 0.8 + 0.15 * min(delta_f2, 12)))
    return round(mean_g, 4)


def _infer_sample_from_hla_filename(path: Path) -> str:
    stem = path.name
    for ext in (".genotype.json", ".genes.json", ".json", ".tsv", ".csv", ".txt"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
            break
    for suf in (
        "-ClassI-class.HLAgenotype4digits",
        ".ClassI-class.HLAgenotype4digits",
        "-ClassII.HLAgenotype4digits",
        ".ClassII.HLAgenotype4digits",
        "-ClassI-class.HLAgenotype2digits",
        ".ClassI-class.HLAgenotype2digits",
        "-ClassII.HLAgenotype2digits",
        ".ClassII.HLAgenotype2digits",
        ".HLAgenotype4digits",
        ".HLAgenotype2digits",
        "_ClassI-class",
        "_ClassII",
        "-ClassI-class",
        "-ClassII",
        "_bestguess_G",
        "_result",
    ):
        if stem.endswith(suf) or stem.lower().endswith(suf.lower()):
            stem = stem[: -len(suf)]
            break
    if stem.lower() in {"r1_bestguess_g", "bestguess_g", "r1", "result", "genotype"} and path.parent.name:
        return path.parent.name
    return stem or path.stem


def _parse_hla_file(path: Path) -> tuple[list[dict[str, Any]], str, list[str]]:
    """Parse one HLA input file into call records: (sample_id, locus, alleles, source_tool)."""
    records: list[dict[str, Any]] = []
    warns: list[str] = []
    fallback_sid = _infer_sample_from_hla_filename(path)

    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(data, dict):
            # Check if flat arcasHLA {A: [...], B: [...]} or nested {sample: {A: [...]}}
            first_val = next(iter(data.values()), None)
            if isinstance(first_val, dict):
                for sid, loc_map in data.items():
                    for loc_raw, alleles_raw in loc_map.items():
                        alleles = [alleles_raw] if isinstance(alleles_raw, str) else list(alleles_raw or [])
                        norm_list = []
                        canon_loc = None
                        for a in alleles:
                            cl, na, is2f = normalize_hla_allele(a, hint_locus=str(loc_raw))
                            if na:
                                canon_loc = canon_loc or cl
                                norm_list.append(na)
                                if not is2f:
                                    warns.append(f"Low-resolution (<2-field) allele '{a}' in {path.name}.")
                        if canon_loc and norm_list:
                            if len(norm_list) == 1:
                                norm_list = [norm_list[0], norm_list[0]]
                            records.append({"sample_id": str(sid), "locus": canon_loc, "alleles": norm_list[:2], "tool": "arcashla"})
            else:
                for loc_raw, alleles_raw in data.items():
                    alleles = [alleles_raw] if isinstance(alleles_raw, str) else list(alleles_raw or [])
                    norm_list = []
                    canon_loc = None
                    for a in alleles:
                        cl, na, is2f = normalize_hla_allele(a, hint_locus=str(loc_raw))
                        if na:
                            canon_loc = canon_loc or cl
                            norm_list.append(na)
                            if not is2f:
                                warns.append(f"Low-resolution (<2-field) allele '{a}' in {path.name}.")
                    if canon_loc and norm_list:
                        if len(norm_list) == 1:
                            norm_list = [norm_list[0], norm_list[0]]
                        records.append({"sample_id": fallback_sid, "locus": canon_loc, "alleles": norm_list[:2], "tool": "arcashla"})
        if not records:
            fail(f"No valid HLA alleles found in JSON file: {path}")
        return records, "arcashla", warns

    df = read_table(path)
    cols_low = {str(c).strip().lower(): str(c).strip() for c in df.columns}

    # Detect seq2HLA
    if "locus" in cols_low and ("allele 1" in cols_low or "allele1" in cols_low):
        l_col = cols_low["locus"]
        a1_col = cols_low.get("allele 1") or cols_low["allele1"]
        a2_col = cols_low.get("allele 2") or cols_low.get("allele2")
        s_col = resolve_column(df, ["sample_id", "sample"], "sample_id", required=False)
        for _, row in df.iterrows():
            sid = str(row[s_col]).strip() if s_col and pd.notna(row[s_col]) else fallback_sid
            loc_hint = str(row[l_col]).strip()
            alleles = []
            for ac in (a1_col, a2_col):
                if ac and pd.notna(row[ac]):
                    cl, na, is2f = normalize_hla_allele(row[ac], hint_locus=loc_hint)
                    if na:
                        alleles.append((cl, na))
                        if not is2f:
                            warns.append(f"Low-resolution (<2-field) allele '{row[ac]}' in {path.name}.")
            if alleles:
                cl = alleles[0][0]
                alist = [a[1] for a in alleles]
                if len(alist) == 1:
                    alist = [alist[0], alist[0]]
                records.append({"sample_id": sid, "locus": cl, "alleles": alist[:2], "tool": "seq2hla"})
        return records, "seq2hla", warns

    # Detect HLA-LA
    if "locus" in cols_low and "chromosome" in cols_low and "allele" in cols_low:
        l_col = cols_low["locus"]
        a_col = cols_low["allele"]
        s_col = resolve_column(df, ["sample_id", "sample"], "sample_id", required=False)
        grouped: dict[tuple[str, str], list[str]] = defaultdict(list)
        for _, row in df.iterrows():
            sid = str(row[s_col]).strip() if s_col and pd.notna(row[s_col]) else fallback_sid
            cl, na, is2f = normalize_hla_allele(row[a_col], hint_locus=str(row[l_col]))
            if cl and na:
                grouped[(sid, cl)].append(na)
                if not is2f:
                    warns.append(f"Low-resolution (<2-field) allele '{row[a_col]}' in {path.name}.")
        for (sid, cl), alist in grouped.items():
            if len(alist) == 1:
                alist = [alist[0], alist[0]]
            records.append({"sample_id": sid, "locus": cl, "alleles": alist[:2], "tool": "hlala"})
        return records, "hlala", warns

    # Detect OptiType
    if {"a1", "a2", "b1", "b2", "c1", "c2"}.issubset(set(cols_low.keys())):
        s_col = resolve_column(df, ["sample_id", "sample"], "sample_id", required=False)
        for _, row in df.iterrows():
            sid = fallback_sid
            if s_col and pd.notna(row[s_col]) and str(row[s_col]).strip():
                sid = str(row[s_col]).strip()
            elif "unnamed: 0" in cols_low:
                v0 = str(row[cols_low["unnamed: 0"]]).strip()
                if v0 and not v0.isdigit():
                    sid = v0
            for loc_short in ("A", "B", "C"):
                c1 = cols_low[f"{loc_short.lower()}1"]
                c2 = cols_low[f"{loc_short.lower()}2"]
                alist = []
                for ac in (c1, c2):
                    if pd.notna(row[ac]):
                        _, na, is2f = normalize_hla_allele(row[ac], hint_locus=loc_short)
                        if na:
                            alist.append(na)
                            if not is2f:
                                warns.append(f"Low-resolution (<2-field) allele '{row[ac]}' in {path.name}.")
                if alist:
                    if len(alist) == 1:
                        alist = [alist[0], alist[0]]
                    records.append({"sample_id": sid, "locus": f"HLA-{loc_short}", "alleles": alist[:2], "tool": "optitype"})
        return records, "optitype", warns

    # Generic clinical cohort table (long or wide)
    s_col = resolve_column(df, ["sample_id", "sample", "patient_id", "subject_id"], "sample_id", required=False)
    if "allele" in cols_low or "hla_allele" in cols_low:
        a_col = cols_low.get("allele") or cols_low["hla_allele"]
        l_col = cols_low.get("locus")
        grouped_long: dict[tuple[str, str], list[str]] = defaultdict(list)
        for _, row in df.iterrows():
            sid = str(row[s_col]).strip() if s_col and pd.notna(row[s_col]) else fallback_sid
            hint = str(row[l_col]).strip() if l_col and pd.notna(row[l_col]) else None
            cl, na, is2f = normalize_hla_allele(row[a_col], hint_locus=hint)
            if cl and na:
                grouped_long[(sid, cl)].append(na)
                if not is2f:
                    warns.append(f"Low-resolution (<2-field) allele '{row[a_col]}' in {path.name}.")
        for (sid, cl), alist in grouped_long.items():
            if len(alist) == 1:
                alist = [alist[0], alist[0]]
            records.append({"sample_id": sid, "locus": cl, "alleles": alist[:2], "tool": "clinical_table"})
        return records, "clinical_table", warns

    # Wide clinical table: look for columns matching HLA loci
    for _, row in df.iterrows():
        sid = str(row[s_col]).strip() if s_col and pd.notna(row[s_col]) else fallback_sid
        loc_alleles: dict[str, list[str]] = defaultdict(list)
        for col in df.columns:
            if col == s_col:
                continue
            m_col = re.match(r"^(?:HLA[-_])?(A|B|C|DRB1|DQB1|DPB1)(?:[_-]?\d)?$", col.strip(), re.I)
            if not m_col:
                continue
            hint = m_col.group(1).upper()
            raw_val = row[col]
            if pd.isna(raw_val):
                continue
            for tok in re.split(r"[,;/\s]+", str(raw_val).strip()):
                if not tok:
                    continue
                cl, na, is2f = normalize_hla_allele(tok, hint_locus=hint)
                if cl and na:
                    loc_alleles[cl].append(na)
                    if not is2f:
                        warns.append(f"Low-resolution (<2-field) allele '{tok}' in {path.name}.")
        for cl, alist in loc_alleles.items():
            if len(alist) == 1:
                alist = [alist[0], alist[0]]
            records.append({"sample_id": sid, "locus": cl, "alleles": alist[:2], "tool": "clinical_table"})

    if not records:
        fail(f"Could not parse HLA alleles from table: {path}")
    return records, "clinical_table", warns


def _consensus_pair(calls: list[dict[str, Any]]) -> tuple[tuple[str, str], str]:
    """Reconcile multiple caller pairs for one (sample_id, locus)."""
    if len(calls) == 1:
        pair = tuple(sorted(calls[0]["alleles"][:2]))
        return (pair[0], pair[1]), "single_caller"

    pair_scores: dict[tuple[str, str], float] = defaultdict(float)
    pair_counts: Counter[tuple[str, str]] = Counter()
    for c in calls:
        pair = tuple(sorted(c["alleles"][:2]))
        p_key = (pair[0], pair[1])
        pair_counts[p_key] += 1
        pair_scores[p_key] += 1.0 + 0.05 * TOOL_PRIORITY.get(c["tool"], 0)

    best_pair = max(pair_scores.keys(), key=lambda k: (pair_counts[k], pair_scores[k]))
    if len(pair_counts) == 1:
        return best_pair, "concordant"
    return best_pair, f"discordant({len(pair_counts)}_calls)"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hla_inputs", nargs="+", required=True)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    input_paths = parse_str_list(args.hla_inputs)
    if not input_paths:
        fail("At least one HLA input file is required.")

    all_records: list[dict[str, Any]] = []
    detected_tools: list[str] = []
    warnings: list[str] = []

    for raw in input_paths:
        recs, tool_name, w_list = _parse_hla_file(Path(raw))
        all_records.extend(recs)
        if tool_name not in detected_tools:
            detected_tools.append(tool_name)
        warnings.extend(w_list)

    warnings = list(dict.fromkeys(warnings))
    by_sample_locus: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    sample_tools: dict[str, set[str]] = defaultdict(set)
    for rec in all_records:
        sid = rec["sample_id"]
        by_sample_locus[sid][rec["locus"]].append(rec)
        sample_tools[sid].add(rec["tool"])

    sample_ids = list(by_sample_locus.keys())
    summary_rows: list[dict[str, Any]] = []
    n_homo_samples = 0

    for sid in sample_ids:
        row: dict[str, Any] = {
            "sample_id": sid,
            "source_tools": ",".join(sorted(sample_tools[sid])),
        }
        disc_loci: list[str] = []
        c1_heds: list[float] = []
        c1_homo = 0
        c1_typed = 0

        for loc in ALL_CANONICAL_LOCI:
            short = loc.replace("HLA-", "")
            calls = by_sample_locus[sid].get(loc, [])
            if not calls:
                row[f"HLA_{short}_1"] = ""
                row[f"HLA_{short}_2"] = ""
                row[f"hed_{short}"] = np.nan
                continue
            (a1, a2), conc = _consensus_pair(calls)
            if conc.startswith("discordant"):
                disc_loci.append(f"{loc}:{conc}")
            row[f"HLA_{short}_1"] = a1
            row[f"HLA_{short}_2"] = a2
            hed_val = approximate_hed(a1, a2)
            row[f"hed_{short}"] = hed_val if hed_val is not None else np.nan

            if loc in CLASS1_LOCI:
                c1_typed += 1
                if a1 == a2:
                    c1_homo += 1
                if hed_val is not None:
                    c1_heds.append(hed_val)

        mean_c1_hed = round(float(np.mean(c1_heds)), 4) if c1_heds else None
        homo_flag = bool(c1_homo >= 1)
        if homo_flag:
            n_homo_samples += 1
        if disc_loci:
            warnings.append(f"Sample '{sid}' has discordant multi-caller HLA calls at {', '.join(disc_loci)}.")

        row["mean_class1_hed"] = mean_c1_hed if mean_c1_hed is not None else np.nan
        row["class1_typed_loci"] = c1_typed
        row["class1_homozygous_loci"] = c1_homo
        row["mhc_i_homozygosity_flag"] = homo_flag
        row["concordance_note"] = "concordant" if not disc_loci else ";".join(disc_loci)
        summary_rows.append(row)

    out_df = pd.DataFrame(summary_rows)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    params = {"detected_tools": detected_tools}
    sci_ctx = {
        "n_samples": len(summary_rows),
        "n_class1_homozygous_samples": n_homo_samples,
        "hed_method": "approximate_pocket_grantham",
    }

    tsv_path = write_tsv(
        out_df,
        out_dir / "hla_consensus_genotypes.tsv",
        tool="imm_hla_summary",
        title="2-Field Consensus HLA Genotypes, Zygosity & Approximate HED",
        inputs=input_paths,
        parameters=params,
        scientific_context=sci_ctx,
    )

    # 2-panel figure: Class I HED + Zygosity heatmap
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.6))

    sids = out_df["sample_id"].tolist()
    x = np.arange(len(sids))
    width = 0.25
    hed_a = out_df["hed_A"].fillna(0.0).to_numpy(dtype=float, copy=True)
    hed_b = out_df["hed_B"].fillna(0.0).to_numpy(dtype=float, copy=True)
    hed_c = out_df["hed_C"].fillna(0.0).to_numpy(dtype=float, copy=True)
    hed_mean = out_df["mean_class1_hed"].to_numpy(dtype=float, copy=True)

    ax1.bar(x - width, hed_a, width=width, color="#0f766e", label="HLA-A HED")
    ax1.bar(x, hed_b, width=width, color="#0284c7", label="HLA-B HED")
    ax1.bar(x + width, hed_c, width=width, color="#f59e0b", label="HLA-C HED")
    valid_m = np.isfinite(hed_mean)
    if np.any(valid_m):
        ax1.plot(x[valid_m], hed_mean[valid_m], "o--", color="#0f172a", linewidth=1.4, label="Mean Class I HED")
    ax1.set_xticks(x)
    ax1.set_xticklabels(sids, rotation=30, ha="right")
    ax1.set_ylabel("Approximate Pocket-Grantham HED")
    ax1.set_title("MHC Class I Evolutionary Divergence (Approx. HED)")
    ax1.legend(frameon=False, fontsize=8)
    ax1.grid(axis="y", alpha=0.25, linestyle="--")

    # Panel B: Zygosity status grid across 6 loci (-1=Not Typed, 0=Homozygous, 1=Heterozygous)
    loci_short = [l.replace("HLA-", "") for l in ALL_CANONICAL_LOCI]
    z_mat = np.full((len(sids), len(loci_short)), -1.0, dtype=float)
    z_txt = [["" for _ in loci_short] for _ in sids]
    for i, r in enumerate(summary_rows):
        for j, sh in enumerate(loci_short):
            a1 = r[f"HLA_{sh}_1"]
            a2 = r[f"HLA_{sh}_2"]
            if not a1 or not a2:
                z_mat[i, j] = -1.0
                z_txt[i][j] = "NT"
            elif a1 == a2:
                z_mat[i, j] = 0.0
                z_txt[i][j] = "HOM"
            else:
                z_mat[i, j] = 1.0
                z_txt[i][j] = "HET"

    from matplotlib.colors import ListedColormap

    cmap = ListedColormap(["#e2e8f0", "#ef4444", "#10b981"])
    ax2.imshow(z_mat, cmap=cmap, vmin=-1.0, vmax=1.0, aspect="auto")
    ax2.set_xticks(np.arange(len(loci_short)))
    ax2.set_xticklabels(loci_short)
    ax2.set_yticks(np.arange(len(sids)))
    ax2.set_yticklabels(sids)
    for i in range(len(sids)):
        for j in range(len(loci_short)):
            ax2.text(j, i, z_txt[i][j], ha="center", va="center", fontsize=8.5, fontweight="bold", color="#0f172a")
    ax2.set_title("HLA Zygosity Status (HET / HOM / NT)")

    fig.tight_layout()
    fig_file = out_dir / "hla_zygosity_hed.png"
    fig.savefig(fig_file, bbox_inches="tight")
    plt.close(fig)
    fig_path = write_artifact(
        fig_file,
        tool="imm_hla_summary",
        title="HLA Class I/II Zygosity & Approximate HED Summary",
        inputs=input_paths,
        parameters=params,
        scientific_context=sci_ctx,
    )

    cohort_hed = out_df["mean_class1_hed"].dropna().to_numpy(dtype=float, copy=True)
    mean_cohort_hed = round(float(np.mean(cohort_hed)), 4) if cohort_hed.size else None

    emit(
        {
            "status": "ok",
            "tool": "imm_hla_summary",
            "n_samples": len(summary_rows),
            "detected_tools": detected_tools,
            "n_class1_homozygous_samples": n_homo_samples,
            "mean_cohort_class1_hed": mean_cohort_hed,
            "samples": out_df.to_dict(orient="records"),
            "outputs": {
                "hla_consensus_genotypes": tsv_path,
                "hla_zygosity_hed": fig_path,
            },
            "artifacts": [tsv_path, fig_path],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
