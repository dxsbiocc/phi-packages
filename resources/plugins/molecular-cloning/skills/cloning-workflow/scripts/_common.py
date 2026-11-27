"""Shared pure-Python + NumPy/Pandas/Matplotlib utilities for molecular-cloning."""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import NoReturn

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
R_GAS = 1.98720425864083  # cal / (mol * K)

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# SantaLucia (1998) PNAS 95:1460-1465 unified NN parameters: (dH kcal/mol, dS cal/(mol*K))
NN_PARAMS: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
INIT_GC = (0.1, -2.8)
INIT_AT = (2.3, 4.1)

# Restriction enzyme table (38 standard Type II + Type IIS enzymes)
# top_cut / bot_cut are 0-based indices relative to start of 5'->3' recognition site on top strand.
RESTRICTION_ENZYMES: dict[str, dict] = {
    # Type IIS
    "BsaI":   {"site": "GGTCTC",  "top_cut": 7,  "bot_cut": 11, "overhang_type": "5prime", "is_type_iis": True},
    "BsmBI":  {"site": "CGTCTC",  "top_cut": 7,  "bot_cut": 11, "overhang_type": "5prime", "is_type_iis": True},
    "Esp3I":  {"site": "CGTCTC",  "top_cut": 7,  "bot_cut": 11, "overhang_type": "5prime", "is_type_iis": True},
    "BbsI":   {"site": "GAAGAC",  "top_cut": 8,  "bot_cut": 12, "overhang_type": "5prime", "is_type_iis": True},
    "BpiI":   {"site": "GAAGAC",  "top_cut": 8,  "bot_cut": 12, "overhang_type": "5prime", "is_type_iis": True},
    "SapI":   {"site": "GCTCTTC", "top_cut": 8,  "bot_cut": 11, "overhang_type": "5prime", "is_type_iis": True},
    "PaqCI":  {"site": "CACCTGC", "top_cut": 11, "bot_cut": 15, "overhang_type": "5prime", "is_type_iis": True},
    # Type II (5' sticky)
    "EcoRI":  {"site": "GAATTC",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "BamHI":  {"site": "GGATCC",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "BglII":  {"site": "AGATCT",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "HindIII":{"site": "AAGCTT",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "NdeI":   {"site": "CATATG",   "top_cut": 2, "bot_cut": 4, "overhang_type": "5prime", "is_type_iis": False},
    "XhoI":   {"site": "CTCGAG",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "SalI":   {"site": "GTCGAC",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "NotI":   {"site": "GCGGCCGC", "top_cut": 2, "bot_cut": 6, "overhang_type": "5prime", "is_type_iis": False},
    "NheI":   {"site": "GCTAGC",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "XbaI":   {"site": "TCTAGA",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "SpeI":   {"site": "ACTAGT",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "AvrII":  {"site": "CCTAGG",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "AgeI":   {"site": "ACCGGT",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "MluI":   {"site": "ACGCGT",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "NcoI":   {"site": "CCATGG",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "ClaI":   {"site": "ATCGAT",   "top_cut": 2, "bot_cut": 4, "overhang_type": "5prime", "is_type_iis": False},
    "AscI":   {"site": "GGCGCGCC", "top_cut": 2, "bot_cut": 6, "overhang_type": "5prime", "is_type_iis": False},
    "BsiWI":  {"site": "CGTACG",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "AflII":  {"site": "CTTAAG",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    "BstBI":  {"site": "TTCGAA",   "top_cut": 2, "bot_cut": 4, "overhang_type": "5prime", "is_type_iis": False},
    "BclI":   {"site": "TGATCA",   "top_cut": 1, "bot_cut": 5, "overhang_type": "5prime", "is_type_iis": False},
    # Type II (3' sticky)
    "KpnI":   {"site": "GGTACC",   "top_cut": 5, "bot_cut": 1, "overhang_type": "3prime", "is_type_iis": False},
    "SacI":   {"site": "GAGCTC",   "top_cut": 5, "bot_cut": 1, "overhang_type": "3prime", "is_type_iis": False},
    "PstI":   {"site": "CTGCAG",   "top_cut": 5, "bot_cut": 1, "overhang_type": "3prime", "is_type_iis": False},
    "SphI":   {"site": "GCATGC",   "top_cut": 5, "bot_cut": 1, "overhang_type": "3prime", "is_type_iis": False},
    "PacI":   {"site": "TTAATTAA", "top_cut": 5, "bot_cut": 3, "overhang_type": "3prime", "is_type_iis": False},
    "FseI":   {"site": "GGCCGGCC", "top_cut": 6, "bot_cut": 2, "overhang_type": "3prime", "is_type_iis": False},
    "SbfI":   {"site": "CCTGCAGG", "top_cut": 6, "bot_cut": 2, "overhang_type": "3prime", "is_type_iis": False},
    "ApaI":   {"site": "GGGCCC",   "top_cut": 5, "bot_cut": 1, "overhang_type": "3prime", "is_type_iis": False},
    # Type II (blunt)
    "EcoRV":  {"site": "GATATC",   "top_cut": 3, "bot_cut": 3, "overhang_type": "blunt", "is_type_iis": False},
    "SmaI":   {"site": "CCCGGG",   "top_cut": 3, "bot_cut": 3, "overhang_type": "blunt", "is_type_iis": False},
    "PvuII":  {"site": "CAGCTG",   "top_cut": 3, "bot_cut": 3, "overhang_type": "blunt", "is_type_iis": False},
    "ScaI":   {"site": "AGTACT",   "top_cut": 3, "bot_cut": 3, "overhang_type": "blunt", "is_type_iis": False},
    "HpaI":   {"site": "GTTAAC",   "top_cut": 3, "bot_cut": 3, "overhang_type": "blunt", "is_type_iis": False},
    "StuI":   {"site": "AGGCCT",   "top_cut": 3, "bot_cut": 3, "overhang_type": "blunt", "is_type_iis": False},
    "SnaBI":  {"site": "TACGTA",   "top_cut": 3, "bot_cut": 3, "overhang_type": "blunt", "is_type_iis": False},
}


def emit(payload: dict) -> None:
    sys.stdout.write(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()


def fail(message: str, code: int = 1) -> NoReturn:
    sys.stdout.write(json.dumps({"error": message}, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()
    raise SystemExit(code)


def complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    u = seq.upper()
    if not u:
        return 0.0
    return (u.count("G") + u.count("C")) / len(u)


def translate_dna(seq: str, *, stop_symbol: str = "*") -> str:
    u = seq.upper().replace("U", "T")
    aa = []
    for i in range(0, len(u) - 2, 3):
        cod = u[i : i + 3]
        if cod not in CODON_TABLE:
            aa.append("X")
        else:
            res = CODON_TABLE[cod]
            aa.append(stop_symbol if res == "*" else res)
    return "".join(aa)


def resolve_enzyme(name: str) -> tuple[str, dict]:
    if name in RESTRICTION_ENZYMES:
        return name, RESTRICTION_ENZYMES[name]
    for k, v in RESTRICTION_ENZYMES.items():
        if k.lower() == name.strip().lower():
            return k, v
    fail(f"Unsupported restriction enzyme: '{name}'")


def parse_list_arg(raw_items: list[str] | None, default: list[str]) -> list[str]:
    if not raw_items:
        return list(default)
    out: list[str] = []
    for item in raw_items:
        for tok in str(item).split(","):
            tok = tok.strip()
            if tok and tok not in out:
                out.append(tok)
    return out or list(default)


def parse_sequences(path: str | Path) -> list[dict]:
    """Parse FASTA or GenBank (.gb/.gbk) file into a list of sequence dicts."""
    p = Path(path)
    if not p.is_file():
        fail(f"Sequence file not found: {path}")
    text = p.read_text(encoding="utf-8", errors="replace")
    stripped = text.lstrip()
    if stripped.startswith(">"):
        records: list[dict] = []
        cur_id: str | None = None
        chunks: list[str] = []
        for raw in text.splitlines():
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seq_s = re.sub(r"[^A-Za-z*]", "", "".join(chunks)).upper()
                    records.append({"id": cur_id, "seq": seq_s, "topology": "linear", "features": []})
                cur_id = line[1:].split()[0] if line[1:].strip() else f"seq_{len(records) + 1}"
                chunks = []
            else:
                chunks.append(line)
        if cur_id is not None:
            seq_s = re.sub(r"[^A-Za-z*]", "", "".join(chunks)).upper()
            records.append({"id": cur_id, "seq": seq_s, "topology": "linear", "features": []})
        if not records:
            fail(f"No FASTA records found in {path}")
        return records

    if "LOCUS" in stripped[:200]:
        return _parse_genbank_text(text, str(path))

    fail(f"Unrecognized sequence format (expected FASTA or GenBank): {path}")


def _parse_genbank_text(text: str, label: str) -> list[dict]:
    records: list[dict] = []
    blocks = re.split(r"(?m)^//\s*$", text)
    for block in blocks:
        if not block.strip():
            continue
        locus_m = re.search(r"(?m)^LOCUS\s+(\S+)(.*)$", block)
        if not locus_m:
            continue
        rec_id = locus_m.group(1).strip()
        locus_rest = locus_m.group(2).lower()
        topology = "circular" if "circular" in locus_rest else "linear"

        features: list[dict] = []
        feat_m = re.search(r"(?m)^FEATURES\s+Location/Qualifiers\s*\n(.*?)(?=^ORIGIN|\Z)", block, flags=re.DOTALL)
        if feat_m:
            cur_feat: dict | None = None
            for fline in feat_m.group(1).splitlines():
                m_head = re.match(r"^     (\S+)\s+(\S.*)$", fline)
                if m_head:
                    if cur_feat:
                        features.append(cur_feat)
                    cur_feat = {"type": m_head.group(1), "location": m_head.group(2).strip(), "qualifiers": {}}
                elif cur_feat is not None:
                    q_m = re.match(r'^\s+/([^=\s]+)=?(.*)$', fline)
                    if q_m:
                        qk = q_m.group(1)
                        qv = q_m.group(2).strip().strip('"')
                        cur_feat["qualifiers"][qk] = qv
            if cur_feat:
                features.append(cur_feat)

        orig_m = re.search(r"(?m)^ORIGIN[^\n]*\n(.*)$", block, flags=re.DOTALL)
        if not orig_m:
            fail(f"GenBank record '{rec_id}' missing ORIGIN sequence in {label}")
        seq = re.sub(r"[^A-Za-z]", "", orig_m.group(1)).upper()
        if not seq:
            fail(f"GenBank record '{rec_id}' has empty ORIGIN sequence in {label}")
        records.append({"id": rec_id, "seq": seq, "topology": topology, "features": features})
    if not records:
        fail(f"No valid GenBank records found in {label}")
    return records


def write_fasta(path: str | Path, records: list[tuple[str, str]]) -> None:
    lines: list[str] = []
    for rec_id, seq in records:
        lines.append(f">{rec_id}")
        for i in range(0, len(seq), 80):
            lines.append(seq[i : i + 80])
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_genbank(
    path: str | Path,
    locus_name: str,
    seq: str,
    *,
    topology: str = "circular",
    definition: str = "Synthetic construct generated by molecular-cloning plugin.",
    features: list[dict] | None = None,
) -> None:
    clean_name = re.sub(r"[^A-Za-z0-9_.-]", "_", locus_name)[:16] or "construct"
    seq_lower = seq.lower()
    topo_str = "circular" if topology == "circular" else "linear"
    date_str = datetime.now(timezone.utc).strftime("%d-%b-%Y").upper()
    lines = [
        f"LOCUS       {clean_name:<16} {len(seq_lower):>7} bp    DNA     {topo_str:<8} SYN {date_str}",
        f"DEFINITION  {definition}",
        f"ACCESSION   {clean_name}",
        f"VERSION     {clean_name}.1",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{len(seq_lower)}",
        '                     /mol_type="other DNA"',
    ]
    for feat in features or []:
        ftype = feat.get("type", "misc_feature")[:15]
        loc = feat.get("location", f"1..{len(seq_lower)}")
        lines.append(f"     {ftype:<15} {loc}")
        for qk, qv in (feat.get("qualifiers") or {}).items():
            val_s = str(qv).replace('"', "'")
            lines.append(f'                     /{qk}="{val_s}"')

    lines.append("ORIGIN")
    for idx in range(0, len(seq_lower), 60):
        chunk = seq_lower[idx : idx + 60]
        blocks = " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{idx + 1:>9} {blocks}")
    lines.append("//")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def calc_tm_santalucia_owczarzy(
    seq: str,
    *,
    na_mm: float = 50.0,
    mg_mm: float = 1.5,
    dntp_mm: float = 0.2,
    primer_nm: float = 250.0,
) -> dict[str, float]:
    """Compute exact SantaLucia (1998) NN + Owczarzy (2008) salt-corrected Tm."""
    s = seq.upper()
    n = len(s)
    if n < 2:
        return {"tm_c": 0.0, "tm_1m_c": 0.0, "dh": 0.0, "ds": 0.0, "dg37": 0.0}

    dh = 0.0
    ds = 0.0
    for i in range(n - 1):
        di = s[i : i + 2]
        if di not in NN_PARAMS:
            # Ambiguous fallback
            dh += -8.0
            ds += -21.9
        else:
            h_i, s_i = NN_PARAMS[di]
            dh += h_i
            ds += s_i

    for term in (s[0], s[-1]):
        if term in ("G", "C"):
            dh += INIT_GC[0]
            ds += INIT_GC[1]
        else:
            dh += INIT_AT[0]
            ds += INIT_AT[1]

    if s == revcomp(s):
        ds -= 1.4

    ct_m = max(1e-12, primer_nm * 1e-9)
    ct_eff = ct_m if s == revcomp(s) else (ct_m / 4.0)
    denom = ds + R_GAS * math.log(ct_eff)
    tm_1m_k = (1000.0 * dh) / denom if abs(denom) > 1e-9 else 273.15

    mon_m = max(0.0, na_mm * 1e-3)
    mg_free_m = max(0.0, (mg_mm - dntp_mm) * 1e-3)
    f_gc = (s.count("G") + s.count("C")) / n

    if mg_free_m <= 0.0:
        if mon_m <= 0.0:
            inv_tm = 1.0 / tm_1m_k
        else:
            ln_mon = math.log(mon_m)
            inv_tm = (
                (1.0 / tm_1m_k)
                + (4.29 * f_gc - 3.95) * 1e-5 * ln_mon
                + 9.40e-6 * (ln_mon ** 2)
            )
    else:
        ratio = math.inf if mon_m <= 0.0 else math.sqrt(mg_free_m) / mon_m
        if mon_m > 0.0 and ratio < 0.22:
            ln_mon = math.log(mon_m)
            inv_tm = (
                (1.0 / tm_1m_k)
                + (4.29 * f_gc - 3.95) * 1e-5 * ln_mon
                + 9.40e-6 * (ln_mon ** 2)
            )
        else:
            ln_mg = math.log(mg_free_m)
            b = -9.11e-6
            c = 6.26e-5
            e = -4.82e-4
            f = 5.25e-4
            if ratio >= 6.0 or mon_m <= 0.0:
                a = 3.92e-5
                d = 1.42e-5
                g = 8.31e-5
            else:
                ln_mon = math.log(mon_m)
                a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(mon_m) * ln_mon)
                d = 1.42e-5 * (1.279 - 4.03e-3 * ln_mon - 8.03e-3 * (ln_mon ** 2))
                g = 8.31e-5 * (0.486 - 0.258 * ln_mon + 5.25e-3 * (ln_mon ** 3))
            inv_tm = (
                (1.0 / tm_1m_k)
                + a
                + b * ln_mg
                + f_gc * (c + d * ln_mg)
                + (1.0 / (2.0 * max(1, n - 1))) * (e + f * ln_mg + g * (ln_mg ** 2))
            )

    tm_c = (1.0 / inv_tm) - 273.15
    dg37 = dh - 310.15 * (ds / 1000.0)
    return {
        "tm_c": round(tm_c, 2),
        "tm_1m_c": round(tm_1m_k - 273.15, 2),
        "dh": round(dh, 2),
        "ds": round(ds, 2),
        "dg37": round(dg37, 2),
    }


def _wc_pair(a: str, b: str) -> bool:
    return (a == "A" and b == "T") or (a == "T" and b == "A") or (a == "G" and b == "C") or (a == "C" and b == "G")


def _stem_dg37(stem_seq: str) -> float:
    if len(stem_seq) < 2:
        return 0.0
    dh = 0.0
    ds = 0.0
    for i in range(len(stem_seq) - 1):
        di = stem_seq[i : i + 2]
        h_i, s_i = NN_PARAMS.get(di, (-8.0, -21.9))
        dh += h_i
        ds += s_i
    for term in (stem_seq[0], stem_seq[-1]):
        h_t, s_t = INIT_GC if term in ("G", "C") else INIT_AT
        dh += h_t
        ds += s_t
    return dh - 310.15 * (ds / 1000.0)


def dimer_dg37(seq1: str, seq2: str) -> float:
    """Minimum (most stable) ungapped duplex dG37 (kcal/mol) between 5'->3' seq1 and 5'->3' seq2."""
    s1 = seq1.upper()
    s2_rev = seq2.upper()[::-1]  # align 5'->3' s1 against 3'->5' s2
    n1, n2 = len(s1), len(s2_rev)
    best = 0.0
    for shift in range(-(n2 - 2), n1 - 1):
        i_start = max(0, shift)
        j_start = max(0, -shift)
        span = min(n1 - i_start, n2 - j_start)
        run_bases: list[str] = []
        for k in range(span):
            a = s1[i_start + k]
            b = s2_rev[j_start + k]
            if _wc_pair(a, b):
                run_bases.append(a)
            else:
                if len(run_bases) >= 3:
                    dg = _stem_dg37("".join(run_bases))
                    # Extra 3' anchor bonus if match touches 3' end of either oligo
                    if (i_start + k == n1) or (j_start + k - len(run_bases) == 0):
                        dg -= 0.8
                    best = min(best, dg)
                run_bases = []
        if len(run_bases) >= 3:
            dg = _stem_dg37("".join(run_bases))
            if (i_start + span == n1) or (j_start + span - len(run_bases) == 0):
                dg -= 0.8
            best = min(best, dg)
    return round(best, 2)


def hairpin_dg37(seq: str) -> float:
    """Estimate strongest intramolecular hairpin dG37 (kcal/mol) with >=3 bp stem and >=3 nt loop."""
    s = seq.upper()
    n = len(s)
    best = 0.0
    for loop_len in range(3, min(12, n - 5)):
        for i in range(0, n - loop_len - 5):
            left_end = i
            right_start = i + loop_len + 1
            # Extend stem outward
            l_idx = left_end
            r_idx = right_start
            stem = []
            while l_idx >= 0 and r_idx < n and _wc_pair(s[l_idx], s[r_idx]):
                stem.append(s[l_idx])
                l_idx -= 1
                r_idx += 1
            if len(stem) >= 3:
                stem_s = "".join(reversed(stem))
                loop_penalty = 3.5 + 0.4 * max(0, loop_len - 4)
                dg = _stem_dg37(stem_s) + loop_penalty
                best = min(best, dg)
    return round(best, 2)


def max_homopolymer(seq: str) -> int:
    s = seq.upper()
    if not s:
        return 0
    best = 1
    cur = 1
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
    return best


def find_restriction_sites(seq: str, enzymes: list[str]) -> list[dict]:
    """Find all forward and reverse-complement recognition sites for given enzymes."""
    u = seq.upper()
    hits: list[dict] = []
    for raw_enz in enzymes:
        enz, spec = resolve_enzyme(raw_enz)
        site = spec["site"]
        rc_site = revcomp(site)
        start = 0
        while True:
            idx = u.find(site, start)
            if idx == -1:
                break
            hits.append({"enzyme": enz, "pos_0": idx, "strand": "+", "site": site})
            start = idx + 1
        if rc_site != site:
            start = 0
            while True:
                idx = u.find(rc_site, start)
                if idx == -1:
                    break
                hits.append({"enzyme": enz, "pos_0": idx, "strand": "-", "site": rc_site})
                start = idx + 1
    return hits


def parse_ab1(path: str | Path) -> dict:
    """Pure-Python binary ABIF (.ab1) chromatogram parser."""
    p = Path(path)
    data = p.read_bytes()
    if len(data) < 128 or data[:4] != b"ABIF":
        fail(f"Invalid ABIF .ab1 file (missing ABIF magic header): {path}")

    version = struct.unpack(">H", data[4:6])[0]
    # Root DirEntry starts at byte 6 (28 bytes)
    _, _, _, _, num_elems, _, dir_offset, _ = struct.unpack(">4sIHHIiIi", data[6:34])
    if dir_offset <= 0 or dir_offset + num_elems * 28 > len(data):
        fail(f"Corrupt ABIF directory offset in {path}")

    tags: dict[tuple[str, int], tuple[int, int, int, bytes]] = {}
    for idx in range(num_elems):
        off = dir_offset + idx * 28
        raw_tag, tag_num, elem_type, elem_size, n_elems, d_size, d_off, _ = struct.unpack(
            ">4sIHHIiIi", data[off : off + 28]
        )
        tag_name = raw_tag.decode("ascii", errors="replace")
        if d_size <= 4:
            payload = data[off + 20 : off + 20 + d_size]
        else:
            if d_off < 0 or d_off + d_size > len(data):
                continue
            payload = data[d_off : d_off + d_size]
        tags[(tag_name, tag_num)] = (elem_type, elem_size, n_elems, payload)

    fwo = b"GATC"
    if ("FWO_", 1) in tags:
        fwo = tags[("FWO_", 1)][3][:4]
    channel_bases = [chr(b).upper() for b in fwo]

    pbas_entry = tags.get(("PBAS", 2)) or tags.get(("PBAS", 1))
    if not pbas_entry:
        fail(f"ABIF file missing PBAS basecalls: {path}")
    seq = pbas_entry[3].decode("ascii", errors="replace").upper()

    pcon_entry = tags.get(("PCON", 2)) or tags.get(("PCON", 1))
    if pcon_entry:
        quals = list(pcon_entry[3])
    else:
        quals = [30] * len(seq)

    ploc_entry = tags.get(("PLOC", 2)) or tags.get(("PLOC", 1))
    ploc: list[int] = []
    if ploc_entry:
        n_loc = len(ploc_entry[3]) // 2
        ploc = list(struct.unpack(f">{n_loc}H", ploc_entry[3][: n_loc * 2]))

    traces: dict[str, list[int]] = {}
    for ch_idx, tag_num in enumerate((9, 10, 11, 12)):
        base = channel_bases[ch_idx] if ch_idx < len(channel_bases) else "ACGT"[ch_idx]
        entry = tags.get(("DATA", tag_num))
        if entry:
            elem_type, _, _, payload = entry
            n_pts = len(payload) // 2
            fmt = f">{n_pts}h" if elem_type == 4 else f">{n_pts}H"
            traces[base] = list(struct.unpack(fmt, payload[: n_pts * 2]))

    return {
        "name": p.stem,
        "format": "ab1",
        "abif_version": version,
        "seq": seq,
        "quals": quals[: len(seq)],
        "ploc": ploc[: len(seq)],
        "channel_order": "".join(channel_bases),
        "traces": traces,
    }


def parse_trace_file(path: str | Path) -> dict:
    """Read .ab1, .fastq, or .fasta Sanger trace into normalized dict."""
    p = Path(path)
    if not p.is_file():
        fail(f"Sanger trace file not found: {path}")
    suf = p.suffix.lower()
    if suf == ".ab1":
        return parse_ab1(p)
    text = p.read_text(encoding="utf-8", errors="replace")
    if suf in (".fastq", ".fq") or text.lstrip().startswith("@"):
        lines = [ln.rstrip("\r\n") for ln in text.splitlines() if ln.strip()]
        if len(lines) < 4 or not lines[0].startswith("@"):
            fail(f"Invalid FASTQ trace file: {path}")
        name = lines[0][1:].split()[0] or p.stem
        seq = lines[1].strip().upper()
        qstr = lines[3].strip()
        quals = [max(0, ord(ch) - 33) for ch in qstr[: len(seq)]]
        if len(quals) < len(seq):
            quals.extend([25] * (len(seq) - len(quals)))
        return {
            "name": name,
            "format": "fastq",
            "seq": seq,
            "quals": quals,
            "ploc": [],
            "channel_order": "GATC",
            "traces": {},
        }
    recs = parse_sequences(p)
    seq = recs[0]["seq"]
    return {
        "name": recs[0]["id"],
        "format": "fasta",
        "seq": seq,
        "quals": [35] * len(seq),
        "ploc": [],
        "channel_order": "GATC",
        "traces": {},
    }


def mott_trim(quals: list[int], error_prob: float = 0.05) -> tuple[int, int]:
    """Return 0-based [start, end) high-quality slice via Mott's trimming algorithm."""
    n = len(quals)
    if n == 0:
        return 0, 0
    best_score = -1.0
    best_end = 0
    best_start = 0
    cur_score = 0.0
    cur_start = 0
    for i, q in enumerate(quals):
        p_err = 10.0 ** (-max(0, min(93, q)) / 10.0)
        cur_score += error_prob - p_err
        if cur_score <= 0.0:
            cur_score = 0.0
            cur_start = i + 1
        elif cur_score > best_score:
            best_score = cur_score
            best_start = cur_start
            best_end = i + 1
    if best_end <= best_start:
        return 0, 0
    return best_start, best_end


def setup_matplotlib():
    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def write_artifact(
    path: str | Path,
    tool: str,
    title: str,
    *,
    skill: str = "cloning-workflow",
    inputs=(),
    parameters: dict | None = None,
    shape: tuple[int, int] | dict | None = None,
    scientific_context: dict | None = None,
) -> str:
    p = Path(path)
    ext = p.suffix.lower()
    figures = {".png": "image/png", ".pdf": "application/pdf", ".svg": "image/svg+xml"}
    tables = {".tsv": "text/tab-separated-values", ".csv": "text/csv"}
    seqs = {".fasta": "text/plain", ".fa": "text/plain", ".gb": "text/plain", ".gbk": "text/plain"}
    if ext in figures:
        kind, media = "figure", figures[ext]
    elif ext in tables:
        kind, media = "table", tables[ext]
    elif ext in seqs:
        kind, media = "sequence", seqs[ext]
    else:
        raise ValueError(f"Unsupported artifact extension: {ext}")

    descriptor: dict = {
        "contractVersion": "1.0.0",
        "kind": kind,
        "file": p.name,
        "mediaType": media,
        "title": title[:200],
        "provenance": {
            "createdAt": datetime.now(timezone.utc).isoformat(),
            "tool": tool,
            "skill": skill,
        },
    }
    if kind == "figure":
        block: dict = {"format": ext[1:]}
        if ext == ".png":
            with open(p, "rb") as fh:
                header = fh.read(24)
            if header[:8] == b"\x89PNG\r\n\x1a\n":
                w, h = struct.unpack(">II", header[16:24])
                block["widthPx"], block["heightPx"] = int(w), int(h)
        descriptor["figure"] = block
    elif kind == "table" and shape is not None:
        if isinstance(shape, dict):
            rows_n = int(shape.get("rows", 0))
            cols_n = int(shape.get("columns", shape.get("cols", 0)))
        else:
            rows_n, cols_n = int(shape[0]), int(shape[1])
        descriptor["table"] = {"rows": rows_n, "columns": cols_n}

    if scientific_context:
        descriptor["scientificContext"] = scientific_context

    cwd = Path.cwd().resolve()
    prov_inputs = []
    for ip in inputs:
        if not ip:
            continue
        ip_p = Path(ip)
        try:
            rel = ip_p.resolve().relative_to(cwd)
        except ValueError:
            continue
        entry = {"path": str(rel)}
        if ip_p.is_file():
            entry["sha256"] = hashlib.sha256(ip_p.read_bytes()).hexdigest()
        prov_inputs.append(entry)
    if prov_inputs:
        descriptor["provenance"]["inputs"] = prov_inputs
    if parameters:
        descriptor["provenance"]["parameters"] = parameters

    target = p.with_name(p.name + ".phi-artifact.json")
    tmp = target.with_name(target.name + ".tmp")
    tmp.write_text(json.dumps(descriptor, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, target)
    return str(p)
