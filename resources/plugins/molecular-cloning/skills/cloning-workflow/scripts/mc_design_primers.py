#!/usr/bin/env python3
"""Design PCR, qPCR, Gibson, GoldenGate, and SiteDirectedMutagenesis primers."""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    calc_tm_santalucia_owczarzy,
    dimer_dg37,
    emit,
    fail,
    find_restriction_sites,
    gc_fraction,
    hairpin_dg37,
    max_homopolymer,
    parse_sequences,
    resolve_enzyme,
    revcomp,
    setup_matplotlib,
    write_artifact,
)


def _load_background_12mers(bg_path: str | None) -> set[str]:
    if not bg_path:
        return set()
    recs = parse_sequences(bg_path)
    kmers: set[str] = set()
    for rec in recs:
        for seq in (rec["seq"], revcomp(rec["seq"])):
            for i in range(0, max(0, len(seq) - 11)):
                kmers.add(seq[i : i + 12])
    return kmers


def _eval_primer(
    core: str,
    full_seq: str,
    *,
    na_mm: float,
    mg_mm: float,
    dntp_mm: float,
    primer_nm: float,
    bg_12mers: set[str],
) -> dict:
    ann = calc_tm_santalucia_owczarzy(
        core, na_mm=na_mm, mg_mm=mg_mm, dntp_mm=dntp_mm, primer_nm=primer_nm
    )
    full = calc_tm_santalucia_owczarzy(
        full_seq, na_mm=na_mm, mg_mm=mg_mm, dntp_mm=dntp_mm, primer_nm=primer_nm
    )
    gc_pct = round(100.0 * gc_fraction(core), 1)
    clamp = core[-5:].count("G") + core[-5:].count("C")
    homo = max_homopolymer(core)
    hp_dg = hairpin_dg37(full_seq)
    sd_dg = dimer_dg37(full_seq, full_seq)
    misprime = 1 if (bg_12mers and len(core) >= 12 and core[-12:] in bg_12mers) else 0
    return {
        "core_seq": core,
        "full_seq": full_seq,
        "core_len": len(core),
        "full_len": len(full_seq),
        "annealing_tm": ann["tm_c"],
        "full_tm": full["tm_c"],
        "gc_pct": gc_pct,
        "gc_clamp_3p": clamp,
        "max_homopolymer": homo,
        "hairpin_dg37": hp_dg,
        "self_dimer_dg37": sd_dg,
        "mispriming_3p_12mer": misprime,
    }


def _score_pair(
    fwd: dict,
    rev: dict,
    *,
    tm_opt: float,
    cross_dg: float,
    internal_gg_hits: int = 0,
) -> float:
    score = 100.0
    for p in (fwd, rev):
        score -= 2.2 * abs(p["annealing_tm"] - tm_opt)
        if p["gc_pct"] < 40.0:
            score -= 1.0 * (40.0 - p["gc_pct"])
        elif p["gc_pct"] > 60.0:
            score -= 1.0 * (p["gc_pct"] - 60.0)
        if p["gc_clamp_3p"] == 0:
            score -= 5.0
        elif p["gc_clamp_3p"] == 3:
            score -= 3.0
        elif p["gc_clamp_3p"] >= 4:
            score -= 12.0
        if p["max_homopolymer"] >= 4:
            score -= 6.0 * (p["max_homopolymer"] - 3)
        if p["hairpin_dg37"] < -3.0:
            score -= 2.5 * (-3.0 - p["hairpin_dg37"])
        if p["self_dimer_dg37"] < -5.0:
            score -= 2.0 * (-5.0 - p["self_dimer_dg37"])
        if p["mispriming_3p_12mer"] > 0:
            score -= 25.0
    dtm = abs(fwd["annealing_tm"] - rev["annealing_tm"])
    if dtm > 2.5:
        score -= 5.0 * (dtm - 2.5)
    else:
        score -= 0.8 * dtm
    if cross_dg < -5.0:
        score -= 2.5 * (-5.0 - cross_dg)
    if internal_gg_hits > 0:
        score -= 20.0 * internal_gg_hits
    return round(max(0.0, min(100.0, score)), 2)


def _gg_adapter(enzyme: str, overhang: str) -> str:
    enz, spec = resolve_enzyme(enzyme)
    site = spec["site"]
    expected_len = spec["bot_cut"] - spec["top_cut"]
    oh = overhang.upper() if overhang else ("AATG" if expected_len == 4 else "ATG")
    if len(oh) != expected_len or not re.fullmatch(r"[ACGT]+", oh):
        fail(f"GoldenGate overhang '{oh}' must be {expected_len} ACGT bases for {enz}")
    spacer_len = spec["top_cut"] - len(site)
    spacer = "A" * max(1, spacer_len)
    return site + spacer + oh


def _pick_primers_pcr_qpcr(
    seq: str,
    *,
    mode: str,
    t_start: int,
    t_end: int,
    amp_min: int,
    amp_max: int,
    tm_opt: float,
    na_mm: float,
    mg_mm: float,
    dntp_mm: float,
    primer_nm: float,
    bg_12mers: set[str],
) -> list[dict]:
    L = len(seq)
    # Candidate forward windows (0-based start of fwd primer)
    fwd_starts = list(range(max(0, t_start - 150), min(L - 25, t_start + 1)))
    if t_end - t_start + 1 <= amp_min and len(fwd_starts) < 5:
        fwd_starts = list(range(0, max(1, L - amp_min)))
    # Subsample if very large template with omitted target_start
    if len(fwd_starts) > 80:
        step = max(1, len(fwd_starts) // 60)
        fwd_starts = fwd_starts[::step]

    fwd_cands: list[dict] = []
    for s_idx in fwd_starts:
        best_f = None
        best_diff = 999.0
        for plen in range(18, 29):
            if s_idx + plen > L:
                break
            core = seq[s_idx : s_idx + plen]
            if not re.fullmatch(r"[ACGT]+", core):
                continue
            ev = _eval_primer(
                core,
                core,
                na_mm=na_mm,
                mg_mm=mg_mm,
                dntp_mm=dntp_mm,
                primer_nm=primer_nm,
                bg_12mers=bg_12mers,
            )
            diff = abs(ev["annealing_tm"] - tm_opt)
            if diff < best_diff:
                best_diff = diff
                best_f = {"start_0": s_idx, "end_0": s_idx + plen, **ev}
        if best_f:
            fwd_cands.append(best_f)

    pairs: list[dict] = []
    for fc in fwd_cands:
        min_end = max(t_end, fc["start_0"] + amp_min)
        max_end = min(L, fc["start_0"] + amp_max)
        if min_end > max_end:
            continue
        rev_ends = list(range(min_end, max_end + 1))
        if len(rev_ends) > 45:
            step = max(1, len(rev_ends) // 35)
            rev_ends = rev_ends[::step]
        for r_end in rev_ends:
            best_r = None
            best_diff = 999.0
            for plen in range(18, 29):
                r_start = r_end - plen
                if r_start <= fc["end_0"]:
                    continue
                core = revcomp(seq[r_start:r_end])
                if not re.fullmatch(r"[ACGT]+", core):
                    continue
                ev = _eval_primer(
                    core,
                    core,
                    na_mm=na_mm,
                    mg_mm=mg_mm,
                    dntp_mm=dntp_mm,
                    primer_nm=primer_nm,
                    bg_12mers=bg_12mers,
                )
                diff = abs(ev["annealing_tm"] - tm_opt)
                if diff < best_diff:
                    best_diff = diff
                    best_r = {"start_0": r_start, "end_0": r_end, **ev}
            if not best_r:
                continue
            cross_dg = dimer_dg37(fc["full_seq"], best_r["full_seq"])
            sc = _score_pair(fc, best_r, tm_opt=tm_opt, cross_dg=cross_dg)
            pairs.append({
                "design_type": mode,
                "amplicon_start_1": fc["start_0"] + 1,
                "amplicon_end_1": best_r["end_0"],
                "amplicon_len": best_r["end_0"] - fc["start_0"],
                "fwd_core": fc["core_seq"],
                "fwd_full": fc["full_seq"],
                "fwd_len": fc["full_len"],
                "fwd_annealing_tm": fc["annealing_tm"],
                "fwd_full_tm": fc["full_tm"],
                "fwd_gc_pct": fc["gc_pct"],
                "fwd_gc_clamp_3p": fc["gc_clamp_3p"],
                "fwd_hairpin_dg37": fc["hairpin_dg37"],
                "fwd_self_dimer_dg37": fc["self_dimer_dg37"],
                "rev_core": best_r["core_seq"],
                "rev_full": best_r["full_seq"],
                "rev_len": best_r["full_len"],
                "rev_annealing_tm": best_r["annealing_tm"],
                "rev_full_tm": best_r["full_tm"],
                "rev_gc_pct": best_r["gc_pct"],
                "rev_gc_clamp_3p": best_r["gc_clamp_3p"],
                "rev_hairpin_dg37": best_r["hairpin_dg37"],
                "rev_self_dimer_dg37": best_r["self_dimer_dg37"],
                "delta_tm": round(abs(fc["annealing_tm"] - best_r["annealing_tm"]), 2),
                "cross_dimer_dg37": cross_dg,
                "mispriming_hits": fc["mispriming_3p_12mer"] + best_r["mispriming_3p_12mer"],
                "score": sc,
            })
    pairs.sort(key=lambda r: (-r["score"], r["delta_tm"], abs(r["fwd_annealing_tm"] - tm_opt)))
    return pairs[:25]


def _pick_cloning_primers(
    seq: str,
    *,
    mode: str,
    t_start: int,
    t_end: int,
    left_tail: str,
    right_tail_rc: str,
    gg_enzyme: str,
    tm_opt: float,
    na_mm: float,
    mg_mm: float,
    dntp_mm: float,
    primer_nm: float,
    bg_12mers: set[str],
) -> list[dict]:
    if t_end - t_start < 36:
        fail("Target region must be at least 36 bp for cloning primer design.")
    fwd_options: list[dict] = []
    for plen in range(18, 31):
        core = seq[t_start : t_start + plen]
        if not re.fullmatch(r"[ACGT]+", core):
            continue
        full = left_tail + core
        ev = _eval_primer(
            core, full, na_mm=na_mm, mg_mm=mg_mm, dntp_mm=dntp_mm, primer_nm=primer_nm, bg_12mers=bg_12mers
        )
        fwd_options.append({"start_0": t_start, "end_0": t_start + plen, **ev})

    rev_options: list[dict] = []
    for plen in range(18, 31):
        core = revcomp(seq[t_end - plen : t_end])
        if not re.fullmatch(r"[ACGT]+", core):
            continue
        full = right_tail_rc + core
        ev = _eval_primer(
            core, full, na_mm=na_mm, mg_mm=mg_mm, dntp_mm=dntp_mm, primer_nm=primer_nm, bg_12mers=bg_12mers
        )
        rev_options.append({"start_0": t_end - plen, "end_0": t_end, **ev})

    pairs: list[dict] = []
    for fc in fwd_options:
        for rc in rev_options:
            cross_dg = dimer_dg37(fc["full_seq"], rc["full_seq"])
            gg_hits = 0
            if mode == "GoldenGate":
                gg_hits = len(find_restriction_sites(fc["core_seq"], [gg_enzyme])) + len(
                    find_restriction_sites(rc["core_seq"], [gg_enzyme])
                )
            sc = _score_pair(fc, rc, tm_opt=tm_opt, cross_dg=cross_dg, internal_gg_hits=gg_hits)
            pairs.append({
                "design_type": mode,
                "amplicon_start_1": t_start + 1,
                "amplicon_end_1": t_end,
                "amplicon_len": t_end - t_start,
                "fwd_core": fc["core_seq"],
                "fwd_full": fc["full_seq"],
                "fwd_len": fc["full_len"],
                "fwd_annealing_tm": fc["annealing_tm"],
                "fwd_full_tm": fc["full_tm"],
                "fwd_gc_pct": fc["gc_pct"],
                "fwd_gc_clamp_3p": fc["gc_clamp_3p"],
                "fwd_hairpin_dg37": fc["hairpin_dg37"],
                "fwd_self_dimer_dg37": fc["self_dimer_dg37"],
                "rev_core": rc["core_seq"],
                "rev_full": rc["full_seq"],
                "rev_len": rc["full_len"],
                "rev_annealing_tm": rc["annealing_tm"],
                "rev_full_tm": rc["full_tm"],
                "rev_gc_pct": rc["gc_pct"],
                "rev_gc_clamp_3p": rc["gc_clamp_3p"],
                "rev_hairpin_dg37": rc["hairpin_dg37"],
                "rev_self_dimer_dg37": rc["self_dimer_dg37"],
                "delta_tm": round(abs(fc["annealing_tm"] - rc["annealing_tm"]), 2),
                "cross_dimer_dg37": cross_dg,
                "mispriming_hits": fc["mispriming_3p_12mer"] + rc["mispriming_3p_12mer"],
                "score": sc,
            })
    pairs.sort(key=lambda r: (-r["score"], r["delta_tm"]))
    return pairs[:20]


def _pick_sdm_primers(
    seq: str,
    *,
    mut_pos_1: int,
    mut_seq: str,
    tm_opt: float,
    na_mm: float,
    mg_mm: float,
    dntp_mm: float,
    primer_nm: float,
    bg_12mers: set[str],
) -> list[dict]:
    L = len(seq)
    m0 = mut_pos_1 - 1
    m_len = len(mut_seq)
    if m0 < 20 or m0 + m_len + 20 > L:
        fail(f"mutation_pos={mut_pos_1} must leave >=20 bp flanking sequence on both sides (template len={L}).")

    pairs: list[dict] = []
    # 1. Q5 back-to-back non-overlapping SDM primers
    f_cands = []
    for arm_len in range(18, 29):
        if m0 + m_len + arm_len > L:
            break
        core = seq[m0 + m_len : m0 + m_len + arm_len]
        full = mut_seq + core
        ev = _eval_primer(
            core, full, na_mm=na_mm, mg_mm=mg_mm, dntp_mm=dntp_mm, primer_nm=primer_nm, bg_12mers=bg_12mers
        )
        f_cands.append(ev)
    r_cands = []
    for arm_len in range(18, 29):
        if m0 - arm_len < 0:
            break
        core = revcomp(seq[m0 - arm_len : m0])
        ev = _eval_primer(
            core, core, na_mm=na_mm, mg_mm=mg_mm, dntp_mm=dntp_mm, primer_nm=primer_nm, bg_12mers=bg_12mers
        )
        r_cands.append(ev)

    for fc in f_cands:
        for rc in r_cands:
            cross_dg = dimer_dg37(fc["full_seq"], rc["full_seq"])
            sc = _score_pair(fc, rc, tm_opt=tm_opt, cross_dg=cross_dg)
            pairs.append({
                "design_type": "SDM_BackToBack",
                "amplicon_start_1": 1,
                "amplicon_end_1": L,
                "amplicon_len": L,
                "fwd_core": fc["core_seq"],
                "fwd_full": fc["full_seq"],
                "fwd_len": fc["full_len"],
                "fwd_annealing_tm": fc["annealing_tm"],
                "fwd_full_tm": fc["full_tm"],
                "fwd_gc_pct": fc["gc_pct"],
                "fwd_gc_clamp_3p": fc["gc_clamp_3p"],
                "fwd_hairpin_dg37": fc["hairpin_dg37"],
                "fwd_self_dimer_dg37": fc["self_dimer_dg37"],
                "rev_core": rc["core_seq"],
                "rev_full": rc["full_seq"],
                "rev_len": rc["full_len"],
                "rev_annealing_tm": rc["annealing_tm"],
                "rev_full_tm": rc["full_tm"],
                "rev_gc_pct": rc["gc_pct"],
                "rev_gc_clamp_3p": rc["gc_clamp_3p"],
                "rev_hairpin_dg37": rc["hairpin_dg37"],
                "rev_self_dimer_dg37": rc["self_dimer_dg37"],
                "delta_tm": round(abs(fc["annealing_tm"] - rc["annealing_tm"]), 2),
                "cross_dimer_dg37": cross_dg,
                "mispriming_hits": fc["mispriming_3p_12mer"] + rc["mispriming_3p_12mer"],
                "score": sc,
            })

    # 2. QuikChange overlapping SDM primers
    for arm_len in (12, 15, 18):
        if m0 - arm_len < 0 or m0 + m_len + arm_len > L:
            continue
        left_arm = seq[m0 - arm_len : m0]
        right_arm = seq[m0 + m_len : m0 + m_len + arm_len]
        fwd_full = left_arm + mut_seq + right_arm
        rev_full = revcomp(fwd_full)
        fev = _eval_primer(
            right_arm, fwd_full, na_mm=na_mm, mg_mm=mg_mm, dntp_mm=dntp_mm, primer_nm=primer_nm, bg_12mers=bg_12mers
        )
        rev = _eval_primer(
            revcomp(left_arm), rev_full, na_mm=na_mm, mg_mm=mg_mm, dntp_mm=dntp_mm, primer_nm=primer_nm, bg_12mers=bg_12mers
        )
        cross_dg = dimer_dg37(fwd_full, rev_full)
        sc = round(max(0.0, _score_pair(fev, rev, tm_opt=tm_opt, cross_dg=cross_dg) - 5.0), 2)
        pairs.append({
            "design_type": "SDM_Overlapping",
            "amplicon_start_1": m0 - arm_len + 1,
            "amplicon_end_1": m0 + m_len + arm_len,
            "amplicon_len": L,
            "fwd_core": fev["core_seq"],
            "fwd_full": fwd_full,
            "fwd_len": len(fwd_full),
            "fwd_annealing_tm": fev["annealing_tm"],
            "fwd_full_tm": fev["full_tm"],
            "fwd_gc_pct": fev["gc_pct"],
            "fwd_gc_clamp_3p": fev["gc_clamp_3p"],
            "fwd_hairpin_dg37": fev["hairpin_dg37"],
            "fwd_self_dimer_dg37": fev["self_dimer_dg37"],
            "rev_core": rev["core_seq"],
            "rev_full": rev_full,
            "rev_len": len(rev_full),
            "rev_annealing_tm": rev["annealing_tm"],
            "rev_full_tm": rev["full_tm"],
            "rev_gc_pct": rev["gc_pct"],
            "rev_gc_clamp_3p": rev["gc_clamp_3p"],
            "rev_hairpin_dg37": rev["hairpin_dg37"],
            "rev_self_dimer_dg37": rev["self_dimer_dg37"],
            "delta_tm": round(abs(fev["annealing_tm"] - rev["annealing_tm"]), 2),
            "cross_dimer_dg37": cross_dg,
            "mispriming_hits": fev["mispriming_3p_12mer"] + rev["mispriming_3p_12mer"],
            "score": sc,
        })

    pairs.sort(key=lambda r: (-r["score"], r["delta_tm"]))
    return pairs[:20]


def _plot_primer_map(
    out_png: Path, seq_len: int, tm_opt: float, df: pd.DataFrame, mode: str, tmpl_id: str
) -> None:
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.5, 6.2), gridspec_kw={"height_ratios": [1.1, 1.3]})

    ax1.axhline(0, color="#334155", linewidth=3.0)
    ax1.set_xlim(0, max(1, seq_len))
    ax1.set_ylim(-1.2, max(2.2, min(8, len(df)) + 1.2))
    top_k = df.head(6)
    for idx, row in enumerate(top_k.itertuples(index=False)):
        y = 0.7 + idx * 0.65
        s = int(row.amplicon_start_1)
        e = int(row.amplicon_end_1)
        color = "#0f766e" if idx == 0 else "#64748b"
        lw = 2.6 if idx == 0 else 1.4
        ax1.plot([s, e], [y, y], color=color, linewidth=lw)
        ax1.plot(s, y, marker=">", color="#2563eb", markersize=6)
        ax1.plot(e, y, marker="<", color="#dc2626", markersize=6)
        ax1.text(
            min(seq_len * 0.98, e + seq_len * 0.01),
            y,
            f"#{idx + 1} ({row.amplicon_len} bp, Tm {row.fwd_annealing_tm}/{row.rev_annealing_tm}°C)",
            va="center",
            fontsize=8,
            color="#0f172a",
        )
    ax1.set_yticks([])
    ax1.set_xlabel("Template coordinate (bp)")
    ax1.set_title(f"{tmpl_id} — Top Primer Pairs ({mode})")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_visible(False)

    x_idx = list(range(1, len(df) + 1))
    ax2.plot(x_idx, df["fwd_annealing_tm"].to_numpy(copy=True), "o-", color="#2563eb", label="Forward annealing Tm")
    ax2.plot(x_idx, df["rev_annealing_tm"].to_numpy(copy=True), "s-", color="#dc2626", label="Reverse annealing Tm")
    ax2.axhline(tm_opt, color="#0f766e", linestyle="--", linewidth=1.2, label=f"Target Tm ({tm_opt:.1f}°C)")
    ax2.set_xlabel("Candidate Pair Rank")
    ax2.set_ylabel("Tm (°C)")
    ax2.set_title("SantaLucia (1998) + Owczarzy (2008) Primer Tm Profile")
    ax2.legend(frameon=False, fontsize=8.5)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Design molecular cloning primers.")
    parser.add_argument("--template", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument(
        "--mode",
        default="PCR",
        choices=["PCR", "qPCR", "Gibson", "GoldenGate", "SiteDirectedMutagenesis"],
    )
    parser.add_argument("--target_start", type=int, default=None)
    parser.add_argument("--target_end", type=int, default=None)
    parser.add_argument("--amplicon_min", type=int, default=None)
    parser.add_argument("--amplicon_max", type=int, default=None)
    parser.add_argument("--tm_opt", type=float, default=60.0)
    parser.add_argument("--na_mm", type=float, default=50.0)
    parser.add_argument("--mg_mm", type=float, default=1.5)
    parser.add_argument("--dntp_mm", type=float, default=0.2)
    parser.add_argument("--primer_nm", type=float, default=250.0)
    parser.add_argument("--left_homology", default="")
    parser.add_argument("--right_homology", default="")
    parser.add_argument("--goldengate_enzyme", default="BsaI", choices=["BsaI", "BsmBI", "BbsI", "SapI"])
    parser.add_argument("--overhang_5p", default="AATG")
    parser.add_argument("--overhang_3p", default="AGGT")
    parser.add_argument("--mutation_pos", type=int, default=None)
    parser.add_argument("--mutation_seq", default="")
    parser.add_argument("--background_fasta", default=None)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records = parse_sequences(args.template)
    tmpl = records[0]
    seq = tmpl["seq"]
    L = len(seq)
    if L < 50:
        fail(f"Template '{tmpl['id']}' is too short ({L} bp; need >=50 bp).")

    amp_min = args.amplicon_min if args.amplicon_min is not None else (80 if args.mode == "qPCR" else 200)
    amp_max = args.amplicon_max if args.amplicon_max is not None else (180 if args.mode == "qPCR" else 1200)

    t_start = (args.target_start - 1) if args.target_start is not None else 0
    t_end = args.target_end if args.target_end is not None else L
    if t_start < 0 or t_end > L or t_start >= t_end:
        fail(f"Invalid target coordinates [{args.target_start}, {args.target_end}] for length {L} bp.")

    bg_12mers = _load_background_12mers(args.background_fasta)
    warnings: list[str] = []

    if args.mode in ("PCR", "qPCR"):
        if args.target_start is None and args.target_end is None:
            t_start = 0
            t_end = min(L, amp_min)
        elif (t_end - t_start) > amp_max:
            amp_max = t_end - t_start + 300
            warnings.append(f"Expanded amplicon_max to {amp_max} bp to enclose requested target span.")
        pairs = _pick_primers_pcr_qpcr(
            seq,
            mode=args.mode,
            t_start=t_start,
            t_end=t_end,
            amp_min=amp_min,
            amp_max=amp_max,
            tm_opt=args.tm_opt,
            na_mm=args.na_mm,
            mg_mm=args.mg_mm,
            dntp_mm=args.dntp_mm,
            primer_nm=args.primer_nm,
            bg_12mers=bg_12mers,
        )
    elif args.mode == "Gibson":
        left_tail = re.sub(r"[^ACGTacgt]", "", args.left_homology).upper()
        right_tail = re.sub(r"[^ACGTacgt]", "", args.right_homology).upper()
        if not left_tail or not right_tail:
            warnings.append("left_homology or right_homology omitted; designing core Gibson primers without vector tails.")
        pairs = _pick_cloning_primers(
            seq,
            mode="Gibson",
            t_start=t_start,
            t_end=t_end,
            left_tail=left_tail,
            right_tail_rc=revcomp(right_tail) if right_tail else "",
            gg_enzyme=args.goldengate_enzyme,
            tm_opt=args.tm_opt,
            na_mm=args.na_mm,
            mg_mm=args.mg_mm,
            dntp_mm=args.dntp_mm,
            primer_nm=args.primer_nm,
            bg_12mers=bg_12mers,
        )
    elif args.mode == "GoldenGate":
        left_tail = _gg_adapter(args.goldengate_enzyme, args.overhang_5p)
        right_tail_rc = _gg_adapter(args.goldengate_enzyme, revcomp(args.overhang_3p.upper()))
        pairs = _pick_cloning_primers(
            seq,
            mode="GoldenGate",
            t_start=t_start,
            t_end=t_end,
            left_tail=left_tail,
            right_tail_rc=right_tail_rc,
            gg_enzyme=args.goldengate_enzyme,
            tm_opt=args.tm_opt,
            na_mm=na_mm if False else args.na_mm,
            mg_mm=args.mg_mm,
            dntp_mm=args.dntp_mm,
            primer_nm=args.primer_nm,
            bg_12mers=bg_12mers,
        )
        # Check internal GoldenGate sites inside amplicon
        amp_sites = find_restriction_sites(seq[t_start:t_end], [args.goldengate_enzyme])
        if amp_sites:
            warnings.append(
                f"Template target span contains {len(amp_sites)} internal {args.goldengate_enzyme} site(s); domesticate before Golden Gate cloning."
            )
    else:  # SiteDirectedMutagenesis
        if args.mutation_pos is None or not args.mutation_seq:
            fail("SiteDirectedMutagenesis mode requires --mutation_pos and --mutation_seq.")
        mut_seq = args.mutation_seq.strip().upper()
        if not re.fullmatch(r"[ACGT]+", mut_seq):
            fail(f"mutation_seq must contain only ACGT bases, got '{args.mutation_seq}'.")
        pairs = _pick_sdm_primers(
            seq,
            mut_pos_1=args.mutation_pos,
            mut_seq=mut_seq,
            tm_opt=args.tm_opt,
            na_mm=args.na_mm,
            mg_mm=args.mg_mm,
            dntp_mm=args.dntp_mm,
            primer_nm=args.primer_nm,
            bg_12mers=bg_12mers,
        )

    if not pairs:
        fail("No valid primer pairs satisfied the template coordinates and amplicon length bounds.")

    df = pd.DataFrame(pairs)
    df.insert(0, "rank", range(1, len(df) + 1))
    tsv_path = out_dir / "primer_candidates.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)

    png_path = out_dir / "primer_design_map.png"
    _plot_primer_map(png_path, L, args.tm_opt, df, args.mode, tmpl["id"])

    top = df.iloc[0].to_dict()
    if float(top["delta_tm"]) > 2.5:
        warnings.append(f"Top primer pair |delta_Tm| is {top['delta_tm']} °C (> 2.5 °C).")
    if int(top["mispriming_hits"]) > 0:
        warnings.append("Top primer pair has a 3'-12mer mispriming hit in background_fasta.")

    tool_name = "mc_design_primers"
    params = {
        "mode": args.mode,
        "tm_opt": args.tm_opt,
        "na_mm": args.na_mm,
        "mg_mm": args.mg_mm,
        "dntp_mm": args.dntp_mm,
        "primer_nm": args.primer_nm,
    }
    artifacts = [
        write_artifact(
            tsv_path,
            tool_name,
            f"Primer design candidates ({args.mode})",
            inputs=[args.template, args.background_fasta],
            parameters=params,
            shape=df.shape,
        ),
        write_artifact(
            png_path,
            tool_name,
            f"Primer design map & Tm profile ({args.mode})",
            inputs=[args.template],
            parameters=params,
        ),
    ]

    emit({
        "status": "ok",
        "tool": tool_name,
        "mode": args.mode,
        "template_id": tmpl["id"],
        "template_length": L,
        "n_candidates": int(len(df)),
        "top_pair": top,
        "candidates_preview": df.head(5).to_dict(orient="records"),
        "outputs": {
            "primer_candidates": str(tsv_path),
            "design_map": str(png_path),
        },
        "artifacts": artifacts,
        "warnings": warnings,
    })


if __name__ == "__main__":
    main()
