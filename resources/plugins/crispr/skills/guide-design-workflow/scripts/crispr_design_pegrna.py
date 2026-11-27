#!/usr/bin/env python3
"""crispr_design_pegrna: Prime Editing (PE2 / PE3 / PE3b) pegRNA and nicking
sgRNA (ngRNA) design with SantaLucia PBS Tm scoring, RTT 5'-C scaffold-pairing
guard, PAM-disrupting synonymous codon discovery, and optional 3' tevopreQ1 epegRNA motif.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, revcomp, setup_matplotlib, write_artifact

TEVOPREQ1 = "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA"
SPCAS9_SCAFFOLD = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

NN_DH_DS = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

CODON_TABLE = {
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


def _pbs_tm_c(seq: str) -> float:
    s = seq.upper()
    if len(s) < 2:
        return 0.0
    dh = sum(NN_DH_DS.get(s[i : i + 2], (-8.0, -21.9))[0] for i in range(len(s) - 1)) + 0.2
    ds = sum(NN_DH_DS.get(s[i : i + 2], (-8.0, -21.9))[1] for i in range(len(s) - 1)) - 5.7
    # 50 mM monovalent equivalent + 250 nM effective strand concentration
    tm_k = (1000.0 * dh) / (ds + 1.987 * math.log(2.5e-7 / 4.0))
    inv_tm = (1.0 / tm_k) + (4.29 * ((s.count("G") + s.count("C")) / len(s)) - 3.95) * 1e-5 * math.log(0.15)
    return round((1.0 / inv_tm) - 273.15, 2)


def _locate_edit(wt: str, ed: str) -> tuple[int, int, int, str, int]:
    if wt == ed:
        fail("wt_seq and edited_seq are identical; provide a sequence carrying the desired edit.")
    n_wt, n_ed = len(wt), len(ed)
    prefix = 0
    while prefix < n_wt and prefix < n_ed and wt[prefix] == ed[prefix]:
        prefix += 1
    suffix = 0
    while (suffix < n_wt - prefix) and (suffix < n_ed - prefix) and wt[n_wt - 1 - suffix] == ed[n_ed - 1 - suffix]:
        suffix += 1
    wt_end = n_wt - suffix
    ed_end = n_ed - suffix
    wt_len = wt_end - prefix
    ed_len = ed_end - prefix
    if wt_len == 0 and ed_len > 0:
        etype = "insertion"
    elif ed_len == 0 and wt_len > 0:
        etype = "deletion"
    elif wt_len == ed_len:
        etype = "substitution"
    else:
        etype = "replacement"
    return prefix, wt_end, ed_end, etype, ed_len - wt_len


def _suggest_pam_silent_mutation(ed_sense: str, pam_sense_start: int, pam_strand: str, frame_offset: int | None) -> str:
    """Suggest an in-frame synonymous mutation disrupting SpCas9 NGG PAM on sense sequence if possible."""
    if frame_offset is None:
        return ""
    # For '+' strand, PAM is ed_sense[pam_sense_start : pam_sense_start+3] (NGG -> disrupt G at +1 or +2)
    # For '-' strand, PAM on '-' corresponds to CCN at ed_sense[pam_sense_start : pam_sense_start+3] (disrupt C at +0 or +1)
    mut_positions = [pam_sense_start + 1, pam_sense_start + 2] if pam_strand == "+" else [pam_sense_start, pam_sense_start + 1]
    for pos in mut_positions:
        if pos < 0 or pos >= len(ed_sense):
            continue
        codon_start = frame_offset + ((pos - frame_offset) // 3) * 3
        if codon_start < 0 or codon_start + 3 > len(ed_sense):
            continue
        wt_codon = ed_sense[codon_start : codon_start + 3]
        wt_aa = CODON_TABLE.get(wt_codon)
        if not wt_aa or wt_aa == "*":
            continue
        in_codon = pos - codon_start
        orig_b = ed_sense[pos]
        for alt_b in "ACGT":
            if alt_b == orig_b:
                continue
            cand_codon = wt_codon[:in_codon] + alt_b + wt_codon[in_codon + 1 :]
            if CODON_TABLE.get(cand_codon) == wt_aa:
                return f"pos{pos + 1}:{orig_b}>{alt_b}({wt_codon}>{cand_codon},p.{wt_aa}=)"
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--wt_seq", required=True)
    ap.add_argument("--edited_seq", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--nuclease", choices=["SpCas9", "SpCas9-NG"], default="SpCas9")
    ap.add_argument("--pbs_min", type=int, default=10)
    ap.add_argument("--pbs_max", type=int, default=15)
    ap.add_argument("--rtt_homology_min", type=int, default=8)
    ap.add_argument("--rtt_max", type=int, default=32)
    ap.add_argument("--add_epgrna_motif", action="store_true", default=False)
    ap.add_argument("--cds_frame_offset", type=int, choices=[0, 1, 2], default=None)
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    wt = re.sub(r"\s+", "", args.wt_seq).upper()
    ed = re.sub(r"\s+", "", args.edited_seq).upper()
    if len(wt) < 60 or len(ed) < 60:
        fail("wt_seq and edited_seq must be >= 60 bp flanking the target edit.")
    if not all(b in "ACGT" for b in wt) or not all(b in "ACGT" for b in ed):
        fail("wt_seq and edited_seq must contain unambiguous A/C/G/T bases only.")

    e_start, wt_end, ed_end, edit_type, delta_bp = _locate_edit(wt, ed)
    pam_re = re.compile(r"[ACGT]GG" if args.nuclease == "SpCas9" else r"[ACGT]G")
    pam_len = 3 if args.nuclease == "SpCas9" else 2
    L = 20

    pegrna_rows: list[dict[str, Any]] = []
    warnings: list[str] = []

    # Scan '+' and '-' strands of WT for Cas9 H840A nickase sites upstream of the edit on that strand
    for strand in ("+", "-"):
        wt_s = wt if strand == "+" else revcomp(wt)
        ed_s = ed if strand == "+" else revcomp(ed)
        if strand == "+":
            es_idx, ee_wt, ee_ed = e_start, wt_end, ed_end
        else:
            es_idx = len(wt) - wt_end
            ee_wt = len(wt) - e_start
            ee_ed = len(ed) - e_start

        for i in range(0, len(wt_s) - L - pam_len + 1):
            spacer = wt_s[i : i + L]
            pam = wt_s[i + L : i + L + pam_len]
            if not pam_re.fullmatch(pam):
                continue
            nick_pos = i + 17  # H840A nicks non-target strand (= same sense as spacer!) 3 bp upstream of PAM
            nick_to_edit = es_idx - nick_pos
            if nick_to_edit < 0 or nick_to_edit > 25:
                continue

            # Check if edit disrupts PAM (i+20..i+20+pam_len) or seed (i+12..i+20)
            disrupts_pam_seed = ee_wt > (i + 12) and es_idx < (i + L + pam_len)
            pam_sense_start = (i + L) if strand == "+" else (len(ed) - (i + L + pam_len))
            silent_hint = (
                "" if disrupts_pam_seed else _suggest_pam_silent_mutation(ed, pam_sense_start, strand, args.cds_frame_offset)
            )

            min_rtt_end = ee_ed + args.rtt_homology_min
            max_rtt_end = min(len(ed_s), nick_pos + args.rtt_max)
            if min_rtt_end > max_rtt_end:
                continue

            for pbs_len in range(args.pbs_min, args.pbs_max + 1):
                if nick_pos - pbs_len < 0:
                    continue
                # Non-target strand 5'->3' left of nick is wt_s[nick_pos - pbs_len : nick_pos]
                # 3' PBS of pegRNA hybridizes to it -> reverse complement!
                pbs_seq = revcomp(wt_s[nick_pos - pbs_len : nick_pos])
                pbs_tm = _pbs_tm_c(pbs_seq)
                pbs_gc = round((pbs_seq.count("G") + pbs_seq.count("C")) / pbs_len, 4)

                for rtt_end in range(min_rtt_end, max_rtt_end + 1, 2):
                    # Newly synthesized DNA 5'->3' from nick is ed_s[nick_pos : rtt_end]
                    # pegRNA RTT templates that synthesis -> reverse complement!
                    rtt_seq = revcomp(ed_s[nick_pos:rtt_end])
                    rtt_len = len(rtt_seq)
                    hom_len = rtt_end - ee_ed
                    ext_3p = rtt_seq + pbs_seq + (f"TTTT{TEVOPREQ1}"[:0] if not args.add_epgrna_motif else f"TTAA{TEVOPREQ1}")
                    full_pegrna = spacer + SPCAS9_SCAFFOLD + ext_3p

                    rtt_5p_c = rtt_seq.startswith("C")
                    has_tttt = "TTTT" in (spacer + ext_3p)

                    # Composite PrimeDesign-inspired heuristic score (0..100)
                    score = 78.0
                    score -= 1.2 * nick_to_edit
                    score -= 1.8 * abs(pbs_tm - 39.0)
                    score += min(10.0, 0.8 * hom_len)
                    if disrupts_pam_seed:
                        score += 12.0
                    elif silent_hint:
                        score += 4.0
                    if rtt_5p_c:
                        score -= 16.0
                    if has_tttt:
                        score -= 35.0
                    if not (0.35 <= pbs_gc <= 0.70):
                        score -= 8.0
                    score = round(max(1.0, min(99.0, score)), 2)

                    # Cloning oligos for Golden Gate pegRNA assembly (spacer CACC/AAAC + 3' extension oligo)
                    sp_fwd = f"CACC{spacer}" if spacer.startswith("G") else f"CACCG{spacer}"
                    sp_rev = f"AAAC{revcomp(spacer)}" if spacer.startswith("G") else f"AAAC{revcomp(spacer)}C"

                    pegrna_rows.append(
                        {
                            "pegrna_id": f"peg_{strand}_d{nick_to_edit}_P{pbs_len}_R{rtt_len}",
                            "strand": strand,
                            "spacer": spacer,
                            "pam": pam,
                            "nick_to_edit_bp": int(nick_to_edit),
                            "pbs_len": pbs_len,
                            "pbs_seq": pbs_seq,
                            "pbs_tm_c": pbs_tm,
                            "pbs_gc": pbs_gc,
                            "rtt_len": rtt_len,
                            "rtt_homology_bp": int(hom_len),
                            "rtt_seq": rtt_seq,
                            "rtt_starts_with_c": rtt_5p_c,
                            "disrupts_pam_or_seed": disrupts_pam_seed,
                            "pam_silent_mutation_hint": silent_hint,
                            "has_pol3_tttt": has_tttt,
                            "extension_3prime": ext_3p,
                            "full_pegrna_seq": full_pegrna,
                            "spacer_oligo_fwd": sp_fwd,
                            "spacer_oligo_rev": sp_rev,
                            "pegrna_score": score,
                        }
                    )

    if not pegrna_rows:
        fail(
            f"No valid {args.nuclease} pegRNA candidates found within 25 bp upstream of the edit; "
            "try --nuclease SpCas9-NG or increase --rtt_max / flanking sequence length."
        )

    pegrna_rows.sort(key=lambda r: (-int(not r["has_pol3_tttt"]), -r["pegrna_score"]))
    for idx, r in enumerate(pegrna_rows, start=1):
        r["rank"] = idx

    # Design PE3 (40..90 bp opposite strand on WT) and PE3b (overlapping edited locus on opposite strand of EDITED seq!)
    best_strand = pegrna_rows[0]["strand"]
    opp_strand = "-" if best_strand == "+" else "+"
    ngrna_rows: list[dict[str, Any]] = []

    best_nick_sense = (
        wt.find(pegrna_rows[0]["spacer"]) + 17
        if best_strand == "+"
        else len(wt) - (revcomp(wt).find(pegrna_rows[0]["spacer"]) + 17)
    )

    # Scan opposite strand on EDITED sequence so PE3b edit-specific nicking guides are captured alongside distal PE3 guides
    ed_opp = ed if opp_strand == "+" else revcomp(ed)
    wt_opp = wt if opp_strand == "+" else revcomp(wt)
    ed_opp_start = e_start if opp_strand == "+" else len(ed) - ed_end
    ed_opp_end = ed_end if opp_strand == "+" else len(ed) - e_start

    for i in range(0, len(ed_opp) - L - pam_len + 1):
        sp_ng = ed_opp[i : i + L]
        pam_ng = ed_opp[i + L : i + L + pam_len]
        if not pam_re.fullmatch(pam_ng) or "TTTT" in sp_ng:
            continue
        nick_opp_sense = (i + 17) if opp_strand == "+" else (len(ed) - (i + 17))
        dist = nick_opp_sense - best_nick_sense

        overlaps_edit = (i + L + pam_len > ed_opp_start) and (i < ed_opp_end) and (sp_ng not in wt_opp)
        if overlaps_edit:
            ng_type = "PE3b"
            ng_score = 92.0 - 0.15 * abs(dist)
        elif 40 <= abs(dist) <= 100:
            ng_type = "PE3"
            ng_score = 78.0 - 0.35 * abs(abs(dist) - 65.0)
        else:
            continue

        ngrna_rows.append(
            {
                "ngrna_id": f"{ng_type}_{opp_strand}_{abs(int(dist))}bp",
                "mode": ng_type,
                "strand": opp_strand,
                "spacer": sp_ng,
                "pam": pam_ng,
                "nick_distance_bp": int(dist),
                "edit_specific_pe3b": overlaps_edit,
                "ngrna_score": round(ng_score, 2),
            }
        )

    ngrna_rows.sort(key=lambda r: (-int(r["edit_specific_pe3b"]), -r["ngrna_score"]))
    n_pe3b = sum(1 for r in ngrna_rows if r["edit_specific_pe3b"])
    if not ngrna_rows:
        warnings.append("No clean PE3/PE3b opposite-strand nicking sgRNA found within 40-100 bp; use PE2 or epegRNA.")
    if not any(r["disrupts_pam_or_seed"] for r in pegrna_rows[:5]):
        warnings.append(
            "Top pegRNA candidates do not inherently disrupt the PAM/seed; consider introducing a silent PAM-disrupting mutation or PE3b nicking guide."
        )

    peg_df = pd.DataFrame(pegrna_rows)
    ng_df = pd.DataFrame(
        ngrna_rows,
        columns=["ngrna_id", "mode", "strand", "spacer", "pam", "nick_distance_bp", "edit_specific_pe3b", "ngrna_score"],
    )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    peg_path = out_dir / "pegrna_candidates.tsv"
    ng_path = out_dir / "ngrna_candidates.tsv"
    fig_path = out_dir / "pegrna_architecture_plot.png"

    peg_df.to_csv(peg_path, sep="\t", index=False)
    ng_df.to_csv(ng_path, sep="\t", index=False)

    write_artifact(
        peg_path,
        "crispr_design_pegrna",
        f"Prime Editing pegRNA Candidates ({args.nuclease}, {edit_type})",
        parameters={"nuclease": args.nuclease, "edit_type": edit_type, "add_epgrna_motif": bool(args.add_epgrna_motif)},
        shape=(int(len(peg_df)), int(peg_df.shape[1])),
    )
    write_artifact(
        ng_path,
        "crispr_design_pegrna",
        "PE3 / PE3b Nicking sgRNA Candidates",
        parameters={"nuclease": args.nuclease},
        shape=(int(len(ng_df)), int(ng_df.shape[1])),
    )

    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.6), dpi=180)
    sc = ax1.scatter(
        peg_df["pbs_tm_c"],
        peg_df["rtt_len"],
        c=peg_df["pegrna_score"],
        cmap="viridis",
        s=40,
        alpha=0.9,
    )
    ax1.axvspan(37.0, 42.0, color="#dcfce7", alpha=0.5, label="Optimal PBS Tm (37-42 C)")
    fig.colorbar(sc, ax=ax1, label="pegRNA Score")
    ax1.set_xlabel("PBS Melting Temperature Tm (C)")
    ax1.set_ylabel("RTT Length (nt)")
    ax1.set_title(f"pegRNA PBS vs RTT Design Space ({edit_type})")
    ax1.legend(frameon=False, fontsize=8, loc="best")

    top10 = peg_df.head(10).iloc[::-1]
    y_pos = np.arange(len(top10))
    cols_b = ["#059669" if not c else "#d97706" for c in top10["rtt_starts_with_c"]]
    ax2.barh(y_pos, top10["pegrna_score"], color=cols_b)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(top10["pegrna_id"], fontsize=8)
    ax2.set_xlabel("Composite pegRNA Score (green = non-C 5' RTT)")
    ax2.set_title("Top Ranked pegRNA Candidates")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "crispr_design_pegrna",
        f"Prime Editing pegRNA Design Space ({edit_type})",
        parameters={"nuclease": args.nuclease, "edit_type": edit_type},
    )

    emit(
        {
            "status": "ok",
            "tool": "crispr_design_pegrna",
            "nuclease": args.nuclease,
            "edit_type": edit_type,
            "edit_start_0based": int(e_start),
            "edit_delta_bp": int(delta_bp),
            "epeg_motif_added": bool(args.add_epgrna_motif),
            "n_pegrna_candidates": int(len(peg_df)),
            "n_ngrna_candidates": int(len(ng_df)),
            "n_pe3b_candidates": int(n_pe3b),
            "top_pegrnas": peg_df.head(10).to_dict(orient="records"),
            "top_ngrnas": ng_df.head(10).to_dict(orient="records"),
            "outputs": {
                "pegrna_candidates": str(peg_path),
                "ngrna_candidates": str(ng_path),
                "architecture_plot": str(fig_path),
            },
            "artifacts": [str(peg_path), str(ng_path), str(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
