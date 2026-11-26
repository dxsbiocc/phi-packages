#!/usr/bin/env python3
"""Design CRISPR guide RNAs with efficiency, Bae MMEJ frameshift, CBE/ABE window, and Golden Gate oligo outputs."""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    find_enzyme_collisions,
    parse_enzyme_list,
    parse_fasta,
    revcomp,
    setup_matplotlib,
    write_artifact,
)

NUCLEASE_SPECS = {
    "SpCas9": {"pam_side": "3prime", "spacer_len": 20, "pam_regex": r"[ACGT]GG", "pam_len": 3, "cut_offset": 17},
    "SpCas9-NG": {"pam_side": "3prime", "spacer_len": 20, "pam_regex": r"[ACGT]G", "pam_len": 2, "cut_offset": 17},
    "SaCas9": {"pam_side": "3prime", "spacer_len": 21, "pam_regex": r"[ACGT]{2}G[AG]{2}T", "pam_len": 6, "cut_offset": 18},
    "Cas12a": {"pam_side": "5prime", "spacer_len": 23, "pam_regex": r"TTT[ACG]", "pam_len": 4, "cut_offset": 20},
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


def _efficiency_score(spacer: str, pam: str, nuclease: str) -> float:
    """Deterministic sequence-context efficiency model (0..100)."""
    L = len(spacer)
    gc = (spacer.count("G") + spacer.count("C")) / max(1, L)
    z = 0.35

    if nuclease == "Cas12a":
        # 5' PAM TTTV -> seed is positions 0..7
        if spacer[0] == "T":
            z -= 0.40
        if spacer[0] in ("G", "C"):
            z += 0.25
        seed = spacer[:8]
        z -= 0.30 * seed.count("TT")
        purines = seed.count("A") + seed.count("G")
        if 3 <= purines <= 6:
            z += 0.25
        if pam[-1] == "C":
            z += 0.15
    else:
        # 3' PAM Cas9 architectures: pos -1 is immediately PAM-proximal
        if spacer[-1] == "G":
            z += 0.60
        elif spacer[-1] == "C":
            z -= 0.55
        elif spacer[-1] == "T":
            z -= 0.18
        if spacer[-2] in ("G", "A"):
            z += 0.25
        elif spacer[-2] == "C":
            z -= 0.25
        if L >= 5 and spacer[-5] == "C":
            z += 0.28
        if L >= 4 and spacer[-4] in ("C", "G"):
            z += 0.20
        seed = spacer[-8:]
        z -= 0.35 * seed.count("TT")
        purines = seed.count("A") + seed.count("G")
        if 3 <= purines <= 6:
            z += 0.25
        if pam[0] in ("G", "C"):
            z += 0.12
        elif pam[0] == "T":
            z -= 0.10

    if 0.40 <= gc <= 0.65:
        z += 0.30
    else:
        dist = 0.40 - gc if gc < 0.40 else gc - 0.65
        z -= 4.0 * dist

    if "TTTT" in spacer:
        z -= 1.60
    elif any(b * 4 in spacer for b in ("A", "C", "G")):
        z -= 0.65

    score = 100.0 / (1.0 + math.exp(-z))
    return round(max(1.0, min(99.0, score)), 2)


def _mmej_out_of_frame_score(seq: str, cut_pos: int) -> float:
    """Bae et al. (2014) microhomology out-of-frame score (0..100) around cut_pos."""
    left_start = max(0, cut_pos - 30)
    right_end = min(len(seq), cut_pos + 30)
    left = seq[left_start:cut_pos]
    right = seq[cut_pos:right_end]
    if len(left) < 2 or len(right) < 2:
        return 66.7

    total_w = 0.0
    oof_w = 0.0
    max_k = min(10, len(left), len(right))
    seen_del: set[tuple[int, int]] = set()

    for k in range(max_k, 1, -1):
        for i in range(0, len(left) - k + 1):
            sub = left[i : i + k]
            start_j = 0
            while True:
                j = right.find(sub, start_j)
                if j == -1:
                    break
                del_len = (len(left) - i) + j
                if del_len >= 2 and (del_len, i) not in seen_del:
                    # Skip strict subsets of already counted longer microhomologies at same deletion length
                    seen_del.add((del_len, i))
                    w = k * math.exp(-del_len / 17.0)
                    total_w += w
                    if del_len % 3 != 0:
                        oof_w += w
                start_j = j + 1

    if total_w <= 1e-12:
        return 66.7
    return round(100.0 * oof_w / total_w, 2)


def _cloning_oligos(spacer: str, backbone: str) -> tuple[str, str]:
    if backbone == "none":
        return "", ""
    # lentiCRISPRv2_BsmBI, pX330_BbsI, pLKO5_BsmBI all use 5'-CACC / 5'-AAAC Golden Gate overhangs
    if spacer.startswith("G"):
        fwd = f"CACC{spacer}"
        rev = f"AAAC{revcomp(spacer)}"
    else:
        fwd = f"CACCG{spacer}"
        rev = f"AAAC{revcomp(spacer)}C"
    return fwd, rev


def _base_editing_annotations(
    seq: str, spacer: str, strand: str, spacer_start_sense: int, mode: str
) -> dict:
    """Audit base-editing window [4..8] (1-based) and simulate in-frame codon changes."""
    win_start_idx = 3  # 0-based index for pos 4
    win_end_idx = min(len(spacer), 8)  # 0-based exclusive for pos 8
    win_seq = spacer[win_start_idx:win_end_idx]
    target_base = "C" if mode == "CBE" else ("A" if mode == "ABE" else "")
    edit_to = "T" if mode == "CBE" else ("G" if mode == "ABE" else "")

    if not target_base:
        return {
            "base_edit_window_seq": win_seq,
            "target_bases_in_window": 0,
            "bystander_count": 0,
            "predicted_aa_edits": "",
            "creates_premature_stop": False,
        }

    editable_positions = [i for i in range(win_start_idx, win_end_idx) if spacer[i] == target_base]
    n_target = len(editable_positions)
    bystander = max(0, n_target - 1)

    aa_edits: list[str] = []
    creates_stop = False
    for pos0 in editable_positions:
        if strand == "+":
            sense_idx = spacer_start_sense + pos0
            new_sense_base = edit_to
        else:
            # Spacer pos0 on '-' strand maps from right end of sense protospacer
            sense_idx = spacer_start_sense + (len(spacer) - 1 - pos0)
            new_sense_base = revcomp(edit_to)
        codon_start = (sense_idx // 3) * 3
        if codon_start + 3 <= len(seq):
            wt_codon = seq[codon_start : codon_start + 3]
            if all(b in "ACGT" for b in wt_codon):
                offset = sense_idx - codon_start
                mut_codon = wt_codon[:offset] + new_sense_base + wt_codon[offset + 1 :]
                wt_aa = CODON_TABLE.get(wt_codon, "X")
                mut_aa = CODON_TABLE.get(mut_codon, "X")
                if wt_aa != mut_aa:
                    aa_pos = (codon_start // 3) + 1
                    aa_edits.append(f"p.{wt_aa}{aa_pos}{mut_aa}@{pos0 + 1}")
                    if mut_aa == "*":
                        creates_stop = True

    return {
        "base_edit_window_seq": win_seq,
        "target_bases_in_window": n_target,
        "bystander_count": bystander,
        "predicted_aa_edits": ";".join(aa_edits),
        "creates_premature_stop": creates_stop,
    }


def _scan_target(
    target_id: str,
    seq: str,
    nuclease: str,
    mode: str,
    backbone: str,
    enzymes: list[str],
    gc_min: float,
    gc_max: float,
) -> list[dict]:
    spec = NUCLEASE_SPECS[nuclease]
    L = spec["spacer_len"]
    pam_len = spec["pam_len"]
    pam_re = re.compile(spec["pam_regex"])
    seq_len = len(seq)
    rc_seq = revcomp(seq)
    candidates: list[dict] = []

    for strand, strand_seq in (("+", seq), ("-", rc_seq)):
        if spec["pam_side"] == "3prime":
            for i in range(0, seq_len - L - pam_len + 1):
                spacer = strand_seq[i : i + L]
                pam = strand_seq[i + L : i + L + pam_len]
                if not pam_re.fullmatch(pam) or not all(b in "ACGT" for b in spacer):
                    continue
                if strand == "+":
                    sense_start = i
                    cut_pos = i + spec["cut_offset"]
                else:
                    sense_start = seq_len - (i + L)
                    cut_pos = seq_len - (i + spec["cut_offset"])
                candidates.append(
                    _build_candidate(
                        target_id, seq, seq_len, spacer, pam, strand, sense_start, cut_pos,
                        nuclease, mode, backbone, enzymes, gc_min, gc_max,
                    )
                )
        else:
            # 5' PAM (Cas12a)
            for i in range(0, seq_len - pam_len - L + 1):
                pam = strand_seq[i : i + pam_len]
                spacer = strand_seq[i + pam_len : i + pam_len + L]
                if not pam_re.fullmatch(pam) or not all(b in "ACGT" for b in spacer):
                    continue
                if strand == "+":
                    sense_start = i + pam_len
                    cut_pos = i + pam_len + spec["cut_offset"]
                else:
                    sense_start = seq_len - (i + pam_len + L)
                    cut_pos = seq_len - (i + pam_len + spec["cut_offset"])
                candidates.append(
                    _build_candidate(
                        target_id, seq, seq_len, spacer, pam, strand, sense_start, cut_pos,
                        nuclease, mode, backbone, enzymes, gc_min, gc_max,
                    )
                )
    return candidates


def _build_candidate(
    target_id: str,
    seq: str,
    seq_len: int,
    spacer: str,
    pam: str,
    strand: str,
    sense_start: int,
    cut_pos: int,
    nuclease: str,
    mode: str,
    backbone: str,
    enzymes: list[str],
    gc_min: float,
    gc_max: float,
) -> dict:
    gc = round((spacer.count("G") + spacer.count("C")) / len(spacer), 4)
    cds_pct = round(100.0 * cut_pos / max(1, seq_len), 2)
    eff = _efficiency_score(spacer, pam, nuclease)
    oof = _mmej_out_of_frame_score(seq, cut_pos)
    oligo_fwd, oligo_rev = _cloning_oligos(spacer, backbone)

    has_tttt = "TTTT" in spacer
    has_homopol = any(b * 4 in spacer for b in ("A", "C", "G", "T"))
    enz_hits = find_enzyme_collisions(spacer, enzymes)
    gc_out = gc < gc_min or gc > gc_max

    pos_flag = "ok"
    if mode == "knockout":
        if cds_pct < 5.0:
            pos_flag = "early_nterm_risk"
        elif cds_pct > 85.0:
            pos_flag = "late_cterm_risk"

    be = _base_editing_annotations(seq, spacer, strand, sense_start, mode)

    flags: list[str] = []
    if has_tttt:
        flags.append("pol3_tttt")
    if enz_hits:
        flags.append("enzyme_site:" + "+".join(enz_hits))
    if gc_out:
        flags.append("gc_out_of_range")
    if pos_flag != "ok":
        flags.append(pos_flag)
    if mode in ("CBE", "ABE") and be["target_bases_in_window"] == 0:
        flags.append("no_editable_base_in_window")

    passes = not (has_tttt or bool(enz_hits) or gc_out or (mode == "knockout" and pos_flag != "ok") or (mode in ("CBE", "ABE") and be["target_bases_in_window"] == 0))

    if mode == "knockout":
        comp = 0.55 * eff + 0.45 * oof
    elif mode in ("CBE", "ABE"):
        stop_bonus = 12.0 if be["creates_premature_stop"] else 0.0
        target_bonus = 10.0 if be["target_bases_in_window"] >= 1 else -40.0
        bystander_pen = 8.0 * be["bystander_count"]
        comp = 0.65 * eff + target_bonus + stop_bonus - bystander_pen
    else:
        comp = eff

    if not passes:
        comp -= 50.0

    return {
        "target_id": target_id,
        "guide_id": f"{target_id}_{nuclease}_{strand}_{cut_pos}",
        "nuclease": nuclease,
        "mode": mode,
        "strand": strand,
        "cut_site_bp": int(cut_pos),
        "cds_pct": cds_pct,
        "spacer": spacer,
        "pam": pam,
        "gc_content": gc,
        "efficiency_score": eff,
        "out_of_frame_score": oof,
        "has_pol3_terminator": has_tttt,
        "has_homopolymer": has_homopol,
        "enzyme_collisions": ",".join(enz_hits),
        "positional_flag": pos_flag,
        "base_edit_window_seq": be["base_edit_window_seq"],
        "target_bases_in_window": be["target_bases_in_window"],
        "bystander_count": be["bystander_count"],
        "predicted_aa_edits": be["predicted_aa_edits"],
        "creates_premature_stop": be["creates_premature_stop"],
        "oligo_fwd": oligo_fwd,
        "oligo_rev": oligo_rev,
        "passes_filters": passes,
        "flags": ";".join(flags) if flags else "none",
        "composite_score": round(comp, 2),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Design CRISPR guide RNAs from target FASTA")
    ap.add_argument("--target_fasta", required=True, help="Target FASTA path")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--nuclease", choices=list(NUCLEASE_SPECS.keys()), default="SpCas9")
    ap.add_argument("--mode", choices=["knockout", "CRISPRi", "CRISPRa", "CBE", "ABE"], default="knockout")
    ap.add_argument(
        "--cloning_backbone",
        choices=["lentiCRISPRv2_BsmBI", "pX330_BbsI", "pLKO5_BsmBI", "none"],
        default="lentiCRISPRv2_BsmBI",
    )
    ap.add_argument("--avoid_enzymes", action="append", default=None)
    ap.add_argument("--gc_min", type=float, default=0.35)
    ap.add_argument("--gc_max", type=float, default=0.75)
    ap.add_argument("--top_n", type=int, default=10)
    args = ap.parse_args()

    if args.gc_min < 0 or args.gc_max > 1 or args.gc_min >= args.gc_max:
        fail(f"Invalid GC bounds: gc_min={args.gc_min}, gc_max={args.gc_max}")

    import pandas as pd

    records = parse_fasta(args.target_fasta)
    enzymes = parse_enzyme_list(args.avoid_enzymes)
    all_rows: list[dict] = []
    warnings: list[str] = []

    for tid, seq in records:
        rows = _scan_target(
            tid, seq, args.nuclease, args.mode, args.cloning_backbone,
            enzymes, args.gc_min, args.gc_max,
        )
        if not rows:
            warnings.append(f"Target '{tid}' ({len(seq)} bp): no valid {args.nuclease} PAM sites found.")
            continue
        rows.sort(key=lambda r: (-int(r["passes_filters"]), -r["composite_score"], -r["efficiency_score"]))
        for idx, r in enumerate(rows, start=1):
            r["rank"] = idx
        all_rows.extend(rows)

    if not all_rows:
        fail(f"No candidate {args.nuclease} guide RNAs found across {len(records)} target sequence(s).")

    df = pd.DataFrame(all_rows)
    n_passing = int(df["passes_filters"].sum())
    if n_passing == 0:
        warnings.append(
            "Zero candidates passed all hard quality/position filters; inspect flagged candidates carefully."
        )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / "designed_guides.tsv"
    fig_path = out_dir / "guide_design_overview.png"
    df.to_csv(tsv_path, sep="\t", index=False)

    write_artifact(
        tsv_path,
        "crispr_design_guides",
        f"Designed {args.nuclease} Guides ({args.mode})",
        inputs=[args.target_fasta],
        parameters={
            "nuclease": args.nuclease,
            "mode": args.mode,
            "cloning_backbone": args.cloning_backbone,
            "avoid_enzymes": enzymes,
            "gc_min": args.gc_min,
            "gc_max": args.gc_max,
        },
        shape=(int(len(df)), int(df.shape[1])),
    )

    # 2-panel visualization: Locus track + Efficiency vs Out-of-Frame scatter
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=180)

    first_tid = records[0][0]
    sub = df[df["target_id"] == first_tid].copy()
    passing = sub[sub["passes_filters"]]
    flagged = sub[~sub["passes_filters"]]

    if args.mode == "knockout":
        ax1.axvspan(0, 5, color="#fee2e2", alpha=0.65, label="<5% N-term risk")
        ax1.axvspan(85, 100, color="#fee2e2", alpha=0.65, label=">85% C-term risk")
    if not flagged.empty:
        ax1.scatter(flagged["cds_pct"], flagged["composite_score"], c="#94a3b8", s=26, alpha=0.7, label="Flagged")
    if not passing.empty:
        ax1.scatter(passing["cds_pct"], passing["composite_score"], c="#0d9488", s=38, alpha=0.9, label="Passing")
        for _, row in passing.head(5).iterrows():
            ax1.annotate(
                f"#{int(row['rank'])}",
                (row["cds_pct"], row["composite_score"]),
                textcoords="offset points",
                xytext=(4, 4),
                fontsize=7.5,
                color="#0f172a",
            )
    ax1.set_xlim(-2, 102)
    ax1.set_xlabel("Cut Position Along Target (% length)")
    ax1.set_ylabel("Composite Guide Score")
    ax1.set_title(f"Target Track: {first_tid} ({args.nuclease})")
    ax1.legend(frameon=False, fontsize=7.5, loc="best")

    df_pass = df[df["passes_filters"]]
    df_flag = df[~df["passes_filters"]]
    if not df_flag.empty:
        ax2.scatter(df_flag["efficiency_score"], df_flag["out_of_frame_score"], c="#cbd5e1", s=24, label="Flagged")
    if not df_pass.empty:
        sc = ax2.scatter(
            df_pass["efficiency_score"],
            df_pass["out_of_frame_score"],
            c=df_pass["composite_score"],
            cmap="viridis",
            s=36,
            alpha=0.9,
            label="Passing",
        )
        fig.colorbar(sc, ax=ax2, label="Composite Score")
    ax2.set_xlabel("Sequence-Context Efficiency Score (0-100)")
    ax2.set_ylabel("Bae MMEJ Out-of-Frame Score (0-100)")
    ax2.set_title(f"Efficiency vs Frameshift Score ({args.mode})")
    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "crispr_design_guides",
        f"Guide Design Overview ({args.nuclease}, {args.mode})",
        inputs=[args.target_fasta],
        parameters={"nuclease": args.nuclease, "mode": args.mode},
    )

    top_guides = (
        df[df["rank"] <= args.top_n]
        .head(args.top_n * min(len(records), 5))
        .to_dict(orient="records")
    )

    emit(
        {
            "status": "ok",
            "tool": "crispr_design_guides",
            "nuclease": args.nuclease,
            "mode": args.mode,
            "cloning_backbone": args.cloning_backbone,
            "n_targets": len(records),
            "n_total_candidates": int(len(df)),
            "n_passing_candidates": n_passing,
            "top_guides": top_guides,
            "outputs": {
                "designed_guides": str(tsv_path),
                "overview_plot": str(fig_path),
            },
            "artifacts": [str(tsv_path), str(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
