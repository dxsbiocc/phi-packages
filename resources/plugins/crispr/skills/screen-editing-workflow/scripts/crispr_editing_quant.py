#!/usr/bin/env python3
"""Post-editing amplicon quantification: NHEJ frameshift vs in-frame, HDR, and CBE/ABE bystander profiling."""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, pick_col, read_table, revcomp, setup_matplotlib, write_artifact


def _locate_guide(amplicon: str, guide: str) -> tuple[int, str, int]:
    """Return (guide_start_0based_on_amp, strand, cut_site_bp)."""
    amp = amplicon.upper()
    g = guide.upper()
    idx = amp.find(g)
    if idx != -1:
        cut = idx + max(1, len(g) - 3)
        return idx, "+", cut
    rc = revcomp(g)
    idx_rc = amp.find(rc)
    if idx_rc != -1:
        cut = idx_rc + 3
        return idx_rc, "-", cut
    fail("guide_seq was not found on either strand of amplicon_seq.")
    return 0, "+", 0


def _global_align_affine(seq_a: str, seq_b: str) -> tuple[str, str]:
    """Affine-gap Needleman-Wunsch global alignment (seq_a=read, seq_b=ref)."""
    n, m = len(seq_a), len(seq_b)
    if seq_a == seq_b:
        return seq_a, seq_b
    NEG_INF = -10**9
    match_s, mis_s = 2, -3
    gap_open, gap_ext = -6, -1

    # M: match/mismatch, X: gap in seq_b (insertion in read), Y: gap in seq_a (deletion in read)
    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    trM = [[0] * (m + 1) for _ in range(n + 1)]
    trX = [[0] * (m + 1) for _ in range(n + 1)]
    trY = [[0] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = gap_open + (i - 1) * gap_ext
        trX[i][0] = 1 if i > 1 else 0
    for j in range(1, m + 1):
        Y[0][j] = gap_open + (j - 1) * gap_ext
        trY[0][j] = 2 if j > 1 else 0

    for i in range(1, n + 1):
        ca = seq_a[i - 1]
        M_i, X_i, Y_i = M[i], X[i], Y[i]
        M_p, X_p, Y_p = M[i - 1], X[i - 1], Y[i - 1]
        for j in range(1, m + 1):
            sc = match_s if ca == seq_b[j - 1] else mis_s
            # Transition into M
            prev_vals = (M_p[j - 1], X_p[j - 1], Y_p[j - 1])
            best_m = max(prev_vals)
            M_i[j] = best_m + sc
            trM[i][j] = prev_vals.index(best_m)

            # Transition into X (gap in ref -> insertion)
            o_x = M_p[j] + gap_open
            e_x = X_p[j] + gap_ext
            if o_x >= e_x:
                X_i[j] = o_x
                trX[i][j] = 0
            else:
                X_i[j] = e_x
                trX[i][j] = 1

            # Transition into Y (gap in read -> deletion)
            o_y = M_i[j - 1] + gap_open
            e_y = Y_i[j - 1] + gap_ext
            if o_y >= e_y:
                Y_i[j] = o_y
                trY[i][j] = 0
            else:
                Y_i[j] = e_y
                trY[i][j] = 2

    end_vals = (M[n][m], X[n][m], Y[n][m])
    state = end_vals.index(max(end_vals))
    aln_a: list[str] = []
    aln_b: list[str] = []
    i, j = n, m
    while i > 0 or j > 0:
        if state == 0:
            prev_st = trM[i][j]
            aln_a.append(seq_a[i - 1])
            aln_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
            state = prev_st
        elif state == 1:
            prev_st = trX[i][j]
            aln_a.append(seq_a[i - 1])
            aln_b.append("-")
            i -= 1
            state = prev_st
        else:
            prev_st = trY[i][j]
            aln_a.append("-")
            aln_b.append(seq_b[j - 1])
            j -= 1
            state = prev_st

    return "".join(reversed(aln_a)), "".join(reversed(aln_b))


def _load_mode_a(path: str | Path) -> list[dict]:
    p = Path(path)
    if not p.is_file():
        fail(f"Alleles input not found: {path}")
    if p.suffix.lower() in (".fa", ".fasta", ".fna"):
        rows: list[dict] = []
        cur_hdr = None
        chunks: list[str] = []
        with open(p, "rt", encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if cur_hdr is not None:
                        rows.append(_fasta_read_entry(cur_hdr, "".join(chunks).upper()))
                    cur_hdr = line[1:]
                    chunks = []
                else:
                    chunks.append(line)
        if cur_hdr is not None:
            rows.append(_fasta_read_entry(cur_hdr, "".join(chunks).upper()))
        if not rows:
            fail(f"No FASTA reads found in {path}")
        return rows

    df = read_table(p).copy()
    cols = list(df.columns)
    seq_c = pick_col(cols, "", ["sequence", "seq", "allele", "aligned_sequence", "read_seq"])
    cnt_c = pick_col(cols, "", ["count", "reads", "#reads", "n_reads"], required=False)
    smp_c = pick_col(cols, "", ["sample", "sample_id", "name"], required=False)

    out: list[dict] = []
    for _, r in df.iterrows():
        seq = str(r[seq_c]).strip().upper().replace("-", "")
        if not seq:
            continue
        cnt = int(float(r[cnt_c])) if cnt_c else 1
        smp = str(r[smp_c]).strip() if smp_c else p.stem
        out.append({"sample": smp, "seq": seq, "count": max(0, cnt)})
    if not out:
        fail("No valid read/allele records found in alleles table.")
    return out


def _fasta_read_entry(hdr: str, seq: str) -> dict:
    m_cnt = re.search(r"(?:count|reads|size)[=:](\d+)", hdr, flags=re.IGNORECASE)
    cnt = int(m_cnt.group(1)) if m_cnt else 1
    m_smp = re.search(r"sample[=:](\S+)", hdr, flags=re.IGNORECASE)
    smp = m_smp.group(1) if m_smp else "Sample1"
    return {"sample": smp, "seq": seq, "count": cnt}


def _load_crispresso_table(path: str | Path, amplicon_seq: str) -> tuple[list[dict], int]:
    df = read_table(path).copy()
    cols = list(df.columns)
    aln_r_c = pick_col(cols, "", ["Aligned_Sequence", "aligned_sequence"])
    aln_ref_c = pick_col(cols, "", ["Reference_Sequence", "reference_sequence"])
    cnt_c = pick_col(cols, "", ["#Reads", "Reads", "reads", "count"])
    smp_c = pick_col(cols, "", ["sample", "Sample"], required=False)

    # Determine offset of ungapped Reference_Sequence inside amplicon_seq
    first_ref = str(df[aln_ref_c].iloc[0]).strip().upper().replace("-", "")
    amp_u = amplicon_seq.upper()
    offset = amp_u.find(first_ref)
    if offset == -1:
        fail("CRISPResso2 Reference_Sequence could not be located inside amplicon_seq.")

    stem = Path(path).stem
    records: list[dict] = []
    for _, r in df.iterrows():
        ar = str(r[aln_r_c]).strip().upper()
        af = str(r[aln_ref_c]).strip().upper()
        cnt = int(float(r[cnt_c]))
        if cnt <= 0:
            continue
        smp = str(r[smp_c]).strip() if smp_c else stem
        records.append({"sample": smp, "aln_read": ar, "aln_ref": af, "count": cnt, "ref_offset": offset})
    if not records:
        fail("No positive-count allele rows found in CRISPResso2 table.")
    return records, offset


def _hdr_signature(amplicon_seq: str, hdr_seq: str) -> dict[int, str]:
    """Return map of {amplicon_0based_pos: hdr_base} for SNVs distinguishing hdr_seq from WT."""
    aln_h, aln_w = _global_align_affine(hdr_seq.upper(), amplicon_seq.upper())
    sig: dict[int, str] = {}
    w_pos = 0
    for ch, cw in zip(aln_h, aln_w):
        if cw != "-":
            if ch != cw and ch != "-":
                sig[w_pos] = ch
            w_pos += 1
    return sig


def _inspect_alignment(
    aln_read: str,
    aln_ref: str,
    ref_offset: int,
    win_start: int,
    win_end: int,
    hdr_sig: dict[int, str] | None,
    protospacer_amp_map: list[tuple[int, int]],
    guide_strand: str,
) -> dict:
    """Inspect aligned read vs ref in [win_start, win_end) and across protospacer positions."""
    n_ins = 0
    n_del = 0
    n_sub = 0
    net_indel = 0
    amp_base_at_pos: dict[int, str] = {}

    ref_pos = ref_offset
    for cr, cf in zip(aln_read, aln_ref):
        if cf == "-":
            # Insertion immediately before ref_pos
            if win_start <= ref_pos <= win_end:
                n_ins += 1
                net_indel += 1
        elif cr == "-":
            if win_start <= ref_pos < win_end:
                n_del += 1
                net_indel -= 1
            amp_base_at_pos[ref_pos] = "-"
            ref_pos += 1
        else:
            amp_base_at_pos[ref_pos] = cr
            if win_start <= ref_pos < win_end and cr != cf:
                n_sub += 1
            ref_pos += 1

    has_indel = (n_ins + n_del) > 0
    net_abs = abs(n_ins - n_del)

    # Check HDR signature
    hdr_match = False
    if hdr_sig:
        matched_sig = sum(1 for p, base in hdr_sig.items() if amp_base_at_pos.get(p) == base)
        if matched_sig == len(hdr_sig) and len(hdr_sig) > 0:
            hdr_match = True

    if hdr_match and not has_indel:
        category = "HDR"
    elif hdr_match and has_indel:
        category = "Imperfect_HDR"
    elif has_indel:
        category = "Frameshift_Indel" if (net_abs % 3 != 0) else "In_Frame_Indel"
    elif n_sub > 0:
        category = "Substitution_Only"
    else:
        category = "Unmodified"

    # Extract bases at protospacer positions 1..L oriented 5'->3' on guide strand
    proto_bases: dict[int, str] = {}
    comp_map = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"}
    for p1, amp_idx in protospacer_amp_map:
        b = amp_base_at_pos.get(amp_idx, "-")
        if guide_strand == "-":
            b = comp_map.get(b, b)
        proto_bases[p1] = b

    return {
        "category": category,
        "net_indel": net_indel if has_indel else 0,
        "has_indel": has_indel,
        "proto_bases": proto_bases,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Quantify CRISPR editing outcomes (NHEJ frameshift, HDR, CBE/ABE)")
    ap.add_argument("--alleles", default="", help="Mode A: allele count TSV/CSV or FASTA")
    ap.add_argument("--crispresso_table", default="", help="Mode B: CRISPResso2 Alleles_frequency_table TSV")
    ap.add_argument("--amplicon_seq", required=True, help="WT reference amplicon sequence")
    ap.add_argument("--guide_seq", required=True, help="Guide RNA spacer sequence")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--hdr_seq", default="", help="Optional HDR donor amplicon sequence")
    ap.add_argument("--editor_mode", choices=["nuclease", "CBE", "ABE"], default="nuclease")
    ap.add_argument("--window_bp", type=int, default=10)
    args = ap.parse_args()

    if not args.alleles and not args.crispresso_table:
        fail("Supply either --alleles (Mode A) or --crispresso_table (Mode B).")

    import numpy as np
    import pandas as pd

    amp = args.amplicon_seq.strip().upper()
    guide = args.guide_seq.strip().upper()
    g_start, g_strand, cut_site = _locate_guide(amp, guide)
    L = len(guide)
    win_start = max(0, cut_site - args.window_bp)
    win_end = min(len(amp), cut_site + args.window_bp)

    if g_strand == "+":
        proto_map = [(p + 1, g_start + p) for p in range(L)]
    else:
        proto_map = [(p + 1, g_start + L - 1 - p) for p in range(L)]

    hdr_sig = _hdr_signature(amp, args.hdr_seq.strip().upper()) if args.hdr_seq.strip() else None

    aligned_cache: dict[str, tuple[str, str]] = {}
    records: list[dict] = []
    if args.crispresso_table:
        input_mode = "crispresso2"
        src_path = args.crispresso_table
        raw_recs, _ = _load_crispresso_table(args.crispresso_table, amp)
        for r in raw_recs:
            res = _inspect_alignment(
                r["aln_read"], r["aln_ref"], r["ref_offset"], win_start, win_end, hdr_sig, proto_map, g_strand
            )
            records.append({"sample": r["sample"], "count": r["count"], **res})
    else:
        input_mode = "alleles"
        src_path = args.alleles
        raw_recs = _load_mode_a(args.alleles)
        for r in raw_recs:
            seq = r["seq"]
            if seq not in aligned_cache:
                aligned_cache[seq] = _global_align_affine(seq, amp)
            ar, af = aligned_cache[seq]
            res = _inspect_alignment(ar, af, 0, win_start, win_end, hdr_sig, proto_map, g_strand)
            records.append({"sample": r["sample"], "count": r["count"], **res})

    categories = ["Unmodified", "Frameshift_Indel", "In_Frame_Indel", "Substitution_Only", "HDR", "Imperfect_HDR"]
    by_sample: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        by_sample[r["sample"]].append(r)

    sample_summaries: list[dict] = []
    indel_rows: list[dict] = []
    base_matrix_rows: list[dict] = []
    warnings: list[str] = []

    target_base = "C" if args.editor_mode == "CBE" else ("A" if args.editor_mode == "ABE" else "")
    edit_to = "T" if args.editor_mode == "CBE" else ("G" if args.editor_mode == "ABE" else "")

    for smp, recs in by_sample.items():
        tot = sum(r["count"] for r in recs)
        cat_counts = {c: 0 for c in categories}
        indel_hist: dict[int, int] = defaultdict(int)
        pos_base_counts = {p: {"A": 0, "C": 0, "G": 0, "T": 0, "-": 0} for p in range(1, L + 1)}

        for r in recs:
            cnt = r["count"]
            cat_counts[r["category"]] += cnt
            if r["has_indel"]:
                sz = int(max(-30, min(30, r["net_indel"])))
                indel_hist[sz] += cnt
            for p1, b in r["proto_bases"].items():
                if b in pos_base_counts[p1]:
                    pos_base_counts[p1][b] += cnt

        pct = {c: round(100.0 * cat_counts[c] / max(1, tot), 2) for c in categories}
        total_indel_pct = round(pct["Frameshift_Indel"] + pct["In_Frame_Indel"] + pct["Imperfect_HDR"], 2)

        # Base-editing metrics across protospacer positions 1..L
        window_editable_rates: list[tuple[int, float]] = []
        for p1 in range(1, L + 1):
            wt_b = guide[p1 - 1]
            c_dict = pos_base_counts[p1]
            row_b = {
                "sample": smp,
                "position_1based": p1,
                "wt_base": wt_b,
                "in_be_window_4_8": 4 <= p1 <= 8,
                "A_pct": round(100.0 * c_dict["A"] / max(1, tot), 2),
                "C_pct": round(100.0 * c_dict["C"] / max(1, tot), 2),
                "G_pct": round(100.0 * c_dict["G"] / max(1, tot), 2),
                "T_pct": round(100.0 * c_dict["T"] / max(1, tot), 2),
                "del_pct": round(100.0 * c_dict["-"] / max(1, tot), 2),
            }
            base_matrix_rows.append(row_b)
            if target_base and 4 <= p1 <= 8 and wt_b == target_base:
                conv_pct = 100.0 * c_dict[edit_to] / max(1, tot)
                window_editable_rates.append((p1, round(conv_pct, 2)))

        if window_editable_rates:
            sorted_rates = sorted(window_editable_rates, key=lambda x: x[1], reverse=True)
            max_target_conv = sorted_rates[0][1]
            bystander_max = sorted_rates[1][1] if len(sorted_rates) > 1 else 0.0
            mean_target_conv = round(sum(v for _, v in window_editable_rates) / len(window_editable_rates), 2)
        else:
            max_target_conv = 0.0
            bystander_max = 0.0
            mean_target_conv = 0.0

        if args.editor_mode in ("CBE", "ABE") and bystander_max >= 5.0:
            warnings.append(f"Sample '{smp}': bystander {target_base}->{edit_to} editing reached {bystander_max:.1f}% inside window [4..8].")
        if args.editor_mode in ("CBE", "ABE") and total_indel_pct >= 5.0:
            warnings.append(f"Sample '{smp}': elevated indel byproduct rate ({total_indel_pct:.1f}%) for {args.editor_mode} base editing.")

        for sz in sorted(indel_hist.keys()):
            indel_rows.append(
                {
                    "sample": smp,
                    "indel_bp": sz,
                    "read_count": indel_hist[sz],
                    "pct_reads": round(100.0 * indel_hist[sz] / max(1, tot), 3),
                }
            )

        sample_summaries.append(
            {
                "sample": smp,
                "total_reads": tot,
                "unmodified_pct": pct["Unmodified"],
                "total_indel_pct": total_indel_pct,
                "frameshift_indel_pct": pct["Frameshift_Indel"],
                "in_frame_indel_pct": pct["In_Frame_Indel"],
                "substitution_only_pct": pct["Substitution_Only"],
                "hdr_pct": pct["HDR"],
                "imperfect_hdr_pct": pct["Imperfect_HDR"],
                "max_be_target_conversion_pct": max_target_conv,
                "mean_be_target_conversion_pct": mean_target_conv,
                "be_bystander_edit_pct": bystander_max,
                "indel_byproduct_pct": total_indel_pct,
            }
        )

    sum_df = pd.DataFrame(sample_summaries)
    indel_df = pd.DataFrame(indel_rows, columns=["sample", "indel_bp", "read_count", "pct_reads"])
    bm_df = pd.DataFrame(base_matrix_rows)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sum_path = out_dir / "editing_allele_summary.tsv"
    ind_path = out_dir / "indel_size_distribution.tsv"
    bm_path = out_dir / "protospacer_base_matrix.tsv"
    fig_path = out_dir / "editing_outcomes.png"

    sum_df.to_csv(sum_path, sep="\t", index=False)
    indel_df.to_csv(ind_path, sep="\t", index=False)
    bm_df.to_csv(bm_path, sep="\t", index=False)

    write_artifact(
        sum_path,
        "crispr_editing_quant",
        "CRISPR Post-Editing Outcome Summary (Frameshift, HDR, Base Editing)",
        inputs=[src_path],
        parameters={"editor_mode": args.editor_mode, "window_bp": args.window_bp},
        shape=(int(len(sum_df)), int(sum_df.shape[1])),
    )
    write_artifact(
        ind_path,
        "crispr_editing_quant",
        "CRISPR Amplicon Indel Size Distribution",
        inputs=[src_path],
        shape=(int(len(indel_df)), int(indel_df.shape[1])),
    )
    write_artifact(
        bm_path,
        "crispr_editing_quant",
        "Protospacer Position-by-Base Conversion Matrix",
        inputs=[src_path],
        shape=(int(len(bm_df)), int(bm_df.shape[1])),
    )

    # 3-panel plot
    plt = setup_matplotlib()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15.0, 4.5), dpi=180)

    # Panel A: Stacked outcome bar
    x = np.arange(len(sum_df))
    outcome_cols = [
        ("unmodified_pct", "Unmodified", "#cbd5e1"),
        ("frameshift_indel_pct", "Frameshift KO", "#2563eb"),
        ("in_frame_indel_pct", "In-frame Indel", "#60a5fa"),
        ("substitution_only_pct", "Substitution", "#f59e0b"),
        ("hdr_pct", "Clean HDR", "#059669"),
        ("imperfect_hdr_pct", "Imperfect HDR", "#34d399"),
    ]
    bottom = np.zeros(len(sum_df), dtype=float)
    for col_k, lbl, col_hex in outcome_cols:
        vals = sum_df[col_k].to_numpy(dtype=float)
        if np.any(vals > 0):
            ax1.bar(x, vals, bottom=bottom, label=lbl, color=col_hex, width=0.55)
            bottom += vals
    ax1.set_xticks(x)
    ax1.set_xticklabels(sum_df["sample"], rotation=20, ha="right", fontsize=8)
    ax1.set_ylim(0, 105)
    ax1.set_ylabel("Reads (%)")
    ax1.set_title("Editing Outcome Composition")
    ax1.legend(frameon=False, fontsize=7, loc="upper right")

    # Panel B: Indel spectrum for first sample
    first_smp = sample_summaries[0]["sample"]
    sub_ind = indel_df[indel_df["sample"] == first_smp]
    if not sub_ind.empty:
        colors_ind = ["#dc2626" if abs(v) % 3 != 0 else "#64748b" for v in sub_ind["indel_bp"]]
        ax2.bar(sub_ind["indel_bp"], sub_ind["pct_reads"], color=colors_ind, width=0.8)
    ax2.axvline(0, color="#0f172a", ls=":", lw=0.9)
    ax2.set_xlabel("Indel Size (bp; red = frameshift)")
    ax2.set_ylabel("Reads (%)")
    ax2.set_title(f"Indel Spectrum ({first_smp})")

    # Panel C: Protospacer base-conversion heatmap (4 x L) for first sample
    sub_bm = bm_df[bm_df["sample"] == first_smp].sort_values("position_1based")
    heat = np.vstack(
        [
            sub_bm["A_pct"].to_numpy(),
            sub_bm["C_pct"].to_numpy(),
            sub_bm["G_pct"].to_numpy(),
            sub_bm["T_pct"].to_numpy(),
        ]
    )
    im = ax3.imshow(heat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=100)
    ax3.set_yticks([0, 1, 2, 3])
    ax3.set_yticklabels(["A", "C", "G", "T"])
    ax3.set_xticks(np.arange(L))
    ax3.set_xticklabels([f"{p}\n{guide[p-1]}" for p in range(1, L + 1)], fontsize=6.2)
    ax3.axvspan(2.5, 7.5, edgecolor="#2563eb", facecolor="none", lw=1.4, ls="--")
    ax3.set_xlabel("Protospacer Position (5'->3'; [4..8] dashed)")
    ax3.set_title(f"Protospacer Base Matrix ({args.editor_mode})")
    fig.colorbar(im, ax=ax3, fraction=0.046, pad=0.04, label="% Reads")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "crispr_editing_quant",
        f"CRISPR Editing Outcomes ({args.editor_mode})",
        inputs=[src_path],
        parameters={"editor_mode": args.editor_mode},
    )

    emit(
        {
            "status": "ok",
            "tool": "crispr_editing_quant",
            "input_mode": input_mode,
            "editor_mode": args.editor_mode,
            "guide_strand": g_strand,
            "cut_site_bp": int(cut_site),
            "n_samples": len(sample_summaries),
            "total_reads": int(sum(r["total_reads"] for r in sample_summaries)),
            "samples": sample_summaries,
            "outputs": {
                "editing_allele_summary": str(sum_path),
                "indel_size_distribution": str(ind_path),
                "protospacer_base_matrix": str(bm_path),
                "outcomes_plot": str(fig_path),
            },
            "artifacts": [str(sum_path), str(ind_path), str(bm_path), str(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
