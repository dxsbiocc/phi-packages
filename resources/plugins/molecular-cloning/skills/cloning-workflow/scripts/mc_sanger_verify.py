#!/usr/bin/env python3
"""Verify plasmid/amplicon constructs from Sanger .ab1, .fastq, or .fasta traces."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    COMPLEMENT,
    emit,
    fail,
    mott_trim,
    parse_list_arg,
    parse_sequences,
    parse_trace_file,
    revcomp,
    setup_matplotlib,
    write_artifact,
)


def _seed_align_window(read_seq: str, ref_seq: str, is_circular: bool) -> tuple[str, int, int]:
    """Find best strand ('+' or '-') and reference window [start, end) via k-mer seeding."""
    search_ref = (ref_seq + ref_seq) if (is_circular and len(read_seq) <= len(ref_seq)) else ref_seq
    k = 11 if len(read_seq) >= 35 else 7
    index: dict[str, list[int]] = {}
    for pos in range(0, max(0, len(search_ref) - k + 1)):
        km = search_ref[pos : pos + k]
        index.setdefault(km, []).append(pos)

    best_strand = "+"
    best_votes = -1
    best_offset = 0
    for strand, q in (("+", read_seq), ("-", revcomp(read_seq))):
        offset_counts: dict[int, int] = {}
        for qi in range(0, max(0, len(q) - k + 1)):
            km = q[qi : qi + k]
            for rpos in index.get(km, ()):
                off = rpos - qi
                # Bin offsets by 4 bp to tolerate small indels
                b = round(off / 4) * 4
                cnt = offset_counts.get(b, 0) + 1
                offset_counts[b] = cnt
                if cnt > best_votes:
                    best_votes = cnt
                    best_strand = strand
                    best_offset = off

    if best_votes <= 0:
        return "+", 0, len(search_ref)
    w_start = max(0, best_offset - 35)
    w_end = min(len(search_ref), best_offset + len(read_seq) + 35)
    return best_strand, w_start, w_end


def _semiglobal_align(ref_win: str, query: str) -> tuple[str, str, int]:
    """Affine-gap semi-global alignment (free leading/trailing gaps in ref_win)."""
    n = len(ref_win)
    m = len(query)
    match_s = 2
    mismatch_s = -4
    gap_open = -6
    gap_ext = -2

    neg_inf = -10**9
    M = np.full((m + 1, n + 1), neg_inf, dtype=np.int32)
    Ix = np.full((m + 1, n + 1), neg_inf, dtype=np.int32)  # insertion in query (gap in ref)
    Iy = np.full((m + 1, n + 1), neg_inf, dtype=np.int32)  # deletion from ref (gap in query)
    tb = np.zeros((m + 1, n + 1), dtype=np.int8)  # 0:M, 1:Ix, 2:Iy

    # Free start column in reference
    M[0, :] = 0
    for i in range(1, m + 1):
        Ix[i, 0] = gap_open + (i - 1) * gap_ext
        tb[i, 0] = 1

    for i in range(1, m + 1):
        qi = query[i - 1]
        for j in range(1, n + 1):
            rj = ref_win[j - 1]
            sc = match_s if qi == rj else mismatch_s
            prev_m = int(M[i - 1, j - 1])
            prev_x = int(Ix[i - 1, j - 1])
            prev_y = int(Iy[i - 1, j - 1])
            best_prev = max(prev_m, prev_x, prev_y)
            M[i, j] = best_prev + sc

            open_x = int(M[i - 1, j]) + gap_open
            ext_x = int(Ix[i - 1, j]) + gap_ext
            Ix[i, j] = open_x if open_x >= ext_x else ext_x

            open_y = int(M[i, j - 1]) + gap_open
            ext_y = int(Iy[i, j - 1]) + gap_ext
            Iy[i, j] = open_y if open_y >= ext_y else ext_y

            cur_best = max(int(M[i, j]), int(Ix[i, j]), int(Iy[i, j]))
            if cur_best == int(M[i, j]):
                tb[i, j] = 0
            elif cur_best == int(Ix[i, j]):
                tb[i, j] = 1
            else:
                tb[i, j] = 2

    # Free end column in reference
    end_scores = np.maximum(M[m, :], np.maximum(Ix[m, :], Iy[m, :]))
    j = int(np.argmax(end_scores))
    i = m
    aln_r: list[str] = []
    aln_q: list[str] = []
    state = int(tb[i, j])
    while i > 0:
        if state == 0:
            aln_r.append(ref_win[j - 1])
            aln_q.append(query[i - 1])
            prev_m, prev_x, prev_y = int(M[i - 1, j - 1]), int(Ix[i - 1, j - 1]), int(Iy[i - 1, j - 1])
            bp = max(prev_m, prev_x, prev_y)
            state = 0 if bp == prev_m else (1 if bp == prev_x else 2)
            i -= 1
            j -= 1
        elif state == 1:
            aln_r.append("-")
            aln_q.append(query[i - 1])
            open_x = int(M[i - 1, j]) + gap_open
            ext_x = int(Ix[i - 1, j]) + gap_ext
            state = 0 if open_x >= ext_x else 1
            i -= 1
        else:
            aln_r.append(ref_win[j - 1])
            aln_q.append("-")
            open_y = int(M[i, j - 1]) + gap_open
            ext_y = int(Iy[i, j - 1]) + gap_ext
            state = 0 if open_y >= ext_y else 2
            j -= 1

    return "".join(reversed(aln_r)), "".join(reversed(aln_q)), j


def _analyze_trace(
    trace: dict,
    ref_seq: str,
    is_circular: bool,
    trim_err: float,
    min_sec_ratio: float,
) -> tuple[dict, list[dict], list[tuple[int, int]], dict | None]:
    raw_len = len(trace["seq"])
    ts, te = mott_trim(trace["quals"], trim_err)
    if te - ts < 20:
        # Fallback if all bases low quality: keep full read and flag
        ts, te = 0, raw_len
        trim_warn = True
    else:
        trim_warn = False

    trimmed_seq = trace["seq"][ts:te]
    trimmed_quals = trace["quals"][ts:te]
    trimmed_ploc = trace["ploc"][ts:te] if trace["ploc"] else []

    strand, w_start, w_end = _seed_align_window(trimmed_seq, ref_seq, is_circular)
    search_ref = (ref_seq + ref_seq) if (is_circular and len(trimmed_seq) <= len(ref_seq)) else ref_seq
    ref_win = search_ref[w_start:w_end]
    q_oriented = trimmed_seq if strand == "+" else revcomp(trimmed_seq)
    q_quals_or = trimmed_quals if strand == "+" else list(reversed(trimmed_quals))
    q_ploc_or = trimmed_ploc if strand == "+" else list(reversed(trimmed_ploc))

    aln_r, aln_q, win_off = _semiglobal_align(ref_win, q_oriented)
    ref_len = len(ref_seq)
    ref_pos_0 = w_start + win_off
    q_idx = 0

    matches = 0
    aligned_cols = 0
    covered_intervals: list[tuple[int, int]] = []
    variants: list[dict] = []
    mixed_peaks: list[dict] = []

    first_ref_1 = (ref_pos_0 % ref_len) + 1
    last_ref_1 = first_ref_1

    i = 0
    while i < len(aln_r):
        rb = aln_r[i]
        qb = aln_q[i]
        if rb != "-" and qb != "-":
            cur_ref_1 = (ref_pos_0 % ref_len) + 1
            covered_intervals.append((cur_ref_1, cur_ref_1))
            last_ref_1 = cur_ref_1
            aligned_cols += 1
            if rb == qb:
                matches += 1
            else:
                variants.append({
                    "trace_name": trace["name"],
                    "strand": strand,
                    "ref_pos_1": cur_ref_1,
                    "variant_type": "SNV",
                    "ref_base": rb,
                    "alt_base": qb,
                    "indel_len": 0,
                    "frameshift": False,
                    "phred_q": q_quals_or[q_idx] if q_idx < len(q_quals_or) else 30,
                    "secondary_base": "",
                    "secondary_ratio": 0.0,
                })

            # Check ABIF secondary fluorescence peak at PLOC
            if trace["traces"] and q_idx < len(q_ploc_or):
                p_loc = q_ploc_or[q_idx]
                intens: list[tuple[str, float]] = []
                for b in ("A", "C", "G", "T"):
                    arr = trace["traces"].get(b, [])
                    if arr:
                        lo = max(0, p_loc - 1)
                        hi = min(len(arr), p_loc + 2)
                        val = float(max(arr[lo:hi])) if hi > lo else 0.0
                    else:
                        val = 0.0
                    out_b = b if strand == "+" else b.translate(COMPLEMENT)
                    intens.append((out_b, val))
                intens.sort(key=lambda x: -x[1])
                b1, i1 = intens[0]
                b2, i2 = intens[1]
                ratio = i2 / (i1 + i2) if (i1 + i2) > 0 else 0.0
                if ratio >= min_sec_ratio:
                    mp = {
                        "trace_name": trace["name"],
                        "strand": strand,
                        "ref_pos_1": cur_ref_1,
                        "variant_type": "MIXED_PEAK",
                        "ref_base": rb,
                        "alt_base": f"{b1}/{b2}",
                        "indel_len": 0,
                        "frameshift": False,
                        "phred_q": q_quals_or[q_idx] if q_idx < len(q_quals_or) else 30,
                        "secondary_base": b2,
                        "secondary_ratio": round(ratio, 3),
                        "ploc": p_loc,
                    }
                    mixed_peaks.append(mp)
                    variants.append({k: v for k, v in mp.items() if k != "ploc"})
            ref_pos_0 += 1
            q_idx += 1
            i += 1
        elif rb == "-":
            # Insertion in read
            j = i
            ins_bases = []
            while j < len(aln_r) and aln_r[j] == "-":
                ins_bases.append(aln_q[j])
                q_idx += 1
                j += 1
            cur_ref_1 = (ref_pos_0 % ref_len) + 1
            ilen = len(ins_bases)
            aligned_cols += ilen
            variants.append({
                "trace_name": trace["name"],
                "strand": strand,
                "ref_pos_1": cur_ref_1,
                "variant_type": "INS",
                "ref_base": "-",
                "alt_base": "".join(ins_bases),
                "indel_len": ilen,
                "frameshift": (ilen % 3 != 0),
                "phred_q": 30,
                "secondary_base": "",
                "secondary_ratio": 0.0,
            })
            i = j
        else:
            # Deletion from reference
            j = i
            del_bases = []
            cur_ref_1 = (ref_pos_0 % ref_len) + 1
            while j < len(aln_r) and aln_q[j] == "-":
                del_bases.append(aln_r[j])
                ref_pos_0 += 1
                j += 1
            dlen = len(del_bases)
            aligned_cols += dlen
            variants.append({
                "trace_name": trace["name"],
                "strand": strand,
                "ref_pos_1": cur_ref_1,
                "variant_type": "DEL",
                "ref_base": "".join(del_bases),
                "alt_base": "-",
                "indel_len": dlen,
                "frameshift": (dlen % 3 != 0),
                "phred_q": 30,
                "secondary_base": "",
                "secondary_ratio": 0.0,
            })
            i = j

    identity_pct = round(100.0 * matches / max(1, aligned_cols), 2)
    n_snv_indel = sum(1 for v in variants if v["variant_type"] != "MIXED_PEAK")
    n_mixed = len(mixed_peaks)
    mean_q = round(float(np.mean(trimmed_quals)) if trimmed_quals else 0.0, 1)
    passed = (n_snv_indel == 0 and n_mixed == 0 and identity_pct >= 99.5 and not trim_warn)

    summary = {
        "trace_name": trace["name"],
        "format": trace["format"],
        "strand": strand,
        "raw_length": raw_len,
        "trimmed_start_1": ts + 1,
        "trimmed_end_1": te,
        "trimmed_length": te - ts,
        "mean_phred_q": mean_q,
        "ref_start_1": first_ref_1,
        "ref_end_1": last_ref_1,
        "identity_pct": identity_pct,
        "n_variants": n_snv_indel,
        "n_mixed_peaks": n_mixed,
        "pass": bool(passed),
    }

    zoom_ctx = None
    if trace["traces"] and trace["ploc"]:
        # Center chromatogram zoom around first mixed peak or variant, else midpoint
        focus_ploc = mixed_peaks[0]["ploc"] if mixed_peaks else trace["ploc"][len(trace["ploc"]) // 2]
        zoom_ctx = {
            "trace_name": trace["name"],
            "traces": trace["traces"],
            "ploc": trace["ploc"],
            "seq": trace["seq"],
            "focus_ploc": int(focus_ploc),
            "focus_label": (
                f"Mixed peak at ref {mixed_peaks[0]['ref_pos_1']} ({mixed_peaks[0]['alt_base']}, ratio={mixed_peaks[0]['secondary_ratio']})"
                if mixed_peaks
                else f"ABIF trace zoom ({trace['name']})"
            ),
        }

    return summary, variants, covered_intervals, zoom_ctx


def _plot_chromatogram(
    out_png: Path,
    ref_id: str,
    ref_len: int,
    cov_vec: np.ndarray,
    summaries: list[dict],
    var_df: pd.DataFrame,
    zoom_ctx: dict | None,
) -> None:
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11.0, 6.2), gridspec_kw={"height_ratios": [1.0, 1.35]})

    x = np.arange(1, ref_len + 1)
    ax1.fill_between(x, 0, cov_vec, color="#99f6e4", alpha=0.8, label="Trace coverage depth")
    ax1.plot(x, cov_vec, color="#0f766e", lw=1.2)
    if not var_df.empty:
        for _, vr in var_df.iterrows():
            col = "#f59e0b" if vr["variant_type"] == "MIXED_PEAK" else "#dc2626"
            ax1.axvline(int(vr["ref_pos_1"]), color=col, linestyle="--", lw=1.3, alpha=0.85)
    ax1.set_xlim(1, max(1, ref_len))
    ax1.set_ylim(0, max(1.5, float(np.max(cov_vec)) + 0.5))
    ax1.set_xlabel("Reference Coordinate (bp)")
    ax1.set_ylabel("Coverage Depth")
    ax1.set_title(f"Sanger Verification Coverage & Variants — {ref_id} ({len(summaries)} traces)")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    base_colors = {"A": "#16a34a", "C": "#2563eb", "G": "#1e293b", "T": "#dc2626"}
    if zoom_ctx is not None:
        center = zoom_ctx["focus_ploc"]
        lo = max(0, center - 120)
        hi = center + 120
        for b, col in base_colors.items():
            arr = zoom_ctx["traces"].get(b, [])
            if arr:
                hi_b = min(len(arr), hi)
                xs = np.arange(lo, hi_b)
                ax2.plot(xs, arr[lo:hi_b], color=col, lw=1.4, label=b)
        max_y = max((max(zoom_ctx["traces"].get(b, [100])[lo:hi] or [100]) for b in base_colors), default=100)
        for idx, pl in enumerate(zoom_ctx["ploc"]):
            if lo <= pl < hi:
                b = zoom_ctx["seq"][idx]
                ax2.text(pl, max_y * 1.04, b, ha="center", va="bottom", fontsize=7.5, color=base_colors.get(b, "#334155"), fontweight="bold")
        ax2.axvline(center, color="#f59e0b", linestyle=":", lw=1.5)
        ax2.set_xlim(lo, hi)
        ax2.set_ylim(0, max_y * 1.18)
        ax2.set_xlabel("ABIF Trace Sample Index")
        ax2.set_ylabel("Fluorescence Intensity")
        ax2.set_title(zoom_ctx["focus_label"])
        ax2.legend(frameon=False, ncol=4, loc="upper right", fontsize=8.5)
    else:
        ax2.text(
            0.5,
            0.5,
            "No binary .ab1 fluorescence trace channels present (FASTQ/FASTA input).",
            ha="center",
            va="center",
            fontsize=10,
            color="#475569",
        )
        ax2.set_xticks([])
        ax2.set_yticks([])
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify construct from Sanger .ab1/FASTQ/FASTA traces.")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--traces", required=True, action="append")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--trim_error_prob", type=float, default=0.05)
    parser.add_argument("--min_secondary_peak_ratio", type=float, default=0.15)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    refs = parse_sequences(args.reference)
    ref = refs[0]
    ref_seq = ref["seq"]
    ref_len = len(ref_seq)
    is_circular = ref.get("topology", "linear") == "circular" or ref_len >= 1200

    trace_paths = parse_list_arg(args.traces, [])
    if not trace_paths:
        fail("At least one Sanger trace path must be supplied via --traces.")

    cov_vec = np.zeros(ref_len, dtype=np.int32)
    summaries: list[dict] = []
    all_variants: list[dict] = []
    chosen_zoom: dict | None = None

    for tp in trace_paths:
        t_obj = parse_trace_file(tp)
        summ, vars_list, cov_pts, z_ctx = _analyze_trace(
            t_obj,
            ref_seq,
            is_circular,
            args.trim_error_prob,
            args.min_secondary_peak_ratio,
        )
        summaries.append(summ)
        all_variants.extend(vars_list)
        for s1, _ in cov_pts:
            cov_vec[s1 - 1] += 1
        if z_ctx is not None and (chosen_zoom is None or summ["n_mixed_peaks"] > 0):
            chosen_zoom = z_ctx

    sum_df = pd.DataFrame(summaries)
    sum_path = out_dir / "sanger_verification_summary.tsv"
    sum_df.to_csv(sum_path, sep="\t", index=False)

    var_cols = [
        "trace_name", "strand", "ref_pos_1", "variant_type", "ref_base",
        "alt_base", "indel_len", "frameshift", "phred_q", "secondary_base", "secondary_ratio",
    ]
    var_df = pd.DataFrame(all_variants, columns=var_cols) if all_variants else pd.DataFrame(columns=var_cols)
    var_path = out_dir / "sanger_variants.tsv"
    var_df.to_csv(var_path, sep="\t", index=False)

    png_path = out_dir / "sanger_alignment_chromatogram.png"
    _plot_chromatogram(png_path, ref["id"], ref_len, cov_vec, summaries, var_df, chosen_zoom)

    cov_pct = round(100.0 * float(np.count_nonzero(cov_vec)) / max(1, ref_len), 2)
    n_passed = int(sum(1 for s in summaries if s["pass"]))
    total_snv_indel = int(sum(s["n_variants"] for s in summaries))
    total_mixed = int(sum(s["n_mixed_peaks"] for s in summaries))
    overall_pass = (n_passed == len(summaries)) and (total_snv_indel == 0) and (total_mixed == 0)

    warnings: list[str] = []
    if total_snv_indel > 0:
        warnings.append(f"Detected {total_snv_indel} SNV/indel variant(s) across Sanger traces.")
    if total_mixed > 0:
        warnings.append(f"Detected {total_mixed} secondary chromatogram peak(s) >= {args.min_secondary_peak_ratio}; inspect for mixed colony or heterozygous edit.")

    tool_name = "mc_sanger_verify"
    params = {"trim_error_prob": args.trim_error_prob, "min_secondary_peak_ratio": args.min_secondary_peak_ratio}
    inputs = [args.reference, *trace_paths]
    artifacts = [
        write_artifact(sum_path, tool_name, "Sanger verification summary", inputs=inputs, parameters=params, shape=sum_df.shape),
        write_artifact(var_path, tool_name, "Sanger variants & mixed peaks", inputs=inputs, parameters=params, shape=var_df.shape),
        write_artifact(png_path, tool_name, "Sanger alignment & chromatogram", inputs=inputs, parameters=params),
    ]

    emit({
        "status": "ok",
        "tool": tool_name,
        "reference_id": ref["id"],
        "reference_length_bp": ref_len,
        "n_traces": len(summaries),
        "n_passed_traces": n_passed,
        "reference_coverage_pct": cov_pct,
        "overall_pass": overall_pass,
        "total_variants": total_snv_indel,
        "total_mixed_peaks": total_mixed,
        "trace_summaries": summaries,
        "outputs": {
            "verification_summary": str(sum_path),
            "variants_table": str(var_path),
            "chromatogram_plot": str(png_path),
        },
        "artifacts": artifacts,
        "warnings": warnings,
    })


if __name__ == "__main__":
    main()
