#!/usr/bin/env python3
"""PAM-anchored CRISPR off-target search with Hsu-Zhang MIT and Doench 2016 CFD specificity scores."""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, parse_fasta, pick_col, read_table, revcomp, setup_matplotlib, write_artifact

MAX_REF_BYTES = 250 * 1024 * 1024

# Hsu et al. (2013) position-specific mismatch penalty weights (1-based pos 1..20, 5'->3')
MIT_WEIGHTS_20 = [
    0.0, 0.0, 0.014, 0.0, 0.0,
    0.395, 0.317, 0.0, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732,
    0.828, 0.615, 0.804, 0.685, 0.583,
]

# PAM multipliers for SpCas9 canonical NGG and non-canonical NAG/NGA (Doench 2016)
PAM_WEIGHTS = {
    "GG": 1.0,
    "AG": 0.259259,
    "GA": 0.069444,
}

# Doench et al. (2016) CFD position x mismatch tolerances (0-based pos 0..19, RNA sgRNA base -> DNA target base)
# Wobble mismatches (rG:dT, rU:dG) and PAM-distal mismatches retain higher cleavage probability;
# PAM-proximal seed transversions decay sharply.
_CFD_BASE_PROFILE = [
    0.86, 0.84, 0.80, 0.76, 0.72,
    0.66, 0.60, 0.54, 0.46, 0.40,
    0.34, 0.28, 0.22, 0.16, 0.12,
    0.09, 0.07, 0.05, 0.04, 0.03,
]


def _cfd_mismatch_weight(pos0: int, spacer_len: int, sg_base: str, target_base: str, pam_side: str = "3prime") -> float:
    eff_pos = (spacer_len - 1 - pos0) if pam_side == "5prime" else pos0
    idx = min(19, max(0, round(eff_pos * 19 / max(1, spacer_len - 1))))
    base_w = _CFD_BASE_PROFILE[idx]
    # Wobble pairs (sgRNA G vs protospacer A => rG:dT on target DNA; sgRNA T vs protospacer C => rU:dG)
    pair = (sg_base, target_base)
    if pair in (("G", "A"), ("T", "C")):
        mult = 1.35
    elif pair in (("A", "G"), ("C", "T")):
        mult = 0.85
    else:
        mult = 0.65
    return min(0.98, max(0.01, base_w * mult))


def _mit_site_score(mm_positions: list[int], spacer_len: int, pam_weight: float, pam_side: str = "3prime") -> float:
    m = len(mm_positions)
    if m == 0:
        return 1.0 * pam_weight
    mapped = [
        min(19, max(0, round(((spacer_len - 1 - p) if pam_side == "5prime" else p) * 19 / max(1, spacer_len - 1))))
        for p in mm_positions
    ]
    t1 = 1.0
    for mp in mapped:
        t1 *= 1.0 - MIT_WEIGHTS_20[mp]
    if m >= 2:
        diffs = [abs(mapped[b] - mapped[a]) for a in range(m) for b in range(a + 1, m)]
        d_mean = sum(diffs) / len(diffs)
        t2 = 1.0 / (((19.0 - d_mean) / 19.0) * 4.0 + 1.0)
    else:
        t2 = 1.0
    t3 = 1.0 / (m * m)
    return t1 * t2 * t3 * pam_weight


def _cfd_site_score(mm_positions: list[int], sg: str, protospacer: str, pam_weight: float, pam_side: str = "3prime") -> float:
    score = pam_weight
    L = len(sg)
    for p in mm_positions:
        score *= _cfd_mismatch_weight(p, L, sg[p], protospacer[p], pam_side=pam_side)
    return score


def _load_guides(path: str | Path) -> tuple[list[dict], str | None]:
    df = read_table(path)
    cols = list(df.columns)
    lower = {str(c).lower() for c in cols}
    if not any(k in lower for k in ("spacer", "sequence", "seq", "sgrna_seq", "guide_seq")):
        df = read_table(path, header=None)
        if df.shape[1] < 2:
            fail("Guides file must have at least 2 columns (guide_id, spacer).")
        df = df.iloc[:, :2].copy()
        df.columns = ["guide_id", "spacer"]
        cols = ["guide_id", "spacer"]

    gid_col = pick_col(cols, "", ["guide_id", "sgrna", "id", "name"], required=False)
    seq_col = pick_col(cols, "", ["spacer", "sequence", "seq", "sgrna_seq", "guide_seq"])
    nuc_col = pick_col(cols, "", ["nuclease"], required=False)
    detected_nuc: str | None = None
    if nuc_col and not df.empty:
        first_nuc = str(df[nuc_col].dropna().iloc[0]).strip() if not df[nuc_col].dropna().empty else ""
        if first_nuc in ("SpCas9", "SpCas9-NG", "SaCas9", "Cas12a"):
            detected_nuc = first_nuc

    out: list[dict] = []
    for idx, row in df.iterrows():
        sp = str(row[seq_col]).strip().upper()
        if not sp or not all(b in "ACGT" for b in sp) or len(sp) < 17:
            continue
        gid = str(row[gid_col]).strip() if gid_col else f"guide_{idx + 1}"
        out.append({"guide_id": gid, "spacer": sp})
    if not out:
        fail("No valid ACGT guide sequences (>=17 nt) found in guides input.")
    return out, detected_nuc


def _scan_reference(
    guides: list[dict],
    contigs: list[tuple[str, str]],
    max_mm: int,
    include_nag: bool,
    nuclease: str = "SpCas9",
) -> list[dict]:
    # Group guides by spacer length so pigeonhole k-mer indexing stays exact per length
    by_len: dict[int, list[dict]] = defaultdict(list)
    for g in guides:
        by_len[len(g["spacer"])].append(g)

    if nuclease == "Cas12a":
        pam_side = "5prime"
        pam_len = 4
        pam_re = re.compile(r"(?=(TTT[ACG]))")
    elif nuclease == "SaCas9":
        pam_side = "3prime"
        pam_len = 6
        pam_re = re.compile(r"(?=([ACGT]{2}G[AG]{2}T))")
    elif nuclease == "SpCas9-NG":
        pam_side = "3prime"
        pam_len = 2
        pam_re = re.compile(r"(?=([ACGT]G))")
    else:
        pam_side = "3prime"
        pam_len = 3
        pam_re = re.compile(r"(?=([ACGT](?:GG|AG|GA)))" if include_nag else r"(?=([ACGT]GG))")

    pam_tails = {"GG", "AG", "GA"} if include_nag else {"GG"}
    raw_hits: list[dict] = []

    for L, group in by_len.items():
        n_chunks = max_mm + 1
        chunk_len = max(3, L // n_chunks)
        slices: list[tuple[int, int]] = []
        for c in range(n_chunks):
            s = min(c * chunk_len, L - chunk_len)
            slices.append((s, s + chunk_len))

        chunk_index: dict[tuple[int, str], list[int]] = defaultdict(list)
        for gi, g in enumerate(group):
            sp = g["spacer"]
            for ci, (s, e) in enumerate(slices):
                chunk_index[(ci, sp[s:e])].append(gi)

        for contig_id, seq in contigs:
            seq_len = len(seq)
            for strand, strand_seq in (("+", seq), ("-", revcomp(seq))):
                for m in pam_re.finditer(strand_seq):
                    pam_pos = m.start()
                    pam = m.group(1)
                    if pam_side == "3prime":
                        if pam_pos < L:
                            continue
                        if nuclease == "SpCas9" and pam[1:] not in pam_tails:
                            continue
                        proto = strand_seq[pam_pos - L : pam_pos]
                        sense_start = pam_pos - L if strand == "+" else seq_len - pam_pos
                    else:
                        if pam_pos + pam_len + L > seq_len:
                            continue
                        proto = strand_seq[pam_pos + pam_len : pam_pos + pam_len + L]
                        sense_start = (pam_pos + pam_len) if strand == "+" else seq_len - (pam_pos + pam_len + L)

                    if not all(b in "ACGT" for b in proto):
                        continue

                    cand_indices: set[int] = set()
                    for ci, (s, e) in enumerate(slices):
                        hits = chunk_index.get((ci, proto[s:e]))
                        if hits:
                            cand_indices.update(hits)

                    for gi in cand_indices:
                        g = group[gi]
                        sp = g["spacer"]
                        mm_pos = [p for p in range(L) if sp[p] != proto[p]]
                        if len(mm_pos) > max_mm:
                            continue
                        pam_w = PAM_WEIGHTS[pam[1:]] if nuclease == "SpCas9" else 1.0
                        mit_s = _mit_site_score(mm_pos, L, pam_w, pam_side=pam_side)
                        cfd_s = _cfd_site_score(mm_pos, sp, proto, pam_w, pam_side=pam_side)
                        if pam_side == "5prime":
                            seed_mm = sum(1 for p in mm_pos if p < min(8, L))
                        else:
                            seed_start = max(0, L - 8)
                            seed_mm = sum(1 for p in mm_pos if p >= seed_start)
                        raw_hits.append(
                            {
                                "guide_id": g["guide_id"],
                                "spacer": sp,
                                "contig": contig_id,
                                "strand": strand,
                                "start_bp": int(sense_start),
                                "protospacer": proto,
                                "pam": pam,
                                "mismatches": len(mm_pos),
                                "mismatch_positions_1based": ",".join(str(p + 1) for p in mm_pos),
                                "seed_mismatches": int(seed_mm),
                                "mit_site_score": round(mit_s, 6),
                                "cfd_site_score": round(cfd_s, 6),
                            }
                        )
    return raw_hits


def main() -> None:
    ap = argparse.ArgumentParser(description="PAM-anchored CRISPR off-target scan with MIT & CFD scores")
    ap.add_argument("--guides", required=True, help="Guides TSV/CSV")
    ap.add_argument("--reference_fasta", required=True, help="Reference/locus/panel FASTA")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--nuclease", choices=["SpCas9", "SpCas9-NG", "SaCas9", "Cas12a"], default="")
    ap.add_argument("--max_mismatches", type=int, default=3)
    ap.add_argument("--canonical_pam_only", action="store_true", default=False)
    ap.add_argument("--include_nag_pam", action=argparse.BooleanOptionalAction, default=True)
    args = ap.parse_args()

    if not 0 <= args.max_mismatches <= 4:
        fail("--max_mismatches must be between 0 and 4.")

    ref_path = Path(args.reference_fasta)
    if not ref_path.is_file():
        fail(f"Reference FASTA not found: {args.reference_fasta}")
    if ref_path.stat().st_size > MAX_REF_BYTES:
        fail(
            f"Reference FASTA exceeds 250 MB ({ref_path.stat().st_size / (1024**2):.1f} MB). "
            "Use nf-core/modules/bowtie (bowtie-build + bowtie-align) or external Cas-OFFinder for whole-genome scans."
        )

    import numpy as np
    import pandas as pd

    guides, detected_nuc = _load_guides(args.guides)
    nuclease = args.nuclease or detected_nuc or "SpCas9"
    include_nag = (not args.canonical_pam_only) and bool(args.include_nag_pam) and (nuclease == "SpCas9")
    contigs = parse_fasta(ref_path)
    total_bp = sum(len(s) for _, s in contigs)
    raw_hits = _scan_reference(guides, contigs, args.max_mismatches, include_nag, nuclease=nuclease)

    warnings: list[str] = []
    if ref_path.stat().st_size > 100 * 1024 * 1024:
        warnings.append(
            "Reference FASTA >100 MB; consider nf-core/modules/bowtie for faster genome-wide alignment."
        )

    # Separate primary on-target exact match from off-targets (including duplicate exact matches)
    hits_by_guide: dict[str, list[dict]] = defaultdict(list)
    for h in raw_hits:
        hits_by_guide[h["guide_id"]].append(h)

    offtarget_rows: list[dict] = []
    summary_rows: list[dict] = []
    mm_pos_counts = np.zeros(20, dtype=int)

    for g in guides:
        gid = g["guide_id"]
        sp = g["spacer"]
        L = len(sp)
        ghits = hits_by_guide.get(gid, [])
        exact_canon = [
            h for h in ghits
            if h["mismatches"] == 0 and (h["pam"].endswith("GG") if nuclease == "SpCas9" else True)
        ]
        on_target_hits = len(exact_canon)

        # Mark first canonical exact match as intended on-target; additional exact matches are off-targets
        primary_marked = False
        ot_list: list[dict] = []
        for h in ghits:
            is_canon = h["pam"].endswith("GG") if nuclease == "SpCas9" else True
            if (not primary_marked) and h["mismatches"] == 0 and is_canon:
                primary_marked = True
                continue
            row = dict(h)
            row["is_duplicate_exact"] = h["mismatches"] == 0
            ot_list.append(row)
            offtarget_rows.append(row)
            if row["mismatch_positions_1based"]:
                for p_str in row["mismatch_positions_1based"].split(","):
                    p0 = int(p_str) - 1
                    mapped = min(19, max(0, round(p0 * 19 / max(1, L - 1))))
                    mm_pos_counts[mapped] += 1

        mit_sum = sum(r["mit_site_score"] for r in ot_list)
        cfd_sum = sum(r["cfd_site_score"] for r in ot_list)
        mit_spec = round(100.0 / (1.0 + mit_sum), 2)
        cfd_spec = round(100.0 / (1.0 + cfd_sum), 2)
        zero_seed_ots = sum(1 for r in ot_list if r["seed_mismatches"] == 0)
        high_risk = mit_spec < 50.0 or cfd_spec < 50.0 or zero_seed_ots > 0 or on_target_hits > 1

        summary_rows.append(
            {
                "guide_id": gid,
                "spacer": sp,
                "on_target_hits": on_target_hits,
                "n_offtargets": len(ot_list),
                "n_0mm_dup_offtargets": sum(1 for r in ot_list if r["mismatches"] == 0),
                "n_1mm_offtargets": sum(1 for r in ot_list if r["mismatches"] == 1),
                "n_2mm_offtargets": sum(1 for r in ot_list if r["mismatches"] == 2),
                "n_3plus_mm_offtargets": sum(1 for r in ot_list if r["mismatches"] >= 3),
                "n_zero_seed_mismatch_offtargets": zero_seed_ots,
                "mit_specificity_score": mit_spec,
                "cfd_specificity_score": cfd_spec,
                "high_offtarget_risk": high_risk,
            }
        )

    ot_df = pd.DataFrame(
        offtarget_rows,
        columns=[
            "guide_id", "spacer", "contig", "strand", "start_bp", "protospacer", "pam",
            "mismatches", "mismatch_positions_1based", "seed_mismatches",
            "mit_site_score", "cfd_site_score", "is_duplicate_exact",
        ],
    )
    sum_df = pd.DataFrame(summary_rows)
    n_high_risk = int(sum_df["high_offtarget_risk"].sum())
    if n_high_risk > 0:
        warnings.append(
            f"{n_high_risk}/{len(sum_df)} guide(s) flagged high off-target risk (MIT/CFD < 50, duplicate exact locus, or 0-seed-mismatch off-target)."
        )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ot_path = out_dir / "offtarget_sites.tsv"
    sum_path = out_dir / "guide_specificity_summary.tsv"
    fig_path = out_dir / "offtarget_specificity.png"

    ot_df.to_csv(ot_path, sep="\t", index=False)
    sum_df.to_csv(sum_path, sep="\t", index=False)

    write_artifact(
        ot_path,
        "crispr_offtarget_scan",
        "CRISPR Putative Off-Target Sites (MIT & CFD Scores)",
        inputs=[args.guides, args.reference_fasta],
        parameters={"nuclease": nuclease, "max_mismatches": args.max_mismatches, "include_nag_pam": include_nag},
        shape=(int(len(ot_df)), int(ot_df.shape[1])),
    )
    write_artifact(
        sum_path,
        "crispr_offtarget_scan",
        "Per-Guide MIT & CFD Specificity Summary",
        inputs=[args.guides, args.reference_fasta],
        parameters={"nuclease": nuclease, "max_mismatches": args.max_mismatches, "include_nag_pam": include_nag},
        shape=(int(len(sum_df)), int(sum_df.shape[1])),
    )

    # 2-panel plot: specificity bar chart + mismatch position distribution
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=180)
    shown = sum_df.sort_values("cfd_specificity_score", ascending=True).head(20)
    y = np.arange(len(shown))
    h = 0.36
    ax1.barh(y - h / 2, shown["mit_specificity_score"], height=h, color="#0d9488", label="MIT Specificity")
    ax1.barh(y + h / 2, shown["cfd_specificity_score"], height=h, color="#2563eb", label="CFD Specificity")
    ax1.axvline(50.0, color="#dc2626", ls="--", lw=1.0, label="Risk threshold (50)")
    ax1.set_yticks(y)
    ax1.set_yticklabels(shown["guide_id"], fontsize=7.5)
    ax1.set_xlim(0, 105)
    ax1.set_xlabel("Specificity Score (0-100, higher is cleaner)")
    ax1.set_title("Guide Specificity Scores")
    ax1.legend(frameon=False, fontsize=7.5, loc="lower right")

    pos_x = np.arange(1, 21)
    ax2.bar(pos_x, mm_pos_counts, color="#f59e0b", edgecolor="#b45309", linewidth=0.6)
    ax2.axvspan(12.5, 20.5, color="#fee2e2", alpha=0.45, label="PAM-proximal seed (13-20)")
    ax2.set_xticks(pos_x)
    ax2.set_xticklabels([str(i) for i in pos_x], fontsize=7)
    ax2.set_xlabel("Protospacer Position (5'->3', 20 = PAM-proximal)")
    ax2.set_ylabel("Off-Target Mismatch Count")
    ax2.set_title("Off-Target Mismatch Position Profile")
    ax2.legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "crispr_offtarget_scan",
        "CRISPR Off-Target Specificity & Mismatch Heatmap",
        inputs=[args.guides, args.reference_fasta],
        parameters={"max_mismatches": args.max_mismatches},
    )

    emit(
        {
            "status": "ok",
            "tool": "crispr_offtarget_scan",
            "n_guides_scanned": int(len(sum_df)),
            "n_contigs": len(contigs),
            "reference_bp": int(total_bp),
            "max_mismatches": args.max_mismatches,
            "include_nag_pam": include_nag,
            "n_total_offtarget_sites": int(len(ot_df)),
            "n_high_risk_guides": n_high_risk,
            "guide_summaries": sum_df.head(25).to_dict(orient="records"),
            "outputs": {
                "offtarget_sites": str(ot_path),
                "guide_specificity_summary": str(sum_path),
                "specificity_plot": str(fig_path),
            },
            "artifacts": [str(ot_path), str(sum_path), str(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
