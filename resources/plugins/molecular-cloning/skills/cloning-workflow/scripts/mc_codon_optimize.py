#!/usr/bin/env python3
"""Host-specific codon optimization preserving 100% AA identity & clearing cloning sites."""
from __future__ import annotations

import argparse
import math
import os
import random
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    CODON_TABLE,
    emit,
    fail,
    find_restriction_sites,
    gc_fraction,
    parse_list_arg,
    parse_sequences,
    setup_matplotlib,
    translate_dna,
    write_artifact,
    write_fasta,
    write_genbank,
)

# Relative codon frequencies per 1000 codons (curated from Kazusa / CoCoPUTs standard tables)
HOST_CODON_FREQ: dict[str, dict[str, float]] = {
    "E_coli_K12": {
        "TTT": 22.3, "TTC": 16.6, "TTA": 13.8, "TTG": 13.7,
        "CTT": 11.0, "CTC": 11.1, "CTA": 3.9,  "CTG": 52.6,
        "ATT": 30.3, "ATC": 25.1, "ATA": 4.4,  "ATG": 27.9,
        "GTT": 18.3, "GTC": 15.3, "GTA": 10.9, "GTG": 26.4,
        "TCT": 8.5,  "TCC": 8.6,  "TCA": 7.2,  "TCG": 8.9,
        "CCT": 7.0,  "CCC": 5.5,  "CCA": 8.4,  "CCG": 23.2,
        "ACT": 9.0,  "ACC": 23.4, "ACA": 7.1,  "ACG": 14.4,
        "GCT": 15.3, "GCC": 25.5, "GCA": 20.1, "GCG": 33.7,
        "TAT": 16.2, "TAC": 12.2, "TAA": 2.0,  "TAG": 0.2,
        "CAT": 12.9, "CAC": 9.7,  "CAA": 15.3, "CAG": 28.8,
        "AAT": 17.7, "AAC": 21.7, "AAA": 33.6, "AAG": 10.3,
        "GAT": 32.1, "GAC": 19.1, "GAA": 39.4, "GAG": 17.8,
        "TGT": 5.2,  "TGC": 6.5,  "TGA": 1.0,  "TGG": 15.2,
        "CGT": 20.9, "CGC": 22.0, "CGA": 3.6,  "CGG": 5.4,
        "AGT": 8.8,  "AGC": 16.1, "AGA": 2.1,  "AGG": 1.2,
        "GGT": 24.7, "GGC": 29.6, "GGA": 8.0,  "GGG": 11.1,
    },
    "HEK293_Human": {
        "TTT": 17.6, "TTC": 20.3, "TTA": 7.7,  "TTG": 12.9,
        "CTT": 13.2, "CTC": 19.6, "CTA": 7.2,  "CTG": 39.6,
        "ATT": 16.0, "ATC": 20.8, "ATA": 7.5,  "ATG": 22.0,
        "GTT": 11.0, "GTC": 14.5, "GTA": 7.1,  "GTG": 28.1,
        "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4,
        "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
        "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1,
        "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
        "TAT": 12.2, "TAC": 15.3, "TAA": 1.0,  "TAG": 0.8,
        "CAT": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2,
        "AAT": 17.0, "AAC": 19.1, "AAA": 24.4, "AAG": 31.9,
        "GAT": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6,
        "TGT": 10.6, "TGC": 12.6, "TGA": 1.6,  "TGG": 13.2,
        "CGT": 4.5,  "CGC": 10.4, "CGA": 6.2,  "CGG": 11.4,
        "AGT": 12.1, "AGC": 19.5, "AGA": 12.2, "AGG": 12.0,
        "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
    },
    "CHO_Cricetulus": {
        "TTT": 19.1, "TTC": 21.2, "TTA": 6.7,  "TTG": 13.6,
        "CTT": 13.4, "CTC": 19.8, "CTA": 6.9,  "CTG": 39.5,
        "ATT": 16.5, "ATC": 22.6, "ATA": 6.8,  "ATG": 22.7,
        "GTT": 11.1, "GTC": 15.3, "GTA": 7.2,  "GTG": 29.1,
        "TCT": 15.4, "TCC": 17.8, "TCA": 11.3, "TCG": 4.2,
        "CCT": 17.9, "CCC": 19.0, "CCA": 16.7, "CCG": 6.1,
        "ACT": 13.4, "ACC": 19.5, "ACA": 15.3, "ACG": 5.5,
        "GCT": 19.9, "GCC": 26.8, "GCA": 15.8, "GCG": 6.2,
        "TAT": 12.4, "TAC": 16.1, "TAA": 0.9,  "TAG": 0.7,
        "CAT": 10.5, "CAC": 14.7, "CAA": 11.6, "CAG": 34.4,
        "AAT": 16.5, "AAC": 19.9, "AAA": 22.9, "AAG": 34.1,
        "GAT": 22.6, "GAC": 26.2, "GAA": 27.6, "GAG": 40.3,
        "TGT": 9.8,  "TGC": 11.8, "TGA": 1.3,  "TGG": 12.9,
        "CGT": 5.1,  "CGC": 10.1, "CGA": 6.5,  "CGG": 10.9,
        "AGT": 11.9, "AGC": 18.6, "AGA": 11.2, "AGG": 10.9,
        "GGT": 11.7, "GGC": 22.2, "GGA": 16.4, "GGG": 15.4,
    },
    "S_cerevisiae": {
        "TTT": 26.1, "TTC": 18.4, "TTA": 26.2, "TTG": 27.2,
        "CTT": 12.3, "CTC": 5.4,  "CTA": 13.4, "CTG": 10.5,
        "ATT": 30.1, "ATC": 17.2, "ATA": 17.8, "ATG": 20.9,
        "GTT": 22.1, "GTC": 11.8, "GTA": 11.8, "GTG": 10.8,
        "TCT": 23.5, "TCC": 14.2, "TCA": 18.7, "TCG": 8.6,
        "CCT": 13.5, "CCC": 6.8,  "CCA": 18.3, "CCG": 5.3,
        "ACT": 20.3, "ACC": 12.7, "ACA": 17.8, "ACG": 8.0,
        "GCT": 21.2, "GCC": 12.6, "GCA": 16.2, "GCG": 6.2,
        "TAT": 18.8, "TAC": 14.8, "TAA": 1.1,  "TAG": 0.5,
        "CAT": 13.6, "CAC": 7.8,  "CAA": 27.3, "CAG": 12.1,
        "AAT": 35.7, "AAC": 24.8, "AAA": 41.9, "AAG": 30.8,
        "GAT": 37.6, "GAC": 20.2, "GAA": 45.6, "GAG": 19.2,
        "TGT": 8.1,  "TGC": 4.8,  "TGA": 0.7,  "TGG": 10.4,
        "CGT": 6.4,  "CGC": 2.6,  "CGA": 3.0,  "CGG": 1.7,
        "AGT": 14.2, "AGC": 9.8,  "AGA": 21.3, "AGG": 9.2,
        "GGT": 23.9, "GGC": 9.8,  "GGA": 10.9, "GGG": 6.0,
    },
}

DEFAULT_AVOID = ["BsaI", "BsmBI", "BbsI", "EcoRI", "BamHI", "HindIII", "NdeI", "XhoI"]
HOMOPOLYMERS = ["TTTT", "AAAAA", "GGGGGG", "CCCCCC"]


def _build_weights(host: str) -> tuple[dict[str, float], dict[str, list[tuple[str, float]]]]:
    freqs = HOST_CODON_FREQ[host]
    aa_to_codons: dict[str, list[tuple[str, float]]] = {}
    for cod, aa in CODON_TABLE.items():
        aa_to_codons.setdefault(aa, []).append((cod, freqs[cod]))
    weights: dict[str, float] = {}
    ranked: dict[str, list[tuple[str, float]]] = {}
    for aa, pairs in aa_to_codons.items():
        max_f = max(f for _, f in pairs)
        w_pairs = sorted(((c, f / max_f) for c, f in pairs), key=lambda x: -x[1])
        ranked[aa] = w_pairs
        for c, w in w_pairs:
            weights[c] = w
    return weights, ranked


def _compute_cai(codons: list[str], weights: dict[str, float]) -> float:
    vals = [max(1e-4, weights[c]) for c in codons if CODON_TABLE.get(c, "*") != "*"]
    if not vals:
        return 1.0
    return round(math.exp(sum(math.log(v) for v in vals) / len(vals)), 4)


def _sliding_gc(seq: str, win: int = 30) -> list[float]:
    if len(seq) <= win:
        return [gc_fraction(seq)]
    return [gc_fraction(seq[i : i + win]) for i in range(0, len(seq) - win + 1, 3)]


def _count_violations(
    seq: str, *, enzymes: list[str], host: str, gc_min: float, gc_max: float
) -> tuple[int, int, int, int]:
    enz_hits = len(find_restriction_sites(seq, enzymes))
    homo_hits = sum(seq.count(hp) for hp in HOMOPOLYMERS)
    sd_hits = seq.count("AGGAGG") if host == "E_coli_K12" else 0
    gc_wins = _sliding_gc(seq, 30)
    gc_oob = sum(1 for g in gc_wins if g < gc_min or g > gc_max)
    return enz_hits, homo_hits, sd_hits, gc_oob


def _objective(
    codons: list[str],
    weights: dict[str, float],
    *,
    enzymes: list[str],
    host: str,
    gc_min: float,
    gc_max: float,
) -> float:
    seq = "".join(codons)
    enz_hits, homo_hits, sd_hits, gc_oob = _count_violations(
        seq, enzymes=enzymes, host=host, gc_min=gc_min, gc_max=gc_max
    )
    cai = _compute_cai(codons, weights)
    # Mild 5' GC relaxation penalty across codons 2..8
    five_prime_seq = "".join(codons[1:8])
    gc_5p = gc_fraction(five_prime_seq) if five_prime_seq else 0.45
    five_p_pen = max(0.0, gc_5p - 0.55) * 2.0
    return (
        1000.0 * enz_hits
        + 500.0 * homo_hits
        + 400.0 * sd_hits
        + 50.0 * gc_oob
        + five_p_pen
        - 10.0 * cai
    )


def _optimize_codons(
    aa_seq: str,
    *,
    host: str,
    enzymes: list[str],
    gc_min: float,
    gc_max: float,
    seed: int,
) -> list[str]:
    rng = random.Random(seed)
    weights, ranked = _build_weights(host)

    codons: list[str] = []
    for idx, aa in enumerate(aa_seq):
        opts = ranked[aa]
        if idx in range(1, 8) and len(opts) > 1:
            # Prefer high-adaptiveness codon with lower GC in 5' leader codons 2..8
            viable = [c for c, w in opts if w >= 0.35] or [opts[0][0]]
            viable.sort(key=lambda c: (c.count("G") + c.count("C"), -weights[c]))
            codons.append(viable[0])
        else:
            # Deterministic weighted choice biased strongly to top codons
            if len(opts) == 1 or rng.random() < 0.85:
                codons.append(opts[0][0])
            else:
                top_two = opts[:2]
                codons.append(top_two[1][0] if top_two[1][1] >= 0.45 else top_two[0][0])

    # Iterative local synonymous search to eliminate all enzyme sites, homopolymers, SD motifs & GC outliers
    best_obj = _objective(codons, weights, enzymes=enzymes, host=host, gc_min=gc_min, gc_max=gc_max)
    for _pass in range(6):
        improved = False
        indices = list(range(len(codons)))
        rng.shuffle(indices)
        for pos in indices:
            aa = aa_seq[pos]
            opts = ranked[aa]
            if len(opts) <= 1:
                continue
            orig = codons[pos]
            best_c = orig
            for cand_c, _w in opts:
                if cand_c == orig:
                    continue
                codons[pos] = cand_c
                obj = _objective(
                    codons, weights, enzymes=enzymes, host=host, gc_min=gc_min, gc_max=gc_max
                )
                if obj < best_obj - 1e-9:
                    best_obj = obj
                    best_c = cand_c
                    improved = True
            codons[pos] = best_c
        if not improved and best_obj < 0.5:
            break

    return codons


def _plot_profile(
    out_png: Path,
    seq_before: str | None,
    seq_after: str,
    *,
    host: str,
    cai_before: float | None,
    cai_after: float,
    gc_min: float,
    gc_max: float,
) -> None:
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.3), gridspec_kw={"width_ratios": [2.2, 1.0]})

    win_after = _sliding_gc(seq_after, 30)
    x_after = np.arange(len(win_after)) * 3 + 15
    if seq_before:
        win_before = _sliding_gc(seq_before, 30)
        x_before = np.arange(len(win_before)) * 3 + 15
        ax1.plot(x_before, np.array(win_before) * 100.0, color="#94a3b8", lw=1.3, label="Input CDS GC%")
    ax1.plot(x_after, np.array(win_after) * 100.0, color="#0f766e", lw=2.0, label="Optimized CDS GC%")
    ax1.axhspan(gc_min * 100.0, gc_max * 100.0, color="#ccfbf1", alpha=0.45, label=f"Target [{gc_min*100:.0f}%, {gc_max*100:.0f}%]")
    ax1.set_xlabel("CDS Position (bp)")
    ax1.set_ylabel("30-bp Sliding GC (%)")
    ax1.set_title(f"Sliding-Window GC Profile ({host})")
    ax1.set_ylim(0, 100)
    ax1.legend(frameon=False, fontsize=8.5)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    labels = ["Before", "After"] if cai_before is not None else ["Optimized"]
    vals = [cai_before, cai_after] if cai_before is not None else [cai_after]
    colors = ["#94a3b8", "#2563eb"] if cai_before is not None else ["#2563eb"]
    bars = ax2.bar(labels, vals, color=colors, width=0.5)
    for b in bars:
        h = b.get_height()
        ax2.text(b.get_x() + b.get_width() / 2, h + 0.02, f"{h:.3f}", ha="center", va="bottom", fontsize=9)
    ax2.set_ylim(0, 1.12)
    ax2.set_ylabel("Codon Adaptation Index (CAI)")
    ax2.set_title("CAI Improvement")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Optimize coding sequence for host expression.")
    parser.add_argument("--sequence_input", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--is_dna", action="store_true")
    parser.add_argument(
        "--host",
        default="E_coli_K12",
        choices=["E_coli_K12", "HEK293_Human", "CHO_Cricetulus", "S_cerevisiae"],
    )
    parser.add_argument("--avoid_enzymes", action="append", default=None)
    parser.add_argument("--gc_window_min", type=float, default=0.30)
    parser.add_argument("--gc_window_max", type=float, default=0.72)
    parser.add_argument("--add_kozak", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    enzymes = parse_list_arg(args.avoid_enzymes, DEFAULT_AVOID)

    recs = parse_sequences(args.sequence_input)
    rec = recs[0]
    raw_seq = rec["seq"].replace(" ", "").upper()

    is_strict_acgt = bool(re.fullmatch(r"[ACGT]+", raw_seq))
    if args.is_dna:
        input_type = "DNA"
        if not is_strict_acgt or len(raw_seq) % 3 != 0:
            fail(f"--is_dna set, but sequence length ({len(raw_seq)}) is not a multiple of 3 or contains non-ACGT bases.")
    else:
        if is_strict_acgt and len(raw_seq) % 3 == 0 and (raw_seq.startswith("ATG") or len(raw_seq) >= 30):
            input_type = "DNA"
        else:
            input_type = "Protein"

    weights, _ = _build_weights(args.host)
    seq_before: str | None = None
    cai_before: float | None = None
    gc_before: float | None = None

    if input_type == "DNA":
        seq_before = raw_seq
        aa_seq = translate_dna(seq_before)
        if "*" in aa_seq[:-1]:
            fail("Input DNA CDS contains an internal premature stop codon.")
        orig_codons = [seq_before[i : i + 3] for i in range(0, len(seq_before), 3)]
        cai_before = _compute_cai(orig_codons, weights)
        gc_before = round(gc_fraction(seq_before), 4)
    else:
        aa_seq = raw_seq
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY*")
        bad = set(aa_seq) - valid_aa
        if bad:
            fail(f"Protein sequence contains unsupported amino acid codes: {sorted(bad)}")
        if "*" in aa_seq[:-1]:
            fail("Input protein sequence contains an internal stop '*' before the C-terminus.")

    opt_codons = _optimize_codons(
        aa_seq,
        host=args.host,
        enzymes=enzymes,
        gc_min=args.gc_window_min,
        gc_max=args.gc_window_max,
        seed=args.seed,
    )
    opt_cds = "".join(opt_codons)
    translated_opt = translate_dna(opt_cds)
    if translated_opt != aa_seq:
        fail("Internal error: optimized CDS translation diverged from target protein sequence.")

    full_seq = ("GCCACC" + opt_cds) if args.add_kozak else opt_cds
    cds_start = 7 if args.add_kozak else 1
    cds_end = len(full_seq)

    enz_hits = find_restriction_sites(full_seq, enzymes)
    homo_hits = sum(opt_cds.count(hp) for hp in HOMOPOLYMERS)
    sd_hits = opt_cds.count("AGGAGG") if args.host == "E_coli_K12" else 0
    gc_wins = _sliding_gc(opt_cds, 30)
    gc_oob = sum(1 for g in gc_wins if g < args.gc_window_min or g > args.gc_window_max)

    cai_after = _compute_cai(opt_codons, weights)
    gc_after = round(gc_fraction(opt_cds), 4)

    fasta_path = out_dir / "optimized_cds.fasta"
    write_fasta(fasta_path, [(f"{rec['id']}_{args.host}_opt", full_seq)])

    gb_features = []
    if args.add_kozak:
        gb_features.append({"type": "misc_feature", "location": "1..6", "qualifiers": {"label": "Kozak_consensus"}})
    gb_features.append({
        "type": "CDS",
        "location": f"{cds_start}..{cds_end}",
        "qualifiers": {
            "label": f"{rec['id']}_opt",
            "host": args.host,
            "translation": aa_seq.rstrip("*"),
        },
    })
    gb_path = out_dir / "optimized_construct.gb"
    write_genbank(gb_path, f"{rec['id']}_opt", full_seq, topology="linear", features=gb_features)

    audit_rows = []
    for idx, (aa, cod) in enumerate(zip(aa_seq, opt_codons), start=1):
        audit_rows.append({
            "codon_index": idx,
            "aa": aa,
            "orig_codon": seq_before[(idx - 1) * 3 : idx * 3] if seq_before else "",
            "opt_codon": cod,
            "adaptiveness_w": round(weights[cod], 4),
        })
    audit_df = pd.DataFrame(audit_rows)
    audit_path = out_dir / "codon_audit.tsv"
    audit_df.to_csv(audit_path, sep="\t", index=False)

    png_path = out_dir / "codon_optimization_profile.png"
    _plot_profile(
        png_path,
        seq_before,
        opt_cds,
        host=args.host,
        cai_before=cai_before,
        cai_after=cai_after,
        gc_min=args.gc_window_min,
        gc_max=args.gc_window_max,
    )

    warnings: list[str] = []
    if enz_hits:
        warnings.append(f"{len(enz_hits)} forbidden restriction site(s) could not be removed cleanly.")
    if gc_oob > 0:
        warnings.append(f"{gc_oob} 30-bp GC windows remain outside [{args.gc_window_min}, {args.gc_window_max}].")

    tool_name = "mc_codon_optimize"
    params = {"host": args.host, "avoid_enzymes": enzymes, "add_kozak": bool(args.add_kozak), "seed": args.seed}
    artifacts = [
        write_artifact(fasta_path, tool_name, "Codon-optimized CDS FASTA", inputs=[args.sequence_input], parameters=params),
        write_artifact(gb_path, tool_name, "Codon-optimized construct GenBank", inputs=[args.sequence_input], parameters=params),
        write_artifact(audit_path, tool_name, "Per-codon optimization audit", inputs=[args.sequence_input], parameters=params, shape=audit_df.shape),
        write_artifact(png_path, tool_name, "Codon optimization GC & CAI profile", inputs=[args.sequence_input], parameters=params),
    ]

    emit({
        "status": "ok",
        "tool": tool_name,
        "host": args.host,
        "input_type": input_type,
        "sequence_id": rec["id"],
        "aa_length": len(aa_seq.rstrip("*")),
        "cds_length_bp": len(opt_cds),
        "cai_before": cai_before,
        "cai_after": cai_after,
        "gc_before": gc_before,
        "gc_after": gc_after,
        "translation_identical": True,
        "remaining_restriction_sites": len(enz_hits),
        "remaining_homopolymers": homo_hits,
        "remaining_shine_dalgarno": sd_hits,
        "out_of_range_gc_windows": gc_oob,
        "kozak_prepended": bool(args.add_kozak),
        "outputs": {
            "optimized_cds_fasta": str(fasta_path),
            "optimized_construct_gb": str(gb_path),
            "codon_audit": str(audit_path),
            "profile_plot": str(png_path),
        },
        "artifacts": artifacts,
        "warnings": warnings,
    })


if __name__ == "__main__":
    main()
