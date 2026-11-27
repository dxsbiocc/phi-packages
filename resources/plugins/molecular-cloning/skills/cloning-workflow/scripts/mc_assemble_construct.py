#!/usr/bin/env python3
"""Virtual DNA construct assembly via GoldenGate, Gibson, or RestrictionLigation."""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    calc_tm_santalucia_owczarzy,
    emit,
    fail,
    find_restriction_sites,
    parse_list_arg,
    parse_sequences,
    resolve_enzyme,
    revcomp,
    setup_matplotlib,
    write_artifact,
    write_fasta,
    write_genbank,
)


def _digest_type_iis(part: dict, enzyme: str) -> list[dict]:
    """Digest a linear or circular part with a Type IIS enzyme into sticky-ended sub-fragments."""
    enz, spec = resolve_enzyme(enzyme)
    site = spec["site"]
    rc_site = revcomp(site)
    top_cut = spec["top_cut"]
    bot_cut = spec["bot_cut"]
    oh_len = bot_cut - top_cut
    site_len = len(site)
    seq = part["seq"]
    L = len(seq)
    is_circ = part.get("topology", "linear") == "circular"

    # Represent cuts as (top_cut_0, bot_cut_0) in [0, L)
    # On '+' strand site at pos: cuts at pos + top_cut, pos + bot_cut
    # On '-' strand rc_site at pos: cuts at pos + site_len - bot_cut, pos + site_len - top_cut
    cuts: list[tuple[int, int, str]] = []
    scan_seq = (seq + seq[: site_len + bot_cut]) if is_circ else seq
    for i in range(L if is_circ else max(0, L - site_len + 1)):
        sub = scan_seq[i : i + site_len]
        if sub == site:
            t_c = (i + top_cut) % L if is_circ else (i + top_cut)
            b_c = (i + bot_cut) % L if is_circ else (i + bot_cut)
            if is_circ or (0 <= t_c and b_c <= L):
                cuts.append((t_c, b_c, "+"))
        elif sub == rc_site:
            t_c = (i + site_len - bot_cut) % L if is_circ else (i + site_len - bot_cut)
            b_c = (i + site_len - top_cut) % L if is_circ else (i + site_len - top_cut)
            if is_circ or (0 <= t_c and b_c <= L):
                cuts.append((t_c, b_c, "-"))

    if not cuts:
        return []
    cuts.sort(key=lambda x: x[0])

    frags: list[dict] = []
    if is_circ:
        doubled = seq + seq
        for idx in range(len(cuts)):
            t1, b1, _ = cuts[idx]
            t2, b2, _ = cuts[(idx + 1) % len(cuts)]
            if t2 <= t1:
                t2 += L
                b2 += L
            if b1 < t1:
                b1 += L
            if b2 < t2:
                b2 += L
            left_oh = doubled[t1:b1]
            right_oh = doubled[t2:b2]
            body = doubled[b1:t2]
            if len(left_oh) == oh_len and len(right_oh) == oh_len:
                frags.append({
                    "name": f"{part['id']}_frag{idx + 1}",
                    "left_oh": left_oh,
                    "right_oh": right_oh,
                    "body": body,
                    "retained_sites": len(find_restriction_sites(left_oh + body + right_oh, [enz])),
                })
    else:
        for idx in range(len(cuts) - 1):
            t1, b1, s1 = cuts[idx]
            t2, b2, s2 = cuts[idx + 1]
            if b1 <= t2:
                left_oh = seq[t1:b1]
                right_oh = seq[t2:b2]
                body = seq[b1:t2]
                # Keep interior fragments free of outward recognition sites (or check retained)
                frags.append({
                    "name": f"{part['id']}_frag{idx + 1}",
                    "left_oh": left_oh,
                    "right_oh": right_oh,
                    "body": body,
                    "retained_sites": len(find_restriction_sites(left_oh + body + right_oh, [enz])),
                })
    return frags


def _flip_frag(frag: dict) -> dict:
    return {
        "name": frag["name"] + "_rc",
        "left_oh": revcomp(frag["right_oh"]),
        "right_oh": revcomp(frag["left_oh"]),
        "body": revcomp(frag["body"]),
        "retained_sites": frag["retained_sites"],
    }


def _assemble_goldengate(parts: list[dict], enzyme: str, linear_product: bool) -> tuple[str, list[dict], list[str], list[dict]]:
    enz, _ = resolve_enzyme(enzyme)
    # Extract the productive cloning fragment from each input part (lowest retained sites, longest body)
    part_frags: list[dict] = []
    for p in parts:
        cands = _digest_type_iis(p, enz)
        if not cands:
            fail(f"Part '{p['id']}' contains no valid {enz} cut pairs for GoldenGate assembly.")
        # Prefer fragments with 0 retained sites and largest body
        cands.sort(key=lambda f: (f["retained_sites"], -len(f["body"])))
        part_frags.append(cands[0])

    # Build directed overhang graph chain starting from part 0 in forward orientation
    ordered: list[dict] = [part_frags[0]]
    remaining = list(range(1, len(part_frags)))
    while remaining:
        cur_right = ordered[-1]["right_oh"]
        matched_idx = None
        matched_frag = None
        for idx in remaining:
            f = part_frags[idx]
            if f["left_oh"] == cur_right:
                matched_idx = idx
                matched_frag = f
                break
            f_rc = _flip_frag(f)
            if f_rc["left_oh"] == cur_right:
                matched_idx = idx
                matched_frag = f_rc
                break
        if matched_idx is None or matched_frag is None:
            avail = [(part_frags[i]["name"], part_frags[i]["left_oh"], part_frags[i]["right_oh"]) for i in remaining]
            fail(f"GoldenGate overhang graph broken after '{ordered[-1]['name']}' (3' overhang {cur_right}); available parts: {avail}")
        ordered.append(matched_frag)
        remaining.remove(matched_idx)

    if not linear_product and ordered[-1]["right_oh"] != ordered[0]["left_oh"]:
        fail(
            f"Circular GoldenGate assembly failed to close: final overhang '{ordered[-1]['right_oh']}' != initial '{ordered[0]['left_oh']}'."
        )

    # Stitch sequence: for circular, each junction contributes one copy of oh + body; for linear, include flanking overhangs
    junctions: list[dict] = []
    features: list[dict] = []
    if linear_product:
        seq_pieces = [ordered[0]["left_oh"]]
        cur_pos = len(ordered[0]["left_oh"])
    else:
        seq_pieces = []
        cur_pos = 0

    for idx, f in enumerate(ordered):
        if not linear_product or idx > 0:
            j_oh = f["left_oh"]
            j_start = cur_pos + 1
            seq_pieces.append(j_oh)
            cur_pos += len(j_oh)
            junctions.append({
                "junction_index": idx + 1,
                "upstream_part": ordered[(idx - 1) % len(ordered)]["name"],
                "downstream_part": f["name"],
                "overhang_or_overlap": j_oh,
                "length_bp": len(j_oh),
                "tm_c": None,
                "position_1": j_start,
                "status": "ok",
            })
        body_start = cur_pos + 1
        seq_pieces.append(f["body"])
        cur_pos += len(f["body"])
        features.append({
            "type": "misc_feature",
            "location": f"{body_start}..{cur_pos}",
            "qualifiers": {"label": f["name"]},
        })

    if linear_product:
        seq_pieces.append(ordered[-1]["right_oh"])
    assembled_seq = "".join(seq_pieces)

    # Overhang fidelity auditing
    fidelity_flags: list[str] = []
    oh_list = [j["overhang_or_overlap"] for j in junctions]
    for oh in oh_list:
        if oh == revcomp(oh):
            fidelity_flags.append(f"PALINDROMIC_OVERHANG:{oh}")
    if len(set(oh_list)) < len(oh_list):
        fidelity_flags.append("DUPLICATE_OVERHANG")
    for i in range(len(oh_list)):
        for j in range(i + 1, len(oh_list)):
            a, b = oh_list[i], oh_list[j]
            if a == revcomp(b):
                fidelity_flags.append(f"REVCOMP_COLLISION:{a}/{b}")
            else:
                mm_direct = sum(x != y for x, y in zip(a, b))
                mm_rc = sum(x != y for x, y in zip(a, revcomp(b)))
                if min(mm_direct, mm_rc) <= 1:
                    fidelity_flags.append(f"1BP_CROSSTALK:{a}/{b}")

    return assembled_seq, junctions, fidelity_flags, features


def _find_overlap(a: str, b: str, min_ov: int = 10, max_ov: int = 80) -> int:
    limit = min(max_ov, len(a), len(b))
    for k in range(limit, min_ov - 1, -1):
        if a[-k:] == b[:k]:
            return k
    return 0


def _assemble_gibson(parts: list[dict], linear_product: bool) -> tuple[str, list[dict], list[str], list[dict]]:
    if len(parts) < 2:
        fail("Gibson assembly requires at least 2 DNA fragments.")
    ordered = [dict(parts[0])]
    remaining = list(parts[1:])
    while remaining:
        cur_seq = ordered[-1]["seq"]
        best_idx = None
        best_ov = 0
        best_part = None
        for idx, p in enumerate(remaining):
            ov_f = _find_overlap(cur_seq, p["seq"])
            if ov_f > best_ov:
                best_ov, best_idx, best_part = ov_f, idx, p
            rc_seq = revcomp(p["seq"])
            ov_r = _find_overlap(cur_seq, rc_seq)
            if ov_r > best_ov:
                best_ov, best_idx = ov_r, idx
                best_part = {"id": p["id"] + "_rc", "seq": rc_seq, "topology": p.get("topology", "linear")}
        if best_idx is None or best_part is None or best_ov < 10:
            fail(f"No >=10 bp terminal homology overlap found after '{ordered[-1]['id']}' in Gibson assembly.")
        ordered.append(best_part)
        remaining.remove(remaining[best_idx])

    junctions: list[dict] = []
    fidelity_flags: list[str] = []
    features: list[dict] = []
    pieces: list[str] = []
    cur_pos = 0

    for idx in range(len(ordered)):
        cur = ordered[idx]
        nxt = ordered[(idx + 1) % len(ordered)]
        if linear_product and idx == len(ordered) - 1:
            start_1 = cur_pos + 1
            pieces.append(cur["seq"])
            cur_pos += len(cur["seq"])
            features.append({"type": "misc_feature", "location": f"{start_1}..{cur_pos}", "qualifiers": {"label": cur["id"]}})
            break
        ov = _find_overlap(cur["seq"], nxt["seq"])
        if ov < 10:
            fail(f"Circular Gibson closure failed between '{cur['id']}' and '{nxt['id']}' (overlap < 10 bp).")
        ov_seq = cur["seq"][-ov:]
        tm = calc_tm_santalucia_owczarzy(ov_seq)["tm_c"]
        status = "ok"
        if ov < 15:
            status = "short_overlap"
            fidelity_flags.append(f"SHORT_OVERLAP:{cur['id']}->{nxt['id']}({ov}bp)")
        if tm < 48.0:
            status = "low_tm"
            fidelity_flags.append(f"LOW_OVERLAP_TM:{cur['id']}->{nxt['id']}({tm}C)")

        frag_core = cur["seq"][:-ov]
        start_1 = cur_pos + 1
        pieces.append(frag_core)
        cur_pos += len(frag_core)
        features.append({"type": "misc_feature", "location": f"{start_1}..{cur_pos}", "qualifiers": {"label": cur["id"]}})
        junctions.append({
            "junction_index": idx + 1,
            "upstream_part": cur["id"],
            "downstream_part": nxt["id"],
            "overhang_or_overlap": ov_seq,
            "length_bp": ov,
            "tm_c": tm,
            "position_1": cur_pos + 1,
            "status": status,
        })

    return "".join(pieces), junctions, fidelity_flags, features


def _digest_type_ii(part: dict, enzymes: list[str]) -> list[dict]:
    seq = part["seq"]
    L = len(seq)
    is_circ = part.get("topology", "linear") == "circular"
    cuts: list[dict] = []
    for raw_e in enzymes:
        enz, spec = resolve_enzyme(raw_e)
        site = spec["site"]
        slen = len(site)
        scan = (seq + seq[: slen - 1]) if is_circ else seq
        for i in range(L if is_circ else max(0, L - slen + 1)):
            if scan[i : i + slen] == site:
                tc = (i + spec["top_cut"]) % L if is_circ else (i + spec["top_cut"])
                bc = (i + spec["bot_cut"]) % L if is_circ else (i + spec["bot_cut"])
                oh_type = spec["overhang_type"]
                if oh_type == "5prime":
                    oh = (seq + seq)[tc:bc] if is_circ and bc < tc else seq[tc:bc]
                    c_start, c_end = tc, bc
                elif oh_type == "3prime":
                    oh = (seq + seq)[bc:tc] if is_circ and tc < bc else seq[bc:tc]
                    c_start, c_end = bc, tc
                else:
                    oh = ""
                    c_start, c_end = tc, tc
                cuts.append({"enzyme": enz, "c_start": c_start, "c_end": c_end, "oh": oh, "oh_type": oh_type})
    if not cuts:
        return []
    cuts.sort(key=lambda c: c["c_start"])
    frags: list[dict] = []
    doubled = seq + seq
    if is_circ:
        for i in range(len(cuts)):
            c1 = cuts[i]
            c2 = cuts[(i + 1) % len(cuts)]
            s1, e1 = c1["c_start"], c1["c_end"]
            s2, e2 = c2["c_start"], c2["c_end"]
            if s2 <= s1:
                s2 += L
                e2 += L
            if e1 < s1:
                e1 += L
            if e2 < s2:
                e2 += L
            body = doubled[e1:s2]
            frags.append({
                "name": f"{part['id']}_{c1['enzyme']}_{c2['enzyme']}",
                "left_oh": c1["oh"],
                "left_type": c1["oh_type"],
                "left_enz": c1["enzyme"],
                "right_oh": c2["oh"],
                "right_type": c2["oh_type"],
                "right_enz": c2["enzyme"],
                "body": body,
            })
    else:
        for i in range(len(cuts) - 1):
            c1 = cuts[i]
            c2 = cuts[i + 1]
            if c1["c_end"] <= c2["c_start"]:
                frags.append({
                    "name": f"{part['id']}_{c1['enzyme']}_{c2['enzyme']}",
                    "left_oh": c1["oh"],
                    "left_type": c1["oh_type"],
                    "left_enz": c1["enzyme"],
                    "right_oh": c2["oh"],
                    "right_type": c2["oh_type"],
                    "right_enz": c2["enzyme"],
                    "body": seq[c1["c_end"] : c2["c_start"]],
                })
    return frags


def _ends_compatible(f1: dict, f2: dict) -> bool:
    return f1["right_type"] == f2["left_type"] and f1["right_oh"] == f2["left_oh"]


def _flip_type_ii(f: dict) -> dict:
    return {
        "name": f["name"] + "_rc",
        "left_oh": revcomp(f["right_oh"]),
        "left_type": f["right_type"],
        "left_enz": f["right_enz"],
        "right_oh": revcomp(f["left_oh"]),
        "right_type": f["left_type"],
        "right_enz": f["left_enz"],
        "body": revcomp(f["body"]),
    }


def _assemble_restriction(parts: list[dict], enzyme_str: str, linear_product: bool) -> tuple[str, list[dict], list[str], list[dict]]:
    if not enzyme_str:
        fail("RestrictionLigation requires --enzyme (e.g. 'EcoRI,BamHI').")
    enzymes = parse_list_arg([enzyme_str], [])
    part_frags: list[list[dict]] = []
    for p in parts:
        cands = _digest_type_ii(p, enzymes)
        if not cands:
            fail(f"Part '{p['id']}' has no cut sites for enzymes {enzymes}.")
        # Keep longest digested fragment per part (vector backbone or released insert)
        cands.sort(key=lambda f: -len(f["body"]))
        part_frags.append(cands)

    # Pick longest fragment for each part and orient into compatible chain
    chosen = [c[0] for c in part_frags]
    ordered = [chosen[0]]
    remaining = chosen[1:]
    while remaining:
        matched = None
        for idx, cand in enumerate(remaining):
            if _ends_compatible(ordered[-1], cand):
                matched = (idx, cand)
                break
            rc = _flip_type_ii(cand)
            if _ends_compatible(ordered[-1], rc):
                matched = (idx, rc)
                break
        if matched is None:
            fail(
                f"Incompatible restriction overhangs after '{ordered[-1]['name']}' ({ordered[-1]['right_type']}:{ordered[-1]['right_oh']})."
            )
        ordered.append(matched[1])
        remaining.pop(matched[0])

    if not linear_product and not _ends_compatible(ordered[-1], ordered[0]):
        fail(
            f"Circular RestrictionLigation failed to close: '{ordered[-1]['right_oh']}' != '{ordered[0]['left_oh']}'."
        )

    junctions: list[dict] = []
    features: list[dict] = []
    pieces: list[str] = []
    cur_pos = 0
    for idx, f in enumerate(ordered):
        if not linear_product or idx > 0:
            oh = f["left_oh"]
            j_pos = cur_pos + 1
            pieces.append(oh)
            cur_pos += len(oh)
            up = ordered[(idx - 1) % len(ordered)]
            junctions.append({
                "junction_index": idx + 1,
                "upstream_part": up["name"],
                "downstream_part": f["name"],
                "overhang_or_overlap": oh or "BLUNT",
                "length_bp": len(oh),
                "tm_c": None,
                "position_1": j_pos,
                "status": "compatible_hybrid" if up["right_enz"] != f["left_enz"] else "ok",
            })
        b_start = cur_pos + 1
        pieces.append(f["body"])
        cur_pos += len(f["body"])
        features.append({"type": "misc_feature", "location": f"{b_start}..{cur_pos}", "qualifiers": {"label": f["name"]}})

    return "".join(pieces), junctions, [], features


def _simulate_digest_bands(seq: str, topology: str, enzyme_names: list[str]) -> list[int]:
    hits = sorted({h["pos_0"] for h in find_restriction_sites(seq, enzyme_names)})
    L = len(seq)
    if not hits:
        return [L]
    if topology == "circular":
        if len(hits) == 1:
            return [L]
        bands = [hits[i + 1] - hits[i] for i in range(len(hits) - 1)]
        bands.append((L - hits[-1]) + hits[0])
        return sorted(bands, reverse=True)
    pts = [0] + hits + [L]
    return sorted([pts[i + 1] - pts[i] for i in range(len(pts) - 1) if pts[i + 1] > pts[i]], reverse=True)


def _plot_map_and_gel(
    out_png: Path,
    seq: str,
    topology: str,
    construct_name: str,
    features: list[dict],
    junctions: list[dict],
    diag_enzymes: list[str],
) -> None:
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 5.0), gridspec_kw={"width_ratios": [1.35, 1.0]})

    L = max(1, len(seq))
    palette = ["#0f766e", "#2563eb", "#d97706", "#7c3aed", "#dc2626", "#059669"]
    if topology == "circular":
        theta_all = np.linspace(0, 2 * np.pi, 360)
        ax1.plot(np.cos(theta_all), np.sin(theta_all), color="#cbd5e1", lw=4.0)
        for idx, feat in enumerate(features):
            s_s, e_s = feat["location"].split("..")
            s, e = int(s_s) - 1, int(e_s)
            t0 = 0.5 * np.pi - 2 * np.pi * (s / L)
            t1 = 0.5 * np.pi - 2 * np.pi * (e / L)
            th = np.linspace(t0, t1, 80)
            col = palette[idx % len(palette)]
            ax1.plot(np.cos(th), np.sin(th), color=col, lw=7.5, solid_capstyle="butt")
            mid_t = 0.5 * (t0 + t1)
            ax1.text(
                1.22 * math.cos(mid_t),
                1.22 * math.sin(mid_t),
                feat["qualifiers"].get("label", f"part_{idx + 1}")[:18],
                ha="center",
                va="center",
                fontsize=8,
                color=col,
                fontweight="bold",
            )
        for j in junctions:
            pos = int(j["position_1"]) - 1
            tj = 0.5 * np.pi - 2 * np.pi * (pos / L)
            ax1.plot([0.88 * math.cos(tj), 1.08 * math.cos(tj)], [0.88 * math.sin(tj), 1.08 * math.sin(tj)], color="#0f172a", lw=1.4)
        ax1.text(0, 0.05, construct_name, ha="center", va="center", fontsize=11, fontweight="bold", color="#0f172a")
        ax1.text(0, -0.13, f"{len(seq)} bp ({topology})", ha="center", va="center", fontsize=9, color="#475569")
        ax1.set_xlim(-1.55, 1.55)
        ax1.set_ylim(-1.45, 1.45)
        ax1.set_aspect("equal")
        ax1.axis("off")
    else:
        ax1.axhline(0, color="#334155", lw=3.0)
        for idx, feat in enumerate(features):
            s_s, e_s = feat["location"].split("..")
            s, e = int(s_s) - 1, int(e_s)
            col = palette[idx % len(palette)]
            ax1.plot([s, e], [0.25, 0.25], color=col, lw=8.0, solid_capstyle="butt")
            ax1.text(0.5 * (s + e), 0.48, feat["qualifiers"].get("label", f"part_{idx + 1}")[:18], ha="center", fontsize=8, color=col)
        ax1.set_xlim(0, L)
        ax1.set_ylim(-0.6, 1.0)
        ax1.set_yticks([])
        ax1.set_xlabel("Position (bp)")
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.spines["left"].set_visible(False)
    ax1.set_title("Assembled Construct Map")

    # Virtual 1% agarose gel electrophoresis (log10 bp migration)
    ladder = [10000, 6000, 3000, 2000, 1500, 1000, 750, 500, 250]
    lanes = [("1kb Ladder", ladder), ("Uncut", [len(seq)])]
    for enz in diag_enzymes[:2]:
        lanes.append((f"{enz} digest", _simulate_digest_bands(seq, topology, [enz])))

    ax2.set_facecolor("#0f172a")
    for lane_idx, (lname, bands) in enumerate(lanes):
        for b in bands:
            y = math.log10(max(100, b))
            ax2.plot([lane_idx - 0.28, lane_idx + 0.28], [y, y], color="#38bdf8" if lane_idx > 0 else "#e2e8f0", lw=3.0)
            if lane_idx == 0 and b in (10000, 3000, 1000, 500, 250):
                ax2.text(lane_idx - 0.34, y, f"{b}", ha="right", va="center", fontsize=7.5, color="#cbd5e1")
    ax2.set_xticks(range(len(lanes)))
    ax2.set_xticklabels([ln[0] for ln in lanes], fontsize=8.5, rotation=15)
    ax2.set_ylim(2.1, 4.25)
    ax2.set_ylabel("log10(Fragment Size bp)")
    ax2.set_title("Virtual 1% Agarose Gel")

    plt.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Virtually assemble DNA constructs.")
    parser.add_argument("--fragments", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--method", required=True, choices=["GoldenGate", "Gibson", "RestrictionLigation"])
    parser.add_argument("--vector", default=None)
    parser.add_argument("--enzyme", default=None)
    parser.add_argument("--linear_product", action="store_true")
    parser.add_argument("--construct_name", default="pAssembled")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    parts: list[dict] = []
    if args.vector:
        parts.extend(parse_sequences(args.vector))
    parts.extend(parse_sequences(args.fragments))
    if len(parts) < 2:
        fail("Construct assembly requires at least 2 parts (fragments + optional vector).")

    topology = "linear" if args.linear_product else "circular"
    if args.method == "GoldenGate":
        enz = args.enzyme or "BsaI"
        assembled_seq, junctions, fidelity_flags, features = _assemble_goldengate(parts, enz, args.linear_product)
        retained = len(find_restriction_sites(assembled_seq, [enz]))
        diag_enzymes = [enz, "EcoRI"]
    elif args.method == "Gibson":
        assembled_seq, junctions, fidelity_flags, features = _assemble_gibson(parts, args.linear_product)
        retained = len(find_restriction_sites(assembled_seq, ["BsaI"]))
        diag_enzymes = ["EcoRI", "BamHI"]
    else:
        assembled_seq, junctions, fidelity_flags, features = _assemble_restriction(
            parts, args.enzyme or "", args.linear_product
        )
        retained = 0
        diag_enzymes = parse_list_arg([args.enzyme] if args.enzyme else None, ["EcoRI"])

    warnings = list(fidelity_flags)
    if args.method == "GoldenGate" and retained > 0:
        warnings.append(f"Assembled construct retains {retained} internal {args.enzyme or 'BsaI'} recognition site(s).")

    fasta_path = out_dir / "assembled_construct.fasta"
    write_fasta(fasta_path, [(args.construct_name, assembled_seq)])

    gb_path = out_dir / "assembled_construct.gb"
    write_genbank(gb_path, args.construct_name, assembled_seq, topology=topology, features=features)

    j_df = pd.DataFrame(junctions)
    j_path = out_dir / "assembly_junctions.tsv"
    j_df.to_csv(j_path, sep="\t", index=False)

    png_path = out_dir / "construct_map_and_gel.png"
    _plot_map_and_gel(png_path, assembled_seq, topology, args.construct_name, features, junctions, diag_enzymes)

    tool_name = "mc_assemble_construct"
    params = {"method": args.method, "enzyme": args.enzyme, "linear_product": bool(args.linear_product)}
    inputs = [args.fragments, args.vector]
    artifacts = [
        write_artifact(fasta_path, tool_name, "Assembled construct FASTA", inputs=inputs, parameters=params),
        write_artifact(gb_path, tool_name, "Assembled construct GenBank", inputs=inputs, parameters=params),
        write_artifact(j_path, tool_name, "Assembly junction audit", inputs=inputs, parameters=params, shape=j_df.shape),
        write_artifact(png_path, tool_name, "Construct map and virtual gel", inputs=inputs, parameters=params),
    ]

    emit({
        "status": "ok",
        "tool": tool_name,
        "method": args.method,
        "construct_name": args.construct_name,
        "topology": topology,
        "assembled_length_bp": len(assembled_seq),
        "n_fragments_assembled": len(parts),
        "n_junctions": len(junctions),
        "fidelity_flags": fidelity_flags,
        "retained_type_iis_sites": retained,
        "junctions": junctions,
        "outputs": {
            "assembled_fasta": str(fasta_path),
            "assembled_gb": str(gb_path),
            "junctions_table": str(j_path),
            "map_and_gel": str(png_path),
        },
        "artifacts": artifacts,
        "warnings": warnings,
    })


if __name__ == "__main__":
    main()
