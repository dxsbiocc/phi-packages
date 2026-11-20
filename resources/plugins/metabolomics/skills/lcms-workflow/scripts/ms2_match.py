#!/usr/bin/env python3
"""mb_ms2_match: cosine matching of query MGF spectra against MSP/MGF libraries."""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from _common import emit, fail, write_artifact

SKILL = "lcms-workflow"


def parse_msp(path):
    entries = []
    cur = {}
    peaks = []
    in_peaks = False
    for line in Path(path).read_text(errors="replace").splitlines():
        line = line.strip()
        if not line:
            if cur and peaks:
                entries.append((cur, peaks))
            cur, peaks, in_peaks = {}, [], False
            continue
        if in_peaks:
            parts = re.split(r"[\s:;]+", line)
            if len(parts) >= 2:
                try:
                    peaks.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
            continue
        m = re.match(r"^([^:=]+)[:=]\s*(.*)$", line)
        if not m:
            continue
        key, val = m.group(1).strip().lower(), m.group(2).strip()
        if key in {"name", "compound", "title"}:
            cur["name"] = val
        elif key in {"precursormz", "precursor_mz", "precursor mz", "parent", "parentmass"}:
            try:
                cur["precursor"] = float(val.split()[0])
            except (ValueError, IndexError):
                pass
        elif key in {"num peaks", "num_peaks", "num peaks:"}:
            in_peaks = True
        elif key == "comments" and "precursor" not in cur:
            pm = re.search(r"(?:parent|precursor)[=:\s]+([0-9.]+)", val, re.I)
            if pm:
                cur["precursor"] = float(pm.group(1))
    if cur and peaks:
        entries.append((cur, peaks))
    return entries


def parse_mgf(path):
    entries = []
    cur = None
    peaks = []
    for line in Path(path).read_text(errors="replace").splitlines():
        line = line.strip()
        if line.upper() == "BEGIN IONS":
            cur, peaks = {}, []
        elif line.upper() == "END IONS":
            if cur is not None and peaks:
                entries.append((cur, peaks))
            cur, peaks = None, []
        elif cur is not None:
            if "=" in line:
                key, val = line.split("=", 1)
                key = key.strip().upper()
                if key == "TITLE":
                    cur["name"] = val.strip()
                elif key == "PEPMASS":
                    try:
                        cur["precursor"] = float(val.split()[0])
                    except (ValueError, IndexError):
                        pass
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        peaks.append((float(parts[0]), float(parts[1])))
                    except ValueError:
                        pass
    return entries


def load_library(path):
    text = Path(path).read_text(errors="replace")[:100000].upper()
    if "BEGIN IONS" in text:
        return parse_mgf(path)
    return parse_msp(path)


def cosine(qm, qi, lm, li, frag_tol):
    qs, ls = np.sqrt(qi), np.sqrt(li)
    nq, nl = np.linalg.norm(qs), np.linalg.norm(ls)
    if nq == 0 or nl == 0:
        return 0.0, 0
    pairs = []
    for a in range(len(qm)):
        for b in range(len(lm)):
            d = abs(qm[a] - lm[b])
            if d <= frag_tol:
                pairs.append((d, a, b))
    pairs.sort()
    used_q, used_l = set(), set()
    dot = 0.0
    for d, a, b in pairs:
        if a in used_q or b in used_l:
            continue
        used_q.add(a)
        used_l.add(b)
        dot += qs[a] * ls[b]
    return float(dot / (nq * nl)), len(used_q)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True)
    ap.add_argument("--library", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--fragment_tol", type=float, default=0.02)
    ap.add_argument("--min_score", type=float, default=0.7)
    ap.add_argument("--top", type=int, default=5)
    args = ap.parse_args()

    queries = parse_mgf(args.query)
    if not queries:
        fail("no spectra parsed from query MGF")
    library = load_library(args.library)
    if not library:
        fail("no spectra parsed from library")
    lib = [(meta.get("name", f"lib_{i}"), meta.get("precursor"),
            np.array([p[0] for p in peaks]), np.array([p[1] for p in peaks]))
           for i, (meta, peaks) in enumerate(library)]
    n_no_precursor = sum(1 for _, prec, _, _ in lib if prec is None)

    rows = []
    n_skipped = 0
    for qi, (qmeta, qpeaks) in enumerate(queries):
        qname = qmeta.get("name", f"query_{qi}")
        qprec = qmeta.get("precursor")
        if qprec is None:
            n_skipped += 1
            continue
        qm = np.array([p[0] for p in qpeaks])
        qv = np.array([p[1] for p in qpeaks])
        tol = args.ppm * qprec / 1e6
        hits = []
        for lname, lprec, lm, lv in lib:
            if lprec is None or abs(lprec - qprec) > tol:
                continue
            score, nmatch = cosine(qm, qv, lm, lv, args.fragment_tol)
            if score >= args.min_score:
                hits.append((score, nmatch, lname, lprec))
        hits.sort(key=lambda h: -h[0])
        for rank, (score, nmatch, lname, lprec) in enumerate(hits[:args.top], 1):
            rows.append({"query": qname, "query_precursor": qprec,
                         "library_name": lname, "library_precursor": lprec,
                         "cosine": score, "matched_fragments": nmatch,
                         "query_fragments": len(qm), "rank": rank})

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(rows, columns=["query", "query_precursor", "library_name",
                                        "library_precursor", "cosine",
                                        "matched_fragments", "query_fragments", "rank"])
    table_path = out / "ms2-matches.tsv"
    table.to_csv(table_path, sep="\t", index=False)

    params = {"ppm": args.ppm, "fragment_tol": args.fragment_tol,
              "min_score": args.min_score, "top": args.top}
    artifacts = [
        write_artifact(table_path, "mb_ms2_match", "MS2 library match table",
                       skill=SKILL, inputs=[args.query, args.library],
                       parameters=params, shape=table.shape),
    ]
    emit({
        "query": args.query,
        "library": args.library,
        "n_query_spectra": len(queries),
        "n_library_spectra": len(lib),
        "n_library_without_precursor": n_no_precursor,
        "n_query_skipped_no_precursor": n_skipped,
        "n_queries_with_hits": int(table["query"].nunique()) if len(table) else 0,
        "n_hits": int(len(table)),
        "ppm": args.ppm,
        "fragment_tol": args.fragment_tol,
        "min_score": args.min_score,
        "caveat": "library spectral match = MSI level 2 (putative); "
                  "level 1 requires an authentic standard (RT + MS2) on the same method",
        "output_dir": str(out),
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()