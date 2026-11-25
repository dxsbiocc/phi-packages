#!/usr/bin/env python3
"""mg_merge: merge per-sample taxonomic profiles into one taxa x sample matrix."""
import argparse
from pathlib import Path

import pandas as pd

from _common import emit, fail, read_table, sample_stem, write_artifact

RANK_MARKER = {
    "D": "d__", "K": "k__", "P": "p__", "C": "c__", "O": "o__",
    "F": "f__", "G": "g__", "S": "s__", "T": "t__",
}


def parse_kraken2(path, rank):
    rows = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != 6 or f[3].strip() != rank:
                continue
            taxid, name = f[4].strip(), f[5].strip()
            rows[f"{name} [{taxid}]"] = float(f[2])
    return rows, "counts", {}


def parse_bracken(path, rank):
    rows = {}
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "taxonomy_lvl" in header:
            idx = {c: i for i, c in enumerate(header)}
            lvl_i, name_i, reads_i = idx["taxonomy_lvl"], idx["name"], idx["new_est_reads"]
            taxid_i = idx.get("taxonomy_id")
        else:
            lvl_i, name_i, reads_i, taxid_i = 2, 0, 3, 1
            fh.seek(0)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(lvl_i, name_i, reads_i) or f[lvl_i].strip() != rank:
                continue
            name = f[name_i].strip()
            taxid = f[taxid_i].strip() if taxid_i is not None else ""
            key = f"{name} [{taxid}]" if taxid else name
            rows[key] = float(f[reads_i])
    return rows, "counts", {}


def parse_metaphlan(path, rank):
    marker = RANK_MARKER.get(rank.upper())
    if marker is None:
        fail(f"metaphlan profiles support ranks D/K/P/C/O/F/G/S/T, got {rank}")
    deeper = [m for m in RANK_MARKER.values() if m != marker and
              list(RANK_MARKER.values()).index(m) > list(RANK_MARKER.values()).index(marker)]
    rows = {}
    sgb_species = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            clade = f[0]
            if marker not in clade or any(m in clade for m in deeper):
                continue
            name = clade.split(marker)[-1]
            rows[name] = float(f[2])
            if marker == "t__":
                for part in clade.split("|t__")[0].split("|"):
                    if part.startswith("s__") and part[3:]:
                        sgb_species[name] = part[3:]
    return rows, "relative", {"sgb_species": sgb_species}


def parse_mpa(path, rank):
    marker = RANK_MARKER.get(rank.upper(), rank)
    rows = {}
    sgb_species = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            clade = f[0]
            parts = clade.split("|")
            if not parts[-1].startswith(marker):
                continue
            name = parts[-1][len(marker):]
            rows[name] = float(f[-1])
            if marker == "t__":
                for part in parts[:-1]:
                    if part.startswith("s__") and part[3:]:
                        sgb_species[name] = part[3:]
    return rows, "counts", {"sgb_species": sgb_species}


PARSERS = {
    "kraken2": parse_kraken2,
    "bracken": parse_bracken,
    "metaphlan": parse_metaphlan,
    "mpa": parse_mpa,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--kind", required=True, choices=list(PARSERS))
    ap.add_argument("--rank", default="S",
                    help="rank code D/K/P/C/O/F/G/S/T (T = SGB/strain level, MetaPhlAn/MPA only)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--metadata")
    ap.add_argument("--min_abundance", type=float)
    ap.add_argument("--top", type=int)
    args = ap.parse_args()

    series = {}
    value_kind = None
    sgb_species = {}
    for path in args.inputs:
        if not Path(path).is_file():
            fail(f"not a file: {path}")
        rows, kind, extra = PARSERS[args.kind](path, args.rank)
        if value_kind is None:
            value_kind = kind
        elif value_kind != kind:
            fail(f"mixed value kinds across inputs: {value_kind} vs {kind}")
        for sgb, sp in extra.get("sgb_species", {}).items():
            sgb_species.setdefault(sgb, sp)
        sample = sample_stem(path)
        if sample in series:
            fail(f"duplicate sample name {sample}; rename files or pass distinct inputs")
        series[sample] = rows
    if not series:
        fail("no inputs parsed")

    matrix = pd.DataFrame(series).fillna(0.0)
    matrix.index.name = "taxon"
    taxa_before = len(matrix)
    if args.min_abundance is not None:
        matrix = matrix[matrix.max(axis=1) >= args.min_abundance]
    if args.top is not None and args.top > 0:
        matrix = matrix.loc[matrix.mean(axis=1).sort_values(ascending=False).head(args.top).index]
    if matrix.empty:
        fail(f"no taxa at rank {args.rank} survived filtering")

    if args.metadata:
        meta = read_table(args.metadata)
        meta.index = meta.index.astype(str)
        missing = [s for s in matrix.columns if s not in meta.index]
        if missing:
            fail(f"metadata missing samples: {missing}")

    matrix.to_csv(args.output, sep="\t")
    artifact = write_artifact(
        args.output, "mg_merge",
        f"Taxa x sample matrix ({args.kind}, rank {args.rank})",
        skill="profiling-workflow", inputs=args.inputs,
        parameters={"kind": args.kind, "rank": args.rank, "value_kind": value_kind},
        shape=matrix.shape,
    )
    emit({
        "kind": args.kind,
        "rank": args.rank,
        "value_kind": value_kind,
        "samples": matrix.columns.tolist(),
        "taxa_before": taxa_before,
        "taxa_after": int(len(matrix)),
        "output": args.output,
        "scientificContext": {
            str(args.output): {
                "valueSemantics": "raw_counts" if value_kind == "counts" else "relative_abundance",
                "sampleIds": [str(s) for s in matrix.columns],
            }
        },
        "artifacts": [artifact],
        **({"sgb_species": sgb_species} if sgb_species else {}),
    })


if __name__ == "__main__":
    main()