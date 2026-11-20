"""ORA: gated over-representation analysis over certified results + explicit GMT.

Revalidates the certificate and exact result bytes on every call. Eligible
only; never computes or creates output for zero/invalid/missing input.
One-sided hypergeometric test, BH across all tested terms. Stdlib only.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deqclib as lib


def blocked(reasons, choices, **extra):
    obj = {"status": "blocked", "reasons": reasons, "choices": choices}
    obj.update(extra)
    lib.emit(obj)


def main():
    ap = argparse.ArgumentParser(description="Gated ORA over certified DESeq2 results")
    ap.add_argument("--results", required=True)
    ap.add_argument("--qc", required=True)
    ap.add_argument("--gmt", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min_size", type=int, default=5)
    ap.add_argument("--max_size", type=int, default=500)
    args = ap.parse_args()

    if args.min_size < 1 or args.max_size < args.min_size:
        blocked(["require 1 <= min_size <= max_size"],
                        ["choose valid gene-set size bounds, e.g. --min_size 5 --max_size 500"])

    verdict = lib.verified_evidence(args.results, args.qc)
    if verdict["disposition"] != "eligible":
        reasons = verdict["reasons"] or ["No genes met the prespecified thresholds; candidate-list ORA was skipped."]
        lib.emit({"status": verdict["disposition"], "disposition": verdict["disposition"], "reasons": reasons,
                  "evaluated": verdict["evaluated"], "hits": verdict["hits"], "outputs": []})

    try:
        with open(args.gmt, "rb") as fh:
            gmt_data = fh.read()
    except OSError:
        blocked(["GMT file not found: %s" % args.gmt],
                        ["provide an explicit GMT gene-set file for this organism/identifier space"])
    try:
        gmt_text = gmt_data.decode("utf-8")
    except UnicodeDecodeError:
        blocked(["GMT is not valid UTF-8"], ["re-export the GMT as UTF-8"])
    terms = {}
    for lineno, raw in enumerate(gmt_text.split("\n"), start=1):
        line = raw[:-1] if raw.endswith("\r") else raw
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3 or not fields[0].strip():
            blocked(["malformed GMT line %d: need name, description, and >= 1 gene" % lineno],
                            ["fix the GMT format: name<TAB>description<TAB>gene1<TAB>gene2..."])
        name = fields[0].strip()
        if name in terms:
            blocked(["the GMT repeats a gene-set name on line %d" % lineno],
                            ["make gene-set names unique"])
        genes = {g.strip() for g in fields[2:] if g.strip() and g.strip() != "NA"}
        terms[name] = (fields[1], genes)
    if not terms:
        blocked(["GMT contains no gene sets"], ["provide a non-empty GMT"])

    annotated = set()
    for _, genes in terms.values():
        annotated |= genes
    universe = verdict["evaluable_ids"] & annotated
    mapped_hits = verdict["hit_ids"] & annotated
    if not mapped_hits:
        lib.emit({"status": "skip",
                  "reasons": ["no significant genes map to the annotated GMT universe; enrichment is not computed"],
                  "choices": ["provide a GMT covering this organism and gene-identifier space",
                              "deliver the gene screen and report the mapping limitation"],
                  "universe_size": len(universe), "mapped_hits": 0,
                  "inputs": [{"id": "results", "path": lib.display_path(args.results)},
                             {"id": "research_qc", "path": lib.display_path(args.qc)},
                             {"id": "gmt", "path": lib.display_path(args.gmt)}],
                  "outputs": []})

    N = len(universe)
    n = len(mapped_hits)
    tested = []
    for name, (description, genes) in terms.items():
        in_universe = genes & universe
        K = len(in_universe)
        if args.min_size <= K <= args.max_size:
            k = len(genes & mapped_hits)
            tested.append((name, description, K, k, genes & mapped_hits))
    if not tested:
        blocked(["no gene sets fall within size bounds [%d, %d] after universe intersection"
                         % (args.min_size, args.max_size)],
                        ["adjust --min_size/--max_size", "provide a GMT with better universe coverage"])

    pvals = [lib.hypergeom_sf_ge(k, N, K, n) for _, _, K, k, _ in tested]
    qvals = lib.bh_adjust(pvals)
    rows = sorted(zip(tested, pvals, qvals), key=lambda t: (t[1], t[0][0]))

    out = args.out
    lines = ["gene_set\tdescription\tset_size\toverlap\tp_value\tp_adjust\toverlap_genes"]
    for (name, description, K, k, overlap), p, q in rows:
        lines.append("\t".join([name, description, str(K), str(k),
                                format(p, ".6g"), format(q, ".6g"),
                                ";".join(sorted(overlap))]))
    try:
        write_state = lib.write_text_atomic(out, "\n".join(lines) + "\n")
    except lib.OutputConflict:
        blocked(
            ["a different ORA output already exists at %s; refusing to silently overwrite another run" % out],
            ["choose a different --out path for this run",
             "keep the existing file"])
    except OSError:
        blocked(["cannot write ORA output"], ["choose a writable --out project path"])

    lib.emit({
        "status": "ok",
        "contrast": verdict["contrast"],
        "universe_size": N,
        "mapped_hits": n,
        "terms_tested": len(rows),
        "terms_q05": sum(1 for _, _, q in rows if q <= 0.05),
        "write": write_state,
        "inputs": [{"id": "results", "path": lib.display_path(args.results)},
                   {"id": "research_qc", "path": lib.display_path(args.qc)},
                   {"id": "gmt", "path": lib.display_path(args.gmt)}],
        "outputs": [{"id": "ora", "type": "ora", "primary": True, "path": lib.display_path(out)}],
    }, 0)


if __name__ == "__main__":
    main()