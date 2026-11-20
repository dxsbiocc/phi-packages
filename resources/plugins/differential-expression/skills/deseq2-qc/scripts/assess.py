"""Assess: read-only verification of DESeq2 results against a QC certificate.

Never writes files. Compact JSON on stdout; exit 0 for every scientific
decision (disposition eligible/skip/invalid).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deqclib as lib

DECISIONS = {
    "eligible": "QC eligible. Run ORA only if the user explicitly authorized enrichment "
                "(write approval) and the original question asks for it; otherwise deliver the QC report.",
    "skip": "No genes pass the thresholds (evaluated > 0, hits = 0). Do not run ORA. Deliver the valid zero-hit screen. Ask only if an unresolved research goal requires choosing a different analysis.",
    "invalid": "Results are invalid or unevaluable. Do not run ORA. Report the reasons and request "
               "corrected results from the wrapper.",
}


def main():
    ap = argparse.ArgumentParser(description="Read-only DESeq2 results QC verification")
    ap.add_argument("--results", required=True)
    ap.add_argument("--qc", required=True)
    args = ap.parse_args()

    ev = lib.verified_evidence(args.results, args.qc)
    lib.emit({
        "status": "ok",
        "disposition": ev["disposition"],
        "evaluated": ev["evaluated"],
        "hits": ev["hits"],
        "geneIdColumn": ev["geneIdColumn"],
        "contrast": ev["contrast"],
        "padjThreshold": ev["padjThreshold"],
        "log2fcThreshold": ev["log2fcThreshold"],
        "resultMd5": ev["resultMd5"],
        "reasons": ev["reasons"],
        "decision": DECISIONS[ev["disposition"]],
        "inputs": [{"id": "results", "path": lib.display_path(args.results)},
                   {"id": "research_qc", "path": lib.display_path(args.qc)}],
        "outputs": [],
    }, 0)


if __name__ == "__main__":
    main()