"""Fixed, independently constructed results/certificate fixtures.

Built from literal text and hashlib only; never uses deqclib writers.
"""

import hashlib

QC_HEADER = ("schemaVersion\tdisposition\tcontrast\tgeneIdColumn\tpadjThreshold\t"
             "log2fcThreshold\tevaluated\thits\tresultMd5")

# At padj < 0.05 and |lfc| >= 1: evaluated g1,g2,g3; hits g1,g2.
ELIGIBLE_RESULTS = (
    "gene\tlog2FoldChange\tpadj\n"
    "g1\t2\t0.001\n"
    "g2\t-1.5\t0.04\n"
    "g3\t0.5\t0.001\n"
    "g4\t1\tNA\n"
    "g5\tNA\t0.5\n"
).encode("utf-8")

# g1 sits exactly on the padj threshold (not a hit); g2 exactly on |lfc| (hit).
BOUNDARY_RESULTS = (
    "gene\tlog2FoldChange\tpadj\n"
    "g1\t1\t0.05\n"
    "g2\t-1\t0.049\n"
).encode("utf-8")

SKIP_RESULTS = (
    "gene\tlog2FoldChange\tpadj\n"
    "g1\t2\t0.5\n"
    "g2\t-2\t0.9\n"
).encode("utf-8")

ALLNA_RESULTS = (
    "gene\tlog2FoldChange\tpadj\n"
    "g1\tNA\tNA\n"
).encode("utf-8")

UNNAMED_RESULTS = (
    "\tlog2FoldChange\tpadj\n"
    "g1\t2\t0.001\n"
).encode("utf-8")


def md5(data):
    return hashlib.md5(data).hexdigest()


def certificate(results, disposition, evaluated, hits,
                padj="0.05", lfc="1", contrast="treated vs control", gene_col="gene"):
    row = "\t".join(["1", disposition, contrast, gene_col, padj, lfc,
                     str(evaluated), str(hits), md5(results)])
    return (QC_HEADER + "\n" + row + "\n").encode("utf-8")


ELIGIBLE_QC = certificate(ELIGIBLE_RESULTS, "eligible", 3, 2)
BOUNDARY_QC = certificate(BOUNDARY_RESULTS, "eligible", 2, 1)
SKIP_QC = certificate(SKIP_RESULTS, "skip", 2, 0)
ALLNA_QC = certificate(ALLNA_RESULTS, "invalid", 0, 0)
UNNAMED_QC = certificate(UNNAMED_RESULTS, "eligible", 1, 1, gene_col="(unnamed-column-1)")