Purpose: write the proposed labels back into obs and map them in space.
Principle: labels come either from an explicit mapping JSON (cluster -> label)
or from the highest-scoring panel type per domain; a domain below min_score stays
"unassigned". Interpretation: edit the mapping by hand if your reading of the dot
plot disagrees with the automatic call — that judgement IS the analysis. Read the
spatial cell-type map: a label should paint a coherent region. Pitfalls: a
reference panel forces a label onto every cluster, so an unassigned domain is
information, not a failure; in spatial data a label names the dominant cell type
or region of a spot/bin, not a single pure cell. A mapping must cover every
cluster or it is rejected.