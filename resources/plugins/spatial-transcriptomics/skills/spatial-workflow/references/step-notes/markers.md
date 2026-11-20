Purpose: turn each spatial domain from a cluster ID into a biological hypothesis.
Principle: Wilcoxon ranks genes by how specifically they distinguish one domain
from the rest, on the pre-scale lognorm snapshot (scaled z-scores have no
meaningful fold change). Interpretation: read the top genes and log2 fold changes
per domain, and map them in space (spatial-plots) to confirm anatomy. Pitfalls:
the p-values are DESCRIPTIVE only — the domains were defined from the same data
(double dipping), so treat the ranking as evidence, not significance; a wrong
identifier type (Ensembl vs symbol) makes every marker disappear. avg_log2FC is
not a biological effect size without the test and comparison it came from.