Purpose: produce the evidence a cell-type call must rest on. Principle:
rank_genes_groups tests each group against the rest (Wilcoxon by default, the
best-calibrated rank test) on the pre-scale lognorm snapshot; pts says how
specific a gene is, not just how different. Interpretation: read top genes with
their pts and the heatmap; the tables are ranked evidence, not significance.
Pitfalls: clusters were defined from the same data, so p-values are double-
dipped and descriptive only; tests must never run on scaled or batch-corrected
values; avg_log2FC is this comparison's fold change, not a portable effect size.
