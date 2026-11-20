# Differential abundance

**Purpose.** Which features (taxa or functional groups) differ between
groups, with effect sizes and multiple-testing control.

**Principle.** Mann-Whitney U (two groups) tests whether one group's
abundances tend to exceed the other's — rank-based, robust to the
zero-inflation and skew of microbiome data. Kruskal-Wallis extends to >2
groups (omnibus only). BH FDR controls the expected false-discovery
proportion across taxa. log2FC uses a pseudocount (0.5 counts / 1e-6
relative) so zeros do not explode; effect r = |Z|/√N.

**Interpretation.** Significance = FDR < threshold AND |log2FC| ≥ threshold
(two-group). A taxon can be significant with tiny effect at large N — read
both columns. The volcano's dashed lines are your thresholds, not truth.

**Pitfalls.** Compositionality: one bloom mechanically depresses everything
else; Wilcoxon on relative abundances cannot distinguish bloom from
depletion — consider count-based tools or ANCOM-BC (project environment) for
strong claims. The same applies to functional tables from mg_function
(pass `values: relative`); regroup sparse gene families to KO before
testing. LEfSe's LDA effect size and cladogram are a different method
(lefse + graphlan wrappers); this tool deliberately does not imitate them.