Purpose: compare spatial organisation between conditions/samples — does a gene's
spatial pattern, or two types' colocalization, differ across groups? Principle:
the same statistic (spatial autocorrelation for SVG, co_occurrence for
colocalization) is computed per sample, then contrasted between groups, so a
difference reflects changed spatial structure rather than changed mean
expression. Interpretation: read the per-gene/per-pair cross-group difference and
its consistency across replicates; map the pattern in each group to see what
changed. Pitfalls: the SAMPLE is the experimental unit, not the spot — with fewer
than ~2 biological replicates per group there is no basis for a condition-level
claim; differences can be driven by composition or depth shifts, so inspect the
maps; this is a screen for follow-up, not a definitive differential test.