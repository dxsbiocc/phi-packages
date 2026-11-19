Purpose: keep the genes that carry structure and drop the noise that dilutes
every distance. Principle: vst ranks genes by variance beyond the mean-variance
trend of technical noise; selection runs on counts, not on scaled values.
Interpretation: read the top-10 list and the mean-variance plot — ribosomal,
mitochondrial or one dominant population monopolizing the list is a warning;
with batches, per-batch selection (nbatches) stops one sample choosing every
gene. Pitfalls: too few HVGs merges rare types, too many re-admits noise; 2000
is a starting point, not a law.
