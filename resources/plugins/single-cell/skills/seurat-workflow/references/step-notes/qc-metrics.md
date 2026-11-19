Purpose: quantify droplet quality before any threshold is chosen. Principle:
cytoplasmic RNA leaks out of a torn membrane while mitochondria stay inside, so
a high mitochondrial fraction marks damage, not biology; genes-per-cell and
UMIs-per-cell separate empty droplets (low) from likely doublets (high).
Interpretation: read the quantiles, not a tutorial's round number — adaptive
rules such as 3 MADs below the median log-depth or above the median logit-mito
fraction (OSCA, sc-best-practices) adapt to each dataset, and mitochondrial caps
are tissue-dependent (blood and tumour samples tolerate more than 10-20%).
Pitfalls: a wrong species regex silently reports 0% and the gate never fires;
ribosomal and hemoglobin fractions are context (cell type, contamination), not
universal gates.
