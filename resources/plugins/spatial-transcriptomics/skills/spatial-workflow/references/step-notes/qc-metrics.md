Purpose: quantify spatial-unit quality before any threshold is chosen.
Principle: cytoplasmic RNA leaks out of a torn membrane while mitochondria stay
inside, so a high mitochondrial fraction marks damage, not biology;
genes-per-unit and UMIs-per-unit separate empty/off-tissue units (low) from
merged bins or dense tissue (high). Interpretation: read the quantiles, not a
tutorial's round number — spatial tissue tolerates more mitochondrial signal
than dissociated cells (mt_max ~20 is a starting point, not a rule), and the low
shoulder is often tissue edge rather than empty droplets. Pitfalls: a wrong
species regex silently matches no gene, pct_counts_mt is absent and the gate is
skipped (use ^mt- for mouse); ribosomal and hemoglobin fractions are context
(region, blood contamination), not universal gates. There is no doublet gate —
a spot/bin is not a droplet.