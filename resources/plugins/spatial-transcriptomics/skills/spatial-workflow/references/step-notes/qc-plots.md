Purpose: see the QC distributions and their spatial layout before choosing
thresholds. Principle: violin panels show the shape of each metric (a low
shoulder of empty units, a high tail of merged bins); the same metrics mapped
onto tissue show whether quality follows anatomy. Interpretation: a QC gradient
that follows a known structure (low counts at the tissue edge, high
mitochondrial in a necrotic core) is biology to interpret, not a threshold to
filter away; a scattered low-quality patch off-tissue is what filtering removes.
Read these next to the tissue image when one exists. Pitfalls: filtering on a
metric that tracks anatomy erases real regions; plot_tag is "before" on the first
pass and "after" once filtering has run, so compare the two passes rather than
trusting one.