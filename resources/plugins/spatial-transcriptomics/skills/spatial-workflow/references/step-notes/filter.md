Purpose: remove low-quality spatial units with thresholds read from the data.
Principle: each gate (min_genes, min_counts, max_counts, mt_max) is applied in
sequence and reports how many units IT removed, so a surprising total traces to
one threshold; a None bound skips that gate. Interpretation: read kept_fraction
— below half usually means a threshold is too strict or the bin is too small,
not that the tissue is half empty. Re-read the distributions and the tissue QC
map before continuing. Pitfalls: an empty result (every unit removed) is a
configuration error, not a finding; the mitochondrial gate is skipped when no
gene matched mt_pattern, so a mouse dataset run with the human ^MT- default
loses that gate silently. Counts, not normalized values, are what QC gates on.