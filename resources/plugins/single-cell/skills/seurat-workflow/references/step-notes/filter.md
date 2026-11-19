Purpose: remove empty, damaged and doublet droplets. Principle: each gate is a
statement about the QC metrics (genes/UMIs low = empty, high = doublet, mito
high = stressed); gates apply in sequence and report their own removed count so
a surprising total traces to one threshold. Interpretation: kept_fraction and
the per-gate counts are the audit record; below half surviving, stop and
reconsider. Pitfalls: thresholds are tissue- and protocol-dependent; aggressive
filtering deletes real rare populations, and filtering must happen on raw
counts, never on normalized or scaled values.
