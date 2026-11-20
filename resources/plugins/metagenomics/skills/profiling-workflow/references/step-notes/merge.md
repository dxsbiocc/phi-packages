# Merge profiles

**Purpose.** One taxa × sample matrix per kind and rank — the single source
for diversity, differential and figures.

**Principle.** Kraken2 column 3 is clade-assigned reads (a read counted at
the lowest common ancestor); Bracken re-estimates species reads Bayesianly
from those assignments — prefer Bracken for abundance questions. MetaPhlAn
reports marker-gene relative abundance (0–100), a fundamentally different
estimator: do not mix with k-mer counts in one matrix.

**Interpretation.** `taxa_before` vs `taxa_after` shows filter impact. Many
`[taxid]`-suffixed duplicate names mean the database contains same-named taxa
across lineages; the suffix keeps them distinct.

**Pitfalls.** Bracken needs the correct read length and a database built for
the same Kraken2 library. MetaPhlAn `t__` (strain) rows are excluded at S
rank by design. Missing taxa become 0 — true absence vs below-detection is
not distinguishable; say so when interpreting.