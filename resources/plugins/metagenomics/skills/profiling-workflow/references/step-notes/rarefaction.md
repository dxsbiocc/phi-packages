# Rarefaction curves

**Purpose.** Check whether sequencing depth is sufficient to capture richness
before comparing alpha diversity across samples. Runs on the merged count
matrix (`mg_merge` kraken2/bracken output), between merge and `mg_diversity`.

**Principle.** Each library is subsampled without replacement to a grid of
depths (sequential hypergeometric draws, equivalent to drawing reads from the
whole community); observed richness is averaged over `--iterations` seeded
draws, with the spread shown as a band. The dashed vertical line marks the
smallest library size — the deepest point at which every sample can still be
compared on equal footing.

**Interpretation.** A curve that flattens means deeper sequencing would add
few new taxa at that depth. Samples listed in `not_plateaued` still gain >5%
richness over the last depth step: their observed richness is depth-limited,
so raw observed counts are not comparable across libraries. Either compare at
a common depth, or report an estimator (`mg_diversity` Chao1) alongside
observed richness.

**Pitfalls.** Integer counts only — MetaPhlAn relative abundances and Bracken
fractional `new_est_reads` are rejected with guidance (round estimates first
if you must). Rarefaction throws data away: never rarefy before differential
testing; `mg_differential` works on the full matrix. The plateau flag is a
coarse last-two-depths heuristic, not a formal coverage estimate, and the
curve metric is observed richness only (not Shannon).