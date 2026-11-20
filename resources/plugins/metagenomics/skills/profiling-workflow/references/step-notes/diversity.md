# Diversity

**Purpose.** Within-sample richness/evenness (alpha), between-sample
dissimilarity (beta), ordination (PCoA) and a group-level test (PERMANOVA).

**Principle.** Shannon = −Σ p ln p (evenness-sensitive); Simpson = 1 − Σ p²
(dominance-sensitive); Chao1 = S_obs + f1²/(2 f2) estimates unseen species
from singletons/doubletons — counts only. Bray-Curtis = Σ|x−y| / Σ(x+y)
ignores joint absences, the standard for community composition. PCoA is
classical MDS on the distance matrix; PERMANOVA partitions distance sums of
squares between/within groups and permutes labels (Anderson 2001).

**Interpretation.** Compare alpha only at comparable sequencing depth.
PCoA axis percentages tell you how much of the structure the 2D plot
captures. PERMANOVA p < 0.05 with visible separation → groups differ in
composition; significant PERMANOVA with overlapping clouds can be a
dispersion effect — check group spread before claiming location shifts.

**Pitfalls.** Relative abundances (MetaPhlAn) make Chao1/observed meaningless
— the tool skips them. PERMANOVA assumes exchangeable samples; nested designs
(patient within hospital) need restricted permutations, which this tool does
not do — report the limitation. Seed is fixed for reproducibility; changing
it changes p slightly.