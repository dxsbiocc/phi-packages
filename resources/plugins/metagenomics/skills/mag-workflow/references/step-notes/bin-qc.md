# Bin quality

**Purpose.** Grade every bin so downstream steps (filter, dereplicate,
classify) operate on a defensible quality floor.

**Principle.** CheckM places bins in a reference genome tree and counts
lineage-specific single-copy marker genes: completeness = fraction of
expected markers present, contamination = fraction present in multiple
copies. CheckM2 replaces marker sets with machine-learning models trained
on genome-wide features — faster, database-lighter, but values are not
interchangeable with CheckM's.

**Interpretation.** MIMAG (Bowers et al. 2017): HQ ≥90% complete, ≤5%
contaminated; MQ ≥50%, ≤10%; below that, LQ. The scatter's shape is
diagnostic: a tight high-completeness/low-contamination cloud means binning
worked; a diagonal trade-off tail means bins merged relatives; a cloud at
low completeness means depth or complexity defeated binning.

**Pitfalls.** Strict MIMAG HQ additionally requires 5S/16S/23S rRNA genes
and ≥18 tRNAs — neither table reports these; the tool's caveat says so and
reports must repeat it. CheckM contamination can exceed 100% (multi-copy
markers); do not clamp it. Small bins (<1 Mbp) have unstable estimates —
few markers, wide error bars.