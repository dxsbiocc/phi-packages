Purpose: screen for ligand-receptor pairs whose expression co-locates in space,
suggesting which cell types may signal to each other globally or within specific
tissue niches. Principle: squidpy ligrec tests whether a ligand in one cluster
and its receptor in a neighbour cluster are expressed together more than a
permutation null, using a curated ligand-receptor database; setting stratify_by
(e.g. niche) repeats the screen inside each tissue microenvironment.
Interpretation: read the significant interactions with their clusters, strata,
and p-values; a hit is a hypothesis about an interface, to be confirmed with the
spatial map and markers. Pitfalls: this is a SCREEN — co-expression is not proof
of signalling, and the curated offline database (~70 pairs) is not exhaustive;
ligrec needs lognorm expression (use_raw=False) and a categorical cluster key;
small niches with too few units or < 2 groups are skipped; omnipath/liana for
richer databases are project-env gaps.
