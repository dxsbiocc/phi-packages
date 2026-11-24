Purpose: map continuous state transitions or developmental trajectories onto
tissue coordinates, and optionally score pathway/signature gene programs in
space. Principle: PAGA quantifies coarse-grained connectivity between clusters
on the expression manifold, diffusion pseudotime (DPT) orders units from a root
cluster along that manifold and projects pseudotime back onto obsm[coord_key],
and sc.tl.score_genes evaluates user-supplied gene sets on the unscaled lognorm
layer. Interpretation: compare PAGA connectivity with spatial adjacency to see
whether transcriptional transitions align with anatomical axes (e.g. germinal
centre dark-to-light zone); inspect spatial gene-program maps for localized
pathway activity. Pitfalls: without an explicit root_cluster the script picks
the first category automatically (flagged in notes) — always set root_cluster
to a biologically justified starting state; DPT assumes a connected manifold
and can mislead across unrelated cell lineages; never score gene sets on scaled
expression.
