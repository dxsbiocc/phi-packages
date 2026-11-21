Purpose: sample-level PCA on log-normalized counts — the first check of whether
the dominant variation is biology (condition) or a nuisance (batch, depth,
antibody lot). Principle: log2(normalized+1) is computed over variable peaks
(variance>0), then prcomp on the transposed matrix (samples as rows); samples are
colored by the sheet's Condition. Interpretation: replicates of one condition
should sit together; a PC1 split that tracks condition is the signal the
differential test will quantify, while a split tracking depth/batch is a confound
to address first. variance_explained tells how much PC1/PC2 carry. Pitfalls: PCA
is visualization, not a differential test — it carries no significance and models
no counts. It needs normalized counts (normalize first); without a Condition
column every sample is labeled "all" and the biology-vs-batch read is lost.