Purpose: give every gene equal weight in PCA. Principle: z-scoring each gene
removes the abundance artifact where a handful of highly expressed genes define
the axes; clipping at 10 stops one outlier cell stretching an axis.
Interpretation: scaling is a modelling choice, not a truth — sc-best-practices
argues against automatic z-scaling because it amplifies noise in low-expression
genes. Pitfalls: scaled values are z-scores: marker tests, fold changes and
pseudobulk must never read them (the lognorm layer exists for that);
regress_vars removes covariates but can remove real biology with them.
