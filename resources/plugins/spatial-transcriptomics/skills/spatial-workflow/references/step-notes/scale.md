Purpose: put every gene on equal footing before PCA. Principle: PCA weights
genes by variance, so without scaling a handful of abundant genes define the
axes; scaling to unit variance (clipped at 10) lets each gene contribute, and the
clip stops one outlier unit from stretching an axis. Interpretation: scaling
changes X into z-scores — this is why the lognorm snapshot was taken first, so
marker tests and SVG still read un-scaled log expression. Pitfalls: scoring
panels or testing markers on scaled X is meaningless (a z-score "fold change" is
not a fold change); the clip value trades outlier suppression against signal in
rare, high-expressing domains.