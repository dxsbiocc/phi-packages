Purpose: make cells comparable despite orders-of-magnitude depth differences.
Principle: LogNormalize rescales each cell to the same total then log1p-stabilizes
variance; SCTransform (the commented alternative) models counts with regularized
negative-binomial residuals and handles depth variation better on some datasets.
Interpretation: the data slot now holds log expression — the mean of logs is not
the log of the mean, so group averages here are geometric-mean-like, not raw
fold changes. Pitfalls: never normalize twice; keep raw counts for QC, HVG and
pseudobulk, and remember log values from different datasets are not directly
comparable.
