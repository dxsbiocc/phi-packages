Purpose: make cells comparable despite orders-of-magnitude depth differences.
Principle: normalize_total rescales each cell to the same total, then log1p
stabilizes variance so a few abundant genes do not dominate every distance.
Interpretation: the lognorm layer snapshots this pre-scale state — marker tests
and panel scoring read it, while counts keeps raw UMIs; the mean of logs is not
the log of the mean, so group averages here are not raw fold changes. Pitfalls:
never normalize twice; log values from different datasets are not directly
comparable, and scaling (next cell) must not leak into any test.
