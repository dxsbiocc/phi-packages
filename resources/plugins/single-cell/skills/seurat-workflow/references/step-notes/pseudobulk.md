Purpose: move the comparison from cells to experimental units. Principle: cells
from one donor are not independent observations; summing RAW counts per
sample x cell type x condition restores the donor as the unit, and edgeR/DESeq2
quasi-likelihood tests on those sums control false discovery. Interpretation:
read the design sheet and the notes — thin groups (<2 samples) cannot support a
contrast, and dropped-gene counts say what the expression filter removed;
filterByExpr-style filtering belongs in the wrapper. Pitfalls: cell-level tests
(pseudoreplication) inflate significance dramatically; power comes from more
donors, not more cells per donor; sex genes and ambient RNA are classic
confounds in the resulting lists.
