Purpose: remove technical between-sample variation so clusters are not
one-per-donor. Principle: Harmony corrects the PCA embedding iteratively;
ComBat fits a linear batch model on log-normalized expression and needs a fresh
PCA; bbknn (project environment) corrects the neighbour graph instead.
Interpretation: judge on the before/after UMAPs — samples should mix within each
cell type while markers and known populations survive; benchmark suites (scIB)
weight batch removal and biology conservation jointly. Pitfalls: over-
correction deletes real biology (patient-specific tumour states are the classic
case); corrected values must never feed differential expression — test on
uncorrected lognorm with sample as the unit; without a batch column the cell
says so and the chain continues uncorrected.
