Purpose: judge whether a batch effect exists and whether integration removed it
without erasing biology. Principle: batch is technical variation (donor,
chemistry, run); colouring the UMAP by sample shows it as one-donor-per-corner
structure. Interpretation: after correction, samples should mix within each cell
type while marker genes and known populations survive; benchmark suites (scIB)
weight batch removal and biology conservation jointly. Pitfalls: over-correction
deletes real biology (patient-specific tumour states are the classic case), and
corrected/integrated values must never feed differential expression — test on
the uncorrected log-normalized counts with sample as the unit.
