Purpose: build the DiffBind DBA object from a sample sheet — the container for
occupancy and affinity analysis of TF/narrow-mark ChIP and CUT&Tag. Principle: the
sheet is a CSV with one row per sample; DiffBind reads the peak files named in the
Peaks column and stores the binding matrix. Standard columns: SampleID, Tissue,
Factor, Condition, Treatment, Replicate, bamReads, ControlID, bamControl, Peaks,
PeakCaller. Interpretation: required columns are SampleID, Condition, Replicate,
Peaks, PeakCaller; a missing one aborts with the expected format. ControlID/
bamControl are optional (input control); bamReads is needed for counting and FRiP,
and its absence is noted, not fatal. It reports the sample count, binding-matrix
peak count and the conditions. Pitfalls: relative Peaks/bamReads paths resolve
against the sheet directory; the Condition column drives contrasts and csaw's
design, so mislabeled conditions silently test the wrong comparison.