# A DiffBind sample sheet builds the DBA object: the container for occupancy and
# affinity analysis of TF / narrow-mark ChIP and CUT&Tag. The sheet is a CSV with
# one row per sample; DiffBind reads the peak files named in the Peaks column.
#
# Standard columns: SampleID, Tissue, Factor, Condition, Treatment, Replicate,
# bamReads, ControlID, bamControl, Peaks, PeakCaller. ControlID/bamControl are
# optional (input control); bamReads is needed for counting and FRiP.
path <- if (!is.null(sample_sheet_path) && nzchar(sample_sheet_path)) sample_sheet_path else data_path
if (!file.exists(path)) stop("sample sheet not found: ", path, call. = FALSE)
sample_sheet_path <- path

sample_sheet <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("SampleID", "Condition", "Replicate", "Peaks", "PeakCaller")
missing <- setdiff(required, colnames(sample_sheet))
if (length(missing) > 0) {
  stop("sample sheet is missing required columns: ", paste(missing, collapse = ", "),
       ". Expected a DiffBind-format CSV (see references/analysis-playbook.md).",
       call. = FALSE)
}
if (!"bamReads" %in% colnames(sample_sheet)) {
  cat("note: no bamReads column; occupancy analysis works, but counting/FRiP need BAM paths\n")
}

dba_object <- DiffBind::dba(sampleSheet = sample_sheet)
cat(sprintf("DBA object: %d samples, %d peaks in the binding matrix\n",
            nrow(dba_object$samples), nrow(dba_object$binding)))
cat("conditions:", paste(unique(dba_object$samples$Condition), collapse = ", "), "\n")
analysis_steps <- c(analysis_steps, "load-sample-sheet")