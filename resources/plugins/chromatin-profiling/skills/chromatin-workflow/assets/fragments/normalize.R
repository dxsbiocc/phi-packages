# Library-size normalization removes sequencing-depth differences so samples are
# comparable. CPM (counts per million) is the default for occupancy visualisation
# and PCA; RPKM additionally divides by peak width, which matters when comparing
# peaks of very different lengths. The differential engines (DESeq2/edgeR/csaw)
# do their own internal normalization — this matrix is for exploration, not the
# test itself.
if (is.null(count_matrix)) stop("count reads before normalizing (count or load-counts)", call. = FALSE)

lib_size <- colSums(count_matrix)
if (any(lib_size == 0)) {
  cat("warning: some samples have zero total counts; they will not normalize meaningfully\n")
}
method <- if (is.null(normalize_method)) "libsize" else normalize_method
if (method == "none") {
  normalized_counts <- count_matrix
} else if (method == "rpkm") {
  if (is.null(consensus_peaks)) stop("RPKM needs consensus peak widths (consensus)", call. = FALSE)
  width_kb <- width(consensus_peaks) / 1000
  cpm <- t(t(count_matrix) / lib_size) * 1e6
  normalized_counts <- cpm / width_kb
} else {
  normalized_counts <- t(t(count_matrix) / lib_size) * 1e6
}

# DiffBind normalizes the DBA object for dba.analyze (library-size by default).
# method= selects the differential engine and normalize= the normalization
# scheme; DBA_NORM_LIB is a normalization token, so it goes to normalize=.
if (!is.null(dba_object)) {
  dba_object <- DiffBind::dba.normalize(dba_object, normalize = DiffBind::DBA_NORM_LIB)
}
cat(sprintf("normalized counts (%s): %d peaks x %d samples\n",
            method, nrow(normalized_counts), ncol(normalized_counts)))
analysis_steps <- c(analysis_steps, "normalize")