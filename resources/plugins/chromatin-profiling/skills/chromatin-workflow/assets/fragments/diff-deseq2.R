# Differential binding with DiffBind + DESeq2: the standard route for TF and
# narrow-mark ChIP, and for CUT&Tag. DiffBind builds a consensus count matrix,
# normalizes for library size, and runs DESeq2's negative-binomial test per peak.
# A contrast compares two conditions; with no explicit contrast, every condition
# pair with enough replicates is tested.
if (is.null(dba_object)) {
  stop("differential binding needs a DBA object (load-sample-sheet, then count)", call. = FALSE)
}

min_members <- if (is.null(min_members)) 2L else as.integer(min_members)
if (!is.null(contrast) && length(contrast) == 3) {
  dba_object <- DiffBind::dba.contrast(dba_object, contrast = contrast, minMembers = min_members)
} else {
  dba_object <- DiffBind::dba.contrast(dba_object, categories = DiffBind::DBA_CONDITION,
                                       minMembers = min_members)
}
dba_object <- DiffBind::dba.analyze(dba_object, method = DiffBind::DBA_DESEQ2)

report <- DiffBind::dba.report(dba_object, method = DiffBind::DBA_DESEQ2,
                               th = fdr_threshold, fold = lfc_threshold,
                               bUsePval = FALSE)
report_df <- as.data.frame(report)
if (nrow(report_df) > 0 && "FDR" %in% colnames(report_df)) {
  report_df <- report_df[order(report_df$FDR), , drop = FALSE]
}
diff_results <- list(
  method = "DESeq2",
  fdr_threshold = fdr_threshold,
  lfc_threshold = lfc_threshold,
  n_significant = nrow(report_df),
  n_tested = nrow(dba_object$binding),
  table = utils::head(report_df, top_n)
)

ma <- DiffBind::dba.plotMA(dba_object, method = DiffBind::DBA_DESEQ2, th = fdr_threshold)
show_plot(ma, "diff-ma-deseq2")
cat(sprintf("DESeq2: %d of %d peaks differential at FDR < %.3g, |log2FC| > %.3g\n",
            nrow(report_df), nrow(dba_object$binding), fdr_threshold, lfc_threshold))
analysis_steps <- c(analysis_steps, "diff-deseq2")