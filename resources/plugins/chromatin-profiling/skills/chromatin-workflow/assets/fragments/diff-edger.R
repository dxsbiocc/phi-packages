# Differential binding with DiffBind + edgeR: the alternative to DESeq2 for TF /
# narrow-mark ChIP and CUT&Tag. edgeR's exact test / QL framework is often more
# sensitive with few replicates; DESeq2 is more conservative. The choice is the
# analyst's — both run on the same DiffBind consensus count matrix, so switching
# method does not change the peaks being tested.
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
dba_object <- DiffBind::dba.analyze(dba_object, method = DiffBind::DBA_EDGER)

report <- DiffBind::dba.report(dba_object, method = DiffBind::DBA_EDGER,
                               th = fdr_threshold, fold = lfc_threshold,
                               bUsePval = FALSE)
report_df <- as.data.frame(report)
if (nrow(report_df) > 0 && "FDR" %in% colnames(report_df)) {
  report_df <- report_df[order(report_df$FDR), , drop = FALSE]
}
diff_results <- list(
  method = "edgeR",
  fdr_threshold = fdr_threshold,
  lfc_threshold = lfc_threshold,
  n_significant = nrow(report_df),
  n_tested = nrow(dba_object$binding),
  table = utils::head(report_df, top_n)
)

ma <- DiffBind::dba.plotMA(dba_object, method = DiffBind::DBA_EDGER, th = fdr_threshold)
show_plot(ma, "diff-ma-edger")
cat(sprintf("edgeR: %d of %d peaks differential at FDR < %.3g, |log2FC| > %.3g\n",
            nrow(report_df), nrow(dba_object$binding), fdr_threshold, lfc_threshold))
analysis_steps <- c(analysis_steps, "diff-edger")