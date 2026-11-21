# Replicate reproducibility: do biological replicates of the same condition call
# the same peaks? The formal test is IDR (irreproducible discovery rate), which
# needs the `idr` pip package this environment does not ship — a reported gap.
# Here we give a transparent proxy: the fraction of each sample's peaks that
# overlap another sample, plus a count correlation when a count matrix exists.
if (is.null(peak_list) || length(peak_list) < 2) {
  cat("reproducibility needs at least two peak sets; skipping\n")
} else {
  samples <- names(peak_list)
  n <- length(samples)
  overlap <- matrix(0, n, n, dimnames = list(samples, samples))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) { overlap[i, j] <- 1; next }
      a <- peak_list[[i]]
      overlap[i, j] <- round(sum(IRanges::overlapsAny(a, peak_list[[j]])) /
                               max(1, length(a)), 4)
    }
  }
  reproducibility <- list(overlap_fraction = overlap)

  # Same-condition pairs are the replicates whose overlap matters most.
  if (!is.null(sample_sheet) && "Condition" %in% colnames(sample_sheet)) {
    pairs <- list()
    for (condition in unique(sample_sheet$Condition)) {
      members <- sample_sheet$SampleID[sample_sheet$Condition == condition]
      members <- intersect(members, samples)
      if (length(members) >= 2) {
        for (k in seq_len(length(members) - 1)) {
          for (l in (k + 1):length(members)) {
            pairs <- c(pairs, list(list(
              condition = condition, a = members[[k]], b = members[[l]],
              overlap_fraction = overlap[members[[k]], members[[l]]])))
          }
        }
      }
    }
    reproducibility$replicate_pairs <- pairs
    for (pair in pairs) {
      cat(sprintf("%s replicates %s / %s: %.1f%% peak overlap\n",
                  pair$condition, pair$a, pair$b, 100 * pair$overlap_fraction))
    }
  }
  if (!is.null(count_matrix)) {
    reproducibility$count_correlation <- round(stats::cor(count_matrix, method = "spearman"), 4)
  }
  qc_metrics$reproducibility <- reproducibility
  cat("note: IDR is the formal reproducibility test; the idr pip package is a gap ",
      "(project environment via env_request). Overlap fraction is a proxy.\n", sep = "")
}
analysis_steps <- c(analysis_steps, "qc-reproducibility")