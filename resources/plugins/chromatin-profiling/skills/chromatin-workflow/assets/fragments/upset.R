# Peak overlap across samples as an UpSet plot: how many sites are private to one
# sample versus shared by two, three, or all. This is the visual companion to the
# reproducibility overlap fractions — a condition whose replicates share most
# peaks is reproducible; samples that share almost nothing differ biologically or
# technically.
if (!is.null(occupancy_matrix) && ncol(occupancy_matrix) >= 2) {
  comb <- ComplexHeatmap::make_comb_mat(occupancy_matrix)
  show_plot(ComplexHeatmap::UpSet(comb), "peak-upset")
  cat(sprintf("UpSet over %d consensus peaks x %d samples\n",
              nrow(occupancy_matrix), ncol(occupancy_matrix)))
} else if (inherits(annotation_result, "csAnno") && have("ggupset")) {
  show_plot(ChIPseeker::upsetplot(annotation_result), "peak-upset")
} else {
  cat("UpSet needs a consensus occupancy matrix (consensus); the single-set annotation ",
      "UpSet fallback also needs ggupset (a ChIPseeker Suggests, not shipped). Skipping\n", sep = "")
}
analysis_steps <- c(analysis_steps, "upset")