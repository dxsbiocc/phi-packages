# Sample-level PCA on log-normalized counts: the first check of whether the
# dominant variation is the biology (condition) or a nuisance (batch, depth,
# antibody lot). Replicates of one condition should sit together; a split along
# PC1 that tracks condition is the signal the differential test will quantify.
if (is.null(normalized_counts)) stop("normalize before PCA (normalize)", call. = FALSE)

log_counts <- log2(normalized_counts + 1)
keep <- apply(log_counts, 1, stats::var) > 0
pca <- stats::prcomp(t(log_counts[keep, , drop = FALSE]), center = TRUE, scale. = FALSE)
variance <- round(summary(pca)$importance[2, ], 4)

scores <- as.data.frame(pca$x)
scores$sample <- rownames(scores)
if (!is.null(sample_sheet) && all(c("SampleID", "Condition") %in% colnames(sample_sheet))) {
  mapped <- sample_sheet$Condition[match(scores$sample, sample_sheet$SampleID)]
  scores$condition <- ifelse(is.na(mapped), "unknown", mapped)
} else {
  scores$condition <- "all"
}
pca_result <- list(scores = pca$x, variance_explained = variance, n_peaks = sum(keep))

axis_label <- function(i) sprintf("PC%d (%.1f%%)", i, 100 * variance[[i]])
plot <- ggplot(scores, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample), vjust = -0.9, size = 3, check_overlap = TRUE) +
  labs(title = "Sample PCA (log-normalized counts)",
       x = axis_label(1), y = axis_label(2), colour = "Condition") +
  theme_minimal()
show_plot(plot, "sample-pca")
cat(sprintf("PCA on %d variable peaks; PC1 %.1f%%, PC2 %.1f%%\n",
            sum(keep), 100 * variance[[1]], 100 * variance[[2]]))
analysis_steps <- c(analysis_steps, "pca")