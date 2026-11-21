# Window-based differential binding with csaw — the preferred route for broad
# marks (H3K27me3, H3K36me3, H3K9me3) whose signal spreads over kilobases rather
# than sharp peaks. csaw counts reads in sliding windows, filters out background
# windows, normalizes for composition bias, fits a quasi-likelihood glm, then
# combines window-level tests into region-level FDR. Peak-based DiffBind can miss
# diffuse changes that this window approach catches.
if (!have("csaw")) {
  csaw_result <- list(error = paste(
    "csaw is not available on this platform (bioconda has no osx-arm64 build).",
    "Use the DiffBind route for broad marks, or run on linux-x64 / darwin-x64.",
    "A project environment (env_request) can provide csaw on Apple Silicon."))
  cat(csaw_result$error, "\n")
} else {
  need("csaw")
  if (is.null(bam_paths) || length(bam_paths) == 0) {
    stop("csaw needs BAM files (load-bam)", call. = FALSE)
  }
  if (is.null(sample_sheet) || !"Condition" %in% colnames(sample_sheet)) {
    stop("csaw needs a sample sheet with a Condition column to build the design", call. = FALSE)
  }
  group <- sample_sheet$Condition[match(tools::file_path_sans_ext(basename(bam_paths)),
                                        sample_sheet$SampleID)]
  if (anyNA(group)) stop("every BAM must match a SampleID in the sample sheet", call. = FALSE)
  group <- factor(group)
  if (nlevels(group) < 2) stop("csaw needs at least two conditions to test", call. = FALSE)

  pe_mode <- if (isTRUE(paired_end)) "both" else "none"
  param <- csaw::readParam(pe = pe_mode, dedup = TRUE, minq = 10)
  ext <- if (isTRUE(paired_end)) NA else fragment_size
  data <- csaw::windowCounts(bam_paths, width = csaw_width, spacing = csaw_spacing,
                             ext = ext, param = param)

  # Filter low-abundance windows against the global background (broad marks).
  keep <- csaw::filterWindows(data, type = csaw_filter)
  filtered <- data[keep$filter, ]

  # TMM normalization factors from large bins correct composition bias.
  binned <- csaw::windowCounts(bam_paths, bin = TRUE, width = 10000, param = param)
  normfacs <- csaw::normFactors(binned, se.out = FALSE)

  design <- stats::model.matrix(~ group)
  y <- csaw::asDGEList(filtered, norm.factors = normfacs)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  res <- edgeR::glmQLFTest(fit, coef = ncol(design))

  # Merge adjacent windows into regions and combine their p-values (Simes).
  merged <- csaw::mergeWindows(rowRanges(filtered), tol = csaw_spacing)
  comb <- csaw::combineTests(merged$id, res$table)
  regions <- merged$region
  S4Vectors::mcols(regions) <- comb

  region_df <- as.data.frame(regions)
  if ("FDR" %in% colnames(region_df)) region_df <- region_df[order(region_df$FDR), , drop = FALSE]
  sig <- if ("FDR" %in% colnames(region_df)) sum(region_df$FDR < fdr_threshold, na.rm = TRUE) else 0L
  csaw_result <- list(
    method = "csaw",
    window_width = csaw_width,
    window_spacing = csaw_spacing,
    filter = csaw_filter,
    n_windows = nrow(filtered),
    n_regions = length(regions),
    n_significant = as.integer(sig),
    fdr_threshold = fdr_threshold,
    table = utils::head(region_df, top_n),
    note = if (nlevels(group) > 2)
      "more than two conditions: coef tests the last level vs the first; set explicit contrasts for others"
      else NULL
  )
  cat(sprintf("csaw: %d windows -> %d regions, %d differential at FDR < %.3g\n",
              nrow(filtered), length(regions), sig, fdr_threshold))
}
analysis_steps <- c(analysis_steps, "csaw")