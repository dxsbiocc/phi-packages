# Doublets are two cells captured in one droplet. They inflate clusters and
# invent cell types that do not exist, so remove them before clustering.
# scDblFinder simulates artificial doublets from the data and scores each real
# cell by how much it looks like one; the expected rate should reflect how many
# cells were LOADED on the instrument, not how many survived QC.
doublet_result <- NULL
if (detect_doublets) {
  if (have("scDblFinder") && have("SingleCellExperiment")) {
    set.seed(123)  # the method samples; the seed makes the call repeatable
    rate <- doublet_expected_rate
    if (is.null(rate)) rate <- max(0.01, min(0.15, 0.008 * (ncol(obj) / 1000)))
    doublet_result <- tryCatch({
      sce <- as.SingleCellExperiment(obj)
      scored <- scDblFinder::scDblFinder(sce, dbr = rate)
      classes <- as.character(scDblFinder::scoreDoublets(scored)$class)
      if (is.null(classes)) classes <- as.character(scored$scDblFinder.class)
      if (is.null(classes)) classes <- as.character(SummarizedExperiment::colData(scored)$scDblFinder.class)
      names(classes) <- colnames(obj)
      obj[["doublet_class"]] <- classes
      counts <- table(classes)
      print(counts)
      list(
        method = "scDblFinder",
        expected_rate = round(rate, 4),
        singlet = as.integer(if ("singlet" %in% names(counts)) counts[["singlet"]] else 0),
        doublet = as.integer(if ("doublet" %in% names(counts)) counts[["doublet"]] else 0)
      )
    }, error = function(error) {
      message("scDblFinder failed: ", conditionMessage(error))
      list(method = "scDblFinder", error = conditionMessage(error))
    })
  } else if (have("DoubletFinder")) {
    # DoubletFinder needs pK and nExp tuned per dataset (paramSweep/find.pK);
    # an untuned call is worse than no call, so it stays a documented fallback.
    doublet_result <- list(
      method = "DoubletFinder",
      error = "DoubletFinder needs a tuned pK and nExp per dataset; run it in a script instead"
    )
    message(doublet_result$error)
  } else {
    doublet_result <- list(
      method = "none",
      error = paste(
        "neither scDblFinder nor DoubletFinder is installed in this environment;",
        "install scDblFinder into a project environment to enable doublet detection"
      )
    )
    message("scDblFinder or SingleCellExperiment unavailable: doublets are not removed.")
  }
} else {
  message("detect_doublets is FALSE: doublet detection skipped.")
}