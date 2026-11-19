# Thresholds come from the config cell. Read the summary() output and the violin
# panels above, then edit those values rather than accepting them blindly:
#   nFeature_RNA > min_genes    removes empty droplets
#   nFeature_RNA < max_genes    removes likely doublets / oversized cells
#   nCount_RNA   > min_counts   removes empty droplets
#   nCount_RNA   < max_counts   removes likely doublets
#   percent.mt   < mt_max       removes stressed cells
# A NULL bound skips that gate. Each gate is applied in sequence and reports how
# many cells IT removed, so a surprising total can be traced to one threshold.
obj_before <- obj
qc_filters <- list()
kept <- rep(TRUE, ncol(obj))
names(kept) <- colnames(obj)

apply_filter <- function(name, condition) {
  removed <- sum(!condition, na.rm = TRUE)
  qc_filters[[length(qc_filters) + 1]] <<- list(filter = name, removed = removed)
  kept[!is.na(condition) & !condition] <<- FALSE
  cat(sprintf("%-28s removed %d cells\n", name, removed))
  invisible(NULL)
}

if (!is.null(min_genes)) apply_filter(
  sprintf("nFeature_RNA >= %s", min_genes),
  obj[["nFeature_RNA", drop = TRUE]] >= min_genes
)
if (!is.null(max_genes)) apply_filter(
  sprintf("nFeature_RNA <= %s", max_genes),
  obj[["nFeature_RNA", drop = TRUE]] <= max_genes
)
if (!is.null(min_counts)) apply_filter(
  sprintf("nCount_RNA >= %s", min_counts),
  obj[["nCount_RNA", drop = TRUE]] >= min_counts
)
if (!is.null(max_counts)) apply_filter(
  sprintf("nCount_RNA <= %s", max_counts),
  obj[["nCount_RNA", drop = TRUE]] <= max_counts
)
if (!is.null(mt_max)) apply_filter(
  sprintf("percent.mt <= %s", mt_max),
  obj[["percent.mt", drop = TRUE]] <= mt_max
)
if (detect_doublets && "doublet_class" %in% colnames(obj[[]])) {
  apply_filter("doublet_class == singlet", obj[["doublet_class", drop = TRUE]] == "singlet")
}

obj <- obj[, kept, drop = FALSE]
cat(sprintf("kept %d of %d cells (%.1f%%)\n",
            ncol(obj), ncol(obj_before), 100 * ncol(obj) / max(1, ncol(obj_before))))

# Keeping under half the cells usually means a threshold is too strict for this
# dataset. Re-read the distributions before continuing.
if (ncol(obj) == 0) {
  warning("Every cell was filtered out; relax the thresholds.", call. = FALSE)
} else if (ncol(obj) < 0.5 * ncol(obj_before)) {
  warning("Fewer than half the cells survived QC; check the thresholds.", call. = FALSE)
}