# Cell-level tests treat every cell as independent, which inflates significance
# because cells from one donor are not independent (pseudoreplication). A
# sample-level comparison sums RAW counts per sample first and tests between
# samples: the experimental unit is the donor, not the cell.
#
# sample_column is set in the config cell; group_column (optional) splits a
# sample across conditions so one sample never becomes one summed column.
if (!sample_column %in% colnames(obj[[]])) {
  message("No '", sample_column, "' column in this object; skipping pseudobulk.")
  message("Add one with obj$", sample_column, " <- ... then re-run this cell.")
  pb_info <- NULL
} else {
  group_provided <- !is.null(group_column) && nzchar(group_column)
  if (group_provided && !group_column %in% colnames(obj[[]])) {
    stop("group column not found: ", group_column, call. = FALSE)
  }
  counts <- LayerData(obj, layer = "counts")
  metadata <- obj[[]]
  sample_values <- as.character(metadata[[sample_column]])
  group_values <- if (group_provided) as.character(metadata[[group_column]]) else rep("all", nrow(metadata))

  unit <- if (group_provided) paste(sample_values, group_values, sep = "|") else sample_values
  names(unit) <- rownames(metadata)

  per_unit_cells <- table(unit)
  small_units <- names(per_unit_cells)[per_unit_cells < min_cells]

  columns <- list()
  design <- list()
  for (name in unique(unit)) {
    cells <- intersect(names(unit)[unit == name], colnames(counts))
    if (length(cells) < min_cells) {
      message("skipping ", name, ": ", length(cells), " cells < ", min_cells)
      next
    }
    summed <- Matrix::rowSums(counts[, cells, drop = FALSE])
    parts <- strsplit(name, "|", fixed = TRUE)[[1]]
    identifier <- gsub("[^A-Za-z0-9._-]", "_", name)
    columns[[identifier]] <- as.numeric(summed)
    design[[length(design) + 1]] <- list(
      column = identifier,
      sample = parts[[1]],
      group = if (group_provided) parts[[2]] else "all",
      cells = length(cells)
    )
  }
  if (length(columns) == 0) {
    stop("no unit has at least ", min_cells, " cells; lower min_cells or check sample_column (",
         sample_column, ")", call. = FALSE)
  }

  matrix_out <- do.call(cbind, columns)
  rownames(matrix_out) <- rownames(counts)
  colnames(matrix_out) <- names(columns)

  # Drop genes that are essentially unexpressed: they add noise and break
  # dispersion estimation. The filter counts UNITS (columns), not cells.
  expressed <- rowSums(matrix_out >= min_count) >= min_cells_expressing
  dropped <- sum(!expressed)
  matrix_out <- matrix_out[expressed, , drop = FALSE]
  if (nrow(matrix_out) == 0) {
    stop("no gene passed the expression filter (min_count ", min_count,
         ", min_cells_expressing ", min_cells_expressing, ")", call. = FALSE)
  }

  # The sample sheet the DE wrappers read: `sample` equals the counts column
  # names; source_sample and the condition column stay separate for contrasts.
  sample_sheet <- data.frame(sample = vapply(design, function(e) e$column, character(1)),
                             stringsAsFactors = FALSE)
  sample_sheet[["source_sample"]] <- vapply(design, function(e) e$sample, character(1))
  condition_column <- if (group_provided) group_column else "group"
  sample_sheet[[condition_column]] <- vapply(design, function(e) e$group, character(1))
  sample_sheet[["cells"]] <- vapply(design, function(e) e$cells, integer(1))
  rownames(sample_sheet) <- NULL

  written <- character(0)
  if (!is.null(counts_output) && nzchar(counts_output)) {
    dir.create(dirname(counts_output), recursive = TRUE, showWarnings = FALSE)
    write.table(data.frame(gene_id = rownames(matrix_out), matrix_out, check.names = FALSE),
                counts_output, sep = "\t", quote = FALSE, row.names = FALSE)
    written <- c(written, counts_output)
  }
  if (!is.null(samples_output) && nzchar(samples_output)) {
    dir.create(dirname(samples_output), recursive = TRUE, showWarnings = FALSE)
    write.table(sample_sheet, samples_output, sep = "\t", quote = FALSE, row.names = FALSE)
    written <- c(written, samples_output)
  }

  pb_notes <- character(0)
  if (group_provided) {
    group_counts <- table(sample_sheet[[condition_column]])
    thin <- names(group_counts)[group_counts < min_samples_per_group]
    if (length(thin) > 0) {
      pb_notes <- c(pb_notes, sprintf(
        "group(s) with fewer than %d samples: %s; a two-group test needs biological replicates",
        min_samples_per_group, paste(thin, collapse = ", ")))
    }
    if (length(group_counts) < 2) {
      pb_notes <- c(pb_notes, "only one group is present, so no comparison can be made")
    }
  } else {
    pb_notes <- c(pb_notes, "no group_column: columns are one per sample and carry no condition")
  }
  if (length(small_units) > 0) {
    pb_notes <- c(pb_notes, sprintf("skipped %d unit(s) under %d cells", length(small_units), min_cells))
  }
  if (dropped > 0) {
    pb_notes <- c(pb_notes, sprintf("dropped %d gene(s) below the expression filter", dropped))
  }

  pb_info <- list(
    matrix = matrix_out, design = sample_sheet, written = written,
    dropped = dropped, notes = pb_notes, condition_column = condition_column,
    group_provided = group_provided
  )
  cat(sprintf("pseudobulk: %d genes x %d units\n", nrow(matrix_out), ncol(matrix_out)))
  cat("Next: run the deseq2 wrapper with counts=", counts_output,
      " and samplesheet=", samples_output, "\n", sep = "")
}