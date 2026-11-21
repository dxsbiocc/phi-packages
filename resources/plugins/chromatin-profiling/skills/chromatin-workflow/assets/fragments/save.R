# Persist the analysis so a later session or script resumes from here. The RDS
# holds the DBA object (or the consensus GRanges when there is no DBA); the tables
# beside it are the reviewable, language-neutral outputs: consensus BED, count
# matrix, differential results and peak annotation.
outdir <- if (!is.null(output_object) && nzchar(output_object)) dirname(output_object) else "."
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!is.null(output_object) && nzchar(output_object)) {
  object <- if (!is.null(dba_object)) dba_object
    else if (!is.null(consensus_peaks)) consensus_peaks
    else NULL
  if (!is.null(object)) {
    saveRDS(object, output_object)
    cat("saved", output_object, "\n")
  } else {
    cat("nothing to save as RDS (no DBA object or consensus peakset)\n")
  }
}

write_table <- function(frame, name) {
  if (is.null(frame)) return(invisible(NULL))
  path <- file.path(outdir, name)
  utils::write.csv(frame, path, row.names = TRUE)
  cat("wrote", path, "\n")
  invisible(path)
}

if (!is.null(consensus_peaks)) {
  bed <- file.path(outdir, "consensus.bed")
  rtracklayer::export(consensus_peaks, bed)
  cat("wrote", bed, "\n")
}
write_table(count_matrix, "counts.csv")
write_table(normalized_counts, "normalized_counts.csv")
write_table(occupancy_matrix, "occupancy.csv")
if (!is.null(diff_results) && !is.null(diff_results$table)) write_table(diff_results$table, "diff_binding.csv")
if (!is.null(csaw_result) && !is.null(csaw_result$table)) write_table(csaw_result$table, "csaw_regions.csv")
write_table(anno_df, "annotation.csv")
analysis_steps <- c(analysis_steps, "save")