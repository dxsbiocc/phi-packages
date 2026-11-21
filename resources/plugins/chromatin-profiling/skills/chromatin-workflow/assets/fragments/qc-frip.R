# FRiP (fraction of reads in peaks): the headline ChIP/CUT&Tag quality metric.
# A low FRiP means the antibody did not enrich — most reads are background.
# ENCODE guidance: TF ChIP ideally > 1%, broad marks lower; CUT&Tag is often much
# higher because background is minimal. This is computed from BAM + peaks, so it
# is independent of any peak caller's own report.
if (is.null(peak_list)) stop("load peaks before computing FRiP (load-peaks)", call. = FALSE)
if (is.null(bam_paths) || length(bam_paths) == 0) stop("load BAM before FRiP (load-bam)", call. = FALSE)

single_end <- !isTRUE(paired_end)
frip <- lapply(seq_along(bam_paths), function(i) {
  bam <- bam_paths[[i]]
  name <- tools::file_path_sans_ext(basename(bam))
  peaks <- peak_list[[name]]
  # Fall back to positional matching when BAM and peak names differ.
  if (is.null(peaks) && i <= length(peak_list)) peaks <- peak_list[[i]]
  if (is.null(peaks)) {
    return(list(sample = name, total_reads = NA_integer_, reads_in_peaks = NA_integer_,
                frip = NA_real_, note = "no matching peak set"))
  }
  total <- Rsamtools::countBam(bam, param = Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)))$records
  overlaps <- GenomicAlignments::summarizeOverlaps(
    peaks, Rsamtools::BamFile(bam), mode = "Union",
    singleEnd = single_end, ignore.strand = TRUE)
  in_peaks <- sum(SummarizedExperiment::assay(overlaps))
  list(sample = name, total_reads = as.integer(total), reads_in_peaks = as.integer(in_peaks),
       frip = round(in_peaks / max(1, total), 4))
})

qc_metrics <- list(frip = frip)
for (entry in frip) {
  if (is.na(entry$frip)) {
    cat(sprintf("%s: FRiP unavailable (%s)\n", entry$sample, entry$note))
  } else {
    cat(sprintf("%s: FRiP %.2f%% (%s reads in peaks / %s mapped)\n",
                entry$sample, 100 * entry$frip,
                format(entry$reads_in_peaks, big.mark = ","),
                format(entry$total_reads, big.mark = ",")))
  }
}
analysis_steps <- c(analysis_steps, "qc-frip")