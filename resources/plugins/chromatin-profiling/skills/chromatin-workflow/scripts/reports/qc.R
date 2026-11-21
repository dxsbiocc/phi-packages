# JSON contract for chip_qc. Reads the bookkeeping the fragments left behind:
# per-sample peak counts, FRiP, replicate reproducibility and blacklist removal.
peaks_summary <- if (!is.null(peak_list)) {
  lapply(names(peak_list), function(name) {
    peaks <- peak_list[[name]]
    list(sample = name, n_peaks = length(peaks),
         median_width = as.integer(stats::median(width(peaks))))
  })
} else NULL

warnings <- character(0)
if (is.null(peak_list) || length(peak_list) < 2) {
  warnings <- c(warnings, "fewer than two peak sets; replicate reproducibility not computed")
}
if (!is.null(qc_metrics$frip)) {
  unavailable <- vapply(qc_metrics$frip, function(e) is.na(e$frip), logical(1))
  if (any(unavailable)) warnings <- c(warnings, "some samples had no matching peak set for FRiP")
}
warnings <- c(warnings, paste(
  "IDR is the formal reproducibility test; the idr pip package is a gap",
  "(project environment via env_request). Overlap fraction is a proxy."))

emit(list(
  input = tool_input,
  output = tool_output,
  samples = chr(names(peak_list)),
  peaks = peaks_summary,
  frip = qc_metrics$frip,
  reproducibility = qc_metrics$reproducibility,
  blacklist_removed = blacklist_removed,
  figures = chr(unlist(figure_paths, use.names = FALSE)),
  steps = chr(analysis_steps),
  warnings = chr(warnings)
))