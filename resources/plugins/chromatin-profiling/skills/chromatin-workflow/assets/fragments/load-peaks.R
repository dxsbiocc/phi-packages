# Called peaks per sample, imported as a GRangesList. This is the pure-downstream
# entry point: peak calling (MACS2 / SEACR) already happened upstream — in the
# nf-core/chipseq wrapper for ChIP, or in the user's SEACR run for CUT&Tag.
#
# data_path may be a single peak file or a directory of peak files; peak_paths
# may also list them explicitly. narrowPeak/broadPeak carry signal columns that
# rtracklayer keeps as metadata; a plain BED is three columns.
paths <- peak_paths
if (is.null(paths) || length(paths) == 0) {
  if (!is.null(sample_sheet) && "Peaks" %in% colnames(sample_sheet)) {
    # DiffBind route: the sample sheet names each sample's peak file.
    paths <- sample_sheet$Peaks
  } else if (dir.exists(data_path)) {
    paths <- sort(list.files(data_path, pattern = "[.](bed|narrowPeak|broadPeak)$",
                             full.names = TRUE))
  } else {
    paths <- data_path
  }
}
# Resolve relative peak paths against the sample sheet's directory.
if (!is.null(sample_sheet_path) && nzchar(sample_sheet_path)) {
  sheet_dir <- dirname(sample_sheet_path)
  paths <- vapply(paths, function(p) {
    if (file.exists(p)) p else file.path(sheet_dir, p)
  }, character(1))
}
if (length(paths) == 0) {
  stop("No peak files found. Pass a directory of BED/narrowPeak/broadPeak files ",
       "as data_path, list them in peak_paths, or use a sample sheet with a Peaks column.",
       call. = FALSE)
}

peak_list <- GRangesList(lapply(paths, function(path) {
  if (!file.exists(path)) stop("peak file not found: ", path, call. = FALSE)
  rtracklayer::import(path)
}))
names(peak_list) <- tools::file_path_sans_ext(basename(paths))

for (name in names(peak_list)) {
  peaks <- peak_list[[name]]
  widths <- width(peaks)
  cat(sprintf("%s: %d peaks, median width %d bp\n",
              name, length(peaks), as.integer(median(widths))))
}
analysis_steps <- c(analysis_steps, "load-peaks")