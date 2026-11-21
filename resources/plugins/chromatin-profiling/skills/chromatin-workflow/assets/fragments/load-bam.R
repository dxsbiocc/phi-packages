# Aligned reads (BAM) for FRiP, read counting and coverage profiles. Alignment
# already happened upstream (bowtie2/BWA in the nf-core/chipseq wrapper); this
# step only opens the BAMs. Counting and FRiP need a .bai index beside each BAM.
paths <- bam_paths
if (is.null(paths) || length(paths) == 0) {
  if (!is.null(sample_sheet) && "bamReads" %in% colnames(sample_sheet)) {
    paths <- sample_sheet$bamReads
  } else {
    directory <- if (!is.null(bam_dir) && nzchar(bam_dir)) bam_dir else data_path
    if (dir.exists(directory)) {
      paths <- sort(list.files(directory, pattern = "[.]bam$", full.names = TRUE))
    } else if (grepl("[.]bam$", data_path)) {
      paths <- data_path
    }
  }
}
if (!is.null(sample_sheet_path) && nzchar(sample_sheet_path)) {
  sheet_dir <- dirname(sample_sheet_path)
  paths <- vapply(paths, function(p) {
    if (file.exists(p)) p else file.path(sheet_dir, p)
  }, character(1))
}
if (length(paths) == 0) {
  stop("No BAM files found. Pass a directory of .bam files as bam_dir, list them ",
       "in bam_paths, or use a sample sheet with a bamReads column.", call. = FALSE)
}

has_index <- function(path) {
  file.exists(paste0(path, ".bai")) || file.exists(sub("[.]bam$", ".bai", path))
}
missing_index <- paths[!vapply(paths, has_index, logical(1))]

bam_files <- Rsamtools::BamFileList(paths)
bam_paths <- paths
for (path in paths) cat("BAM:", basename(path), "\n")
if (length(missing_index) > 0) {
  cat("warning: no .bai index for ", paste(basename(missing_index), collapse = ", "),
      "; counting and FRiP need indexed BAM (samtools index upstream)\n", sep = "")
}
analysis_steps <- c(analysis_steps, "load-bam")