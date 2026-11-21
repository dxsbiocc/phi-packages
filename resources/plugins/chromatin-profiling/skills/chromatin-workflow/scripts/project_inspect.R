#!/usr/bin/env Rscript
# Inspect a chromatin-profiling input and report what it contains: called peaks,
# a DiffBind sample sheet, a count matrix, aligned BAMs, a bigWig, or a saved
# object. Pure downstream — it never expects FASTQ. With --output and a sample
# sheet it builds and saves a DiffBind DBA object.
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
source(file.path(dirname(if (length(here)) here[[1]] else "."), "_common.R"))

flags <- parse_flags()
input <- flag_req(flags, "input")
output <- flag_str(flags, "output")

detect_kind <- function(path) {
  if (dir.exists(path)) {
    entries <- list.files(path)
    if (any(grepl("[.](bed|narrowPeak|broadPeak)$", entries))) return("peaks-dir")
    if (any(grepl("[.]bam$", entries))) return("bam-dir")
    if (any(grepl("[.](bw|bigWig)$", entries))) return("bigwig-dir")
    return("unknown-directory")
  }
  lower <- tolower(path)
  if (endsWith(lower, ".rds")) return("saved-object")
  if (endsWith(lower, ".narrowpeak")) return("narrowpeak")
  if (endsWith(lower, ".broadpeak")) return("broadpeak")
  if (endsWith(lower, ".bed")) return("bed")
  if (endsWith(lower, ".bam")) return("bam-file")
  if (endsWith(lower, ".bw") || endsWith(lower, ".bigwig")) return("bigwig")
  if (endsWith(lower, ".gtf") || endsWith(lower, ".gff") || endsWith(lower, ".gff3")) return("annotation")
  if (endsWith(lower, ".csv") || endsWith(lower, ".tsv") || endsWith(lower, ".txt")) {
    header <- tolower(readLines(path, n = 1))
    if (grepl("sampleid", header) && grepl("peaks", header)) return("sample-sheet")
    return("count-table")
  }
  "unknown-file"
}

kind <- detect_kind(input)
info("input kind: ", kind)
require_bioc("GenomicRanges")
require_bioc("rtracklayer")

shape <- list()
samples <- character(0)
candidate <- character(0)
written <- NULL

if (kind %in% c("peaks-dir", "narrowpeak", "broadpeak", "bed")) {
  files <- if (kind == "peaks-dir")
    sort(list.files(input, pattern = "[.](bed|narrowPeak|broadPeak)$", full.names = TRUE))
    else input
  per_sample <- lapply(files, function(path) {
    peaks <- import_peaks(path)
    list(sample = tools::file_path_sans_ext(basename(path)), peaks = length(peaks),
         median_width = as.integer(median(width(peaks))))
  })
  samples <- vapply(per_sample, function(entry) entry$sample, character(1))
  shape <- list(n_samples = length(files), per_sample = per_sample)
  candidate <- c("qc", "consensus", "count", "annotate")
} else if (kind == "sample-sheet") {
  sheet <- utils::read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  samples <- if ("SampleID" %in% colnames(sheet)) as.character(sheet$SampleID) else character(0)
  conditions <- if ("Condition" %in% colnames(sheet)) unique(as.character(sheet$Condition)) else character(0)
  shape <- list(n_samples = nrow(sheet), columns = chr(colnames(sheet)),
                conditions = chr(conditions),
                has_bam = "bamReads" %in% colnames(sheet),
                has_control = "bamControl" %in% colnames(sheet))
  candidate <- c("qc", "consensus", "count", "diff", "annotate")
  if (!is.null(output)) {
    require_diffbind()
    dba <- DiffBind::dba(sampleSheet = sheet)
    write_object(dba, output)
    written <- output
    shape$n_peaks <- nrow(dba$binding)
  }
} else if (kind == "count-table") {
  counts <- read_count_table(input)
  samples <- colnames(counts)
  shape <- list(n_peaks = nrow(counts), n_samples = ncol(counts),
                library_sizes = as.list(colSums(counts)))
  candidate <- c("normalize", "pca", "diff")
} else if (kind %in% c("bam-dir", "bam-file")) {
  files <- if (kind == "bam-dir") sort(list.files(input, pattern = "[.]bam$", full.names = TRUE)) else input
  samples <- tools::file_path_sans_ext(basename(files))
  indexed <- vapply(files, function(path)
    file.exists(paste0(path, ".bai")) || file.exists(sub("[.]bam$", ".bai", path)), logical(1))
  shape <- list(n_bam = length(files), indexed = as.list(indexed))
  candidate <- c("qc", "count", "csaw")
} else if (kind %in% c("bigwig", "bigwig-dir")) {
  shape <- list(note = "bigWig coverage track; used for profiles and visualization")
  candidate <- c("annotate", "coverage-profile")
} else if (kind == "saved-object") {
  object <- read_object(input)
  shape <- list(class = chr(class(object)))
  if (inherits(object, "DBA")) candidate <- c("diff", "annotate")
  else if (inherits(object, "GRanges")) candidate <- c("annotate", "count")
} else if (kind == "annotation") {
  shape <- list(note = "GTF/GFF3 for building a TxDb (annotation input, not peaks)")
  candidate <- c("annotate")
} else {
  fail("unsupported --input; expected peaks (BED/narrowPeak/broadPeak), a DiffBind ",
       "sample sheet CSV, a count table, BAM, bigWig, or a saved .rds")
}

emit(list(
  input = input,
  kind = kind,
  output = written,
  samples = chr(samples),
  shape = shape,
  candidate_steps = chr(candidate),
  notes = chr(c(
    "pure downstream: alignment and peak calling belong to the nf-core/chipseq wrapper",
    if (kind == "sample-sheet" && !isTRUE(shape$has_bam))
      "no bamReads column: counting and FRiP need BAM paths",
    if (kind %in% c("bam-dir", "bam-file") && any(!unlist(shape$indexed)))
      "some BAMs lack a .bai index; counting and FRiP need indexed BAM"
  ))
))