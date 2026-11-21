#!/usr/bin/env Rscript
# Shared helpers for the chromatin-workflow script toolkit.
#
# Source this file from every CLI script:
#
#     source(file.path(script_dir(), "_common.R"))
#
# Rules for this toolkit:
#   * stdout carries exactly one JSON object; every progress line goes to stderr.
#   * a missing optional package is reported in the JSON result, never fatal.
#   * arguments arrive as `--flag value` pairs built from the declared schema.
#   * pure downstream: inputs are called peaks / counts / BAM / bigWig, never FASTQ.

# --- location ---------------------------------------------------------------

# The directory holding the running script, so `source()` works when the engine
# invokes the script by absolute path from the project directory.
script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file <- sub("^--file=", "", args[startsWith(args, "--file=")])
  if (length(file) == 0) return(getwd())
  normalizePath(dirname(file[[1]]), mustWork = FALSE)
}

# --- output -----------------------------------------------------------------

info <- function(...) message(sprintf("[chip] %s", paste0(..., collapse = "")))

fail <- function(...) {
  message(sprintf("[chip] error: %s", paste0(..., collapse = "")))
  cat(jsonlite::toJSON(list(error = paste0(..., collapse = "")), auto_unbox = TRUE, null = "null"))
  cat("\n")
  quit(status = 1, save = "no")
}

emit <- function(result) {
  cat(jsonlite::toJSON(result, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE))
  cat("\n")
  invisible(result)
}

# A character vector safe to place in the JSON result: accessors return NULL
# rather than character(0) when a slot is empty, and `I(NULL)` is an error.
chr <- function(values) {
  if (is.null(values)) return(I(character(0)))
  I(as.character(unname(values)))
}

# --- argument parsing -------------------------------------------------------

# `--flag value` pairs; a repeated flag collects into a character vector; a flag
# with no following value is a boolean TRUE. Underscored schema names match the
# dashes the engine emits, so both spellings work.
parse_flags <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  index <- 1L
  while (index <= length(argv)) {
    token <- argv[[index]]
    if (!startsWith(token, "--")) {
      index <- index + 1L
      next
    }
    key <- gsub("-", "_", sub("^--", "", token))
    has_value <- index + 1L <= length(argv) && !startsWith(argv[[index + 1L]], "--")
    if (!has_value) {
      out[[key]] <- TRUE
      index <- index + 1L
    } else {
      value <- argv[[index + 1L]]
      if (is.null(out[[key]])) out[[key]] <- value else out[[key]] <- c(out[[key]], value)
      index <- index + 2L
    }
  }
  out
}

flag_str <- function(flags, name, default = NULL) {
  value <- flags[[name]]
  if (is.null(value) || identical(value, TRUE)) return(default)
  as.character(value[[1]])
}

flag_req <- function(flags, name) {
  value <- flag_str(flags, name)
  if (is.null(value) || !nzchar(value)) fail("missing required argument --", gsub("_", "-", name))
  value
}

flag_int <- function(flags, name, default = NULL) {
  value <- flag_str(flags, name)
  if (is.null(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed)) fail("--", gsub("_", "-", name), " must be an integer, got ", value)
  parsed
}

flag_num <- function(flags, name, default = NULL) {
  value <- flag_str(flags, name)
  if (is.null(value)) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed)) fail("--", gsub("_", "-", name), " must be a number, got ", value)
  parsed
}

flag_bool <- function(flags, name, default = FALSE) {
  value <- flags[[name]]
  if (is.null(value)) return(default)
  if (identical(value, TRUE)) return(TRUE)
  tolower(as.character(value[[1]])) %in% c("true", "1", "yes", "on")
}

flag_vec <- function(flags, name, default = character(0)) {
  value <- flags[[name]]
  if (is.null(value) || identical(value, TRUE)) return(default)
  as.character(value)
}

flag_nums <- function(flags, name, default = numeric(0)) {
  value <- flag_vec(flags, name)
  if (length(value) == 0) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (anyNA(parsed)) fail("--", gsub("_", "-", name), " takes numbers only")
  parsed
}

# --- packages ---------------------------------------------------------------

# Report a package the tool can work without, instead of dying.
optional_package <- function(name) {
  if (requireNamespace(name, quietly = TRUE)) return(TRUE)
  info("optional package not installed: ", name)
  FALSE
}

# A required Bioconductor package; the hint names the private environment so a
# gap is actionable (project environment via env_request), never an install here.
require_bioc <- function(name) {
  if (!requireNamespace(name, quietly = TRUE)) {
    fail(name, " is not installed in this environment; use the plugin:chip environment")
  }
  suppressPackageStartupMessages(library(name, character.only = TRUE))
  invisible(TRUE)
}

require_diffbind <- function() require_bioc("DiffBind")
require_chipseeker <- function() require_bioc("ChIPseeker")
require_csaw <- function() require_bioc("csaw")

# --- objects ----------------------------------------------------------------

# Read an RDS saved by a previous step (a DiffBind DBA, a GRanges, or a list).
# `expect` optionally names a class to check; NULL accepts any object.
read_object <- function(path, expect = NULL) {
  if (!file.exists(path)) fail("input not found: ", path)
  info("reading ", path)
  object <- tryCatch(
    readRDS(path),
    error = function(error) fail("cannot read ", path, ": ", conditionMessage(error))
  )
  if (!is.null(expect) && !inherits(object, expect)) {
    fail("expected ", expect, " saved with saveRDS, got ", paste(class(object), collapse = "/"))
  }
  object
}

write_object <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  info("wrote ", path)
  invisible(path)
}

# A directory the engine already validated; create it and return the path.
ensure_dir <- function(path, fallback) {
  target <- if (is.null(path) || !nzchar(path)) fallback else path
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  target
}

# --- peaks / counts ---------------------------------------------------------

# Import called peaks as a GRanges. narrowPeak/broadPeak carry signal columns;
# a plain BED is three columns. rtracklayer picks the format from the extension.
import_peaks <- function(path) {
  if (!file.exists(path)) fail("peak file not found: ", path)
  lower <- tolower(path)
  format <- if (grepl("[.]narrowpeak$", lower)) "narrowPeak"
    else if (grepl("[.]broadpeak$", lower)) "broadPeak"
    else if (grepl("[.]bed$", lower)) "BED"
    else if (grepl("[.]gff[23]?$", lower)) "GFF"
    else "BED"
  info("importing peaks ", path, " as ", format)
  peaks <- tryCatch(
    rtracklayer::import(path, format = format),
    error = function(error) fail("cannot import peaks ", path, ": ", conditionMessage(error))
  )
  if (length(peaks) == 0) info("warning: no peaks imported from ", path)
  peaks
}

# Read a peak x sample count table (first column peak id, header row of samples).
read_count_table <- function(path) {
  if (!file.exists(path)) fail("count table not found: ", path)
  separator <- if (grepl("[.]csv$", tolower(path))) "," else "\t"
  frame <- utils::read.delim(path, sep = separator, row.names = 1, check.names = FALSE)
  counts <- as.matrix(frame)
  storage.mode(counts) <- "numeric"
  info("count table ", nrow(counts), " peaks x ", ncol(counts), " samples")
  counts
}

# --- figures ----------------------------------------------------------------

save_plot <- function(plot, directory, filename, width = 8, height = 5) {
  if (is.null(directory)) return(NULL)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, filename)
  device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  if (identical(device, "png")) {
    grDevices::png(path, width = width, height = height, units = "in", res = 150)
  } else {
    device(path, width = width, height = height, units = "in", res = 150)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  info("figure ", path)
  path
}

# Levels sorted numerically when they are numbers, otherwise as text.
sorted_levels <- function(values) {
  unique_values <- unique(as.character(values))
  numbers <- suppressWarnings(as.numeric(unique_values))
  if (!anyNA(numbers) && length(numbers) > 0) {
    unique_values[order(numbers)]
  } else {
    sort(unique_values)
  }
}