#!/usr/bin/env Rscript
# Shared helpers for the seurat-workflow script toolkit.
#
# Source this file from every CLI script:
#
#     source(file.path(dirname(sys.frame(1)$ofile), "_common.R"))
#
# or, when run through `skill_run`, which executes the absolute script path:
#
#     source(file.path(script_dir(), "_common.R"))
#
# Rules for this toolkit:
#   * stdout carries exactly one JSON object; every progress line goes to stderr.
#   * a missing optional package is reported in the JSON result, never fatal.
#   * arguments arrive as `--flag value` pairs built from the declared schema.

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

info <- function(...) message(sprintf("[seurat] %s", paste0(..., collapse = "")))

fail <- function(...) {
  message(sprintf("[seurat] error: %s", paste0(..., collapse = "")))
  cat(jsonlite::toJSON(list(error = paste0(..., collapse = "")), auto_unbox = TRUE, null = "null"))
  cat("\n")
  quit(status = 1, save = "no")
}

.seurat_tool_map <- c(
  "project_inspect.R" = "sc_inspect",
  "qc_filter.R" = "sc_qc",
  "preprocess_cluster.R" = "sc_preprocess",
  "find_markers.R" = "sc_markers",
  "annotate_clusters.R" = "sc_annotate",
  "pseudobulk.R" = "sc_pseudobulk",
  "env_probe.R" = "sc_env_probe"
)

.running_tool_name <- function(tool_map, prefix) {
  args <- commandArgs(trailingOnly = FALSE)
  file <- sub("^--file=", "", args[startsWith(args, "--file=")])
  if (length(file) == 0) return(paste0(prefix, "_tool"))
  base <- basename(file[[1]])
  if (base %in% names(tool_map)) unname(tool_map[[base]]) else paste0(prefix, "_", tools::file_path_sans_ext(base))
}

.png_dims <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  hdr <- readBin(con, what = "raw", n = 24L)
  if (length(hdr) < 24L || !identical(as.integer(hdr[1:8]), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))) return(NULL)
  b <- as.numeric(hdr[17:24])
  w <- ((b[1] * 256 + b[2]) * 256 + b[3]) * 256 + b[4]
  h <- ((b[5] * 256 + b[6]) * 256 + b[7]) * 256 + b[8]
  list(widthPx = as.integer(w), heightPx = as.integer(h))
}

write_artifact_sidecar <- function(path, tool, skill = "seurat-workflow") {
  if (is.null(path) || !nzchar(path) || !file.exists(path) || dir.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("png", "pdf", "svg")) {
    kind <- "figure"
    media <- switch(ext, png = "image/png", pdf = "application/pdf", svg = "image/svg+xml")
  } else if (ext %in% c("tsv", "csv")) {
    kind <- "table"
    media <- if (ext == "csv") "text/csv" else "text/tab-separated-values"
  } else {
    return(NULL)
  }
  stem <- gsub("[-_]+", " ", tools::file_path_sans_ext(basename(path)))
  title <- substr(if (nzchar(trimws(stem))) paste0(tool, ": ", trimws(stem)) else tool, 1L, 200L)
  desc <- list(
    contractVersion = "1.0.0",
    kind = kind,
    file = basename(path),
    mediaType = media,
    title = title,
    provenance = list(
      createdAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      tool = tool,
      skill = skill
    )
  )
  if (kind == "figure") {
    fig_block <- list(format = ext)
    if (ext == "png") {
      dims <- tryCatch(.png_dims(path), error = function(e) NULL)
      if (!is.null(dims)) {
        fig_block$widthPx <- dims$widthPx
        fig_block$heightPx <- dims$heightPx
      }
    }
    desc$figure <- fig_block
  } else if (kind == "table") {
    sep <- if (ext == "csv") "," else "\t"
    lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
    if (length(lines) > 0) {
      n_cols <- length(strsplit(lines[[1]], sep, fixed = TRUE)[[1]])
      desc$table <- list(rows = as.integer(max(0L, length(lines) - 1L)), columns = as.integer(n_cols))
    }
  }
  sidecar <- paste0(path, ".phi-artifact.json")
  tmp <- paste0(sidecar, ".tmp")
  writeLines(jsonlite::toJSON(desc, auto_unbox = TRUE, null = "null", pretty = TRUE), tmp, useBytes = TRUE)
  file.rename(tmp, sidecar)
  path
}

attach_artifacts <- function(result, tool_map = .seurat_tool_map, prefix = "sc", skill = "seurat-workflow") {
  if (!is.list(result) || !is.null(result$artifacts) || !is.null(result$error)) return(result)
  cands <- character(0)
  for (key in c("figures", "table", "tables", "outputs")) {
    val <- result[[key]]
    if (is.character(val)) cands <- c(cands, as.character(val))
  }
  cands <- unique(cands[nzchar(cands)])
  if (length(cands) == 0) return(result)
  tool <- .running_tool_name(tool_map, prefix)
  written <- character(0)
  for (p in cands) {
    art <- write_artifact_sidecar(p, tool = tool, skill = skill)
    if (!is.null(art)) written <- c(written, art)
  }
  if (length(written) > 0) {
    result$artifacts <- chr(head(written, 8L))
  }
  result
}

emit <- function(result) {
  result <- attach_artifacts(result)
  cat(jsonlite::toJSON(result, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE))
  cat("\n")
  invisible(result)
}

# A character vector safe to place in the JSON result: Seurat's accessors return
# NULL rather than character(0) when a slot is empty, and `I(NULL)` is an error.
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

require_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    fail("Seurat is not installed in this environment; use the plugin:seurat environment")
  }
  suppressPackageStartupMessages(library(Seurat))
  invisible(TRUE)
}

# --- objects ----------------------------------------------------------------

read_object <- function(path) {
  if (!file.exists(path)) fail("input not found: ", path)
  info("reading ", path)
  object <- tryCatch(
    readRDS(path),
    error = function(error) fail("cannot read ", path, ": ", conditionMessage(error))
  )
  if (!inherits(object, "Seurat")) {
    fail("expected a Seurat object saved with saveRDS, got ", paste(class(object), collapse = "/"))
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

# Multi-panel QC violins.
#
# Seurat's own `VlnPlot(..., ncol = 3)` combines panels with patchwork's `&`
# operator, which fails on the ggplot2 4.x / patchwork pair shipped in this environment
# ("Can't find method for generic `&`"). Building one plot per feature and
# assembling them with wrap_plots avoids that operator entirely. If patchwork is
# unavailable, fall back to Seurat's combined plot so a figure is still written.
violin_panels <- function(object, features, ncol = 3) {
  features <- features[features %in% colnames(object[[]])]
  if (length(features) == 0) return(NULL)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    panels <- lapply(features, function(feature) {
      Seurat::VlnPlot(object, features = feature, pt.size = 0) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = NULL, title = feature)
    })
    return(patchwork::wrap_plots(panels, ncol = ncol))
  }
  Seurat::VlnPlot(object, features = features, ncol = ncol, pt.size = 0)
}

# Group of a Seurat object as a data frame, in a stable column name.
metadata_of <- function(object) {
  frame <- object[[]]
  frame$cell <- rownames(frame)
  frame
}

counts_of <- function(object, assay = NULL) {
  SeuratObject::LayerData(assay_or_default(object, assay), layer = "counts")
}

assay_or_default <- function(object, assay = NULL) {
  if (is.null(assay) || !nzchar(assay)) return(object[[Seurat::DefaultAssay(object)]])
  if (!assay %in% names(object@assays)) {
    fail("assay not found: ", assay, "; available: ", paste(names(object@assays), collapse = ", "))
  }
  object[[assay]]
}

resolve_group <- function(object, group_by) {
  columns <- colnames(object[[]])
  if (!group_by %in% columns) {
    fail("group column not found: ", group_by, "; available: ", paste(columns, collapse = ", "))
  }
  invisible(group_by)
}

# Cluster labels sorted numerically when they are numbers, otherwise as text.
sorted_levels <- function(values) {
  unique_values <- unique(as.character(values))
  numbers <- suppressWarnings(as.numeric(unique_values))
  if (!anyNA(numbers) && length(numbers) > 0) {
    unique_values[order(numbers)]
  } else {
    sort(unique_values)
  }
}
