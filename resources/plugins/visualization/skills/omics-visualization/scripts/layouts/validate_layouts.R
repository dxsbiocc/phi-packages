#!/usr/bin/env Rscript

# Validate the layout catalog and optionally run every canonical layout.
# Usage: Rscript validate_layouts.R [--smoke]

suppressPackageStartupMessages({
    library(jsonlite)
    library(yaml)
})

`%||%` <- function(value, fallback) {
    if (is.null(value) || !length(value)) fallback else value
}

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this file with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
skill_root <- normalizePath(file.path(dirname(script_path), "../.."))

args <- commandArgs(trailingOnly = TRUE)
unknown_args <- setdiff(args, "--smoke")
if (length(unknown_args)) {
    stop("Unknown argument(s): ", paste(unknown_args, collapse = ", "), call. = FALSE)
}
run_smoke <- "--smoke" %in% args

catalog_path <- file.path(skill_root, "references", "layouts.yaml")
catalog <- yaml::read_yaml(catalog_path)
layouts <- catalog$layouts
if (!is.list(layouts) || !length(layouts)) {
    stop("Layout catalog contains no layouts.", call. = FALSE)
}

required <- c(
    "id", "title", "aliases", "source", "preview", "panel_count",
    "use_when", "avoid_when"
)
errors <- character()
ids <- vapply(layouts, function(entry) {
    if (is.null(entry$id) || length(entry$id) != 1L) "" else as.character(entry$id)
}, character(1))
if (any(!nzchar(ids))) {
    errors <- c(errors, "Every layout must have one non-empty id.")
}
duplicate_ids <- unique(ids[nzchar(ids) & duplicated(ids)])
if (length(duplicate_ids)) {
    errors <- c(errors, paste("Duplicate layout ids:", paste(duplicate_ids, collapse = ", ")))
}

for (entry in layouts) {
    missing <- setdiff(required, names(entry))
    if (length(missing)) {
        errors <- c(errors, sprintf(
            "%s missing field(s): %s",
            entry$id %||% "<unknown>",
            paste(missing, collapse = ", ")
        ))
        next
    }

    source_path <- file.path(skill_root, entry$source)
    preview_path <- file.path(skill_root, entry$preview)
    if (!file.exists(source_path)) {
        errors <- c(errors, paste(entry$id, "missing source", entry$source))
    }
    if (!file.exists(preview_path)) {
        errors <- c(errors, paste(entry$id, "missing preview", entry$preview))
    }

    if (run_smoke && file.exists(source_path)) {
        smoke_dir <- tempfile(paste0("layout-", entry$id, "-"))
        dir.create(smoke_dir, recursive = TRUE)
        output_path <- file.path(smoke_dir, "preview.png")
        result <- system2(
            "Rscript",
            c(shQuote(source_path), shQuote(output_path)),
            stdout = TRUE,
            stderr = TRUE
        )
        status <- attr(result, "status") %||% 0L
        audit_path <- file.path(smoke_dir, "preview.layout-audit.json")
        if (status != 0L || !file.exists(output_path) || !file.exists(audit_path)) {
            errors <- c(errors, paste(
                entry$id,
                "smoke test failed:",
                paste(result, collapse = " | ")
            ))
        } else {
            audit <- jsonlite::read_json(audit_path, simplifyVector = TRUE)
            if (!identical(audit$status, "PASS")) {
                errors <- c(errors, paste(entry$id, "audit status", audit$status))
            }
        }
    }
}

if (length(errors)) {
    stop(paste(errors, collapse = "\n"), call. = FALSE)
}

cat(sprintf(
    "Validated %d layout(s)%s.\n",
    length(layouts),
    if (run_smoke) " with smoke tests" else ""
))
