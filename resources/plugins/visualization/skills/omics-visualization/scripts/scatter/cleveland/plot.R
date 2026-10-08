#!/usr/bin/env Rscript

# Template-ID: scatter-cleveland
#
# Purpose:
#   Draw a Cleveland dot plot of a numeric value across categories, with
#   a stem from zero to each point.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - Tissue: category label
#     - Count: numeric value
#
# Output:
#   A PDF, PNG, or SVG Cleveland dot plot.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category with a comparable numeric value.
#   Category order is a display choice, not a statistical transform.
#   Values are treated as already summarized, not as raw replicates.

local({
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (!length(file_arg)) {
        stop("Run this file with Rscript.", call. = FALSE)
    }
    dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
    for (i in seq_len(8)) {
        candidate <- file.path(dir, "lib", "common.R")
        if (file.exists(candidate)) {
            source(normalizePath(candidate), chdir = FALSE)
            return(invisible())
        }
        parent <- dirname(dir)
        if (identical(parent, dir)) break
        dir <- parent
    }
    stop("Cannot find scripts/lib/common.R", call. = FALSE)
})

io <- parse_io_args()

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        x = "Count",
        y = "Tissue"
    ),
    labels = list(
        title = "Cleveland Dot Plot",
        x = "Count",
        y = ""
    )
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y

df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
df[[x_col]] <- as.numeric(df[[x_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]]
)) +
    geom_segment(
        aes(x = 0, xend = .data[[config$columns$x]], yend = .data[[config$columns$y]]),
        colour = "#caccd1"
    ) +
    geom_point(size = 5, colour = "#1aafd0", fill = "#ffffffff", shape = 21) +
    geom_point(size = 3, colour = "#1aafd0") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
