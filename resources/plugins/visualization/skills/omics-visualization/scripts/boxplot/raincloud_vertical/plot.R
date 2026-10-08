#!/usr/bin/env Rscript

# Template-ID: boxplot-raincloud-vertical
#
# Purpose:
#   Draw a vertical raincloud plot: half-boxplot, jittered points, and
#   a half-violin on the right of each category.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - gene: category label
#     - expression: numeric value
#
# Output:
#   A PDF, PNG, or SVG vertical raincloud plot.
#
# Dependencies:
#   ggplot2, readr, gghalves, ggprism
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
#   Each row is one observation of a numeric value in a category.
#   The half-violin is a kernel density estimate, not a statistical test.
#   Half-boxplot and jitter show the same observations as the density.

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
        x = "gene",
        y = "expression"
    ),
    labels = list(
        title = "",
        x = "",
        y = "Gene expression level"
    )
)

load_packages(c("ggplot2", "readr", "gghalves", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c(
    "#009392", "#39b185", "#9ccb86", "#e9e29c",
    "#eeb479", "#e88471", "#cf597e"
)
shape_values <- c(21, 22, 23, 24)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[x_col]]
)) +
    geom_half_boxplot(nudge = 0.2, errorbar.length = 0.3) +
    geom_jitter(
        aes(shape = .data[[x_col]]),
        width = 0.1,
        size = 2
    ) +
    geom_half_violin(side = "r", nudge = 0.2) +
    scale_fill_manual(values = fill_values) +
    scale_shape_manual(values = shape_values) +
    guides(
        x = guide_prism_bracket(),
        fill = "none",
        shape = "none"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
