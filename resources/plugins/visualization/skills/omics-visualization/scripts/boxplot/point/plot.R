#!/usr/bin/env Rscript

# Template-ID: boxplot-point
#
# Purpose:
#   Draw a notched boxplot of a numeric value across categories, with
#   jittered points coloured and shaped by category.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - gene: category label
#     - expression: numeric value
#
# Output:
#   A PDF, PNG, or SVG boxplot with points.
#
# Dependencies:
#   ggplot2, readr, ggprism
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
#   Boxes summarize the distribution; jittered points are the raw values.
#   Notches are a display of the median interval, not a formal test.

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

load_packages(c("ggplot2", "readr", "ggprism"))

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
colour_values <- c(
    "#d7689d", "#f3d337", "#a0c443", "#66b5ae", "#5da4dc", "#40769e"
)
shape_values <- c(16, 17, 15, 18)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    colour = .data[[x_col]]
)) +
    geom_boxplot(notch = TRUE) +
    geom_jitter(
        aes(shape = .data[[x_col]]),
        width = 0.2,
        size = 4,
        alpha = 0.8
    ) +
    scale_colour_manual(values = colour_values) +
    scale_shape_manual(values = shape_values) +
    guides(colour = "none", x = guide_prism_bracket()) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
