#!/usr/bin/env Rscript

# Template-ID: scatter-matrix
#
# Purpose:
#   Draw a matrix (dot) plot of a numeric value across two categorical
#   axes, with point size and colour mapped to the value.
#
# Inputs:
#   A table with one row per cell-by-gene. Default example:
#     - cell: column category
#     - gene: row category
#     - count: numeric value
#
# Output:
#   A PDF, PNG, or SVG matrix scatter plot.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change axis category order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category-by-category pair with a comparable numeric value.
#   Values are treated as already summarized, not as raw replicates.
#   Size and colour encode the same value as a display choice.

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
        x = "cell",
        y = "gene",
        colour = "count",
        size = "count"
    ),
    labels = list(
        title = "",
        x = "",
        y = ""
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
colour_col <- config$columns$colour
size_col <- config$columns$size
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
df[[colour_col]] <- as.numeric(df[[colour_col]])
if (!identical(size_col, colour_col)) {
    df[[size_col]] <- as.numeric(df[[size_col]])
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    colour = .data[[colour_col]]
)) +
    geom_point(aes(size = .data[[size_col]])) +
    scale_size_continuous(range = c(2, 12), breaks = c(1000, 3000, 5000)) +
    scale_colour_gradient2(low = "#E8F2D1", high = "#437F79") +
    guides(size = guide_legend(override.aes = list(size = c(2, 4, 6)))) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
