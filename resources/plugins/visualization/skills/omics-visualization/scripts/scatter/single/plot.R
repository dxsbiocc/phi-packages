#!/usr/bin/env Rscript

# Template-ID: scatter-single
#
# Purpose:
#   Draw a single-axis-style scatter of hourly values, one row per day.
#
# Inputs:
#   A table with one row per hour-by-day. Default example:
#     - hours: category on the shared axis
#     - value: numeric magnitude
#     - days: row (series)
#
# Output:
#   A PDF, PNG, or SVG single-axis scatter plot.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or reshape the table.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one already summarized hour-by-day value.
#   Point size encodes the value; zero values remain as rows but are small.
#   Hour and day order in the file are display order.

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
        x = "hours",
        y = "days",
        size = "value"
    ),
    labels = list(
        title = "Single-Axis Scatter",
        x = "Hour",
        y = "Day"
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
size_col <- config$columns$size
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
df[[size_col]] <- as.numeric(df[[size_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    size = .data[[size_col]],
    colour = .data[[y_col]]
)) +
    geom_point(alpha = 0.85) +
    scale_size_area(max_size = 10) +
    scale_colour_manual(
        values = expand_palette(
            c("#1cc7d0",
            "#2dde98",
            "#ffc168",
            "#ff6c5f",
            "#ff4f81",
            "#b84592",
            "#8e43e7"),
            nlevels(df[[y_col]])
        ),
        guide = "none"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        size = "Value"
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        legend.position = "top"
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 14, height = 6)
