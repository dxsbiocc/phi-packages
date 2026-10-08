#!/usr/bin/env Rscript

# Template-ID: boxplot-basic
#
# Purpose:
#   Draw a basic boxplot comparing a numeric value across categories,
#   with jittered points overlaid.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - Country: category label
#     - Income: numeric value
#
# Output:
#   A PDF, PNG, or SVG boxplot.
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
        x = "Country",
        y = "Income"
    ),
    labels = list(
        title = "Basic_Boxplot",
        x = "Country",
        y = "Income"
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
fill_values <- c("steelblue", "orange", "skyblue", "red", "green")

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]]
)) +
    geom_boxplot(aes(fill = .data[[x_col]]), alpha = 0.7) +
    geom_jitter(
        colour = "grey50",
        size = 2,
        alpha = 0.5,
        width = 0.2,
        height = 0
    ) +
    scale_fill_manual(values = fill_values) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism(axis_text_angle = 45) +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
