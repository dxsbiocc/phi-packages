#!/usr/bin/env Rscript

# Template-ID: boxplot-group
#
# Purpose:
#   Draw grouped boxplots comparing a numeric value across categories,
#   with a fill aesthetic for a second grouping variable and jittered points.
#
# Inputs:
#   A table with one row per observation. Default example:
#     - Country: category on the x-axis
#     - Income: numeric value
#     - Year: fill group
#
# Output:
#   A PDF, PNG, or SVG grouped boxplot.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category or group order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one observation of a numeric value in a category and group.
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
        y = "Income",
        fill = "Year"
    ),
    labels = list(
        title = "Group_Boxplot",
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
fill_col <- config$columns$fill
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[fill_col]] <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[fill_col]]
)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_jitter(
        position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
        alpha = 0.5
    ) +
    scale_fill_brewer(palette = "Set3") +
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
