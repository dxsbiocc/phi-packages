#!/usr/bin/env Rscript

# Template-ID: bar-group
#
# Purpose:
#   Draw a grouped bar chart comparing a numeric value across years by land-cover group.
#
# Inputs:
#   A table with one row per year-by-group. Default example:
#     - year: category label
#     - value: numeric value
#     - category: group
#
# Output:
#   A PDF, PNG, or SVG grouped bar chart.
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
#   Each row is one year-by-group pair with a comparable numeric value.
#   Dodging is a display choice, not a statistical transform.
#   Year order in the file is treated as the x-axis order.

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
        x = "year",
        y = "value",
        fill = "category"
    ),
    labels = list(
        title = "Grouped Bar Chart",
        x = "Year",
        y = "Value"
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
    geom_col(position = position_dodge(width = 0.85), width = 0.8, colour = "white") +
    geom_text(
        aes(label = .data[[y_col]]),
        position = position_dodge(width = 0.85),
        vjust = -0.35,
        size = 3
    ) +
    scale_fill_manual(
        values = expand_palette(
            c("#1cc7d0",
            "#2dde98",
            "#ffc168",
            "#ff6c5f"),
            nlevels(df[[fill_col]])
        )
    ) +
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
