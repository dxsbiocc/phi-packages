#!/usr/bin/env Rscript

# Template-ID: line-stack
#
# Purpose:
#   Draw a stacked area chart of a numeric value across ordered categories by group.
#
# Inputs:
#   A table with one row per category-by-group. Default example:
#     - day: ordered category
#     - value: numeric value
#     - category: group
#
# Output:
#   A PDF, PNG, or SVG stacked area chart.
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
#   Each row is one sequential observation in a group.
#   Stacking areas is a display of already summarized values.
#   Category order in the file is treated as the x-axis order.

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
        x = "day",
        y = "value",
        fill = "category"
    ),
    labels = list(
        title = "Stacked Area Chart",
        x = "Day",
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
    fill = .data[[fill_col]],
    group = .data[[fill_col]]
)) +
    geom_area(position = "stack", colour = "white", linewidth = 0.2) +
    scale_fill_manual(
        values = expand_palette(
            c("#3be8b0",
            "#1aafd0",
            "#6a67ce",
            "#ffb900",
            "#fc636b"),
            nlevels(df[[fill_col]])
        )
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        fill = NULL
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal"
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
