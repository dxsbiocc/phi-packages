#!/usr/bin/env Rscript

# Template-ID: line-group
#
# Purpose:
#   Draw a grouped line chart of a numeric value across ordered
#   categories, with a line and points per group.
#
# Inputs:
#   A table with one row per category-by-group. Default example:
#     - day: ordered category
#     - value: numeric value
#     - category: group
#
# Output:
#   A PDF, PNG, or SVG grouped line chart.
#
# Dependencies:
#   ggplot2, readr
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
#   Each row is one sequential observation in a group.
#   Category order in the file is treated as the x-axis order.
#   Connecting points with a line assumes an ordered relationship.

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
        colour = "category"
    ),
    labels = list(
        title = "Lines by Category",
        x = "Country",
        y = "Income"
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
colour_col <- config$columns$colour
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[colour_col]] <- factor(df[[colour_col]], levels = unique(df[[colour_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    colour = .data[[colour_col]]
)) +
    geom_line(aes(group = .data[[colour_col]])) +
    geom_point(colour = "white", size = 3) +
    geom_point(size = 2) +
    scale_colour_viridis_d() +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom",
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
