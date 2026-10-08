#!/usr/bin/env Rscript

# Template-ID: scatter-dumbbell-groups
#
# Purpose:
#   Draw a horizontal multi-group dumbbell plot of a numeric value
#   across categories, with a point per year.
#
# Inputs:
#   A table with one row per category-by-year. Default example:
#     - mean_prev: numeric value
#     - age_group: category label
#     - year_id: group (year)
#
# Output:
#   A PDF, PNG, or SVG multi-group dumbbell plot.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category or year order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category-by-year pair with a comparable numeric value.
#   Values are treated as already summarized, not as raw replicates.
#   Connecting years within a category is a display choice, not a test.

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
        x = "mean_prev",
        y = "age_group",
        fill = "year_id"
    ),
    labels = list(
        title = "",
        x = "Mean Prevalence of obesity",
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
fill_col <- config$columns$fill
df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
df[[fill_col]] <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))
df[[x_col]] <- as.numeric(df[[x_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#f784b6", "#ffdf8b", "#79ceb8")
shape_values <- c(21, 22, 23)
x_breaks <- c(0.06, 0.07, 0.08, 0.09, 0.10, 0.11)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]]
)) +
    geom_line(
        aes(group = .data[[y_col]]),
        colour = "#cba7a2",
        linewidth = 1
    ) +
    geom_vline(
        xintercept = x_breaks,
        linetype = "dashed",
        colour = "grey90"
    ) +
    geom_point(
        aes(
            fill = .data[[fill_col]],
            shape = .data[[fill_col]]
        ),
        size = 4,
        stroke = 1,
        colour = "#ffffff"
    ) +
    scale_fill_manual(values = fill_values) +
    scale_shape_manual(values = shape_values) +
    scale_x_continuous(breaks = x_breaks) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.87),
        legend.direction = "horizontal",
        axis.line.x = element_blank(),
        axis.line.y = element_line(colour = "grey90", linewidth = 0.8),
        axis.ticks = element_blank(),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
