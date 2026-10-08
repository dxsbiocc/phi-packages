#!/usr/bin/env Rscript

# Template-ID: scatter-dumbbell
#
# Purpose:
#   Draw a dumbbell plot comparing a numeric value across two groups
#   on a shared categorical axis.
#
# Inputs:
#   A table with one row per category-by-group. Default example:
#     - age_group: category label
#     - mean_prev: numeric value
#     - Sex: group (two levels)
#
# Output:
#   A PDF, PNG, or SVG dumbbell plot.
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
#   Each row is one category-by-group pair with a comparable numeric value.
#   Values are treated as already summarized, not as raw replicates.
#   Connecting the two groups is a display choice, not a statistical test.

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
        x = "age_group",
        y = "mean_prev",
        colour = "Sex"
    ),
    labels = list(
        title = "",
        x = "",
        y = "Mean Prevalence of obesity"
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
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[colour_col]] <- factor(df[[colour_col]], levels = unique(df[[colour_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c("#f784b6", "#79ceb8")

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    colour = .data[[colour_col]]
)) +
    geom_line(aes(group = .data[[x_col]]), colour = "#cba7a2") +
    geom_point(size = 5, shape = 21, fill = "#ffffff") +
    geom_point(size = 3) +
    scale_colour_manual(values = colour_values) +
    guides(x = guide_prism_bracket()) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.text.x = element_text(angle = 90),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
