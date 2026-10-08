#!/usr/bin/env Rscript

# Template-ID: scatter-jitter
#
# Purpose:
#   Draw a jittered scatter plot of a numeric value across categories.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - day: category label
#     - value: numeric value
#
# Output:
#   A PDF, PNG, or SVG jitter plot.
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
#   Jitter is a display transform to reduce overplotting, not a
#   statistical transform of the values.

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
        y = "value"
    ),
    labels = list(
        title = "Jitter Plot of Values by Day",
        x = "Value",
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
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c(
    "#8F499C", "#4185BE", "#6DC067", "#F6DB35",
    "#F78822", "#232121", "#CF202E"
)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]]
)) +
    geom_jitter(
        aes(colour = .data[[x_col]]),
        size = 3,
        width = 0.18,
        alpha = 0.8,
        show.legend = FALSE
    ) +
    scale_colour_manual(values = colour_values) +
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
