#!/usr/bin/env Rscript

# Template-ID: line-basic
#
# Purpose:
#   Draw a smoothed area chart of a numeric value across an ordered
#   category or sequence, with a spline curve, fill, and points.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - day: ordered category
#     - temperature: numeric value
#
# Output:
#   A PDF, PNG, or SVG smoothed area chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or the spline density.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one sequential observation.
#   Category order in the file is treated as the x-axis order.
#   The curve is a cubic spline through the observed points for display,
#   not a statistical smoother or a model fit.

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
        y = "temperature"
    ),
    labels = list(
        title = "Basic Line",
        x = "Day",
        y = "Temperature"
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
df$x_num <- as.numeric(df[[x_col]])
n_smooth <- max(200L, 40L * max(nrow(df) - 1L, 1L))
spl <- stats::spline(df$x_num, df[[y_col]], n = n_smooth)
smooth <- data.frame(x = spl$x, y = spl$y)
last <- df[nrow(df), , drop = FALSE]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_area(
        data = smooth,
        aes(x = x, y = y),
        fill = "#1aafd0",
        alpha = 0.3
    ) +
    geom_line(
        data = smooth,
        aes(x = x, y = y),
        colour = "#5070dd",
        linewidth = 1
    ) +
    geom_point(
        data = df,
        aes(x = x_num, y = .data[[y_col]]),
        shape = 21,
        fill = "white",
        colour = "#5070dd",
        size = 3,
        stroke = 1
    ) +
    geom_text(
        data = last,
        aes(x = x_num, y = .data[[y_col]], label = .data[[y_col]]),
        hjust = -0.4,
        size = 3.5
    ) +
    scale_x_continuous(
        breaks = df$x_num,
        labels = as.character(df[[x_col]])
    ) +
    scale_y_continuous(
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.08))
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
