#!/usr/bin/env Rscript

# Template-ID: line-highlight
#
# Purpose:
#   Draw a smoothed area line with a highlighted window around the peak.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - date: ordered date
#     - value: numeric value
#
# Output:
#   A PDF, PNG, or SVG line chart with a highlight band.
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
#   The highlight band marks the peak-adjacent dates for display.
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
        x = "date",
        y = "value"
    ),
    labels = list(
        title = "Highlighted Line Chart",
        x = "Date",
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
df[[x_col]] <- as.Date(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df <- df[order(df[[x_col]]), , drop = FALSE]
df$x_num <- as.numeric(df[[x_col]])
n_smooth <- max(200L, 40L * max(nrow(df) - 1L, 1L))
spl <- stats::spline(df$x_num, df[[y_col]], n = n_smooth)
smooth <- data.frame(
    x = as.Date(spl$x, origin = "1970-01-01"),
    y = spl$y
)
peak_i <- which.max(df[[y_col]])
peak_lo <- max(1L, peak_i - 1L)
peak_hi <- min(nrow(df), peak_i + 1L)
xmin <- df[[x_col]][[peak_lo]]
xmax <- df[[x_col]][[peak_hi]]
hl <- smooth[smooth$x >= xmin & smooth$x <= xmax, , drop = FALSE]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_area(
        data = hl,
        aes(x = x, y = y),
        fill = "#3BE8B0",
        alpha = 0.7
    ) +
    geom_line(
        data = smooth,
        aes(x = x, y = y),
        colour = "#30c39e",
        linewidth = 1.2
    ) +
    geom_point(
        data = df,
        aes(x = .data[[x_col]], y = .data[[y_col]]),
        shape = 21,
        fill = "white",
        colour = "#30c39e",
        size = 2.6,
        stroke = 1
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
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
