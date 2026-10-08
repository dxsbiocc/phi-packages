#!/usr/bin/env Rscript

# Template-ID: line-double
#
# Purpose:
#   Draw a dual-axis hydrograph: flow rises from the bottom, rainfall hangs
#   from the top on a reversed secondary axis.
#
# Inputs:
#   A table with one row per time point. Default example:
#     - date: timestamp
#     - Flow(m³/s): primary numeric series (bottom axis)
#     - Rainfall(mm): secondary numeric series (top, inverted)
#
# Output:
#   A PDF, PNG, or SVG dual-axis time-series chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change time parsing or axis limits.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one time point with two already measured series.
#   The secondary axis is an inverted linear rescaling for display,
#   not a converted unit.
#   Connecting points with a line assumes an ordered time relationship.

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
        y1 = "Flow(m³/s)",
        y2 = "Rainfall(mm)"
    ),
    labels = list(
        title = "Rainfall and Flow Relationship",
        x = NULL,
        y1 = "Flow(m³/s)",
        y2 = "Rainfall(mm)",
        series1 = "Flow",
        series2 = "Rainfall"
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y1_col <- config$columns$y1
y2_col <- config$columns$y2

parse_time <- function(x) {
    parsed <- as.POSIXct(x, format = "%Y/%m/%d %H:%M", tz = "UTC")
    if (all(is.na(parsed))) {
        parsed <- as.POSIXct(x, tz = "UTC")
    }
    parsed
}

df[[x_col]] <- parse_time(df[[x_col]])
df[[y1_col]] <- as.numeric(df[[y1_col]])
df[[y2_col]] <- as.numeric(df[[y2_col]])
keep <- !is.na(df[[x_col]]) & !is.na(df[[y1_col]]) & !is.na(df[[y2_col]])
df <- df[keep, , drop = FALSE]
df <- df[order(df[[x_col]]), , drop = FALSE]

flow_max <- max(pretty(c(0, max(df[[y1_col]]))), na.rm = TRUE)
rain_max <- max(pretty(c(0, max(df[[y2_col]]))), na.rm = TRUE)
if (!is.finite(flow_max) || flow_max <= 0) {
    flow_max <- 1
}
if (!is.finite(rain_max) || rain_max <= 0) {
    rain_max <- 1
}
coeff <- flow_max / rain_max
df$flow_y <- df[[y1_col]]
df$rain_y <- flow_max - df[[y2_col]] * coeff

x_breaks <- pretty(range(df[[x_col]]), n = 8L)
nearest <- function(target) {
    which.min(abs(as.numeric(df[[x_col]]) - as.numeric(target)))
}
pts <- df[unique(vapply(x_breaks, nearest, integer(1))), , drop = FALSE]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
col_flow <- "#037ef3"
col_rain <- "#f85a40"
series1 <- config$labels$series1
series2 <- config$labels$series2

p <- ggplot(df, aes(x = .data[[x_col]])) +
    geom_ribbon(
        aes(ymin = 0, ymax = flow_y, fill = series1),
        alpha = 0.7,
        colour = NA
    ) +
    geom_ribbon(
        aes(ymin = rain_y, ymax = flow_max, fill = series2),
        alpha = 0.7,
        colour = NA
    ) +
    geom_line(aes(y = flow_y, colour = series1), linewidth = 0.55) +
    geom_line(aes(y = rain_y, colour = series2), linewidth = 0.55) +
    geom_point(
        data = pts,
        aes(y = flow_y, colour = series1),
        shape = 21,
        fill = "white",
        size = 2.2,
        stroke = 0.9
    ) +
    geom_point(
        data = pts,
        aes(y = rain_y, colour = series2),
        shape = 21,
        fill = "white",
        size = 2.2,
        stroke = 0.9
    ) +
    scale_x_datetime(breaks = x_breaks, date_labels = "%Y/%m/%d\n%H:%M") +
    scale_y_continuous(
        name = config$labels$y1,
        limits = c(0, flow_max),
        expand = c(0, 0),
        sec.axis = sec_axis(
            ~ (flow_max - .) / coeff,
            name = config$labels$y2
        )
    ) +
    scale_fill_manual(
        values = setNames(c(col_flow, col_rain), c(series1, series2))
    ) +
    scale_colour_manual(
        values = setNames(c(col_flow, col_rain), c(series1, series2))
    ) +
    labs(title = config$labels$title, x = config$labels$x) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom"
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
