#!/usr/bin/env Rscript

# Template-ID: line-independent
#
# Purpose:
#   Draw two independent-scale time series as stacked panels sharing the time
#   axis: the top series rises from zero, the bottom series hangs downward.
#
# Inputs:
#   A table with one row per time point. Default example:
#     - date: timestamp
#     - Flow(m³/s): first numeric series
#     - Rainfall(mm): second numeric series
#
# Output:
#   A PDF, PNG, or SVG two-panel time-series chart.
#
# Dependencies:
#   ggplot2, readr, patchwork, ggprism
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
#   Each row is one time point with two already measured series.
#   Panels have independent y-scales; no dual-axis conversion is applied.
#   The bottom panel reverses its y-axis so zero sits at the shared midline
#   and larger values point downward.
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
        title = "Independent Axis Time Series",
        x = "Date",
        y1 = "Flow(m³/s)",
        y2 = "Rainfall(mm)"
    )
)

load_packages(c("ggplot2", "readr", "patchwork", "ggprism"))

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

x_range <- range(df[[x_col]])
x_breaks <- pretty(x_range, n = 8L)
nearest <- function(target) {
    which.min(abs(as.numeric(df[[x_col]]) - as.numeric(target)))
}
pts <- df[unique(vapply(x_breaks, nearest, integer(1))), , drop = FALSE]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
col_flow <- "#037ef3"
col_rain <- "#f85a40"

shared_x <- scale_x_datetime(
    limits = x_range,
    breaks = x_breaks,
    date_labels = "%m/%d\n%H:%M",
    expand = expansion(mult = c(0.01, 0.01))
)

p1 <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y1_col]])) +
    geom_area(fill = col_flow, alpha = 0.35, colour = NA) +
    geom_line(colour = col_flow, linewidth = 0.55) +
    geom_point(
        data = pts,
        shape = 21,
        fill = "white",
        colour = col_flow,
        size = 2.2,
        stroke = 0.9
    ) +
    shared_x +
    scale_y_continuous(
        limits = c(0, flow_max),
        expand = c(0, 0)
    ) +
    labs(x = NULL, y = config$labels$y1) +
    theme_prism() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.background = element_blank()
    )

p2 <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y2_col]])) +
    geom_area(fill = col_rain, alpha = 0.35, colour = NA) +
    geom_line(colour = col_rain, linewidth = 0.55) +
    geom_point(
        data = pts,
        shape = 21,
        fill = "white",
        colour = col_rain,
        size = 2.2,
        stroke = 0.9
    ) +
    shared_x +
    scale_y_reverse(
        limits = c(0, rain_max),
        expand = c(0, 0)
    ) +
    labs(x = config$labels$x, y = config$labels$y2) +
    theme_prism() +
    theme(plot.background = element_blank())

p <- p1 / p2 + plot_annotation(title = config$labels$title)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
