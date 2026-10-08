#!/usr/bin/env Rscript

# Template-ID: line-bump
#
# Purpose:
#   Draw a smoothed bump chart of ranks over time for multiple series.
#
# Inputs:
#   A wide table with one row per time point. Default example:
#     - year: ordered time
#     - remaining columns: one rank series per item (lower is better)
#
# Output:
#   A PDF, PNG, or SVG bump chart.
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
#   Each cell is an already computed rank at one time point.
#   The y-axis is reversed so rank 1 appears at the top.
#   Connecting ranks with a spline is a display interpolation, not a
#   model of rank trajectories.
#   The y-axis is reversed so rank 1 appears at the top.

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
        x = "year"
    ),
    labels = list(
        title = "Bump Chart",
        x = "Year",
        y = "Rank"
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
series_cols <- setdiff(names(df), x_col)
if (!length(series_cols)) {
    stop("Bump chart input needs at least one series column besides x.",
        call. = FALSE
    )
}
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
long <- data.frame(
    x = rep(df[[x_col]], times = length(series_cols)),
    series = factor(
        rep(series_cols, each = nrow(df)),
        levels = series_cols
    ),
    rank = as.numeric(as.matrix(df[series_cols])),
    stringsAsFactors = FALSE
)
long$x_num <- as.numeric(long$x)
smooth <- do.call(rbind, lapply(split(long, long$series), function(d) {
    d <- d[order(d$x_num), ]
    n_smooth <- max(80L, 30L * max(nrow(d) - 1L, 1L))
    spl <- stats::spline(d$x_num, d$rank, n = n_smooth)
    data.frame(
        x = spl$x,
        rank = spl$y,
        series = d$series[[1]],
        stringsAsFactors = FALSE
    )
}))
smooth$series <- factor(smooth$series, levels = levels(long$series))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_line(
        data = smooth,
        aes(x = x, y = rank, colour = series, group = series),
        linewidth = 1.2
    ) +
    geom_point(
        data = long,
        aes(x = x_num, y = rank, colour = series),
        shape = 21,
        fill = "white",
        size = 3.8,
        stroke = 1.3
    ) +
    scale_colour_manual(
        values = expand_palette(
            c("#5070dd",
            "#b6d634",
            "#505372",
            "#ff994d",
            "#0ca8df",
            "#ffd10a",
            "#fb628b",
            "#785db0",
            "#3fbe95"),
            nlevels(long$series)
        )
    ) +
    scale_x_continuous(
        breaks = sort(unique(long$x_num)),
        labels = levels(long$x)
    ) +
    scale_y_reverse(breaks = sort(unique(long$rank))) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        colour = ""
    ) +
    theme_prism() +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
