#!/usr/bin/env Rscript

# Template-ID: radar-group
#
# Purpose:
#   Draw one radar chart per series in a faceted grid (not overlaid).
#
# Inputs:
#   A wide table with one row per axis. Default example:
#     - name: axis label
#     - remaining columns: one numeric series per model
#
# Output:
#   A PDF, PNG, or SVG grouped radar chart.
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
#   Each series is drawn in its own panel. All panels share the same axes.
#   Closing the polygon repeats the first axis and is a display device.
#   Axes are projected to Cartesian coordinates so edges stay straight.
#   Scores are treated as already computed, not as raw replicates.

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
        axis = "name"
    ),
    labels = list(
        title = "Multiple Radar",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
axis_col <- config$columns$axis
series_cols <- setdiff(names(df), axis_col)
df[[axis_col]] <- factor(df[[axis_col]], levels = unique(df[[axis_col]]))
n_axes <- nlevels(df[[axis_col]])
long <- data.frame(
    axis = rep(df[[axis_col]], times = length(series_cols)),
    series = factor(
        rep(series_cols, each = nrow(df)),
        levels = series_cols
    ),
    value = as.numeric(as.matrix(df[series_cols])),
    stringsAsFactors = FALSE
)
long$idx <- as.numeric(long$axis)
n_split <- 5L
y_max <- max(pretty(c(0, max(long$value, na.rm = TRUE))))
if (!is.finite(y_max) || y_max <= 0) {
    y_max <- 1
}

to_xy <- function(idx, r) {
    theta <- pi / 2 + 2 * pi * (idx - 1) / n_axes
    data.frame(
        px = r * cos(theta),
        py = r * sin(theta),
        stringsAsFactors = FALSE
    )
}

xy <- to_xy(long$idx, long$value)
long$px <- xy$px
long$py <- xy$py
close_i <- long$idx == 1L
close <- long[close_i, , drop = FALSE]
long <- rbind(long, close)

ring_cols <- c("#E8F8F2", "#FDF6E6")
rings <- do.call(rbind, lapply(seq_len(n_split), function(i) {
    d <- to_xy(c(seq_len(n_axes), 1L), y_max * i / n_split)
    d$ring <- i
    d
}))
ring_geoms <- lapply(rev(seq_len(n_split)), function(i) {
    geom_polygon(
        data = rings[rings$ring == i, , drop = FALSE],
        aes(x = px, y = py),
        fill = ring_cols[[((i - 1L) %% 2L) + 1L]],
        colour = NA,
        inherit.aes = FALSE
    )
})
spokes <- to_xy(seq_len(n_axes), y_max)
axis_lab <- to_xy(seq_len(n_axes), y_max * 1.16)
axis_lab$label <- levels(df[[axis_col]])
axis_lab$hjust <- ifelse(abs(axis_lab$px) < y_max * 0.05, 0.5,
    ifelse(axis_lab$px < 0, 1.08, -0.08)
)
axis_lab$vjust <- ifelse(abs(axis_lab$py) < y_max * 0.05, 0.5,
    ifelse(axis_lab$py < 0, 1.15, -0.15)
)
lim <- y_max * 1.38

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(long, aes(
    x = px,
    y = py,
    colour = series,
    fill = series,
    group = series
)) +
    ring_geoms +
    geom_path(
        data = rings,
        aes(x = px, y = py, group = ring),
        colour = "#dbdee4",
        linewidth = 0.35,
        inherit.aes = FALSE
    ) +
    geom_segment(
        data = spokes,
        aes(x = 0, y = 0, xend = px, yend = py),
        colour = "#cfd2d7",
        linewidth = 0.4,
        inherit.aes = FALSE
    ) +
    geom_polygon(alpha = 0.12, linewidth = 0.7) +
    geom_point(size = 2) +
    geom_text(
        data = axis_lab,
        aes(x = px, y = py, label = label, hjust = hjust, vjust = vjust),
        inherit.aes = FALSE,
        size = 3.2
    ) +
    scale_colour_manual(
        values = expand_palette(
            c(
                "#1cc7d0",
                "#2dde98",
                "#ffc168",
                "#ff6c5f"
            ),
            nlevels(long$series)
        )
    ) +
    scale_fill_manual(
        values = expand_palette(
            c(
                "#1cc7d0",
                "#2dde98",
                "#ffc168",
                "#ff6c5f"
            ),
            nlevels(long$series)
        )
    ) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    facet_wrap(~series, ncol = min(2L, nlevels(long$series))) +
    labs(
        title = config$labels$title,
        colour = "",
        fill = "",
        x = NULL,
        y = NULL
    ) +
    theme_prism() +
    theme(
        plot.margin = margin(16, 16, 16, 16),
        plot.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 10, height = 8)
