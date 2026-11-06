#!/usr/bin/env Rscript

# Template-ID: line-spiral
#
# Purpose:
#   Archimedean spiral for periodic series: one turn per cycle (e.g.
#   year), angle = step within the cycle (month / day), radius rises
#   continuously so turns join end to start. Geometry is a line,
#   points, rising parallelogram bars, a full-height spiral heatmap,
#   or an area band (solid or value-coloured).
#
# Inputs:
#   One row per cycle x step (already at the plotted resolution).
#   Default example:
#     - year: cycle
#     - month: step within the cycle (1..steps_per_cycle)
#     - samples: value (colour; also bar / area height)
#     - range: optional bar height column (bar mode)
#
# Output:
#   A PDF, PNG, or SVG spiral plot.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, geom, steps per cycle, step labels,
#   heights, inner hole, palette, cycle labels). Line / point / bar /
#   heatmap / area are variants of one id. Do not aggregate (e.g.
#   monthly means) in this script; supply one row per cycle x step.
#
# Scientific assumptions:
#   Values are supplied. Radius encodes elapsed steps only; bar and
#   area heights are a linear display scaling of a supplied column
#   into one turn's spacing, not a second axis.

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
    columns = list(cycle = "year", step = "month", value = "samples",
        height = "range"),           # height: NULL = use value
    geom = "line",                   # "line", "point", "bar", "heatmap", "area"
    gradient = TRUE,                 # area: colour each step by value
    steps_per_cycle = 12,
    step_breaks = 1:12,              # angles labelled (NULL = none)
    step_labels = month.abb,
    height_range = c(0.25, 0.95),    # bar / area height, fraction of one turn
    inner = 0.6,                     # hole radius, in turns
    line_width = 1.4,
    area_line_width = 1.5,           # gradient area: one segment per step
    point_size = 2.6,
    palette = "Diverging.Spectral",
    reverse_palette = TRUE,
    area_fill = "Qualitative.Vivid:1",
    show_cycle_labels = TRUE,
    labels = list(fill = "Samples"),
    size = list(width = 7, height = 7)
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}
cc <- config$columns
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, cc[c("cycle", "step", "value")])
S <- as.numeric(config$steps_per_cycle)
d <- data.frame(cycle = as.numeric(df[[cc$cycle]]),
    step = as.numeric(df[[cc$step]]), value = as.numeric(df[[cc$value]]))
has_h <- !is.null(cc$height) && cc$height %in% names(df)
d$h_src <- if (has_h) as.numeric(df[[cc$height]]) else d$value
d <- d[stats::complete.cases(d), , drop = FALSE]
if (any(d$step < 1 | d$step > S)) {
    stop("step must lie in 1..steps_per_cycle.", call. = FALSE)
}
d <- d[order(d$cycle, d$step), , drop = FALSE]
c0 <- min(d$cycle)
d$base <- (d$cycle - c0) * S + d$step - 1
hr <- config$height_range
rng <- range(d$h_src)
d$h <- if (diff(rng) > 0) {
    (hr[[1]] + (d$h_src - rng[[1]]) / diff(rng) * diff(hr)) * S
} else {
    hr[[2]] * S
}

stops <- palette_colors(config$palette)
if (isTRUE(config$reverse_palette)) stops <- rev(stops)
n_cycle <- length(unique(d$cycle))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
g <- config$geom
p <- ggplot()
if (g %in% c("line", "point")) {
    d$x <- d$step
    nxt <- d[d$step == 1 & d$cycle > c0, , drop = FALSE]
    if (identical(g, "line") && nrow(nxt)) {
        nxt$cycle <- nxt$cycle - 1
        nxt$x <- S + 1
        d2 <- rbind(d, nxt)
    } else {
        d2 <- d
    }
    d2 <- d2[order(d2$cycle, d2$x), , drop = FALSE]
    if (identical(g, "line")) {
        p <- p + geom_path(data = d2, aes(x = x, y = base, group = cycle,
            colour = value), linewidth = config$line_width,
            lineend = "round")
    } else {
        p <- p + geom_point(data = d, aes(x = x, y = base, colour = value),
            size = config$point_size)
    }
    p <- p + scale_colour_gradientn(colours = stops,
        name = config$labels$fill)
} else if (g %in% c("bar", "heatmap")) {
    hh <- if (identical(g, "heatmap")) rep(S, nrow(d)) else d$h
    poly <- data.frame(
        id = rep(seq_len(nrow(d)), each = 4),
        x = as.vector(rbind(d$step, d$step + 1, d$step + 1, d$step)),
        y = as.vector(rbind(d$base, d$base + 1, d$base + 1 + hh, d$base + hh)),
        value = rep(d$value, each = 4)
    )
    p <- p + geom_polygon(data = poly, aes(x = x, y = y, group = id,
        fill = value), colour = if (identical(g, "heatmap")) "white" else
        "grey25", linewidth = 0.25) +
        scale_fill_gradientn(colours = stops, name = config$labels$fill)
} else if (identical(g, "area")) {
    d$x <- d$step
    d$ymax <- d$base + d$h
    if (isTRUE(config$gradient)) {
        p <- p + geom_linerange(data = d, aes(x = x, ymin = base, ymax = ymax,
            colour = value), linewidth = config$area_line_width) +
            scale_colour_gradientn(colours = stops, name = config$labels$fill)
    } else {
        p <- p + geom_ribbon(data = d, aes(x = x, ymin = base, ymax = ymax,
            group = cycle), fill = pref(config$area_fill))
    }
    p <- p +
        geom_line(data = d, aes(x = x, y = base, group = cycle),
            colour = "grey20", linewidth = 0.3) +
        geom_line(data = d, aes(x = x, y = ymax, group = cycle),
            colour = "grey40", linewidth = 0.25)
} else {
    stop("geom must be line, point, bar, heatmap or area.", call. = FALSE)
}

if (isTRUE(config$show_cycle_labels)) {
    cl <- data.frame(x = 1, y = (sort(unique(d$cycle)) - c0) * S +
        if (g %in% c("bar", "heatmap", "area")) S * 0.5 else 0,
        lab = sort(unique(d$cycle)))
    p <- p + geom_label(data = cl, aes(x = x, y = y, label = lab),
        size = 2.8, linewidth = 0, fill = grDevices::adjustcolor("white",
            0.75), label.padding = unit(0.12, "lines"))
}

top <- n_cycle * S + if (g %in% c("bar", "heatmap", "area")) S else 1
brk <- if (is.null(config$step_breaks)) NULL else
    as.numeric(config$step_breaks) + if (g %in% c("bar", "heatmap")) 0.5 else 0
p <- p +
    scale_x_continuous(limits = c(1, S + 1), breaks = brk,
        labels = config$step_labels, expand = c(0, 0)) +
    scale_y_continuous(limits = c(-config$inner * S, top), expand = c(0, 0)) +
    coord_polar(clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, colour = "grey20"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.key.height = unit(12, "mm")
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width,
    height = config$size$height)
