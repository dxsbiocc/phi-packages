#!/usr/bin/env Rscript

# Template-ID: scatter-xyerror
#
# Purpose:
#   Scatter of two numeric axes with supplied bilateral error bars
#   (x and y), point size from a third magnitude, colour by a discrete
#   group, and optional supplied reference curve(s). Matches resonance
#   / diameter style figures: experimental points with SEM size bars
#   plus a simulation line.
#
# Inputs:
#   1. Points, one row per observation. Default columns:
#        - d, lambda, xerr, yerr, frequency, group
#   2. Curves, one row per polyline vertex. Default columns:
#        - d, lambda, series
#      Pass a curve table with only the header when no reference line
#      is needed.
#
# Output:
#   A PDF, PNG, or SVG scatter with capped x/y error bars and curves.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv curve.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, colours, axis expressions, size range,
#   limits, which groups appear in the legend). Edit DATA PREPARATION
#   when the error columns are ymin/ymax instead of half-widths. Do not
#   fit Gaussians, SEM histograms, or the simulation curve here.
#
# Scientific assumptions:
#   x, y, xerr, yerr, size, and group are supplied. xerr and yerr are
#   half-widths drawn as [x - xerr, x + xerr] and [y - yerr, y + yerr].
#   Size encodes the supplied frequency / magnitude column (area scale).
#   Reference curves are supplied vertices, not a model fit in this
#   script. This script does not estimate errors from replicates.

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

io <- parse_io_args(n_input = 2L, input_names = c("points", "curves"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        points = list(
            x = "d",
            y = "lambda",
            xerr = "xerr",
            yerr = "yerr",
            size = "frequency",
            group = "group"
        ),
        curves = list(
            x = "d",
            y = "lambda",
            series = "series"
        )
    ),
    # NULL = first-appearance order.
    group_order = c("experimental", "low_frequency"),
    # Hex copied from colors.json. Unlisted groups fall back to
    # Qualitative.Safe.
    group_fill = c(
        experimental = "#e73f74",   # Qualitative.Bold
        low_frequency = "#b3b3b3"   # Set2 / Pastel extras
    ),
    group_error = c(
        experimental = "#af6458",   # Qualitative.Antique
        low_frequency = "#b3b3b3"   # Set2
    ),
    # Groups drawn but omitted from the fill legend (figure shows only
    # the secondary group as a square). NULL keeps every group.
    legend_groups = c("low_frequency"),
    series_order = c("Simulation"),
    series_colours = c(
        Simulation = "#f2b701"      # Qualitative.Bold
    ),
    palettes = list(group = "Qualitative.Safe"),
    size_range = c(1.6, 9.5),
    point_stroke = 0.45,
    error_linewidth = 0.45,
    error_cap = 0.25,
    curve_linewidth = 1.55,
    # plotmath strings; plain labels$x / labels$y used when expr is NULL.
    labels = list(
        title = NULL,
        x = "Nano-object diameter d (nm)",
        y = "lambda_res (nm)",
        x_expr = "Nano-object~diameter~italic(d)~(nm)",
        y_expr = "lambda[res]~(nm)",
        group = NULL,
        series = NULL
    ),
    xlim = c(1, 13),
    ylim = c(600, 860),
    x_breaks = c(1, 5, 9, 13),
    y_breaks = c(600, 700, 800),
    show_size_legend = FALSE,
    size = list(width = 6.4, height = 5.4)
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pt_cols <- config$columns$points
cu_cols <- config$columns$curves

pts <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
curves <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
require_columns(pts, pt_cols)

for (nm in c("x", "y", "xerr", "yerr", "size")) {
    pts[[pt_cols[[nm]]]] <- as.numeric(pts[[pt_cols[[nm]]]])
}
pts[[pt_cols$group]] <- as.character(pts[[pt_cols$group]])

ok <- is.finite(pts[[pt_cols$x]]) &
    is.finite(pts[[pt_cols$y]]) &
    is.finite(pts[[pt_cols$xerr]]) & pts[[pt_cols$xerr]] >= 0 &
    is.finite(pts[[pt_cols$yerr]]) & pts[[pt_cols$yerr]] >= 0 &
    is.finite(pts[[pt_cols$size]]) & pts[[pt_cols$size]] > 0 &
    !is.na(pts[[pt_cols$group]]) & nzchar(pts[[pt_cols$group]])
if (any(!ok)) {
    message("Dropped ", sum(!ok), " point row(s) with non-finite values.")
}
pts <- pts[ok, , drop = FALSE]
if (!nrow(pts)) {
    stop("No finite point rows to plot.", call. = FALSE)
}

pts$`.xmin` <- pts[[pt_cols$x]] - pts[[pt_cols$xerr]]
pts$`.xmax` <- pts[[pt_cols$x]] + pts[[pt_cols$xerr]]
pts$`.ymin` <- pts[[pt_cols$y]] - pts[[pt_cols$yerr]]
pts$`.ymax` <- pts[[pt_cols$y]] + pts[[pt_cols$yerr]]

group_levels <- if (!is.null(config$group_order) && length(config$group_order)) {
    c(
        intersect(as.character(config$group_order), unique(pts[[pt_cols$group]])),
        setdiff(unique(pts[[pt_cols$group]]), as.character(config$group_order))
    )
} else {
    unique(pts[[pt_cols$group]])
}
pts[[pt_cols$group]] <- factor(pts[[pt_cols$group]], levels = group_levels)

fill_cols <- config$group_fill
fill_cols <- fill_cols[names(fill_cols) %in% group_levels]
miss <- setdiff(group_levels, names(fill_cols))
if (length(miss)) {
    extra <- expand_palette(palette_colors(config$palettes$group), length(miss))
    names(extra) <- miss
    fill_cols <- c(fill_cols, extra)
}
fill_cols <- fill_cols[group_levels]

err_cols <- config$group_error
err_cols <- err_cols[names(err_cols) %in% group_levels]
miss_err <- setdiff(group_levels, names(err_cols))
if (length(miss_err)) {
    err_cols <- c(err_cols, fill_cols[miss_err])
}
err_cols <- err_cols[group_levels]
legend_groups <- config$legend_groups
if (is.null(legend_groups)) {
    legend_groups <- as.character(group_levels)
} else {
    legend_groups <- intersect(as.character(legend_groups), as.character(group_levels))
}

has_curves <- nrow(curves) > 0L
series_levels <- character(0)
series_cols <- character(0)
if (has_curves) {
    require_columns(curves, cu_cols)
    curves[[cu_cols$x]] <- as.numeric(curves[[cu_cols$x]])
    curves[[cu_cols$y]] <- as.numeric(curves[[cu_cols$y]])
    curves[[cu_cols$series]] <- as.character(curves[[cu_cols$series]])
    cu_ok <- is.finite(curves[[cu_cols$x]]) &
        is.finite(curves[[cu_cols$y]]) &
        !is.na(curves[[cu_cols$series]]) & nzchar(curves[[cu_cols$series]])
    curves <- curves[cu_ok, , drop = FALSE]
    has_curves <- nrow(curves) > 0L
}
if (has_curves) {
    series_levels <- if (!is.null(config$series_order) && length(config$series_order)) {
        c(
            intersect(as.character(config$series_order), unique(curves[[cu_cols$series]])),
            setdiff(unique(curves[[cu_cols$series]]), as.character(config$series_order))
        )
    } else {
        unique(curves[[cu_cols$series]])
    }
    curves[[cu_cols$series]] <- factor(curves[[cu_cols$series]], levels = series_levels)
    series_cols <- config$series_colours
    series_cols <- series_cols[names(series_cols) %in% series_levels]
    miss_s <- setdiff(series_levels, names(series_cols))
    if (length(miss_s)) {
        extra_s <- expand_palette(palette_colors("Qualitative.Bold"), length(miss_s))
        names(extra_s) <- miss_s
        series_cols <- c(series_cols, extra_s)
    }
    series_cols <- series_cols[series_levels]
}

axis_label <- function(expr_str, plain) {
    if (!is.null(expr_str) && nzchar(as.character(expr_str))) {
        return(parse(text = expr_str)[[1]])
    }
    plain
}
x_lab <- axis_label(config$labels$x_expr, config$labels$x)
y_lab <- axis_label(config$labels$y_expr, config$labels$y)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot()

if (has_curves) {
    for (s in series_levels) {
        sub <- curves[curves[[cu_cols$series]] == s, , drop = FALSE]
        p <- p +
            geom_line(
                data = sub,
                aes(
                    x = .data[[cu_cols$x]],
                    y = .data[[cu_cols$y]],
                    group = 1
                ),
                colour = unname(series_cols[[s]]),
                linewidth = as.numeric(config$curve_linewidth),
                lineend = "round"
            )
    }
}

p <- p +
    geom_errorbar(
        data = pts,
        aes(
            x = .data[[pt_cols$x]],
            ymin = .data$`.ymin`,
            ymax = .data$`.ymax`,
            colour = .data[[pt_cols$group]]
        ),
        width = as.numeric(config$error_cap),
        linewidth = as.numeric(config$error_linewidth),
        show.legend = FALSE
    ) +
    geom_errorbar(
        data = pts,
        aes(
            y = .data[[pt_cols$y]],
            xmin = .data$`.xmin`,
            xmax = .data$`.xmax`,
            colour = .data[[pt_cols$group]]
        ),
        orientation = "y",
        width = as.numeric(config$error_cap),
        linewidth = as.numeric(config$error_linewidth),
        show.legend = FALSE
    ) +
    geom_point(
        data = pts,
        aes(
            x = .data[[pt_cols$x]],
            y = .data[[pt_cols$y]],
            size = .data[[pt_cols$size]],
            fill = .data[[pt_cols$group]]
        ),
        shape = 21,
        colour = "grey20",
        stroke = as.numeric(config$point_stroke)
    )

# Square keys drawn in the panel (as in the reference figure), not a
# ggplot legend. Coordinates are in data units; edit legend_marks.
legend_marks <- list()
if (length(legend_groups)) {
    for (i in seq_along(legend_groups)) {
        g <- legend_groups[[i]]
        legend_marks[[length(legend_marks) + 1L]] <- data.frame(
            x = as.numeric(config$xlim[[1]]) + 0.55,
            y = as.numeric(config$ylim[[2]]) - 12 - 28 * (i - 1L),
            fill = unname(fill_cols[[g]]),
            stringsAsFactors = FALSE
        )
    }
}
if (length(series_levels)) {
    for (i in seq_along(series_levels)) {
        s <- series_levels[[i]]
        legend_marks[[length(legend_marks) + 1L]] <- data.frame(
            x = as.numeric(config$xlim[[1]]) + 0.55 + 1.35 * i,
            y = as.numeric(config$ylim[[2]]) - 12,
            fill = unname(series_cols[[s]]),
            stringsAsFactors = FALSE
        )
    }
}
if (length(legend_marks)) {
    marks <- do.call(rbind, legend_marks)
    p <- p +
        geom_point(
            data = marks,
            aes(x = x, y = y),
            shape = 15,
            size = 4.8,
            colour = marks$fill,
            fill = marks$fill,
            show.legend = FALSE
        )
}

p <- p +
    scale_colour_manual(values = err_cols, guide = "none") +
    scale_fill_manual(values = fill_cols, guide = "none") +
    scale_size_area(
        max_size = max(as.numeric(config$size_range)),
        guide = if (isTRUE(config$show_size_legend)) "legend" else "none"
    ) +
    scale_x_continuous(
        breaks = as.numeric(config$x_breaks),
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        breaks = as.numeric(config$y_breaks),
        expand = c(0, 0)
    ) +
    coord_cartesian(
        xlim = as.numeric(config$xlim),
        ylim = as.numeric(config$ylim),
        clip = "on"
    ) +
    labs(
        title = config$labels$title,
        x = x_lab,
        y = y_lab
    )

p <- p +
    theme_bw(base_size = 12) +
    theme(
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.ticks.length = unit(-3.2, "pt"),
        axis.text = element_text(colour = "black", size = 11),
        axis.text.x = element_text(margin = margin(t = 6)),
        axis.text.y = element_text(margin = margin(r = 6)),
        axis.title = element_text(colour = "black", size = 12),
        legend.position = "none",
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 12, 8, 10)
    )

if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
    p <- p + theme(plot.title = element_blank())
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = config$size$width,
    height = config$size$height
)
