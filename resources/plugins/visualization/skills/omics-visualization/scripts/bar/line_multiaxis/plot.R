#!/usr/bin/env Rscript

# Template-ID: bar-line-multiaxis
#
# Purpose:
#   One category axis with any number of series, each on its own
#   coloured y axis: bar series side by side, line series with open
#   points on top, axes stacked on the left and / or right with
#   their titles on top, shared tick count so gridlines align.
#
# Inputs:
#   One row per category (wide table). Default example:
#     - month: category (file order)
#     - Evaporation, Precipitation, Temperature: numeric series
#
# Output:
#   A PDF, PNG, or SVG multi-axis bar + line chart.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG: one `series` entry per column (type bar / line,
#   side left / right, colour, limits, unit, axis title). Extra
#   series are extra entries; axis count is not a new id. The first
#   left-side series sets the gridlines. Two mirrored hydrograph
#   series (flow up, rain down) are line-double.
#
# Scientific assumptions:
#   Series values are supplied. Each series is linearly mapped onto
#   the gridline scale for drawing only; every axis shows its own
#   supplied units. Axis limits are display choices.

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
    category = "month",
    # order = bar dodge / legend order; colours are "Family.Name:index"
    series = list(
        list(column = "Evaporation", type = "bar", side = "right",
            color = "Brand.Algolia:8", unit = " ml", title = "Evaporation",
            limits = NULL),
        list(column = "Precipitation", type = "bar", side = "right",
            color = "Qualitative.Bold:6", unit = " ml",
            title = "Precipitation", limits = NULL),
        list(column = "Temperature", type = "line", side = "left",
            color = "Qualitative.Bold:11", unit = " \u00b0C",
            title = "Temperature", limits = NULL)
    ),
    n_intervals = 6,             # shared tick intervals on every axis
    axis_gap = 1.8,              # spacing between stacked axes, category units
    bar_width = 0.7,             # total width of dodged bars per category
    show_grid = TRUE,
    legend_position = "bottom",
    size = list(width = 10, height = 7)
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pref <- function(ref) {
    if (!grepl(".", ref, fixed = TRUE)) return(ref)
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}
nice_step <- function(x) {
    if (!is.finite(x) || x <= 0) return(1)
    e <- 10^floor(log10(x))
    m <- c(1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10)
    m[which(m * e >= x - 1e-9)[[1]]] * e
}

df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
ser <- config$series
cols_needed <- c(config$category, vapply(ser, `[[`, "", "column"))
require_columns(df, as.list(stats::setNames(cols_needed, cols_needed)))
cats <- as.character(df[[config$category]])
n <- length(cats)
ni <- as.integer(config$n_intervals)

for (i in seq_along(ser)) {
    v <- as.numeric(df[[ser[[i]]$column]])
    if (any(!is.finite(v))) {
        stop(sprintf("Non-finite value in %s.", ser[[i]]$column),
            call. = FALSE)
    }
    lim <- ser[[i]]$limits
    if (is.null(lim)) {
        lo <- min(0, v)
        step <- nice_step((max(v) - lo) / ni)
        lo <- floor(lo / step) * step
        lim <- c(lo, lo + step * ni)
    }
    ser[[i]]$v <- v
    ser[[i]]$lim <- as.numeric(lim)
    ser[[i]]$col <- pref(ser[[i]]$color)
}
sides <- vapply(ser, `[[`, "", "side")
prim <- which(sides == "left")
prim <- if (length(prim)) prim[[1]] else 1L
P <- ser[[prim]]$lim
to_p <- function(v, lim) (v - lim[[1]]) / diff(lim) * diff(P) + P[[1]]

left_i <- which(sides == "left")
right_i <- which(sides == "right")
gap <- config$axis_gap
ax_x <- numeric(length(ser))
ax_x[left_i] <- 0.5 - (seq_along(left_i) - 1) * gap
ax_x[right_i] <- n + 0.5 + (seq_along(right_i) - 1) * gap

bar_i <- which(vapply(ser, `[[`, "", "type") == "bar")
nb <- length(bar_i)
bw <- config$bar_width / max(nb, 1)
bars <- do.call(rbind, lapply(seq_along(bar_i), function(j) {
    s <- ser[[bar_i[[j]]]]
    off <- (j - (nb + 1) / 2) * bw
    data.frame(series = s$title, xmin = seq_len(n) + off - bw * 0.45,
        xmax = seq_len(n) + off + bw * 0.45,
        ymin = to_p(pmax(0, s$lim[[1]]), s$lim), ymax = to_p(s$v, s$lim))
}))
line_i <- which(vapply(ser, `[[`, "", "type") == "line")
lines <- do.call(rbind, lapply(line_i, function(i) {
    s <- ser[[i]]
    data.frame(series = s$title, x = seq_len(n), y = to_p(s$v, s$lim))
}))

pr <- diff(P)
tick_len <- 0.1
axes <- do.call(rbind, lapply(seq_along(ser), function(i) {
    s <- ser[[i]]
    br <- seq(s$lim[[1]], s$lim[[2]], length.out = ni + 1)
    data.frame(i = i, x = ax_x[[i]], y = to_p(br, s$lim),
        lab = paste0(format(br, trim = TRUE, drop0trailing = TRUE), s$unit),
        dir = if (s$side == "left") -1 else 1, col = s$col)
}))
fills <- stats::setNames(vapply(ser[bar_i], `[[`, "", "col"),
    vapply(ser[bar_i], `[[`, "", "title"))
lcols <- stats::setNames(vapply(ser[line_i], `[[`, "", "col"),
    vapply(ser[line_i], `[[`, "", "title"))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot()
if (isTRUE(config$show_grid)) {
    gy <- seq(P[[1]], P[[2]], length.out = ni + 1)[-1]
    p <- p + annotate("segment", x = 0.5, xend = n + 0.5, y = gy, yend = gy,
        colour = "grey88", linewidth = 0.35)
}
if (!is.null(bars)) {
    p <- p + geom_rect(data = bars, aes(xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax, fill = series)) +
        scale_fill_manual(values = fills, breaks = names(fills), name = NULL)
}
if (!is.null(lines)) {
    p <- p + geom_line(data = lines, aes(x = x, y = y, colour = series),
        linewidth = 0.8) +
        geom_point(data = lines, aes(x = x, y = y, colour = series),
            shape = 21, fill = "white", size = 2.3, stroke = 1) +
        scale_colour_manual(values = lcols, breaks = names(lcols),
            name = NULL)
}
p <- p +
    annotate("segment", x = 0.5, xend = n + 0.5, y = P[[1]], yend = P[[1]],
        colour = "grey20", linewidth = 0.4) +
    annotate("segment", x = seq_len(n + 1) - 0.5, xend = seq_len(n + 1) - 0.5,
        y = P[[1]], yend = P[[1]] - pr * 0.012, colour = "grey20",
        linewidth = 0.4)
for (i in seq_along(ser)) {
    a <- axes[axes$i == i, , drop = FALSE]
    s <- ser[[i]]
    d <- a$dir[[1]]
    p <- p +
        annotate("segment", x = a$x[[1]], xend = a$x[[1]], y = P[[1]],
            yend = P[[2]], colour = s$col, linewidth = 0.6) +
        annotate("segment", x = a$x, xend = a$x + d * tick_len, y = a$y,
            yend = a$y, colour = s$col, linewidth = 0.5) +
        annotate("text", x = a$x + d * (tick_len + 0.06), y = a$y,
            label = a$lab, hjust = if (d < 0) 1 else 0, size = 3.4,
            colour = "grey20") +
        annotate("text", x = a$x[[1]], y = P[[2]] + pr * 0.045,
            label = s$title, hjust = 0.5, vjust = 0,
            size = 3.8, colour = s$col)
}
xlim <- c(min(ax_x) - gap * 0.75, max(ax_x) + gap * 0.75)
p <- p +
    scale_x_continuous(breaks = seq_len(n), labels = cats, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = xlim, ylim = c(P[[1]], P[[2]] + pr * 0.1),
        clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "grey20", size = 10),
        legend.position = config$legend_position,
        legend.text = element_text(size = 10),
        legend.key.width = unit(9, "mm"),
        plot.background = element_blank(),
        panel.background = element_blank()
    ) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2))

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width,
    height = config$size$height)
