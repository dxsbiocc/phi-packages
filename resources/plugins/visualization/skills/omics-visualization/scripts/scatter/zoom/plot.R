#!/usr/bin/env Rscript

# Template-ID: scatter-zoom
#
# Purpose:
#   Overview scatter plus magnified panel(s) of a supplied x and/or y
#   window, joined by shaded trapezoid connectors (facet_zoom layout).
#   Optional rounded convex hulls with leader-line labels per group.
#
# Inputs:
#   A long table (tsv/csv). Default example (iris):
#     - Petal.Length: x
#     - Petal.Width: y
#     - Species: group (colour, hull)
#
# Output:
#   A PDF, PNG, or SVG overview + zoom figure.
#
# Dependencies:
#   ggplot2, readr, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG: columns, zoom_x / zoom_y (one or both), split,
#   mark ("none" or "hull"), zoom colour, palette, zoom_size.
#   Layouts: x only -> zoom below; y only -> zoom left; both with
#   split = FALSE -> one zoom panel below; both with split = TRUE ->
#   2 x 2 grid (overview, x zoom, y zoom, x+y zoom).
#
# Scientific assumptions:
#   Groups and the zoom window are supplied; the script does not
#   cluster or choose the window. Zoom panels keep the original axis
#   units. Hulls outline supplied groups and are display only (a
#   rounded convex hull, not a density or confidence region).
#   The layout is drawn directly instead of ggforce::facet_zoom, which
#   fails with ggplot2 >= 4 when ggforce < 0.5.

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
        x = "Petal.Length",
        y = "Petal.Width",
        group = "Species"
    ),
    labels = list(
        x = "Petal.Length",
        y = "Petal.Width",
        group = "Species"
    ),
    zoom_x = c(4, 5),           # NULL = no x zoom
    zoom_y = c(1, 2),           # NULL = no y zoom
    split = TRUE,               # both zooms: 2 x 2 grid instead of one panel
    zoom_size = 2,              # zoom panel size relative to the overview
    connector_size = 0.22,
    zoom_fill = "#F5A623",
    zoom_alpha = 0.45,
    mark = "none",              # "none" or "hull"
    hull_expand = 0.025,        # hull padding as a fraction of the panel
    hull_alpha = 0.3,
    hull_labels = TRUE,
    palettes = list(group = "Qualitative.Safe"),
    group_colours = NULL,
    point_size = 1.8,
    pad = 0.05
)

load_packages(c("ggplot2", "readr", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
cc <- config$columns
has_group <- !is.null(cc$group) && cc$group %in% names(df)
require_columns(df, cc[c("x", "y")])
df$.x <- as.numeric(df[[cc$x]])
df$.y <- as.numeric(df[[cc$y]])
df$.g <- if (has_group) as.character(df[[cc$group]]) else "all"
df <- df[is.finite(df$.x) & is.finite(df$.y) & !is.na(df$.g), , drop = FALSE]
if (!nrow(df)) {
    stop("No finite x-y rows to plot.", call. = FALSE)
}
df$.g <- factor(df$.g, levels = unique(df$.g))

zx <- config$zoom_x
zy <- config$zoom_y
if (is.null(zx) && is.null(zy)) {
    stop("Set zoom_x and/or zoom_y in CONFIG.", call. = FALSE)
}
check_win <- function(w, nm) {
    if (!is.null(w) && (length(w) != 2L || !all(is.finite(w)) || w[[1]] >= w[[2]])) {
        stop(nm, " must be c(low, high) with low < high.", call. = FALSE)
    }
}
check_win(zx, "zoom_x")
check_win(zy, "zoom_y")

g_levels <- levels(df$.g)
gcols <- if (!is.null(config$group_colours)) {
    config$group_colours[g_levels]
} else {
    palette_colors(config$palettes$group, n = length(g_levels))
}
names(gcols) <- g_levels

padded <- function(r) r + c(-1, 1) * diff(r) * config$pad
full_x <- padded(range(df$.x))
full_y <- padded(range(df$.y))
win <- function(w, full) if (is.null(w)) full else w + c(-1, 1) * diff(w) * 0.02

# Rounded convex hull in display units, mapped back to data units.
hull_poly <- function(d, xr, yr) {
    if (nrow(d) < 1) return(NULL)
    u <- (d$.x - xr[[1]]) / diff(xr)
    v <- (d$.y - yr[[1]]) / diff(yr)
    h <- if (length(u) >= 3) grDevices::chull(u, v) else seq_along(u)
    a <- seq(0, 2 * pi, length.out = 25)[-25]
    e <- config$hull_expand
    uu <- as.vector(outer(u[h], e * cos(a), `+`))
    vv <- as.vector(outer(v[h], e * sin(a), `+`))
    h2 <- grDevices::chull(uu, vv)
    data.frame(x = xr[[1]] + uu[h2] * diff(xr), y = yr[[1]] + vv[h2] * diff(yr))
}

panel <- function(xr, yr, bands = list(), show_x = TRUE, show_y = TRUE) {
    p <- ggplot2::ggplot()
    for (b in bands) {
        p <- p + ggplot2::annotate(
            "rect", xmin = b$xmin, xmax = b$xmax, ymin = b$ymin, ymax = b$ymax,
            fill = config$zoom_fill, alpha = config$zoom_alpha, colour = NA
        )
    }
    if (identical(config$mark, "hull")) {
        for (g in g_levels) {
            hp <- hull_poly(df[df$.g == g, , drop = FALSE], xr, yr)
            if (is.null(hp)) next
            p <- p + ggplot2::geom_polygon(
                data = hp, ggplot2::aes(x = x, y = y),
                fill = gcols[[g]], colour = gcols[[g]],
                alpha = config$hull_alpha, linewidth = 0.5
            )
        }
    }
    p <- p +
        ggplot2::geom_point(
            data = df, ggplot2::aes(x = .x, y = .y, colour = .g),
            size = config$point_size
        )
    if (identical(config$mark, "hull") && isTRUE(config$hull_labels)) {
        # Label each group whose points are visible, above its topmost visible point.
        lab <- do.call(rbind, lapply(g_levels, function(g) {
            d <- df[df$.g == g & df$.x >= xr[[1]] & df$.x <= xr[[2]] &
                df$.y >= yr[[1]] & df$.y <= yr[[2]], , drop = FALSE]
            if (!nrow(d)) return(NULL)
            top <- d[which.max(d$.y), ]
            up <- top$.y + 0.16 * diff(yr) <= yr[[2]]
            dy <- if (up) 0.16 else -0.16
            data.frame(
                g = g, x0 = top$.x, y0 = top$.y + sign(dy) * config$hull_expand * diff(yr),
                x1 = top$.x - 0.04 * diff(xr), y1 = top$.y + dy * diff(yr),
                stringsAsFactors = FALSE
            )
        }))
        if (!is.null(lab) && nrow(lab)) {
            lab$x2 <- lab$x1 - 0.10 * diff(xr)
            p <- p +
                ggplot2::geom_segment(data = lab,
                    ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1), linewidth = 0.4) +
                ggplot2::geom_segment(data = lab,
                    ggplot2::aes(x = x1, y = y1, xend = x2, yend = y1), linewidth = 0.4) +
                ggplot2::geom_text(data = lab,
                    ggplot2::aes(x = (x1 + x2) / 2, y = y1, label = g),
                    vjust = -0.4, fontface = "bold", size = 3.6)
        }
    }
    p +
        ggplot2::scale_colour_manual(name = config$labels$group, values = gcols) +
        ggplot2::coord_cartesian(xlim = xr, ylim = yr, expand = FALSE, clip = "on") +
        ggplot2::labs(
            x = if (show_x) config$labels$x else NULL,
            y = if (show_y) config$labels$y else NULL
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.background = ggplot2::element_rect(fill = "white", colour = NA),
            plot.margin = ggplot2::margin(2, 2, 2, 2),
            legend.position = if (has_group) "right" else "none"
        )
}

# Connector trapezoid. vertical: source above (span sr in source display
# range sd), target below; horizontal: source on the right, target left.
connector <- function(win, sd, td, vertical = TRUE) {
    f <- function(v, d) (v - d[[1]]) / diff(d)
    s <- f(win, sd)
    t <- f(win, td)
    poly <- if (vertical) {
        data.frame(x = c(s[[1]], s[[2]], t[[2]], t[[1]]), y = c(1, 1, 0, 0))
    } else {
        data.frame(x = c(1, 1, 0, 0), y = c(s[[1]], s[[2]], t[[2]], t[[1]]))
    }
    ggplot2::ggplot(poly, ggplot2::aes(x, y)) +
        ggplot2::geom_polygon(fill = config$zoom_fill, alpha = config$zoom_alpha) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 2, 0, 2))
}

band_x <- function(yr) list(xmin = zx[[1]], xmax = zx[[2]], ymin = yr[[1]], ymax = yr[[2]])
band_y <- function(xr) list(xmin = xr[[1]], xmax = xr[[2]], ymin = zy[[1]], ymax = zy[[2]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
zs <- config$zoom_size
cs <- config$connector_size
both <- !is.null(zx) && !is.null(zy)

if (both && isTRUE(config$split)) {
    wx <- win(zx, full_x)
    wy <- win(zy, full_y)
    p_o <- panel(full_x, full_y, list(band_x(full_y), band_y(full_x)), show_x = FALSE, show_y = FALSE)
    p_zy <- panel(full_x, wy, list(band_x(wy)), show_x = FALSE)
    p_zx <- panel(wx, full_y, list(band_y(wx)), show_y = FALSE)
    p_zxy <- panel(wx, wy)
    c_top <- connector(zy, full_y, wy, vertical = FALSE)
    c_left <- connector(zx, full_x, wx, vertical = TRUE)
    c_right <- connector(zx, full_x, wx, vertical = TRUE)
    c_bot <- connector(zy, full_y, wy, vertical = FALSE)
    design <- "
ABC
DEF
GHI
"
    p <- patchwork::wrap_plots(
        A = p_zy, B = c_top, C = p_o,
        D = c_left, E = patchwork::plot_spacer(), F = c_right,
        G = p_zxy, H = c_bot, I = p_zx,
        design = design
    ) + patchwork::plot_layout(widths = c(zs, cs, 1), heights = c(1, cs, zs),
        guides = "collect")
    w_in <- 9.5
    h_in <- 7.5
} else if (!is.null(zx)) {
    wx <- win(zx, full_x)
    wy <- if (both) win(zy, full_y) else full_y
    bands_o <- list(band_x(full_y))
    if (both) bands_o <- c(bands_o, list(band_y(full_x)))
    p_o <- panel(full_x, full_y, bands_o, show_x = FALSE)
    p_z <- panel(wx, wy)
    p <- patchwork::wrap_plots(p_o, connector(zx, full_x, wx, TRUE), p_z, ncol = 1) +
        patchwork::plot_layout(heights = c(1, cs, zs), guides = "collect")
    w_in <- 8.5
    h_in <- 7.5
} else {
    wy <- win(zy, full_y)
    p_o <- panel(full_x, full_y, list(band_y(full_x)), show_y = FALSE)
    p_z <- panel(full_x, wy)
    p <- patchwork::wrap_plots(p_z, connector(zy, full_y, wy, FALSE), p_o, nrow = 1) +
        patchwork::plot_layout(widths = c(zs, cs, 1), guides = "collect")
    w_in <- 10
    h_in <- 5.5
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = w_in, height = h_in)
