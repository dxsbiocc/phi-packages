#!/usr/bin/env Rscript

# Template-ID: scatter-density2d
#
# Purpose:
#   Two-dimensional kernel density of supplied x-y points: filled
#   density bands, contour lines, both, a smooth raster, or labeled
#   contours with points coloured by local density, with optional
#   group-coloured contours, per-group facets, and top / right
#   marginal density or histogram panels.
#
# Inputs:
#   One row per observation. Default example:
#     - x, y: numeric coordinates
#     - group: optional label (colour or facet; ignored unless used)
#
# Output:
#   A PDF, PNG, or SVG 2D density plot.
#
# Dependencies:
#   ggplot2, readr, patchwork, MASS
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, mode, bins, bandwidth adjust,
#   contour_var, colour_by_group, facet_by_group, show_points,
#   marginal, palettes). Group colouring, facets, raster, labeled
#   contours and marginal panels are variants, not new ids.
#   Hexagon / square counts are scatter-hex; group hull outlines
#   are scatter-contour.
#
# Scientific assumptions:
#   Points are supplied. The 2D kernel density (MASS::kde2d via
#   ggplot2) is a display smoothing, like a violin's density; it is
#   not a test or a cluster call. Band levels are display bins.

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
    columns = list(x = "x", y = "y", group = "group"),
    # "labeled": points coloured by local density + black contours
    # with density value labels (faithful-style).
    # "points": density evaluated on a coarse grid as sized points.
    mode = "both",               # "filled", "lines", "both", "raster", "labeled", "points"
    bins = 9,
    drop_lowest = TRUE,          # lowest filled band transparent
    adjust = 1,                  # bandwidth multiplier
    contour_var = "density",     # "density", "ndensity", "count"
    colour_by_group = FALSE,     # contour lines coloured by group
    facet_by_group = FALSE,
    show_points = TRUE,
    point_size = 0.7,
    point_alpha = 0.45,
    marginal = "none",           # "none", "density", "histogram"
    fill_palette = "Gradient.YlGnBu",
    group_palette = "Qualitative.Bold",
    labeled = list(
        point_size = 2.2,
        point_alpha = 0.9,
        contour_colour = "grey20",
        contour_linewidth = 0.45,
        label_size = 3.0,
        label_colour = "Qualitative.Bold:1",
        n_grid = 100,
        show_legend = TRUE
    ),
    labels = list(title = "", x = "Eruption (min)", y = "Waiting (min)"),
    size = list(width = 7.5, height = 6.5)
)

load_packages(c("ggplot2", "readr", "patchwork", "MASS"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
cc <- config$columns
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
use_group <- isTRUE(config$colour_by_group) || isTRUE(config$facet_by_group)
need <- if (use_group) cc else cc[c("x", "y")]
require_columns(df, need)
df <- data.frame(x = as.numeric(df[[cc$x]]), y = as.numeric(df[[cc$y]]),
    group = if (use_group) as.character(df[[cc$group]]) else "all",
    stringsAsFactors = FALSE)
df <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
if (nrow(df) < 10L) stop("Need at least 10 finite points.", call. = FALSE)
df$group <- factor(df$group, levels = unique(df$group))

fill_stops <- palette_colors(config$fill_palette)
grp_cols <- stats::setNames(palette_colors(config$group_palette,
    n = nlevels(df$group)), levels(df$group))
ink <- "grey15"

pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

interp_kde <- function(kd, x, y) {
    ix <- findInterval(x, kd$x, all.inside = TRUE)
    iy <- findInterval(y, kd$y, all.inside = TRUE)
    # bilinear on the local cell
    x1 <- kd$x[ix]
    x2 <- kd$x[pmin(ix + 1L, length(kd$x))]
    y1 <- kd$y[iy]
    y2 <- kd$y[pmin(iy + 1L, length(kd$y))]
    z11 <- kd$z[cbind(ix, iy)]
    z21 <- kd$z[cbind(pmin(ix + 1L, length(kd$x)), iy)]
    z12 <- kd$z[cbind(ix, pmin(iy + 1L, length(kd$y)))]
    z22 <- kd$z[cbind(pmin(ix + 1L, length(kd$x)), pmin(iy + 1L, length(kd$y)))]
    tx <- ifelse(x2 > x1, (x - x1) / (x2 - x1), 0)
    ty <- ifelse(y2 > y1, (y - y1) / (y2 - y1), 0)
    (1 - tx) * (1 - ty) * z11 + tx * (1 - ty) * z21 +
        (1 - tx) * ty * z12 + tx * ty * z22
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
mode <- config$mode

if (identical(mode, "labeled")) {
    lc <- config$labeled
    hx <- MASS::bandwidth.nrd(df$x) * config$adjust
    hy <- MASS::bandwidth.nrd(df$y) * config$adjust
    kd <- MASS::kde2d(df$x, df$y, n = lc$n_grid, h = c(hx, hy))
    df$density <- interp_kde(kd, df$x, df$y)
    grid <- expand.grid(x = kd$x, y = kd$y)
    grid$z <- as.vector(kd$z)
    z_rng <- range(kd$z, finite = TRUE)
    levels <- pretty(z_rng, n = config$bins)
    levels <- levels[levels > z_rng[[1]] & levels < z_rng[[2]]]
    if (!length(levels)) {
        levels <- seq(z_rng[[1]], z_rng[[2]], length.out = config$bins + 2L)[-c(1, config$bins + 2L)]
    }
    cl <- grDevices::contourLines(kd$x, kd$y, kd$z, levels = levels)
    lab_df <- NULL
    if (length(cl)) {
        lab_df <- do.call(rbind, lapply(cl, function(seg) {
            n <- length(seg$x)
            i <- max(1L, floor(n * 0.38))
            data.frame(
                x = seg$x[[i]], y = seg$y[[i]],
                label = format(seg$level, digits = 2, trim = TRUE),
                level = seg$level,
                stringsAsFactors = FALSE
            )
        }))
        # One label per level (first segment wins).
        lab_df <- lab_df[!duplicated(lab_df$level), , drop = FALSE]
    }
    p <- ggplot(df, aes(x = x, y = y)) +
        geom_point(
            aes(colour = density),
            size = lc$point_size, alpha = lc$point_alpha
        ) +
        geom_contour(
            data = grid, aes(x = x, y = y, z = z),
            breaks = levels,
            colour = lc$contour_colour,
            linewidth = lc$contour_linewidth
        ) +
        scale_colour_gradientn(
            colours = fill_stops, name = "Density",
            guide = if (isTRUE(lc$show_legend)) "colourbar" else "none"
        ) +
        labs(
            title = config$labels$title,
            x = config$labels$x,
            y = config$labels$y
        ) +
        theme_bw(base_size = 13) +
        theme(
            panel.grid.major = element_line(
                colour = "grey80", linetype = "dashed", linewidth = 0.3
            ),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
            axis.text = element_text(colour = "black"),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            legend.background = element_rect(fill = "white", colour = NA)
        )
    if (!is.null(lab_df) && nrow(lab_df)) {
        p <- p + geom_text(
            data = lab_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = lc$label_size,
            colour = pal_ref(lc$label_colour),
            fontface = "plain"
        )
    }
    if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
        p <- p + theme(plot.title = element_blank())
    }
    save_ggplot(p, io$output, width = config$size$width, height = config$size$height)
    quit(save = "no")
}

if (identical(mode, "points")) {
    # Coarse grid of density values drawn as points (tutorial
    # stat_density_2d(geom = "point", n = 20) style).
    n_grid <- max(8L, as.integer(config$bins) * 2L)
    hx <- MASS::bandwidth.nrd(df$x) * config$adjust
    hy <- MASS::bandwidth.nrd(df$y) * config$adjust
    kd <- MASS::kde2d(df$x, df$y, n = n_grid, h = c(hx, hy))
    grid <- expand.grid(x = kd$x, y = kd$y)
    grid$density <- as.vector(kd$z)
    p <- ggplot(grid, aes(x = x, y = y)) +
        geom_point(
            aes(size = density, colour = density),
            alpha = 0.85
        ) +
        scale_colour_gradientn(colours = fill_stops, name = "Density") +
        scale_size_continuous(name = "Density", range = c(0.4, 6)) +
        labs(
            title = config$labels$title,
            x = config$labels$x,
            y = config$labels$y
        ) +
        theme_bw(base_size = 13) +
        theme(
            panel.grid.major = element_line(
                colour = "grey85", linetype = "dashed", linewidth = 0.3
            ),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            axis.text = element_text(colour = "black"),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA)
        )
    if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
        p <- p + theme(plot.title = element_blank())
    }
    save_ggplot(p, io$output, width = config$size$width, height = config$size$height)
    quit(save = "no")
}

band_cols <- grDevices::colorRampPalette(fill_stops)(config$bins)
if (isTRUE(config$drop_lowest)) band_cols[[1]] <- "transparent"
pad <- if (identical(mode, "raster")) 0 else 0.02
p <- ggplot(df, aes(x = x, y = y))
if (identical(mode, "raster")) {
    p <- p + stat_density_2d(aes(fill = after_stat(density)), geom = "raster",
        contour = FALSE, adjust = config$adjust, n = 200) +
        scale_fill_gradientn(colours = fill_stops, name = "Density")
} else if (mode %in% c("filled", "both")) {
    p <- p + geom_density_2d_filled(bins = config$bins,
        adjust = config$adjust, contour_var = config$contour_var,
        alpha = 0.9) +
        scale_fill_manual(values = band_cols,
            name = if (identical(config$contour_var, "count")) "Count"
            else "Density")
}
if (mode %in% c("lines", "both")) {
    if (isTRUE(config$colour_by_group)) {
        p <- p + geom_density_2d(aes(colour = group), bins = config$bins,
            adjust = config$adjust, contour_var = config$contour_var,
            linewidth = 0.5) +
            scale_colour_manual(values = grp_cols, name = NULL)
    } else {
        p <- p + geom_density_2d(bins = config$bins, adjust = config$adjust,
            contour_var = config$contour_var, colour = ink,
            linewidth = 0.25)
    }
}
if (isTRUE(config$show_points)) {
    if (isTRUE(config$colour_by_group)) {
        p <- p + geom_point(aes(colour = group), size = config$point_size,
            alpha = config$point_alpha)
    } else {
        p <- p + geom_point(size = config$point_size,
            alpha = config$point_alpha, colour = ink)
    }
}
if (isTRUE(config$facet_by_group)) {
    p <- p + facet_wrap(~group)
}
p <- p +
    scale_x_continuous(expand = expansion(mult = pad)) +
    scale_y_continuous(expand = expansion(mult = pad)) +
    labs(title = config$labels$title, x = config$labels$x,
        y = config$labels$y) +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.major = element_line(
            colour = "grey85", linetype = "dashed", linewidth = 0.3
        ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.key.size = unit(4, "mm")
    )

if (config$marginal %in% c("density", "histogram") &&
    !isTRUE(config$facet_by_group)) {
    mfill <- fill_stops[[max(1, length(fill_stops) - 2)]]
    void <- theme_void() + theme(plot.background = element_blank())
    marg <- function(aesx) {
        g <- ggplot(df, aesx)
        if (identical(config$marginal, "density")) {
            g + geom_density(fill = mfill, colour = ink, linewidth = 0.3,
                alpha = 0.85) + void
        } else {
            g + geom_histogram(bins = 35, fill = mfill, colour = "white",
                linewidth = 0.2) + void
        }
    }
    top <- marg(aes(x = x)) +
        scale_x_continuous(limits = range(df$x),
            expand = expansion(mult = 0.02))
    right <- marg(aes(y = y)) +
        scale_y_continuous(limits = range(df$y),
            expand = expansion(mult = 0.02))
    p <- top + patchwork::plot_spacer() + p + right +
        patchwork::plot_layout(ncol = 2, widths = c(4, 1),
            heights = c(1, 4), guides = "collect")
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width,
    height = config$size$height)
