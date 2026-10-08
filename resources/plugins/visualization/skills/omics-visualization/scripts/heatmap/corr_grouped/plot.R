#!/usr/bin/env Rscript

# Template-ID: heatmap-corr-grouped
#
# Purpose:
#   Draw a grouped lower-triangle correlation heatmap: variables are
#   ordered by a supplied group map, blocks are separated by gaps,
#   cell fill is a supplied rho, and stars mark supplied p cuts.
#   Group braces sit on the diagonal hypotenuse; the correlation
#   colorbar is slanted in the empty upper-right wedge.
#
# Inputs:
#   1) Long correlation table (one row per pair):
#        x, y, rho, p
#   2) Variable-to-group map:
#        var, group
#
# Output:
#   A PDF, PNG, or SVG grouped correlation heatmap.
#
# Dependencies:
#   ggplot2, readr, ggprism, scales, grid
#
# Example:
#   Rscript plot.R example.tsv groups.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, p_stars, gap, palette, labels).
#   Upper vs lower display is CONFIG — not a new id. Do not compute
#   Spearman here; supply rho and p. An ungrouped square matrix is
#   heatmap-signif; a 45° rotated triangle is heatmap-corr-rotate;
#   a two-set X-vs-Y grid from sample rows is heatmap-two.
#
# Scientific assumptions:
#   rho and p are supplied. Stars are display cuts of p, not an
#   FDR adjustment performed in this script.

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

io <- parse_io_args(
    n_input = 2L,
    input_names = c("corr", "groups")
)

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        x = "x",
        y = "y",
        value = "rho",
        p = "p",
        var = "var",
        group = "group"
    ),
    p_stars = c(
        "***" = 1e-3,
        "**" = 1e-2,
        "*" = 5e-2
    ),
    triangle = "lower",
    gap = 0.85,
    palettes = list(
        # High rho → purple/blue; low → red (reference Spectral look).
        fill = "Diverging.Spectral / 11"
    ),
    fill_limits = c(-1, 1),
    labels = list(
        title = "",
        fill = "Correlation",
        x_size = 9,
        y_size = 9
    ),
    size = list(
        width = 9.2,
        height = 8.4
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "scales"))

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
# Curly brace along the main diagonal, offset into the empty upper wedge.
# Coordinates are in (i, i) data space before y-reverse.
diag_brace <- function(i0, i1, offset = 0.85, n = 48L) {
    # Outward into iy < ix (empty upper triangle): (+offset, -offset).
    nrm <- c(1, -1) / sqrt(2)
    t <- seq(0, 1, length.out = n)
    # Depth profile of a brace: deep at ends & centre flare.
    depth <- 0.22 * sin(pi * t)^0.65 + 0.55 * pmax(0, 1 - abs(2 * t - 1) * 2.2)^2
    depth[t < 0.04 | t > 0.96] <- 0.05
    base <- i0 + t * (i1 - i0)
    ox <- offset * nrm[[1]]
    oy <- offset * nrm[[2]]
    # Flare perpendicular to diagonal into the empty wedge.
    fx <- depth * nrm[[1]]
    fy <- depth * nrm[[2]]
    data.frame(
        x = base + ox + fx,
        y = base + oy + fy
    )
}

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
corr <- read_table_auto(io$input[[1]])
grp <- read_table_auto(io$input[[2]])

cx <- config$columns$x
cy <- config$columns$y
cv <- config$columns$value
cp <- config$columns$p
gv <- config$columns$var
gg <- config$columns$group

require_columns(corr, list(x = cx, y = cy, value = cv, p = cp))
require_columns(grp, list(var = gv, group = gg))

corr[[cx]] <- as.character(corr[[cx]])
corr[[cy]] <- as.character(corr[[cy]])
corr[[cv]] <- as.numeric(corr[[cv]])
corr[[cp]] <- as.numeric(corr[[cp]])
grp[[gv]] <- as.character(grp[[gv]])
grp[[gg]] <- as.character(grp[[gg]])

if (anyDuplicated(grp[[gv]])) {
    stop("groups table must have unique var rows.", call. = FALSE)
}

# Preserve file order of groups / vars.
grp[[gg]] <- factor(grp[[gg]], levels = unique(grp[[gg]]))
grp <- grp[order(grp[[gg]]), , drop = FALSE]
vars <- grp[[gv]]
g_levels <- levels(grp[[gg]])

pos <- data.frame(
    var = vars,
    group = as.character(grp[[gg]]),
    stringsAsFactors = FALSE
)
gap <- as.numeric(config$gap)
if (!is.finite(gap) || gap < 0) gap <- 0.85
cursor <- 1
pos$i <- NA_real_
for (g in g_levels) {
    hit <- which(pos$group == g)
    pos$i[hit] <- cursor + seq_along(hit) - 1
    cursor <- cursor + length(hit) + gap
}

map_i <- setNames(pos$i, pos$var)
map_g <- setNames(pos$group, pos$var)

ok <- corr[[cx]] %in% vars & corr[[cy]] %in% vars &
    is.finite(corr[[cv]]) & is.finite(corr[[cp]]) & corr[[cp]] >= 0
if (any(!ok)) {
    message("Dropped ", sum(!ok), " corr row(s) outside groups or invalid.")
    corr <- corr[ok, , drop = FALSE]
}
if (!nrow(corr)) stop("No correlation rows left to plot.", call. = FALSE)

corr$ix <- unname(map_i[corr[[cx]]])
corr$iy <- unname(map_i[corr[[cy]]])
corr$gx <- unname(map_g[corr[[cx]]])
corr$gy <- unname(map_g[corr[[cy]]])

tri <- tolower(as.character(config$triangle[[1]]))
if (!tri %in% c("lower", "upper", "full")) {
    stop("config$triangle must be lower, upper, or full.", call. = FALSE)
}
# Overall matrix triangle (not just within-group): matches the reference
# lower-left filled / upper-right empty layout.
if (identical(tri, "lower")) {
    corr <- corr[corr$iy >= corr$ix, , drop = FALSE]
} else if (identical(tri, "upper")) {
    corr <- corr[corr$iy <= corr$ix, , drop = FALSE]
}

star_cuts <- sort(unlist(config$p_stars), decreasing = TRUE)
star_lab <- function(p) {
    out <- character(length(p))
    for (i in seq_along(star_cuts)) {
        hit <- is.finite(p) & p < star_cuts[[i]]
        out[hit] <- names(star_cuts)[[i]]
    }
    out
}
corr$star <- star_lab(corr[[cp]])
corr$star[corr[[cx]] == corr[[cy]]] <- ""

fill_lim <- as.numeric(config$fill_limits)
if (length(fill_lim) != 2L || any(!is.finite(fill_lim))) {
    fill_lim <- c(-1, 1)
}
# Spectral / 11 resolves red→…→blue; high rho should be blue/purple.
fill_cols <- palette_colors(config$palettes$fill)

# Diagonal braces + labels (one per group).
brace_paths <- list()
brace_labs <- list()
for (gi in seq_along(g_levels)) {
    g <- g_levels[[gi]]
    hit <- pos$i[pos$group == g]
    i0 <- min(hit)
    i1 <- max(hit)
    bp <- diag_brace(i0 - 0.15, i1 + 0.15, offset = 0.95)
    bp$id <- g
    brace_paths[[gi]] <- bp
    mid <- (i0 + i1) / 2
    nrm <- c(1, -1) / sqrt(2)
    brace_labs[[gi]] <- data.frame(
        group = g,
        x = mid + 1.85 * nrm[[1]],
        y = mid + 1.85 * nrm[[2]],
        stringsAsFactors = FALSE
    )
}
brace_path_df <- do.call(rbind, brace_paths)
brace_lab_df <- do.call(rbind, brace_labs)

# Slanted colorbar in the empty upper-right wedge, parallel to the
# matrix diagonal (screen: top-left → bottom-right under y-reverse).
i_min <- min(pos$i)
i_max <- max(pos$i)
span <- i_max - i_min
bar_dir <- c(1, 1) / sqrt(2) # along diagonal
bar_n <- c(1, -1) / sqrt(2) # into empty wedge (screen up-right)
diag_mid <- (i_min + i_max) / 2
bar_mid <- c(diag_mid, diag_mid) + 3.6 * bar_n
bar_half_len <- span * 0.38
bar_half_w <- 0.24
n_seg <- 100L
tt <- seq(-1, 1, length.out = n_seg + 1L)
t_mid <- (tt[-length(tt)] + tt[-1]) / 2
# t = −1 (upper-left) → rho = +1; t = +1 (lower-right) → rho = −1.
rho_mid <- fill_lim[[2]] - (t_mid + 1) / 2 * (fill_lim[[2]] - fill_lim[[1]])
pal_fun <- scales::col_numeric(
    palette = fill_cols,
    domain = fill_lim,
    na.color = "grey90"
)
cols <- pal_fun(rho_mid)
cbar_poly <- do.call(rbind, lapply(seq_len(n_seg), function(k) {
    a <- bar_mid + tt[[k]] * bar_half_len * bar_dir
    b <- bar_mid + tt[[k + 1L]] * bar_half_len * bar_dir
    corners <- rbind(
        a + bar_half_w * bar_n,
        b + bar_half_w * bar_n,
        b - bar_half_w * bar_n,
        a - bar_half_w * bar_n
    )
    data.frame(
        id = k,
        fill = cols[[k]],
        x = corners[, 1],
        y = corners[, 2],
        stringsAsFactors = FALSE
    )
}))
brks <- seq(fill_lim[[1]], fill_lim[[2]], by = 0.2)
tick_t <- 1 - 2 * (brks - fill_lim[[1]]) /
    max(fill_lim[[2]] - fill_lim[[1]], 1e-8)
# Ticks then title, both on the outer (+bar_n) side of the bar.
tick_off <- bar_half_w + 0.48
cbar_ticks <- data.frame(
    lab = format(brks, trim = TRUE),
    x = bar_mid[[1]] + tick_t * bar_half_len * bar_dir[[1]] +
        tick_off * bar_n[[1]],
    y = bar_mid[[2]] + tick_t * bar_half_len * bar_dir[[2]] +
        tick_off * bar_n[[2]],
    stringsAsFactors = FALSE
)
# Title further out than the tick lane, parallel to the bar.
title_off <- bar_half_w + 1.55
cbar_title <- data.frame(
    lab = config$labels$fill,
    x = bar_mid[[1]] + title_off * bar_n[[1]],
    y = bar_mid[[2]] + title_off * bar_n[[2]],
    stringsAsFactors = FALSE
)
# Screen angle: y-reverse flips the visual slope of bar_dir.
bar_ang <- atan2(-bar_dir[[2]], bar_dir[[1]]) * 180 / pi
if (bar_ang > 90) bar_ang <- bar_ang - 180
if (bar_ang < -90) bar_ang <- bar_ang + 180

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
# Star colour: light on dark fills, dark otherwise (approx by |rho|).
corr$star_col <- ifelse(abs(corr[[cv]]) > 0.55, "white", "grey15")

p <- ggplot(corr, aes(x = ix, y = iy)) +
    geom_tile(
        aes(fill = .data[[cv]]),
        width = 1,
        height = 1,
        colour = "white",
        linewidth = 0.25
    ) +
    geom_text(
        aes(label = star, colour = I(star_col)),
        size = 3.0,
        na.rm = TRUE,
        show.legend = FALSE
    ) +
    scale_fill_gradientn(
        name = config$labels$fill,
        colours = fill_cols,
        limits = fill_lim,
        oob = scales::squish,
        guide = "none"
    ) +
    scale_x_continuous(
        breaks = pos$i,
        labels = pos$var,
        # Lock panel to the tile grid — colorbar/braces must not train
        # the scale or axis ticks float away from the cells.
        limits = c(min(pos$i) - 0.5, max(pos$i) + 0.5),
        expand = expansion(mult = 0, add = 0),
        oob = scales::oob_keep,
        position = "bottom"
    ) +
    scale_y_reverse(
        breaks = pos$i,
        labels = pos$var,
        limits = c(min(pos$i) - 0.5, max(pos$i) + 0.5),
        expand = expansion(mult = 0, add = 0),
        oob = scales::oob_keep
    ) +
    # Keep square cells via theme aspect.ratio — not coord_fixed.
    # coord_fixed dumps leftover device space into axis 1null units and
    # floats ticks/labels away from the tiles.
    coord_cartesian(clip = "off") +
    geom_path(
        data = brace_path_df,
        aes(x = x, y = y, group = id),
        inherit.aes = FALSE,
        colour = "grey25",
        linewidth = 0.7,
        lineend = "round",
        linejoin = "round"
    ) +
    geom_text(
        data = brace_lab_df,
        aes(x = x, y = y, label = group),
        inherit.aes = FALSE,
        angle = -45,
        hjust = 0,
        vjust = 0.5,
        size = 3.5,
        fontface = "bold",
        colour = "grey15"
    ) +
    geom_polygon(
        data = cbar_poly,
        aes(x = x, y = y, group = id),
        inherit.aes = FALSE,
        fill = cbar_poly$fill,
        colour = NA
    ) +
    geom_text(
        data = cbar_ticks,
        aes(x = x, y = y, label = lab),
        inherit.aes = FALSE,
        angle = bar_ang,
        hjust = 0,
        vjust = 0.5,
        size = 2.7,
        colour = "grey10"
    ) +
    geom_text(
        data = cbar_title,
        aes(x = x, y = y, label = lab),
        inherit.aes = FALSE,
        angle = bar_ang,
        hjust = 0.5,
        vjust = 0.5,
        size = 3.3,
        fontface = "bold",
        colour = "grey10"
    ) +
    labs(title = config$labels$title, x = NULL, y = NULL) +
    theme_prism(base_size = 12) +
    theme(
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "grey20", linewidth = 0.55),
        axis.ticks.length = unit(0.12, "cm"),
        axis.text.x = element_text(
            angle = 90, hjust = 1, vjust = 0.5,
            size = config$labels$x_size,
            margin = margin(t = 1)
        ),
        axis.text.y = element_text(
            size = config$labels$y_size,
            margin = margin(r = 1),
            hjust = 1
        ),
        aspect.ratio = 1,
        # Top/right margin hosts the slanted colorbar (clip = "off").
        plot.margin = margin(36, 80, 8, 8)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
# Two layout traps push ticks away from tiles:
# 1) axis columns/rows embed 1null → leftover device space inflates them
# 2) aspect.ratio letterboxes inside a non-square panel cell
# Resolve both: pinch axes to absolute size, lock panel to a square.
gt <- ggplotGrob(p)
panel_i <- which(gt$layout$name == "panel")[[1]]
panel_col <- gt$layout$l[[panel_i]]
panel_row <- gt$layout$t[[panel_i]]
panel_in <- 6.4

grDevices::pdf(nullfile())
for (ax in c("axis-l", "axis-b")) {
    ii <- which(gt$layout$name == ax)
    if (!length(ii)) next
    if (identical(ax, "axis-l")) {
        col <- gt$layout$l[[ii]]
        gt$widths[col] <- unit(
            grid::convertWidth(gt$widths[col], "cm", valueOnly = TRUE),
            "cm"
        )
    } else {
        row <- gt$layout$t[[ii]]
        gt$heights[row] <- unit(
            grid::convertHeight(gt$heights[row], "cm", valueOnly = TRUE),
            "cm"
        )
    }
}
invisible(grDevices::dev.off())

gt$widths[panel_col] <- unit(panel_in, "in")
gt$heights[panel_row] <- unit(panel_in, "in")

# Figure size = panel + axes + plot.margin (top/right host the colorbar).
grDevices::pdf(nullfile())
fig_w <- grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE)
fig_h <- grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)
invisible(grDevices::dev.off())

save_ggplot(
    gt, io$output,
    width = fig_w,
    height = fig_h
)
