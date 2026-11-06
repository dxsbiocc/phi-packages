#!/usr/bin/env Rscript

# Template-ID: bar-radial-line
#
# Purpose:
#   Circular bar–line composite: one equal-angle slot per item, radial
#   bars of one supplied measure in the inner band, a polar line with
#   points and value labels of a second supplied measure in the outer
#   band, group-tinted background sectors, and an outer ring of
#   per-item label tiles.
#
# Inputs:
#   One row per item. Default example:
#     - species: item label
#     - type: group (contiguous run of slots)
#     - number_of_sites: non-negative bar value
#     - threshold_kpa: line value; empty = no value (the line breaks)
#
# Output:
#   A PDF, PNG, or SVG circular bar–line composite.
#
# Dependencies:
#   ggplot2, readr, grid
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, group order, palettes, radii, scale
#   maxima, neutral groups, axis titles, label sizes). Extra items or
#   groups are extra rows, not extra ids. Groups whose line values are
#   all missing get an arc on the bar band edge; neutral groups get a
#   dot instead. Do not fit thresholds, breakpoints, or regressions in
#   this script.
#
# Scientific assumptions:
#   Both measures are supplied. Bars share one radial scale (0 to
#   bar_max); the line shares one scale (0 to line_max). Missing line
#   values are not imputed; the line is broken at them and at group
#   boundaries.

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
        item = "species",
        group = "type",
        bar = "number_of_sites",
        line = "threshold_kpa"
    ),
    # NULL keeps first-appearance order
    groups = NULL,
    # groups drawn grey with a dot on the bar band edge (no arc)
    neutral_groups = c("T6"),
    palettes = list(
        group = "Qualitative.Bold",
        neutral = "Qualitative.Safe"
    ),
    neutral_index = 12L,
    # white-mix fractions for light tints derived from group colours
    tint = list(
        sector = 0.88,
        label = 0.80
    ),
    darken_text = 0.18,
    bar_max = 100,               # NULL = 1.05 x max(bar)
    line_max = 1.5,              # NULL = 1.2 x max(line)
    bar_reference = NULL,        # dashed circle inside the bar band; NULL = bar_max / 2
    radius = list(
        center = 0.058,
        inner = 0.073,
        site = 0.550,
        line_outer = 0.825,
        label_inner = 0.835,
        label_outer = 1.000
    ),
    value_label = list(
        radius = 0.775,
        min_gap = 0.06,          # move label inside the line when closer
        inset = 0.075,
        size_pt = 8.5,
        digits = 2
    ),
    item_label_pt = 10,
    axis_titles = c(bar = "Number\nof Sites", line = "Threshold\n(kPa)"),
    axis_title_pt = 9,
    size_mm = 183
)

load_packages(c("ggplot2", "readr", "grid"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
i_col <- config$columns$item
g_col <- config$columns$group
b_col <- config$columns$bar
l_col <- config$columns$line
df[[i_col]] <- as.character(df[[i_col]])
df[[g_col]] <- as.character(df[[g_col]])
df[[b_col]] <- as.numeric(df[[b_col]])
df[[l_col]] <- suppressWarnings(as.numeric(df[[l_col]]))
if (any(!nzchar(df[[i_col]]) | !nzchar(df[[g_col]]))) {
    stop("item and group must be non-empty.", call. = FALSE)
}
if (any(!is.finite(df[[b_col]]) | df[[b_col]] < 0)) {
    stop("bar values must be finite and non-negative.", call. = FALSE)
}
if (any(df[[l_col]] < 0, na.rm = TRUE)) {
    stop("line values must be non-negative.", call. = FALSE)
}

groups <- if (is.null(config$groups)) unique(df[[g_col]]) else
    as.character(config$groups)
df <- df[df[[g_col]] %in% groups, , drop = FALSE]
df <- df[order(match(df[[g_col]], groups), seq_len(nrow(df))), , drop = FALSE]
groups <- groups[groups %in% df[[g_col]]]
n <- nrow(df)
if (n < 3L) stop("Need at least three items.", call. = FALSE)

r <- config$radius
bar_max <- if (is.null(config$bar_max)) 1.05 * max(df[[b_col]]) else
    config$bar_max
line_max <- if (is.null(config$line_max)) {
    1.2 * max(df[[l_col]], na.rm = TRUE)
} else {
    config$line_max
}
bar_ref <- if (is.null(config$bar_reference)) bar_max / 2 else
    config$bar_reference

mix <- function(col, with, p) {
    a <- grDevices::col2rgb(col) / 255
    b <- as.vector(grDevices::col2rgb(with) / 255)
    m <- (1 - p) * a + p * b
    grDevices::rgb(m[1, ], m[2, ], m[3, ])
}
neutral <- intersect(config$neutral_groups, groups)
active <- setdiff(groups, neutral)
grp_col <- setNames(character(length(groups)), groups)
grp_col[active] <- palette_colors(config$palettes$group, n = length(active))
grp_col[neutral] <- palette_colors(config$palettes$neutral)[[config$neutral_index]]
line_col <- setNames(mix(grp_col, "#000000", config$darken_text), groups)
sector_fill <- setNames(mix(grp_col, "#ffffff", config$tint$sector), groups)
label_fill <- setNames(mix(grp_col, "#ffffff", config$tint$label), groups)

df$slot <- seq_len(n) - 1
df$bar_r <- r$inner + pmin(df[[b_col]], bar_max) / bar_max * (r$site - r$inner)
df$line_r <- r$site + pmin(df[[l_col]], line_max) / line_max *
    (r$line_outer - r$site)
math_angle <- function(slot) 180 - slot * 360 / n
upright <- function(a) {
    a <- ((a + 180) %% 360) - 180
    ifelse(a > 90, a - 180, ifelse(a < -90, a + 180, a))
}
df$rotation <- upright(math_angle(df$slot))
runs <- rle(df[[g_col]])
df$run <- rep(seq_along(runs$lengths), runs$lengths)
df$seg <- cumsum(c(TRUE, diff(df$run) != 0 |
    diff(is.na(df[[l_col]])) != 0))

bg <- do.call(rbind, lapply(split(df, df$run), function(d) data.frame(
    group = d[[g_col]][[1]],
    xmin = min(d$slot) - 0.5,
    xmax = max(d$slot) + 0.5,
    no_line = all(is.na(d[[l_col]])),
    stringsAsFactors = FALSE
)))
bg$fill <- unname(sector_fill[bg$group])
bg$col <- unname(line_col[bg$group])
arc <- bg[bg$no_line & !bg$group %in% neutral, , drop = FALSE]
arc$xmin <- arc$xmin + 0.14
arc$xmax <- arc$xmax - 0.14

ln <- df[!is.na(df[[l_col]]), , drop = FALSE]
ln$col <- unname(line_col[ln[[g_col]]])
vl <- config$value_label
ln$text_r <- ifelse(ln$line_r >= vl$radius - vl$min_gap,
    ln$line_r - vl$inset, vl$radius)
ln$label <- sub("\\.$", "", sub("0+$", "",
    formatC(ln[[l_col]], format = "f", digits = vl$digits)))

# Item label size: fit the radial length and the tangential slot width.
fig_mm <- config$size_mm
r_scale <- 0.48 * fig_mm
radial_mm <- (r$label_outer - r$label_inner) * r_scale - 0.6
slot_mm <- 2 * pi * r$label_inner * r_scale / n * 0.93
text_w <- vapply(df[[i_col]], function(s) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(s,
        gp = grid::gpar(fontsize = config$item_label_pt, fontface = "bold"))),
        "mm", valueOnly = TRUE)
}, numeric(1))
pt_mm <- 25.4 / 72
df$label_pt <- pmin(
    config$item_label_pt * pmin(1, radial_mm / text_w),
    slot_mm / pt_mm * 0.9
)
df$text_col <- unname(line_col[df[[g_col]]])
df$tile_fill <- unname(label_fill[df[[g_col]]])
df$bar_fill <- unname(grp_col[df[[g_col]]])

neutral_df <- df[df[[g_col]] %in% neutral, , drop = FALSE]
no_line_runs <- bg[bg$no_line, , drop = FALSE]
line_axis_slot <- if (nrow(neutral_df)) {
    neutral_df$slot[[1]] - 0.5
} else if (nrow(no_line_runs)) {
    no_line_runs$xmin[[1]]
} else {
    -0.5
}
arrows <- data.frame(
    slot = c(-0.5, line_axis_slot),
    y0 = c(r$center, r$site + 0.010),
    y1 = c(r$site - 0.003, r$line_outer - 0.006)
)
axis_lab <- data.frame(
    slot = c(n - 0.026 * n, line_axis_slot + 0.076 * n / (2 * pi)),
    radius = c(0.40, 0.70),
    label = unname(config$axis_titles[c("bar", "line")]),
    rotation = c(0, upright(math_angle(line_axis_slot)))
)
bar_ref_r <- r$inner + bar_ref / bar_max * (r$site - r$inner)
circles <- data.frame(
    radius = c(bar_ref_r, r$site, r$line_outer),
    colour = c("#929292", "#111111", "#8b8b8b"),
    linewidth = c(0.21, 0.63, 0.23),
    linetype = c("dashed", "22", "dashed")
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_rect(
        data = bg,
        aes(xmin = xmin, xmax = xmax, ymin = r$inner, ymax = r$site,
            fill = fill),
        colour = NA
    ) +
    geom_segment(
        data = df,
        aes(x = slot, xend = slot, y = 0, yend = r$line_outer),
        colour = "#b8b8b8", linewidth = 0.16
    ) +
    geom_rect(
        data = df,
        aes(xmin = slot - 0.37, xmax = slot + 0.37, ymin = r$inner,
            ymax = bar_r, fill = bar_fill, colour = text_col),
        linewidth = 0.12
    ) +
    geom_hline(
        data = circles,
        aes(yintercept = radius, colour = colour, linewidth = linewidth,
            linetype = linetype)
    ) +
    geom_segment(
        data = arc,
        aes(x = xmin, xend = xmax, y = r$site, yend = r$site, colour = col),
        linewidth = 0.9
    ) +
    geom_line(
        data = ln,
        aes(x = slot, y = line_r, group = seg, colour = col),
        linewidth = 0.85, linejoin = "round"
    ) +
    geom_point(
        data = ln,
        aes(x = slot, y = line_r, colour = col),
        shape = 21, fill = "white", size = 1.2, stroke = 0.7
    ) +
    geom_text(
        data = ln,
        aes(x = slot, y = text_r, label = label, angle = rotation,
            colour = col),
        fontface = "bold", size = vl$size_pt, size.unit = "pt"
    ) +
    geom_rect(
        data = df,
        aes(xmin = slot - 0.465, xmax = slot + 0.465, ymin = r$label_inner,
            ymax = r$label_outer, fill = tile_fill),
        colour = NA
    ) +
    geom_text(
        data = df,
        aes(x = slot, y = (r$label_inner + r$label_outer) / 2,
            label = .data[[i_col]], angle = rotation, colour = text_col,
            size = label_pt),
        fontface = "bold", size.unit = "pt"
    ) +
    geom_hline(yintercept = r$center, colour = "#242424", linewidth = 0.37) +
    geom_point(
        data = neutral_df, aes(x = slot, y = r$site),
        shape = 16, colour = "#969696", size = 1.6
    ) +
    geom_segment(
        data = arrows, aes(x = slot, xend = slot, y = y0, yend = y1),
        arrow = grid::arrow(type = "closed", length = grid::unit(1.8, "mm")),
        arrow.fill = "#111111", colour = "#111111", linewidth = 0.44
    ) +
    geom_text(
        data = axis_lab,
        aes(x = slot, y = radius, label = label, angle = rotation),
        size = config$axis_title_pt, size.unit = "pt", lineheight = 0.95
    ) +
    scale_x_continuous(limits = c(-0.5, n - 0.5), expand = expansion(0)) +
    scale_y_continuous(limits = c(0, r$label_outer), expand = expansion(0)) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_size_identity() +
    scale_linewidth_identity() +
    scale_linetype_identity() +
    coord_polar(theta = "x", start = -pi / 2 - pi / n, direction = 1,
        clip = "off") +
    theme_void(base_size = 9) +
    theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(-0.10 * fig_mm, -0.10 * fig_mm,
            -0.10 * fig_mm, -0.10 * fig_mm, unit = "mm"),
        legend.position = "none"
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = fig_mm / 25.4, height = fig_mm / 25.4)
