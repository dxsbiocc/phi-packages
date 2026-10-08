#!/usr/bin/env Rscript

# Template-ID: scatter-svg
#
# Purpose:
#   Draw an x-y scatter whose marks are a supplied SVG glyph instead of
#   ggplot point shapes. Fill of one SVG path is mapped to a group.
#
# Inputs:
#   One row per observation. Default example:
#     - x, y: coordinates
#     - type: group used to colour the SVG fill
#     - count: optional label drawn in the pin head
#   Bundled glyph: pin.svg in this template directory.
#
# Output:
#   A PDF, PNG, or SVG scatter of SVG marks.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggsvg
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, svg file, css selector,
#   labels). Point to another SVG in this directory or beside the copied
#   plot.R. Edit DATA PREPARATION to change group order. Edit PLOT only
#   when geometry must change. Discrete vs continuous SVG fill is
#   CONFIG (scale_svg_fill_manual vs scale_svg_fill_gradient), not a
#   second template. Do not add an id that only swaps the SVG file.
#
# Scientific assumptions:
#   Each row is one observation with two coordinates and a supplied group.
#   The SVG is a display glyph, not a statistical encoding beyond fill.
#   Use scatter-group for ordinary pch marks. Use bar-svg-icon when the
#   SVG identifies a bar category rather than an x-y point.

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
        x = "x",
        y = "y",
        fill = "type",
        label = "count"
    ),
    svg = "pin.svg",
    css_selector = "path[p-id=\"11477\"]",
    vjust = 0,
    size = 8,
    # Fraction of glyph height from the tip (vjust = 0) to the label.
    # pin.svg hole centre is 640/1024; 0.50 sits optically in the head.
    label_hole = 0.50,
    label_size = 2.6,
    labels = list(
        title = "SVG scatter",
        x = "x",
        y = "y",
        fill = "Group"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggsvg"))

script_dir <- dirname(normalizePath(sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]]
)))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns[c("x", "y", "fill")])
x_col <- config$columns$x
y_col <- config$columns$y
fill_col <- config$columns$fill
label_col <- config$columns$label
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
if (any(!is.finite(df[[x_col]]) | !is.finite(df[[y_col]]))) {
    stop("x and y must be finite numeric values.", call. = FALSE)
}
df[[fill_col]] <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))
df$.svg_fill <- df[[fill_col]]
has_label <- !is.null(label_col) && nzchar(label_col) && label_col %in% names(df)

svg_path <- file.path(script_dir, config$svg)
if (!file.exists(svg_path)) {
    stop(paste("SVG file not found:", svg_path), call. = FALSE)
}
svg_txt <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
selector <- config$css_selector
fill_values <- palette_colors("Qualitative.Safe", nlevels(df[[fill_col]]))
names(fill_values) <- levels(df[[fill_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
# geom_text vjust is in line heights, so it cannot sit in the pin hole
# across panel sizes. Offset the already-transformed y by millimetres.
GeomTextGlyph <- ggplot2::ggproto(
    "GeomTextGlyph",
    ggplot2::GeomText,
    extra_params = c("na.rm", "pin_size", "hole_frac"),
    draw_panel = function(data, panel_params, coord, parse = FALSE,
                          na.rm = FALSE, check_overlap = FALSE,
                          pin_size = 8, hole_frac = 0.625) {
        lab <- data$label
        data <- coord$transform(data, panel_params)
        grid::textGrob(
            lab,
            x = grid::unit(data$x, "native"),
            y = grid::unit(data$y, "native") +
                grid::unit(hole_frac * pin_size, "mm"),
            hjust = data$hjust,
            vjust = data$vjust,
            rot = data$angle,
            gp = grid::gpar(
                col = ggplot2::alpha(data$colour, data$alpha),
                fontsize = data$size * ggplot2::.pt,
                fontfamily = data$family,
                fontface = data$fontface,
                lineheight = data$lineheight
            ),
            check.overlap = check_overlap
        )
    }
)

geom_text_glyph <- function(mapping = NULL, data = NULL,
                            stat = "identity", position = "identity",
                            ..., pin_size = 8, hole_frac = 0.625,
                            na.rm = FALSE, show.legend = NA,
                            inherit.aes = TRUE) {
    ggplot2::layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomTextGlyph,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            na.rm = na.rm,
            pin_size = pin_size,
            hole_frac = hole_frac,
            ...
        )
    )
}

p <- ggplot(df) +
    ggsvg::geom_point_svg(
        ggplot2::aes(
            x = .data[[x_col]],
            y = .data[[y_col]],
            ggsvg::css(selector, fill = .svg_fill)
        ),
        svg = svg_txt,
        vjust = config$vjust,
        size = config$size
    )

if (has_label) {
    p <- p +
        geom_text_glyph(
            ggplot2::aes(
                x = .data[[x_col]],
                y = .data[[y_col]],
                label = .data[[label_col]]
            ),
            pin_size = config$size,
            hole_frac = config$label_hole,
            hjust = 0.5,
            vjust = 0.5,
            size = config$label_size,
            colour = "#333333",
            alpha = 1
        )
}

p <- p +
    ggsvg::scale_svg_fill_manual(
        aesthetics = ggsvg::css(selector, fill = .svg_fill),
        values = fill_values,
        name = config$labels$fill
    ) +
    ggplot2::scale_y_continuous(expand = expansion(c(0.1, 0.12))) +
    ggplot2::labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 7, height = 6)
