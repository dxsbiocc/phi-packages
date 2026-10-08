#!/usr/bin/env Rscript

# Template-ID: bar-svg-icon
#
# Purpose:
#   Draw a horizontal bar for each category, with a supplied SVG glyph
#   as the category mark instead of y-axis text.
#
# Inputs:
#   One row per bar. Default example:
#     - year: category (displayed as the y position)
#     - birth_population: bar length
#     - zodiac: basename of an SVG in config$svg_dir (no extension)
#   Bundled glyphs: svg/鼠.svg … svg/猪.svg in this template directory.
#
# Output:
#   A PDF, PNG, or SVG horizontal bar chart with SVG category marks.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggsvg
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, svg_dir, labels). Place
#   one SVG per category name in svg_dir beside the copied plot.R.
#   Edit DATA PREPARATION to change bar order. Edit PLOT only when
#   geometry must change. A unit/percent icon row (N glyphs = 100%) is
#   not this template; use bar-waffle for square units or scatter-svg
#   for x-y marks. Do not add an id that only swaps the SVG set.
#
# Scientific assumptions:
#   Each row is one already-summarized category with a comparable value.
#   The SVG identifies the category; it is not a second numeric encoding.
#   Fill of the bar is the same value as bar length. Example counts are
#   a supplied series, not computed in this script.

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
        y = "year",
        x = "birth_population",
        icon = "zodiac"
    ),
    svg_dir = "svg",
    labels = list(
        title = "",
        x = "Number of newborns in China (2012-2023)",
        fill = "Population (W)"
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
require_columns(df, config$columns)
y_col <- config$columns$y
x_col <- config$columns$x
icon_col <- config$columns$icon
df[[x_col]] <- as.numeric(df[[x_col]])
if (any(!is.finite(df[[x_col]]) | df[[x_col]] < 0)) {
    stop("Bar values must be finite and non-negative.", call. = FALSE)
}
df[[icon_col]] <- as.character(df[[icon_col]])
if (any(!nzchar(df[[icon_col]]))) {
    stop("Icon names must be non-empty.", call. = FALSE)
}
df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
svg_dir <- file.path(script_dir, config$svg_dir)
icon_files <- file.path(svg_dir, paste0(df[[icon_col]], ".svg"))
missing <- icon_files[!file.exists(icon_files)]
if (length(missing)) {
    stop(
        paste("SVG file(s) not found:", paste(basename(missing), collapse = ", ")),
        call. = FALSE
    )
}
ramp <- palette_colors("Gradient.YlOrRd")
fill_low <- ramp[[2]]
fill_high <- ramp[[length(ramp)]]
icon_width <- max(df[[x_col]]) * 0.12

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
    ggplot2::geom_col(
        ggplot2::aes(fill = .data[[x_col]]),
        width = 0.72,
        colour = NA
    )

n <- nrow(df)
for (i in seq_len(n)) {
    svg_txt <- paste(readLines(icon_files[[i]], warn = FALSE), collapse = "\n")
    grob <- ggsvg::svg_to_rasterGrob(svg_txt)
    y_i <- as.numeric(df[[y_col]][[i]])
    p <- p +
        ggplot2::annotation_custom(
            grob = grob,
            xmin = -icon_width,
            xmax = 0,
            ymin = y_i - 0.45,
            ymax = y_i + 0.45
        )
}

p <- p +
    ggplot2::scale_fill_gradient(
        name = config$labels$fill,
        low = fill_low,
        high = fill_high
    ) +
    ggplot2::scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    ggplot2::coord_cartesian(xlim = c(-icon_width, NA), clip = "off") +
    ggplot2::labs(
        title = config$labels$title,
        x = config$labels$x,
        y = NULL
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.title = element_text()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 6)
