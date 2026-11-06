#!/usr/bin/env Rscript

# Template-ID: scatter-hist-rug
#
# Purpose:
#   Draw overlapping density histograms for a numeric value across
#   categories, with one rug row per category under the x axis.
#
# Inputs:
#   One row per observation. Default example:
#     - continent: group
#     - value: numeric (e.g. life expectancy)
#
# Output:
#   A PDF, PNG, or SVG histogram-with-rugs figure.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, bins, palette, and rug spacing.
#   Scatter axis rugs stay scatter-basic style = "rug". Marginal
#   boxplots stay scatter-marginal. Do not run density tests here.
#
# Scientific assumptions:
#   Each row is one observation. Histogram density and rugs are
#   display summaries of the supplied values.

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
        x = "value",
        group = "continent"
    ),
    labels = list(
        title = NULL,
        x = NULL,
        y = NULL,
        group = NULL
    ),
    palettes = list(
        group = "Qualitative.Bold"
    ),
    bins = 28,
    hist_alpha = 0.55,
    hist_colour = "grey20",
    hist_linewidth = 0.25,
    rug_size = 3.2,
    rug_linewidth = 0.55,
    rug_row_height = 0.012,  # density units per rug row
    rug_gap = 0.004,
    legend_inside = c(0.02, 0.98),
    base_size = 13,
    size = list(width = 7.2, height = 5.2)
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
x_col <- config$columns$x
g_col <- config$columns$group
df[[x_col]] <- as.numeric(df[[x_col]])
df[[g_col]] <- factor(df[[g_col]], levels = unique(df[[g_col]]))
ok <- is.finite(df[[x_col]]) & !is.na(df[[g_col]])
df <- df[ok, , drop = FALSE]
if (!nrow(df)) stop("No finite rows to plot.", call. = FALSE)

n_g <- nlevels(df[[g_col]])
cols <- expand_palette(palette_colors(config$palettes$group), n_g)
names(cols) <- levels(df[[g_col]])

# Rug y positions sit just below 0, one row per group.
df$rug_y <- -(config$rug_gap +
    (as.integer(df[[g_col]]) - 0.5) * config$rug_row_height)
y_floor <- -(config$rug_gap + n_g * config$rug_row_height + config$rug_gap)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = .data[[x_col]], fill = .data[[g_col]])) +
    geom_histogram(
        aes(y = after_stat(density)),
        bins = config$bins,
        alpha = config$hist_alpha,
        colour = config$hist_colour,
        linewidth = config$hist_linewidth,
        position = "identity"
    ) +
    geom_segment(
        aes(
            x = .data[[x_col]], xend = .data[[x_col]],
            y = rug_y - config$rug_row_height * 0.35,
            yend = rug_y + config$rug_row_height * 0.35,
            colour = .data[[g_col]]
        ),
        linewidth = config$rug_linewidth,
        show.legend = FALSE
    ) +
    scale_fill_manual(values = cols, name = config$labels$group) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.05)),
        limits = c(y_floor, NA)
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_bw(base_size = config$base_size) +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = config$legend_inside,
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.text = element_text(colour = "black")
    )

if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
    p <- p + theme(plot.title = element_blank())
}
if (is.null(config$labels$x) || !nzchar(as.character(config$labels$x))) {
    p <- p + theme(axis.title.x = element_blank())
}
if (is.null(config$labels$y) || !nzchar(as.character(config$labels$y))) {
    p <- p + theme(axis.title.y = element_blank())
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p, io$output,
    width = config$size$width,
    height = config$size$height
)
