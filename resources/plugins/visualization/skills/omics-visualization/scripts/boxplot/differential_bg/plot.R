#!/usr/bin/env Rscript

# Template-ID: boxplot-differential-bg
#
# Purpose:
#   Draw a two-level grouped boxplot: categories on x, a within-
#   category contrast as colour, alternating background bands, and
#   per-category significance labels (pancancer Normal vs Tumor style).
#
# Inputs:
#   One row per observation. Default example:
#     - type: category on the x-axis (cancer type)
#     - exprs: numeric value
#     - group: within-category contrast (Normal / Tumor)
#
# Output:
#   A PDF, PNG, or SVG banded grouped boxplot.
#
# Dependencies:
#   ggplot2, readr, ggpubr, ggprism, gground
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Extra categories are extra x levels, not extra ids. show_ns,
#   legend position, and x-label angle stay in CONFIG. Do not filter
#   matched Normal/Tumor patients or call DE in this script.
#
# Scientific assumptions:
#   Each row is one unpaired observation in a category and group.
#   Background bands are a display aid. Significance labels are
#   Wilcoxon (ggpubr default) within each category; they are not
#   multiple-testing corrected unless supplied in the table.

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
        x = "type",
        y = "exprs",
        group = "group"
    ),
    palettes = list(
        group = "Diverging.Temps"
    ),
    show_ns = TRUE,
    labels = list(
        title = "",
        x = "",
        y = "Expression Level (log2 TPM)",
        group = NULL
    ),
    size = list(
        width = 12,
        height = 5.2
    )
)

load_packages(c("ggplot2", "readr", "ggpubr", "ggprism", "gground"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group

df[[x_col]] <- factor(as.character(df[[x_col]]), levels = unique(as.character(df[[x_col]])))
df[[group_col]] <- factor(
    as.character(df[[group_col]]),
    levels = unique(as.character(df[[group_col]]))
)
df[[y_col]] <- as.numeric(df[[y_col]])
ok <- !is.na(df[[x_col]]) & is.finite(df[[y_col]]) & !is.na(df[[group_col]])
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with missing category, group, or value.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

x_unique <- levels(df[[x_col]])
n_x <- length(x_unique)
y_rng <- range(df[[y_col]], na.rm = TRUE)
y_pad <- diff(y_rng)
if (!is.finite(y_pad) || y_pad <= 0) {
    y_pad <- 1
}
fill_data <- data.frame(
    xmin = seq_len(n_x) - 0.5,
    xmax = seq_len(n_x) + 0.5,
    ymin = y_rng[[1]] - 0.08 * y_pad,
    ymax = y_rng[[2]] + 0.22 * y_pad,
    band = factor(((seq_len(n_x) - 1L) %% 2L) + 1L),
    stringsAsFactors = FALSE
)

temps <- palette_colors(config$palettes$group)
group_cols <- expand_palette(c(temps[[1]], temps[[length(temps)]]), nlevels(df[[group_col]]))
names(group_cols) <- levels(df[[group_col]])
# Band fills match the TCGA single-gene pancancer figure (Brand.new balance).
nb <- palette_colors("Brand.new balance")
band_cols <- c(nb[[5]], nb[[4]])
names(band_cols) <- levels(fill_data$band)
outline <- palette_colors("Brand.Aiesec")[[8]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
# roundrectGrob legend keys stroke past corners on PNG/Cairo; PDF is fine.
# Use the standard boxplot glyph for the colour legend only.
round_box <- gground::geom_round_boxplot(outliers = FALSE)
round_box$geom$draw_key <- ggplot2::draw_key_boxplot

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    colour = .data[[group_col]]
)) +
    geom_rect(
        data = fill_data,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            fill = band
        ),
        linetype = "dashed",
        alpha = 0.2,
        colour = outline,
        linewidth = 0.3,
        show.legend = FALSE,
        inherit.aes = FALSE
    ) +
    round_box +
    geom_point(
        position = position_jitterdodge(
            jitter.width = 0.15,
            dodge.width = 0.7
        ),
        alpha = 0.4,
        size = 2,
        stroke = 0,
        show.legend = FALSE
    ) +
    ggpubr::stat_compare_means(
        aes(group = .data[[group_col]]),
        label = "p.signif",
        label.y.npc = 0.95,
        hide.ns = !isTRUE(config$show_ns),
        show.legend = FALSE
    ) +
    scale_colour_manual(
        values = group_cols,
        name = if (is.null(config$labels$group)) waiver() else config$labels$group
    ) +
    scale_fill_manual(values = band_cols, guide = "none") +
    scale_x_discrete(expand = c(0, 0), guide = "prism_bracket") +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(
        ylim = c(fill_data$ymin[[1]], fill_data$ymax[[1]]),
        clip = "off"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = as.numeric(config$size$width),
    height = as.numeric(config$size$height)
)
