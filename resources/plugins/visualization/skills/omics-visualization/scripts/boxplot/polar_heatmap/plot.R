#!/usr/bin/env Rscript

# Template-ID: boxplot-polar-heatmap
#
# Purpose:
#   Draw a radial fan: inner heatmap of group means and outer violins
#   of the raw values, with curved category labels.
#
# Inputs:
#   A table with one row per observation. Default example:
#     - gene: angular category
#     - exprs: numeric value
#     - grade: inner heatmap row
#
# Output:
#   A PDF, PNG, or SVG polar heatmap-violin plot.
#
# Dependencies:
#   ggplot2, readr, ggnewscale, geomtextpath
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change mean aggregation or radial
#   stacking of the heatmap versus violins.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one observation of a numeric value in a category
#   and grade.
#   Heatmap cells are means of the numeric column by gene and grade.
#   Polar coordinates are a display transform.
#   Violins are shifted radially so they sit outside the heatmap.

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
        x = "gene",
        y = "exprs",
        tile = "grade"
    ),
    labels = list(
        title = "",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggnewscale", "geomtextpath"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
tile_col <- config$columns$tile

df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[tile_col]] <- factor(df[[tile_col]], levels = unique(df[[tile_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

x_size <- nlevels(df[[x_col]])
cate_names <- levels(df[[tile_col]])
add_high <- length(cate_names) + 1
y_max <- ceiling(max(df[[y_col]], na.rm = TRUE)) + add_high

agg <- aggregate(
    df[[y_col]],
    by = list(x = df[[x_col]], tile = df[[tile_col]]),
    FUN = mean,
    na.rm = TRUE
)
names(agg) <- c(x_col, tile_col, "mu")
agg[[x_col]] <- factor(agg[[x_col]], levels = levels(df[[x_col]]))
agg[[tile_col]] <- factor(agg[[tile_col]], levels = cate_names)
agg$index <- as.numeric(agg[[tile_col]])

label_data <- data.frame(y = rep(y_max + 0.5, x_size))
label_data[[x_col]] <- levels(df[[x_col]])

y_labels <- as.integer(seq(add_high, y_max, length.out = 4))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_base <- c(
    "#d7689d", "#f3d337", "#a0c443", "#66b5ae", "#5da4dc", "#40769e"
)
fill_values <- if (length(fill_base) < x_size) {
    colorRampPalette(fill_base)(x_size)
} else {
    fill_base
}
tile_colours <- colorRampPalette(c("#E9E29C", "#d7689d"))(50)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[x_col]]
)) +
    geom_violin(
        aes(y = .data[[y_col]] + add_high),
        show.legend = FALSE
    ) +
    scale_fill_manual(values = fill_values) +
    new_scale_fill() +
    geom_tile(
        data = agg,
        aes(
            x = .data[[x_col]],
            y = index,
            fill = mu
        ),
        colour = "#fffffe",
        linewidth = 0.5,
        inherit.aes = FALSE
    ) +
    scale_fill_gradientn(name = NULL, colours = tile_colours) +
    geom_textpath(
        data = label_data,
        aes(x = .data[[x_col]], y = y, label = .data[[x_col]]),
        vjust = 1.5,
        inherit.aes = FALSE
    ) +
    annotate(
        "segment",
        x = c(0.5, 0.5),
        xend = c(x_size + 0.5, 0.5),
        y = c(y_max + 0.5, add_high),
        yend = c(y_max + 0.5, y_max)
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(-1, y_max + 0.5),
        breaks = c(seq_along(cate_names), y_labels),
        labels = c(cate_names, y_labels - add_high),
        guide = guide_axis(angle = -5)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    coord_radial(start = -1.047198, end = 1.047198) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme(
        plot.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        panel.grid.major.x = element_line(
            linetype = "dashed",
            linewidth = 0.5,
            colour = "grey50"
        ),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.517, 0.01),
        legend.key.height = unit(4, "mm"),
        legend.title = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
