#!/usr/bin/env Rscript

# Template-ID: heatmap-obs-dot
#
# Purpose:
#   Draw a category-by-category observation-dot matrix: equal-size
#   circles coloured by a supplied mean, orange stars for supplied
#   significance, top bars for column n, left bars for row n, and
#   right-hand row-category strips (preview.png).
#
# Inputs:
#   One row per cell (example.tsv):
#     - disease, symptom, category, mean_obs, significant
#   Sidecars next to the input:
#     - col_n.tsv: disease, n
#     - row_n.tsv: symptom, n  (category may also appear here)
#
# Output:
#   A PDF, PNG, or SVG observation-dot heatmap with marginal bars.
#
# Dependencies:
#   ggplot2, readr, ggprism, scales, aplot, patchwork, grid
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, fill limits / palette, star colour,
#   panel fill, and labels. Extra diseases or symptoms are extra
#   levels, not a new id. Do not add an id for the marginal bars
#   alone. Correlation bubble grids stay heatmap-corr-bubble.
#
# Scientific assumptions:
#   mean_obs, significant, and the n sidecars are supplied display
#   values. This script does not aggregate raw case reports, test
#   enrichment, or invent significance.

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
        x = "disease",
        y = "symptom",
        category = "category",
        value = "mean_obs",
        significant = "significant"
    ),
    col_n_file = "col_n.tsv",
    row_n_file = "row_n.tsv",
    fill_limits = c(0, 5),
    fill_palette = "Gradient.YlGnBu",
    star_colour = "#E36C2C",
    panel_fill = "#EDE7F6",
    category_colours = c(
        General = "#7E57C2",
        Motor = "#42A5F5",
        `Sensory/autonomic` = "#66BB6A"
    ),
    point_size = 4.2,
    star_size = 2.4,
    labels = list(
        fill = "Mean observations",
        top_n = NULL,
        left_n = NULL
    ),
    size = list(width = 11.5, height = 10.5)
)

load_packages(c(
    "ggplot2", "readr", "ggprism", "scales", "aplot", "patchwork", "grid"
))

sidecar_path <- function(name) {
    if (is.null(name) || !nzchar(as.character(name)[[1]])) return(NULL)
    name <- as.character(name)[[1]]
    if (file.exists(name)) return(normalizePath(name, mustWork = TRUE))
    cand <- file.path(dirname(io$input), name)
    if (file.exists(cand)) return(normalizePath(cand, mustWork = TRUE))
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
        cand2 <- file.path(script_dir, name)
        if (file.exists(cand2)) return(normalizePath(cand2, mustWork = TRUE))
    }
    NULL
}

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
cat_col <- config$columns$category
v_col <- config$columns$value
sig_col <- config$columns$significant

df[[v_col]] <- as.numeric(df[[v_col]])
df[[sig_col]] <- as.logical(as.integer(df[[sig_col]]))
df <- df[is.finite(df[[v_col]]), , drop = FALSE]
if (!nrow(df)) stop("No finite observation rows.", call. = FALSE)

x_lv <- unique(as.character(df[[x_col]]))
y_lv <- unique(as.character(df[[y_col]]))
cat_lv <- unique(as.character(df[[cat_col]]))
df[[x_col]] <- factor(df[[x_col]], levels = x_lv)
# ggplot draws first factor level at the bottom; keep file order top-down.
df[[y_col]] <- factor(df[[y_col]], levels = rev(y_lv))
df[[cat_col]] <- factor(df[[cat_col]], levels = cat_lv)

col_n_path <- sidecar_path(config$col_n_file)
row_n_path <- sidecar_path(config$row_n_file)
if (is.null(col_n_path) || is.null(row_n_path)) {
    stop("Need col_n.tsv and row_n.tsv beside the input.", call. = FALSE)
}
col_n <- as.data.frame(read_table_auto(col_n_path), stringsAsFactors = FALSE)
row_n <- as.data.frame(read_table_auto(row_n_path), stringsAsFactors = FALSE)
require_columns(col_n, list(x_col, "n"))
require_columns(row_n, list(y_col, "n"))
col_n[[x_col]] <- factor(as.character(col_n[[x_col]]), levels = x_lv)
row_n[[y_col]] <- factor(as.character(row_n[[y_col]]), levels = levels(df[[y_col]]))
col_n$n <- as.numeric(col_n$n)
row_n$n <- as.numeric(row_n$n)

fill_cols <- palette_colors(config$fill_palette)
cat_cols <- config$category_colours
miss_cat <- setdiff(levels(df[[cat_col]]), names(cat_cols))
if (length(miss_cat)) {
    extra <- expand_palette(palette_colors("Qualitative.Bold"), length(miss_cat))
    names(extra) <- miss_cat
    cat_cols <- c(cat_cols, extra)
}
cat_cols <- cat_cols[levels(df[[cat_col]])]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
blank_bg <- theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
)

p_main <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_tile(
        fill = config$panel_fill,
        colour = NA,
        width = 1,
        height = 1
    ) +
    geom_point(
        aes(fill = .data[[v_col]]),
        shape = 21,
        colour = "grey35",
        stroke = 0.25,
        size = config$point_size
    ) +
    geom_point(
        data = df[which(as.logical(df[[sig_col]])), , drop = FALSE],
        shape = 8,
        colour = config$star_colour,
        size = config$star_size,
        stroke = 0.55
    ) +
    facet_grid(
        rows = vars(.data[[cat_col]]),
        scales = "free_y",
        space = "free_y"
    ) +
    scale_fill_gradientn(
        colours = fill_cols,
        limits = config$fill_limits,
        oob = scales::squish,
        name = config$labels$fill,
        guide = guide_colorbar(
            title.position = "top",
            barwidth = unit(2.8, "lines"),
            barheight = unit(0.55, "lines"),
            frame.colour = "grey40",
            ticks.colour = "grey40"
        )
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0.15)) +
    scale_y_discrete(expand = c(0, 0.15)) +
    labs(x = NULL, y = NULL) +
    theme_prism(base_size = 10) +
    blank_bg +
    theme(
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(
            angle = 90, vjust = 0.5, hjust = 0, size = rel(0.85)
        ),
        axis.text.y = element_text(size = rel(0.85), hjust = 1),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text.y.right = element_text(
            angle = 0, face = "bold", size = rel(0.85), hjust = 0
        ),
        legend.position = "top",
        legend.justification = "right",
        legend.direction = "horizontal",
        legend.title = element_text(size = rel(0.85), face = "bold"),
        legend.text = element_text(size = rel(0.75)),
        panel.spacing.y = unit(0.35, "lines"),
        plot.margin = margin(2, 6, 2, 2)
    )

# Colour the facet strip text by category via labeller + manual strips:
# ggplot2 strip text colour is global; paint category labels with a
# right-hand annotation strip instead when strip colour is needed.
# Keep facet labels; tint them via a small right panel below.

p_top <- ggplot(col_n, aes(x = .data[[x_col]], y = n)) +
    geom_col(fill = "#5C6BC0", width = 0.7, colour = NA) +
    scale_x_discrete(expand = c(0, 0.15)) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.08)),
        breaks = scales::pretty_breaks(4)
    ) +
    labs(x = NULL, y = config$labels$top_n) +
    theme_prism(base_size = 9) +
    blank_bg +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size = rel(0.8)),
        axis.text.y = element_text(size = rel(0.7)),
        plot.margin = margin(2, 6, 0, 2)
    )

p_left <- ggplot(row_n, aes(x = n, y = .data[[y_col]])) +
    geom_col(fill = "#5C6BC0", width = 0.7, colour = NA) +
    scale_y_discrete(expand = c(0, 0.15)) +
    scale_x_reverse(
        expand = expansion(mult = c(0.08, 0)),
        breaks = scales::pretty_breaks(3)
    ) +
    labs(x = config$labels$left_n, y = NULL) +
    theme_prism(base_size = 9) +
    blank_bg +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = rel(0.75)),
        axis.text.x = element_text(size = rel(0.65), angle = 90, vjust = 0.5),
        plot.margin = margin(2, 0, 2, 2)
    )

# Right category colour rail aligned to facet rows.
rail <- unique(df[, c(y_col, cat_col), drop = FALSE])
p_rail <- ggplot(rail, aes(x = 1, y = .data[[y_col]], fill = .data[[cat_col]])) +
    geom_tile(width = 1, height = 1, colour = NA) +
    facet_grid(
        rows = vars(.data[[cat_col]]),
        scales = "free_y",
        space = "free_y"
    ) +
    scale_fill_manual(values = cat_cols, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0.15)) +
    theme_void() +
    theme(
        strip.text = element_blank(),
        panel.spacing.y = unit(0.35, "lines"),
        plot.margin = margin(2, 2, 2, 0)
    )

p <- p_main %>%
    insert_top(p_top, height = 0.14) %>%
    insert_left(p_left, width = 0.12) %>%
    insert_right(p_rail, width = 0.035)

p <- as.patchwork(p) +
    plot_annotation(
        theme = theme(
            plot.background = element_rect(fill = "white", colour = NA)
        )
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
suppressWarnings(save_ggplot(
    p, io$output,
    width = config$size$width,
    height = config$size$height
))
