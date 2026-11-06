#!/usr/bin/env Rscript

# Template-ID: bar-enrichment-lollipop
#
# Purpose:
#   Draw a horizontal enrichment lollipop chart: terms on y ranked by a
#   supplied score (e.g. EnrichR combined score), stem to a point sized
#   by -log10(adjusted p) and coloured by cluster, with an optional
#   shaded low-score band and vertical guide (preview.png).
#
# Inputs:
#   One row per term. Default example:
#     - term, score, p.adjust, cluster
#
# Output:
#   A PDF, PNG, or SVG enrichment lollipop plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, scales
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, cluster order/colours, size range, shade
#   xmax, line, labels, and size. Extra terms or clusters are extra
#   rows/levels, not a new id. A classic gene-ratio enrichment dotplot
#   stays bar-enrichment-dot; faceted cluster bars stay
#   bar-enrichment-cluster.
#
# Scientific assumptions:
#   score, p.adjust, and cluster are supplied. This script does not run
#   EnrichR, GSEA, or LOESS clustering. -log10(p.adjust) is a display
#   transform. The shaded band is a supplied x cutoff for display.

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
        term = "term",
        score = "score",
        pvalue = "p.adjust",
        cluster = "cluster"
    ),
    cluster_order = c("Cluster 1", "Cluster 3", "Cluster 4"),
    palettes = list(cluster = "Qualitative.Bold"),
    # Optional fixed colours keyed by cluster label (overrides palette).
    cluster_colours = c(
        `Cluster 1` = "#e73f74",
        `Cluster 3` = "#f2b701",
        `Cluster 4` = "#80ba5a"
    ),
    size_range = c(2.5, 8.5),
    stem_linewidth = 0.55,
    stem_colour = "grey55",
    # Shade from x = 0 to shade_xmax; NULL = no band.
    shade_xmax = 8,
    shade_fill = "#F8D7DA",
    shade_alpha = 0.55,
    # Vertical line at this score; NULL = use shade_xmax when set.
    score_line = NULL,
    score_line_colour = "grey40",
    score_line_linetype = "solid",
    score_line_linewidth = 0.45,
    labels = list(
        title = "EnrichR MSigDB analysis: LOESS clusters",
        x = "Combined score",
        y = NULL,
        size = expression(-log[10](P[adj])),
        colour = "Clusters with P adj < 0.05 pathways",
        caption = NULL
    ),
    legend_position = "right",
    size = list(width = 9.5, height = 6.2)
)

load_packages(c("ggplot2", "readr", "ggprism", "scales"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
term_col <- config$columns$term
score_col <- config$columns$score
p_col <- config$columns$pvalue
c_col <- config$columns$cluster

df[[score_col]] <- as.numeric(df[[score_col]])
df[[p_col]] <- as.numeric(df[[p_col]])
df <- df[
    is.finite(df[[score_col]]) & is.finite(df[[p_col]]) & df[[p_col]] > 0,
    ,
    drop = FALSE
]
if (!nrow(df)) stop("No finite enrichment rows.", call. = FALSE)

df$.neglogp <- -log10(df[[p_col]])
df[[term_col]] <- as.character(df[[term_col]])
df[[c_col]] <- as.character(df[[c_col]])

c_lv <- config$cluster_order
if (is.null(c_lv)) {
    c_lv <- unique(df[[c_col]])
} else {
    c_lv <- c(
        intersect(c_lv, unique(df[[c_col]])),
        setdiff(unique(df[[c_col]]), c_lv)
    )
}
df[[c_col]] <- factor(df[[c_col]], levels = c_lv)

# Rank by score descending (top = highest combined score).
df <- df[order(df[[score_col]], decreasing = TRUE), , drop = FALSE]
df[[term_col]] <- factor(df[[term_col]], levels = rev(unique(df[[term_col]])))

if (!is.null(config$cluster_colours) && length(config$cluster_colours)) {
    fill_cols <- unlist(config$cluster_colours)
    miss <- setdiff(c_lv, names(fill_cols))
    if (length(miss)) {
        extra <- expand_palette(palette_colors(config$palettes$cluster), length(miss))
        names(extra) <- miss
        fill_cols <- c(fill_cols, extra)
    }
    fill_cols <- fill_cols[c_lv]
} else {
    fill_cols <- expand_palette(
        palette_colors(config$palettes$cluster),
        length(c_lv)
    )
    names(fill_cols) <- c_lv
}

line_x <- config$score_line
if (is.null(line_x) && !is.null(config$shade_xmax)) {
    line_x <- config$shade_xmax
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = .data[[score_col]], y = .data[[term_col]]))
if (!is.null(config$shade_xmax)) {
    p <- p + annotate(
        "rect",
        xmin = -Inf,
        xmax = config$shade_xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = config$shade_fill,
        alpha = config$shade_alpha,
        colour = NA
    )
}
if (!is.null(line_x)) {
    p <- p + geom_vline(
        xintercept = line_x,
        colour = config$score_line_colour,
        linetype = config$score_line_linetype,
        linewidth = config$score_line_linewidth
    )
}
p <- p +
    geom_segment(
        aes(x = 0, xend = .data[[score_col]], yend = .data[[term_col]]),
        colour = config$stem_colour,
        linewidth = config$stem_linewidth
    ) +
    geom_point(
        aes(size = .neglogp, fill = .data[[c_col]]),
        shape = 21,
        colour = "grey25",
        stroke = 0.35
    ) +
    scale_fill_manual(
        values = fill_cols,
        name = config$labels$colour,
        guide = guide_legend(
            order = 1,
            override.aes = list(size = 4, shape = 21, colour = "grey25")
        )
    ) +
    scale_size_continuous(
        name = config$labels$size,
        range = config$size_range,
        breaks = pretty_breaks(n = 3),
        guide = guide_legend(order = 2, override.aes = list(fill = "grey70"))
    ) +
    scale_x_continuous(
        expand = expansion(mult = c(0, 0.06)),
        breaks = pretty_breaks(n = 6)
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        caption = config$labels$caption
    ) +
    theme_prism(base_size = 11) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = rel(0.85), hjust = 1),
        legend.position = config$legend_position,
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.75)),
        plot.title = element_text(face = "bold", size = rel(1.0), hjust = 0),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 10, 6, 8)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = config$size$width,
    height = config$size$height
)
