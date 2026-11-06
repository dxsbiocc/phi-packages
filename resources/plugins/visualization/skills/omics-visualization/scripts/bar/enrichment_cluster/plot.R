#!/usr/bin/env Rscript

# Template-ID: bar-enrichment-cluster
#
# Purpose:
#   Draw horizontal enrichment bars faceted by cluster, with matching
#   panel / strip colours, rotated strip labels (cluster + dose), and
#   an optional dashed significance threshold.
#
# Inputs:
#   One row per term (example.tsv):
#     - cluster: facet / strip group (e.g. Cluster 1)
#     - dose: optional second strip line (e.g. Low Klf4 dose)
#     - Description: GO / pathway term
#     - p.adjust: adjusted p-value (display as -log10)
#
# Output:
#   A PDF, PNG, or SVG clustered enrichment bar chart.
#
# Dependencies:
#   ggplot2, readr, ggh4x
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, cluster order, palette refs (gold / sky /
#   teal), threshold, and strip label layout. Do not run enricher here;
#   p.adjust is supplied. Same terms across contrasts stay
#   bar-enrichment-groups. Category headers on one y-axis stay
#   bar-enrichment-expand.
#
# Scientific assumptions:
#   Each row is one enrichment term with a supplied adjusted p-value.
#   -log10(p.adjust) is a display transform. The threshold line marks
#   a supplied alpha (default 0.05), not a new test.

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
        description = "Description",
        pvalue = "p.adjust",
        cluster = "cluster",
        dose = "dose"            # NULL to omit the second strip line
    ),
    labels = list(
        title = NULL,
        x = expression(-log[10](P[adj])),
        y = NULL
    ),
    cluster_order = NULL,        # e.g. c("Cluster 1", "Cluster 2", "Cluster 3")
    # One catalog ref per cluster (gold, sky blue, teal).
    palettes = list(
        cluster = c(
            "Qualitative.Bold:4",
            "Qualitative.Safe:1",
            "Qualitative.Bold:2"
        )
    ),
    panel_alpha = 0.18,          # tinted panel fill from cluster colour
    strip_alpha = 0.22,
    bar_width = 0.72,
    threshold = 0.05,            # NULL to hide; drawn at -log10(threshold)
    threshold_colour = "grey55",
    threshold_linetype = "dashed",
    threshold_linewidth = 0.45,
    base_size = 12,
    size = list(width = 8.5, height = 9.5)
)

load_packages(c("ggplot2", "readr", "ggh4x"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols[[1]]
}

lighten <- function(col, amount = 0.72) {
    rgb <- grDevices::col2rgb(col) / 255
    rgb <- rgb + (1 - rgb) * amount
    grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ])
}

df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
desc_col <- config$columns$description
p_col <- config$columns$pvalue
c_col <- config$columns$cluster
d_col <- config$columns$dose
require_columns(df, list(desc_col, p_col, c_col))

df[[p_col]] <- as.numeric(df[[p_col]])
df[[desc_col]] <- as.character(df[[desc_col]])
df[[c_col]] <- as.character(df[[c_col]])
df <- df[is.finite(df[[p_col]]) & df[[p_col]] > 0, , drop = FALSE]
if (!nrow(df)) stop("No finite positive p-values.", call. = FALSE)

df$neglog <- -log10(df[[p_col]])

cl_lv <- unique(df[[c_col]])
if (!is.null(config$cluster_order)) {
    cl_lv <- c(
        intersect(as.character(config$cluster_order), cl_lv),
        setdiff(cl_lv, as.character(config$cluster_order))
    )
}
df[[c_col]] <- factor(df[[c_col]], levels = cl_lv)

# Strip label: "Cluster 1\nLow Klf4 dose"
if (!is.null(d_col) && d_col %in% names(df)) {
    dose_map <- vapply(cl_lv, function(cl) {
        vals <- unique(as.character(df[[d_col]][df[[c_col]] == cl]))
        vals <- vals[!is.na(vals) & nzchar(vals)]
        if (!length(vals)) "" else vals[[1]]
    }, character(1))
    strip_labs <- ifelse(
        nzchar(dose_map),
        paste0(cl_lv, "\n", dose_map),
        cl_lv
    )
} else {
    strip_labs <- cl_lv
}
names(strip_labs) <- cl_lv
df$strip <- factor(
    strip_labs[as.character(df[[c_col]])],
    levels = strip_labs
)

# Order terms within each cluster by decreasing -log10(p).
df <- df[order(df$strip, -df$neglog, df[[desc_col]]), , drop = FALSE]
# Unique y keys so repeated term names across clusters stay distinct.
df$y_key <- paste(as.character(df$strip), df[[desc_col]], sep = " | ")
df$y_key <- factor(df$y_key, levels = rev(unique(df$y_key)))
df$y_label <- df[[desc_col]]

n_cl <- nlevels(df$strip)
bar_cols <- vapply(
    seq_len(n_cl),
    function(i) {
        refs <- config$palettes$cluster
        pal_ref(refs[[((i - 1L) %% length(refs)) + 1L]])
    },
    character(1)
)
names(bar_cols) <- levels(df$strip)
panel_cols <- vapply(bar_cols, lighten, character(1), amount = 0.78)
strip_bg <- vapply(bar_cols, lighten, character(1), amount = 0.55)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
bg <- data.frame(
    strip = factor(levels(df$strip), levels = levels(df$strip)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    stringsAsFactors = FALSE
)

strip <- ggh4x::strip_themed(
    background_y = lapply(levels(df$strip), function(lv) {
        ggplot2::element_rect(
            fill = grDevices::adjustcolor(strip_bg[[lv]], alpha.f = 1),
            colour = NA
        )
    }),
    text_y = lapply(levels(df$strip), function(lv) {
        ggplot2::element_text(
            colour = bar_cols[[lv]],
            face = "bold",
            size = 10,
            angle = 90,
            hjust = 0.5,
            vjust = 0.5,
            lineheight = 0.95
        )
    })
)

p <- ggplot(df, aes(x = .data$neglog, y = .data$y_key, fill = .data$strip)) +
    geom_rect(
        data = bg,
        aes(xmin = .data$xmin, xmax = .data$xmax,
            ymin = .data$ymin, ymax = .data$ymax, fill = .data$strip),
        alpha = config$panel_alpha,
        inherit.aes = FALSE,
        show.legend = FALSE
    ) +
    geom_col(width = config$bar_width, show.legend = FALSE) +
    scale_fill_manual(values = bar_cols, guide = "none") +
    scale_y_discrete(labels = function(x) {
        map <- setNames(as.character(df$y_label), as.character(df$y_key))
        unname(map[x])
    }) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    ggh4x::facet_grid2(
        strip ~ .,
        scales = "free_y",
        space = "free_y",
        switch = "y",
        strip = strip
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_classic(base_size = config$base_size) +
    theme(
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = rel(0.9), hjust = 1),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing = grid::unit(0.55, "lines"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 12, 8, 8)
    )

if (!is.null(config$threshold) && is.finite(config$threshold) &&
    config$threshold > 0 && config$threshold < 1) {
    p <- p + geom_vline(
        xintercept = -log10(config$threshold),
        colour = config$threshold_colour,
        linetype = config$threshold_linetype,
        linewidth = config$threshold_linewidth
    )
}

if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
    p <- p + theme(plot.title = element_blank())
}
if (is.null(config$labels$y) ||
    (is.character(config$labels$y) && !nzchar(config$labels$y))) {
    p <- p + theme(axis.title.y = element_blank())
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width, height = config$size$height)
