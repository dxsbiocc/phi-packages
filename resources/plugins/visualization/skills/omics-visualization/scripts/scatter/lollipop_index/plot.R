#!/usr/bin/env Rscript

# Template-ID: scatter-lollipop-index
#
# Purpose:
#   Draw a faceted signed lollipop / bubble grid of supplied indices
#   (e.g. Predictive Index × Therapeutic Index across treatments):
#   stems to y = 0, points coloured by group and sized by |value|,
#   with Favourable / Unfavourable cues (preview.png).
#
# Inputs:
#   One row per group-by-treatment-by-metric. Default example:
#     - cell, treatment, metric, value
#
# Output:
#   A PDF, PNG, or SVG faceted index lollipop plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggh4x, scales
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, facet orders, colours, size range, and
#   labels. Extra cells or treatments are extra rows/levels, not a new
#   id. A single signed delta lollipop stays scatter-lollipop-delta; a
#   dual GDP bubble panel stays scatter-bubble-panel.
#
# Scientific assumptions:
#   value, cell, treatment, and metric are supplied. This script does
#   not compute predictive or therapeutic indices. Point size defaults
#   to abs(value) as a display encoding.

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
        group = "cell",
        treatment = "treatment",
        metric = "metric",
        value = "value"
    ),
    treatment_order = c("Nab-PTX", "Nab-PTX + ATZ", "PTX + ATZ", "PTX"),
    metric_order = c("Predictive Index", "Therapeutic Index"),
    group_order = c(
        "Myeloid", "CD56+ NK", "pDC", "Mast", "CD4+ T",
        "CD16+ NK", "Plasma", "CD20+ B", "Bfoc", "CD8+ T"
    ),
    palettes = list(group = "Qualitative.Safe"),
    # Named colours for groups (overrides palette picks when set).
    group_colours = c(
        Myeloid = "#7EB6D9",
        `CD56+ NK` = "#F4A6C0",
        pDC = "#C5B0D5",
        Mast = "#2C5F8A",
        `CD4+ T` = "#E88B8B",
        `CD16+ NK` = "#E8C547",
        Plasma = "#C44D8A",
        `CD20+ B` = "#A8D5A2",
        Bfoc = "#9B7EBD",
        `CD8+ T` = "#5BBDBD"
    ),
    # Strip (column header) fills by treatment.
    strip_fill = c(
        `Nab-PTX` = "#A6CEE3",
        `Nab-PTX + ATZ` = "#FDBF6F",
        `PTX + ATZ` = "#CAB2D6",
        `PTX` = "#B2DF8A"
    ),
    size_range = c(1.6, 9.0),
    stem_linewidth = 0.35,
    stem_colour = "grey70",
    stem_linetype = "dashed",
    zero_linewidth = 0.55,
    x_jitter = 0.18,
    seed = 1L,
    show_direction = TRUE,
    labels = list(
        title = NULL,
        x = NULL,
        y = NULL,
        colour = NULL,
        size = NULL,
        favorable = "Favorable",
        unfavorable = "Unfavorable"
    ),
    size = list(width = 11, height = 6.5)
)

load_packages(c("ggplot2", "readr", "ggprism", "ggh4x", "scales"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
g_col <- config$columns$group
t_col <- config$columns$treatment
m_col <- config$columns$metric
v_col <- config$columns$value

df[[v_col]] <- as.numeric(df[[v_col]])
df <- df[is.finite(df[[v_col]]), , drop = FALSE]
if (!nrow(df)) stop("No finite index values.", call. = FALSE)

ord_or_unique <- function(pref, vals) {
    vals <- as.character(vals)
    if (is.null(pref)) {
        return(unique(vals))
    }
    c(intersect(as.character(pref), vals), setdiff(vals, as.character(pref)))
}

t_lv <- ord_or_unique(config$treatment_order, df[[t_col]])
m_lv <- ord_or_unique(config$metric_order, df[[m_col]])
g_lv <- ord_or_unique(config$group_order, df[[g_col]])

df[[t_col]] <- factor(as.character(df[[t_col]]), levels = t_lv)
df[[m_col]] <- factor(as.character(df[[m_col]]), levels = m_lv)
df[[g_col]] <- factor(as.character(df[[g_col]]), levels = g_lv)
df$.abs <- abs(df[[v_col]])

set.seed(as.integer(config$seed))
df$.x <- as.integer(df[[g_col]])
jit <- as.numeric(config$x_jitter)
if (is.finite(jit) && jit > 0) {
    # One offset per group so the same cell lines up across facets.
    off <- stats::setNames(stats::runif(length(g_lv), -jit, jit), g_lv)
    df$.x <- df$.x + unname(off[as.character(df[[g_col]])])
}

if (!is.null(config$group_colours) && length(config$group_colours)) {
    g_cols <- unlist(config$group_colours)
    miss <- setdiff(g_lv, names(g_cols))
    if (length(miss)) {
        extra <- expand_palette(palette_colors(config$palettes$group), length(miss))
        names(extra) <- miss
        g_cols <- c(g_cols, extra)
    }
    g_cols <- g_cols[g_lv]
} else {
    g_cols <- expand_palette(palette_colors(config$palettes$group), length(g_lv))
    names(g_cols) <- g_lv
}

strip_cols <- config$strip_fill
miss_st <- setdiff(t_lv, names(strip_cols))
if (length(miss_st)) {
    extra <- expand_palette(palette_colors("Qualitative.Safe"), length(miss_st))
    names(extra) <- miss_st
    strip_cols <- c(strip_cols, extra)
}
strip_cols <- strip_cols[t_lv]

strip <- ggh4x::strip_themed(
    background_x = lapply(t_lv, function(lv) {
        ggplot2::element_rect(fill = strip_cols[[lv]], colour = NA)
    }),
    text_x = lapply(t_lv, function(lv) {
        ggplot2::element_text(face = "bold", size = 10, colour = "grey15")
    }),
    background_y = lapply(m_lv, function(lv) {
        ggplot2::element_rect(fill = "grey92", colour = NA)
    }),
    text_y = lapply(m_lv, function(lv) {
        ggplot2::element_text(face = "bold", size = 10, colour = "grey20")
    })
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = .x, y = .data[[v_col]])) +
    geom_hline(
        yintercept = 0,
        colour = "grey20",
        linewidth = config$zero_linewidth
    ) +
    geom_segment(
        aes(x = .x, xend = .x, y = 0, yend = .data[[v_col]]),
        colour = config$stem_colour,
        linewidth = config$stem_linewidth,
        linetype = config$stem_linetype
    ) +
    geom_point(
        aes(size = .abs, fill = .data[[g_col]]),
        shape = 21,
        colour = "grey25",
        stroke = 0.3
    ) +
    ggh4x::facet_grid2(
        rows = vars(.data[[m_col]]),
        cols = vars(.data[[t_col]]),
        scales = "free_y",
        strip = strip
    ) +
    scale_fill_manual(values = g_cols, name = config$labels$colour) +
    scale_size_continuous(range = config$size_range, guide = "none") +
    scale_x_continuous(
        breaks = seq_along(g_lv),
        labels = g_lv,
        expand = expansion(add = 0.7)
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism(base_size = 11) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.85)),
        legend.key.width = unit(0.9, "lines"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.spacing = unit(0.45, "lines"),
        plot.margin = margin(8, 10, 6, 8)
    ) +
    guides(
        fill = guide_legend(
            nrow = 2, byrow = TRUE,
            override.aes = list(size = 4)
        )
    )

if (isTRUE(config$show_direction)) {
    left_tr <- t_lv[[1]]
    ann_list <- lapply(m_lv, function(m) {
        sub <- df[df[[m_col]] == m & df[[t_col]] == left_tr, , drop = FALSE]
        if (!nrow(sub)) return(NULL)
        ymax <- max(c(sub[[v_col]], 0.05), na.rm = TRUE)
        ymin <- min(c(sub[[v_col]], -0.05), na.rm = TRUE)
        data.frame(
            treatment = factor(left_tr, levels = t_lv),
            metric = factor(m, levels = m_lv),
            x = 0.2,
            y_fav = ymax * 0.90,
            y_unf = ymin * 0.90,
            stringsAsFactors = FALSE
        )
    })
    ann <- do.call(rbind, ann_list)
    if (!is.null(ann) && nrow(ann)) {
        names(ann)[names(ann) == "treatment"] <- t_col
        names(ann)[names(ann) == "metric"] <- m_col
        p <- p +
            geom_text(
                data = ann,
                aes(
                    x = x, y = y_fav,
                    label = paste0("^ ", config$labels$favorable)
                ),
                inherit.aes = FALSE,
                hjust = 0, size = 2.55, colour = "#C0392B", fontface = "bold"
            ) +
            geom_text(
                data = ann,
                aes(
                    x = x, y = y_unf,
                    label = paste0("v ", config$labels$unfavorable)
                ),
                inherit.aes = FALSE,
                hjust = 0, size = 2.55, colour = "#1E8449", fontface = "bold"
            )
    }
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = config$size$width,
    height = config$size$height
)
