#!/usr/bin/env Rscript

# Template-ID: scatter-pathway-facet
#
# Purpose:
#   Faceted pathway-enrichment scatter: x = -log10(FDR), y = enrichment
#   score, one panel per group (e.g. spatial ecotype). Coloured strip
#   headers, sized points, and leader labels for highlighted pathways.
#
# Inputs:
#   A long table (tsv/csv). Default example:
#     - ecotype: panel / group (e.g. C-I … C-IV)
#     - pathway: pathway name
#     - neg_log10_fdr: x
#     - score: y (enrichment score)
#     - size: point size (often same as -log10(FDR) or a hit ratio)
#     - show_label: 1 to label, 0 to skip (optional; else top-N by size)
#
# Output:
#   A PDF, PNG, or SVG multi-panel enrichment scatter.
#
# Dependencies:
#   ggplot2, readr, ggrepel, ggh4x
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG for columns, labels, palette, facet nrow, and
#   labelling. Edit PLOT only when the geometry must change.
#
# Scientific assumptions:
#   Enrichment statistics are supplied. This script does not run
#   pathway enrichment, spatial clustering, or multiple-testing
#   correction. Point size is a display aesthetic of a supplied column.

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
        group = "ecotype",
        pathway = "pathway",
        x = "neg_log10_fdr",
        y = "score",
        size = "size",
        show_label = "show_label"
    ),
    labels = list(
        x = expression(-log[10](FDR)),
        y = "Pathway enrichment score"
    ),
    palettes = list(
        group = "Qualitative.Safe"
    ),
    # Fixed colours for common ecotype labels; NULL → catalog palette.
    group_colours = c(
        "C-I" = "#d7ae2b",
        "C-II" = "#b85b5b",
        "C-III" = "#5b83b5",
        "C-IV" = "#93b15a"
    ),
    facet_nrow = 2L,
    max_labels_per_panel = 5L,
    point_stroke = 0.45,
    size_range = c(2.0, 5.0),
    label_size = 2.2
)

load_packages(c("ggplot2", "readr", "ggrepel", "ggh4x"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns[c("group", "pathway", "x", "y")])
g_col <- config$columns$group
p_col <- config$columns$pathway
x_col <- config$columns$x
y_col <- config$columns$y
s_col <- config$columns$size
lab_col <- config$columns$show_label

df[[g_col]] <- as.character(df[[g_col]])
df[[p_col]] <- as.character(df[[p_col]])
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
if (!is.null(s_col) && s_col %in% names(df)) {
    df[[s_col]] <- as.numeric(df[[s_col]])
} else {
    s_col <- ".size"
    df[[s_col]] <- df[[x_col]]
}
ok <- is.finite(df[[x_col]]) & is.finite(df[[y_col]]) &
    is.finite(df[[s_col]]) & !is.na(df[[g_col]]) & nzchar(df[[g_col]]) &
    !is.na(df[[p_col]]) & nzchar(df[[p_col]])
df <- df[ok, , drop = FALSE]
if (!nrow(df)) {
    stop("No finite pathway rows to plot.", call. = FALSE)
}

group_levels <- unique(df[[g_col]])
if (length(group_levels) < 2L) {
    stop("Need at least two groups/panels.", call. = FALSE)
}
df[[g_col]] <- factor(df[[g_col]], levels = group_levels)

# Which points get text labels.
df$`.label` <- ""
if (!is.null(lab_col) && lab_col %in% names(df)) {
    show <- as.integer(df[[lab_col]])
    show[is.na(show)] <- 0L
    df$`.label`[show > 0L] <- df[[p_col]][show > 0L]
} else {
    n_lab <- as.integer(config$max_labels_per_panel)
    if (!is.finite(n_lab) || n_lab < 1L) {
        n_lab <- 5L
    }
    for (g in group_levels) {
        idx <- which(df[[g_col]] == g)
        ord <- order(df[[s_col]][idx], decreasing = TRUE)
        take <- idx[ord[seq_len(min(n_lab, length(ord)))]]
        df$`.label`[take] <- df[[p_col]][take]
    }
}

pal <- config$group_colours
if (is.null(pal) || !length(pal)) {
    pal <- palette_colors(config$palettes$group, n = length(group_levels))
    names(pal) <- group_levels
} else {
    miss <- setdiff(group_levels, names(pal))
    if (length(miss)) {
        extra <- palette_colors(config$palettes$group, n = length(miss))
        names(extra) <- miss
        pal <- c(pal, extra)
    }
    pal <- pal[group_levels]
}

strip_fills <- unname(pal[group_levels])
strip_theme <- ggh4x::strip_themed(
    background_x = lapply(strip_fills, function(col) {
        element_rect(fill = col, colour = NA)
    }),
    text_x = element_text(face = "bold", colour = "#2b2b2b", size = 10)
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
lab_df <- df[nzchar(df$`.label`), , drop = FALSE]

p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(
        aes(size = .data[[s_col]], fill = .data[[g_col]]),
        shape = 21,
        colour = "white",
        stroke = as.numeric(config$point_stroke),
        alpha = 0.92
    ) +
    ggh4x::facet_wrap2(
        vars(.data[[g_col]]),
        nrow = as.integer(config$facet_nrow),
        scales = "fixed",
        strip = strip_theme
    ) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_size_continuous(
        range = as.numeric(config$size_range),
        guide = "none"
    ) +
    labs(x = config$labels$x, y = config$labels$y) +
    theme_classic(base_size = 10) +
    theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "#d0d0d0", linewidth = 0.4),
        strip.placement = "outside",
        axis.title = element_text(size = 10.5, colour = "#303030"),
        axis.text = element_text(size = 8.2, colour = "#4a4a4a"),
        plot.margin = margin(10, 12, 10, 12)
    )

if (nrow(lab_df)) {
    p <- p +
        geom_text_repel(
            data = lab_df,
            aes(label = .data$`.label`),
            size = as.numeric(config$label_size),
            colour = "#363636",
            lineheight = 0.9,
            box.padding = 0.22,
            point.padding = 0.15,
            min.segment.length = 0,
            segment.colour = "grey50",
            segment.size = 0.3,
            max.overlaps = Inf,
            seed = 1,
            show.legend = FALSE
        )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
n_panel <- length(group_levels)
nrow_f <- as.integer(config$facet_nrow)
ncol_f <- ceiling(n_panel / nrow_f)
save_ggplot(
    p,
    io$output,
    width = max(7.5, 3.6 * ncol_f + 1.2),
    height = max(5.5, 3.2 * nrow_f + 0.8)
)
