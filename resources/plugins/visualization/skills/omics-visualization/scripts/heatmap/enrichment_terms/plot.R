#!/usr/bin/env Rscript

# Template-ID: heatmap-enrichment-terms
#
# Purpose:
#   Draw a term-by-group enrichment heatmap. Fill is signed
#   -log10(adjusted p); missing cells (term not reported in that
#   group) are marked rather than imputed.
#
# Inputs:
#   One row per term-by-group cell. Default example:
#     - Description: term name
#     - group: contrast or cohort
#     - p.adjust: adjusted p-value (NA if not enriched)
#     - direction: up / down (maps the sign of the fill)
#
# Output:
#   A PDF, PNG, or SVG enrichment-term heatmap.
#
# Dependencies:
#   ggplot2, ggprism, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, palettes, marks).
#   Extra contrasts are extra group levels, not extra ids. Do not
#   run enricher, convert species, or cluster terms in this script.
#
# Scientific assumptions:
#   p.adjust and direction are supplied enrichment results. Signed
#   -log10(p.adjust) is a display transform. Expanding to the full
#   term × group grid only shows combinations absent from the table;
#   it does not test those combinations.

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
        term = "Description",
        group = "group",
        pvalue = "p.adjust",
        direction = "direction"
    ),
    palettes = list(
        fill = "Diverging.TealRose"
    ),
    mark_missing = TRUE,
    cell_ratio = 0.72,
    labels = list(
        fill = "signed -log10(p.adjust)",
        y_size = 13,
        x_size = 10
    ),
    size = list(
        width = 7.2,
        height = 5.4
    )
)

load_packages(c("ggplot2", "ggprism", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
df <- as.data.frame(df, stringsAsFactors = FALSE)

term_col <- config$columns$term
group_col <- config$columns$group
p_col <- config$columns$pvalue
dir_col <- config$columns$direction

df[[term_col]] <- as.character(df[[term_col]])
df[[group_col]] <- as.character(df[[group_col]])
df[[p_col]] <- as.numeric(df[[p_col]])
df[[dir_col]] <- tolower(as.character(df[[dir_col]]))
if (any(!nzchar(df[[term_col]]) | !nzchar(df[[group_col]]))) {
    stop("term and group must be non-empty.", call. = FALSE)
}

terms <- unique(df[[term_col]])
groups <- unique(df[[group_col]])
grid <- expand.grid(
    term = terms,
    group = groups,
    stringsAsFactors = FALSE
)
names(grid) <- c(term_col, group_col)
df <- merge(grid, df, by = c(term_col, group_col), all.x = TRUE)

sign_val <- ifelse(
    df[[dir_col]] %in% c("up", "positive", "1", "+"),
    1,
    ifelse(df[[dir_col]] %in% c("down", "negative", "-1", "-"), -1, NA_real_)
)
finite_p <- is.finite(df[[p_col]]) & df[[p_col]] > 0
df$score <- ifelse(finite_p, -log10(df[[p_col]]) * sign_val, NA_real_)
df$missing <- !is.finite(df$score)

ord <- tapply(abs(df$score), df[[term_col]], function(x) {
    m <- max(x, na.rm = TRUE)
    if (!is.finite(m)) 0 else m
})
term_lv <- names(sort(ord, decreasing = TRUE))
df[[term_col]] <- factor(df[[term_col]], levels = rev(term_lv))
df[[group_col]] <- factor(df[[group_col]], levels = groups)

fill_cols <- palette_colors(config$palettes$fill)
# TealRose: teal (down) → cream → rose (up). No reverse.
miss_col <- palette_colors("Brand.Emma")[[1]]
outline <- palette_colors("Qualitative.Safe")[[12]]
lim <- max(abs(df$score), na.rm = TRUE)
if (!is.finite(lim) || lim == 0) {
    lim <- 1
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = .data[[group_col]], y = .data[[term_col]])) +
    geom_tile(
        aes(fill = score),
        colour = outline,
        linewidth = 0.35
    ) +
    scale_fill_gradientn(
        name = config$labels$fill,
        colours = fill_cols,
        limits = c(-lim, lim),
        na.value = fill_cols[[ceiling(length(fill_cols) / 2)]]
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed(ratio = config$cell_ratio, clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        panel.border = element_rect(
            colour = outline,
            fill = NA,
            linewidth = 0.4
        ),
        axis.text.x = element_text(
            angle = 45, hjust = 1, vjust = 1,
            size = config$labels$x_size
        ),
        axis.text.y = element_text(size = config$labels$y_size),
        plot.margin = margin(6, 8, 6, 6),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank()
    )

if (isTRUE(config$mark_missing) && any(df$missing)) {
    miss <- df[df$missing, , drop = FALSE]
    pad <- 0.08
    miss$x0 <- as.numeric(miss[[group_col]]) - 0.5 + pad
    miss$x1 <- as.numeric(miss[[group_col]]) + 0.5 - pad
    miss$y0 <- as.numeric(miss[[term_col]]) - 0.5 + pad
    miss$y1 <- as.numeric(miss[[term_col]]) + 0.5 - pad
    p <- p +
        geom_segment(
            data = miss,
            aes(x = x0, xend = x1, y = y0, yend = y1),
            inherit.aes = FALSE,
            colour = miss_col,
            linewidth = 0.45,
            lineend = "butt"
        ) +
        geom_segment(
            data = miss,
            aes(x = x0, xend = x1, y = y1, yend = y0),
            inherit.aes = FALSE,
            colour = miss_col,
            linewidth = 0.45,
            lineend = "butt"
        )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p, io$output,
    width = config$size$width,
    height = config$size$height
)
