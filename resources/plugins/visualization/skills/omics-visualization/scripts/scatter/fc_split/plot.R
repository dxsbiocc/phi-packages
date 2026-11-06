#!/usr/bin/env Rscript

# Template-ID: scatter-fc-split
#
# Purpose:
#   Dual-contrast fold-change scatter: each gene sits at (log2FC in
#   contrast X, log2FC in contrast Y) as a disc split into two equal
#   halves, right half filled by -log10(padj) of X and left half by
#   -log10(padj) of Y, with zero and identity reference lines, a glyph
#   key, and repelled labels (optionally highlighted).
#
# Inputs:
#   Long table, one row per gene x contrast. Default example:
#     - gene: feature id
#     - contrast: one of the two contrasts named in CONFIG
#     - log2FC: effect size
#     - padj: adjusted p-value
#     - label: optional display label (empty = no label)
#     - highlight: optional TRUE/FALSE (bold accent label)
#
# Output:
#   A PDF, PNG, or SVG dual-contrast split-disc scatter.
#
# Dependencies:
#   ggplot2, readr, ggforce, ggrepel, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, which contrast is x / y, disc radius,
#   palette, label mode, axis titles). Do not run FindAllMarkers, DE
#   tests or p-value adjustment in this script. A single-contrast
#   volcano is scatter-volcano; several contrasts side by side are
#   scatter-volcano-multi; raw expression vs expression is
#   scatter-diagonal.
#
# Scientific assumptions:
#   Both contrasts are supplied. Only genes present in both are drawn
#   (the count dropped is reported). Halves are always equal; sector
#   angle never encodes p-values. padj = 0 is drawn at the largest
#   finite -log10(padj) of the table.

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
        gene = "gene",
        contrast = "contrast",
        effect = "log2FC",
        padj = "padj",
        label = "label",         # NULL when absent
        highlight = "highlight"  # NULL when absent
    ),
    x_contrast = "1st trimester",  # right half
    y_contrast = "term",           # left half
    radius = NULL,               # disc radius in log2FC units; NULL = auto
    palette = "Gradient.Reds",
    fill_limits = NULL,          # NULL = data range of -log10(padj)
    label_mode = "column",       # "column" or "top" (by |log2FC| sum)
    top_n = 20,
    accent = "#d7191c",
    label_size = 3,
    labels = list(
        x = "Log2(Lm/control) in 1st trimester",
        y = "Log2(Lm/control) in term",
        fill = "-log10(P value)"
    ),
    size = list(
        width = 8,
        height = 8
    )
)

load_packages(c("ggplot2", "readr", "ggforce", "ggrepel", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
cols <- config$columns
require_columns(df, cols[!vapply(cols, is.null, logical(1))])
g_col <- cols$gene
c_col <- cols$contrast
df[[g_col]] <- as.character(df[[g_col]])
df[[c_col]] <- as.character(df[[c_col]])
df[[cols$effect]] <- as.numeric(df[[cols$effect]])
df[[cols$padj]] <- as.numeric(df[[cols$padj]])
if (any(df[[cols$padj]] < 0 | df[[cols$padj]] > 1, na.rm = TRUE)) {
    stop("padj must be in [0, 1].", call. = FALSE)
}
cx <- config$x_contrast
cy <- config$y_contrast
miss <- setdiff(c(cx, cy), df[[c_col]])
if (length(miss)) {
    stop(sprintf("Contrast(s) not in table: %s", paste(miss, collapse = ", ")),
        call. = FALSE)
}
df$nlp <- -log10(df[[cols$padj]])
fin <- is.finite(df$nlp)
df$nlp[!fin & !is.na(df$nlp)] <- max(df$nlp[fin])

pick <- function(ct) {
    d <- df[df[[c_col]] == ct, , drop = FALSE]
    d[!duplicated(d[[g_col]]), , drop = FALSE]
}
dx <- pick(cx)
dy <- pick(cy)
genes <- intersect(dx[[g_col]], dy[[g_col]])
n_drop <- length(union(dx[[g_col]], dy[[g_col]])) - length(genes)
if (n_drop > 0) {
    message(sprintf("Dropped %d gene(s) present in only one contrast.", n_drop))
}
if (!length(genes)) stop("No genes shared by both contrasts.", call. = FALSE)
dx <- dx[match(genes, dx[[g_col]]), , drop = FALSE]
dy <- dy[match(genes, dy[[g_col]]), , drop = FALSE]
pts <- data.frame(gene = genes, x = dx[[cols$effect]], y = dy[[cols$effect]],
    nlp_x = dx$nlp, nlp_y = dy$nlp, stringsAsFactors = FALSE)
pts <- pts[is.finite(pts$x) & is.finite(pts$y), , drop = FALSE]

span <- max(diff(range(pts$x)), diff(range(pts$y)))
rad <- if (is.null(config$radius)) 0.022 * span else config$radius
halves <- rbind(
    data.frame(gene = pts$gene, x0 = pts$x, y0 = pts$y, start = 0,
        end = pi, nlp = pts$nlp_x),
    data.frame(gene = pts$gene, x0 = pts$x, y0 = pts$y, start = pi,
        end = 2 * pi, nlp = pts$nlp_y)
)

if (identical(config$label_mode, "top")) {
    lab <- pts[utils::head(order(-(abs(pts$x) + abs(pts$y))), config$top_n), ,
        drop = FALSE]
    lab$text <- lab$gene
} else {
    if (is.null(cols$label)) stop("label_mode = \"column\" needs a label column.",
        call. = FALSE)
    lv <- as.character(dx[[cols$label]][match(pts$gene, dx[[g_col]])])
    keep <- !is.na(lv) & nzchar(lv)
    lab <- pts[keep, , drop = FALSE]
    lab$text <- lv[keep]
}
lab$hl <- FALSE
if (!is.null(cols$highlight)) {
    hv <- dx[[cols$highlight]][match(lab$gene, dx[[g_col]])]
    lab$hl <- toupper(as.character(hv)) %in% c("TRUE", "T", "1", "YES")
}

xr <- range(pts$x) + c(-1, 1) * (rad * 3)
yr <- range(pts$y) + c(-1, 1) * (rad * 3)
key_x <- xr[[2]] - 0.27 * diff(xr)
key_y <- yr[[1]] + 0.2 * diff(yr)
key_r <- rad * 1.15

# disc radius in mm on the page, so labels are pushed clear of the glyph
panel_mm <- min(config$size$width, config$size$height) * 25.4 * 0.8
rad_mm <- rad / max(diff(xr), diff(yr)) * panel_mm

# every disc is a repel obstacle; unlabelled ones carry an empty label
rep_df <- pts
rep_df$text <- ""
rep_df$hl <- FALSE
rep_df$text[match(lab$gene, rep_df$gene)] <- lab$text
rep_df$hl[match(lab$gene, rep_df$gene)] <- lab$hl

stops <- palette_colors(config$palette)
fill_lims <- if (is.null(config$fill_limits)) range(halves$nlp) else
    config$fill_limits

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4,
        colour = "grey30") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4,
        colour = "grey30") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
        linewidth = 0.4, colour = "grey30") +
    ggforce::geom_arc_bar(
        data = halves,
        aes(x0 = x0, y0 = y0, r0 = 0, r = rad, start = start, end = end,
            fill = nlp),
        colour = "grey25", linewidth = 0.3
    ) +
    ggrepel::geom_text_repel(
        data = rep_df,
        aes(x = x, y = y, label = text,
            colour = ifelse(hl, config$accent, "grey10"),
            fontface = ifelse(hl, "bold.italic", "italic")),
        size = ifelse(rep_df$hl, config$label_size * 1.15, config$label_size),
        point.size = rad_mm * 2 * 1.35, point.padding = grid::unit(0.3, "mm"),
        box.padding = 0.25, force = 2, force_pull = 0.5,
        min.segment.length = 0, segment.size = 0.3, segment.colour = "grey30",
        max.overlaps = Inf, seed = 1
    ) +
    ggforce::geom_arc_bar(
        data = data.frame(start = c(0, pi), end = c(pi, 2 * pi)),
        aes(x0 = key_x, y0 = key_y, r0 = 0, r = key_r, start = start,
            end = end),
        fill = "white", colour = "grey10", linewidth = 0.4
    ) +
    annotate("segment",
        x = key_x + c(-1, 1) * key_r * 1.1, xend = key_x + c(-1, 1) * key_r * 2.6,
        y = key_y, yend = key_y, linewidth = 0.4,
        arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed")) +
    annotate("text", x = key_x - key_r * 2.8, y = key_y, label = cy,
        hjust = 1, size = 3.6) +
    annotate("text", x = key_x + key_r * 2.8, y = key_y, label = cx,
        hjust = 0, size = 3.6) +
    scale_colour_identity() +
    scale_fill_gradientn(colours = stops, limits = fill_lims,
        oob = scales::squish, name = config$labels$fill) +
    coord_fixed(xlim = xr, ylim = yr, expand = FALSE, clip = "off") +
    scale_x_continuous(breaks = scales::breaks_width(2)) +
    scale_y_continuous(breaks = scales::breaks_width(2)) +
    labs(x = config$labels$x, y = config$labels$y) +
    guides(fill = guide_colourbar(title.position = "bottom",
        title.hjust = 0.5, frame.colour = "grey20",
        barwidth = grid::unit(42, "mm"), barheight = grid::unit(4, "mm"))) +
    theme_prism(base_size = 11, border = TRUE) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.03),
        legend.justification = c(1, 0),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill = "white", colour = NA),
        axis.title = element_text(size = 12, face = "plain"),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width,
    height = config$size$height)
