#!/usr/bin/env Rscript

# Template-ID: scatter-volcano-multi
#
# Purpose:
#   Multi-contrast volcano: one panel per contrast side by side with a
#   shared log2FC y axis and a signed -log10(padj) x axis (up-regulated
#   to the right, down-regulated to the left), contrast-coloured points
#   beyond a fold-change threshold, dashed threshold lines, and
#   repelled gene labels pulled toward the panel interior.
#
# Inputs:
#   One row per gene within a contrast. Default example:
#     - contrast: panel (e.g. developmental stage)
#     - gene: feature id
#     - log2FC: effect size
#     - padj: adjusted p-value
#     - label: optional display label (empty = no label)
#
# Output:
#   A PDF, PNG, or SVG multi-panel volcano.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggrepel
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, contrast order, layout, threshold, label
#   mode, palette, sizes). layout = "polar" draws one sector per
#   contrast (angle = log2FC, radius = -log10 P); layout = "strip" draws
#   framed top-strip panels with unsigned -log10(padj) on a shared x and
#   only labelled genes coloured. Neither is a new id.
#   Extra contrasts are extra rows, not extra ids. A single-contrast
#   volcano is scatter-volcano. Do not run a differential test or
#   p-value adjustment in this script.
#
# Scientific assumptions:
#   log2FC and padj are supplied. x = sign(log2FC) x -log10(padj) is a
#   display transform. padj = 0 is drawn at the panel's largest finite
#   -log10(padj) (not a measured value). The fold-change threshold only
#   greys points and places the dashed lines; it is not a test. In the
#   polar layout, rows outside fc_range are dropped with a message.

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
        contrast = "contrast",
        gene = "gene",
        effect = "log2FC",
        padj = "padj",
        label = "label"          # set NULL when there is no label column
    ),
    contrasts = NULL,            # NULL keeps first-appearance order
    layout = "facet",            # "facet", "strip" or "polar"
    # strip: framed panels with top strips, x = unsigned -log10(padj) on a
    # shared range, grey points with only labelled genes in panel colour.
    # polar: one sector per contrast; angle within sector = log2FC in
    # fc_range, radius = -log10(padj). Rows outside fc_range are dropped
    # (counted in a message), not squished.
    polar = list(
        fc_range = c(0, 14),
        r_max = NULL,            # NULL = pretty ceiling of the data
        colour = NULL,           # continuous colour column, e.g. "satt_index"
        colour_palette = "Gradient.YlGnBu",
        colour_reverse = TRUE,
        colour_limits = NULL,    # NULL = data range
        colour_title = NULL,
        sector_alpha = 0.14,
        inner_radius = 0.22,     # empty hole as a fraction of the radius
        sector_pad = 0.05,       # angular margin inside each sector (0-0.5)
        point_alpha = 0.85,
        tag = NULL               # corner label, e.g. "SUMO2\nSATTs"
    ),
    fc_threshold = 0.25,
    label_mode = "column",       # "column" (non-empty label) or "top"
    top_n = 5,                   # per side per panel when label_mode = "top"
    palette = "Qualitative.Set2 / 7",
    neutral = "#dcdcdc",
    point_size = 1.1,
    label_size = 2.7,
    label_italic = TRUE,
    labels = list(
        x = "-log10 (adjusted P value)",
        y = "Atrium vs ventricle (average log2FC)"
    ),
    size = list(
        width = 9,
        height = 6.2
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggrepel"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
cols <- config$columns
require_columns(df, cols[!vapply(cols, is.null, logical(1))])
ct_col <- cols$contrast
g_col <- cols$gene
e_col <- cols$effect
p_col <- cols$padj

df[[ct_col]] <- as.character(df[[ct_col]])
df[[g_col]] <- as.character(df[[g_col]])
df[[e_col]] <- as.numeric(df[[e_col]])
df[[p_col]] <- as.numeric(df[[p_col]])
df <- df[is.finite(df[[e_col]]) & is.finite(df[[p_col]]), , drop = FALSE]
if (any(df[[p_col]] < 0 | df[[p_col]] > 1)) {
    stop("padj must be in [0, 1].", call. = FALSE)
}
contrasts <- if (is.null(config$contrasts)) unique(df[[ct_col]]) else
    as.character(config$contrasts)
df <- df[df[[ct_col]] %in% contrasts, , drop = FALSE]
contrasts <- contrasts[contrasts %in% df[[ct_col]]]
if (!nrow(df)) stop("No rows left.", call. = FALSE)

df$nlp <- -log10(df[[p_col]])
cap <- tapply(df$nlp[is.finite(df$nlp)], df[[ct_col]][is.finite(df$nlp)], max)
inf <- !is.finite(df$nlp)
df$nlp[inf] <- cap[df[[ct_col]][inf]]
df$x <- sign(df[[e_col]]) * df$nlp
df$panel <- factor(df[[ct_col]], levels = contrasts)

pal <- setNames(palette_colors(config$palette, n = length(contrasts)), contrasts)
df$col <- ifelse(abs(df[[e_col]]) >= config$fc_threshold,
    unname(pal[df[[ct_col]]]), config$neutral)
df <- df[order(df$col != config$neutral), , drop = FALSE]

if (identical(config$label_mode, "top")) {
    lab <- do.call(rbind, lapply(split(df, list(df$panel, df[[e_col]] > 0),
        drop = TRUE), function(d) {
        d[utils::head(order(-abs(d[[e_col]]), -d$nlp), config$top_n), ,
            drop = FALSE]
    }))
    lab$text <- lab[[g_col]]
} else {
    l_col <- cols$label
    if (is.null(l_col)) stop("label_mode = \"column\" needs a label column.",
        call. = FALSE)
    lv <- as.character(df[[l_col]])
    keep <- !is.na(lv) & nzchar(lv)
    lab <- df[keep, , drop = FALSE]
    lab$text <- lv[keep]
}

# symmetric x range per panel; labels pulled toward the panel interior
span <- tapply(abs(df$x), df$panel, max)
span <- span[contrasts]
edge <- data.frame(
    panel = factor(rep(contrasts, each = 2), levels = contrasts),
    x = as.vector(rbind(-span, span)),
    y = 0
)
lab$nudge <- -sign(lab$x) * 0.55 * unname(span[as.character(lab$panel)])
lab$hjust <- ifelse(lab$x > 0, 1, 0)

sym_breaks <- function(lims) {
    m <- max(abs(lims))
    c(-m, 0, m)
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_polar <- function() {
    pc <- config$polar
    lo <- pc$fc_range[[1]]
    hi <- pc$fc_range[[2]]
    inside <- df[[e_col]] >= lo & df[[e_col]] <= hi
    if (any(!inside)) {
        message(sprintf("Dropped %d row(s) with %s outside [%g, %g].",
            sum(!inside), e_col, lo, hi))
    }
    d <- df[inside, , drop = FALSE]
    pad <- pc$sector_pad
    to_theta <- function(k, fc) (k - 1) + pad + (fc - lo) / (hi - lo) * (1 - 2 * pad)
    d$theta <- to_theta(as.integer(d$panel), d[[e_col]])
    r_max <- if (is.null(pc$r_max)) {
        pt <- pretty(c(0, max(d$nlp)))
        pt[[which(pt >= max(d$nlp))[[1]]]]
    } else {
        pc$r_max
    }
    d$nlp <- pmin(d$nlp, r_max)
    ld <- d[d[[g_col]] %in% lab[[g_col]] & paste(d$panel, d[[g_col]]) %in%
        paste(lab$panel, lab[[g_col]]), , drop = FALSE]
    ld$text <- lab$text[match(paste(ld$panel, ld[[g_col]]),
        paste(lab$panel, lab[[g_col]]))]
    n_ct <- length(contrasts)
    sectors <- data.frame(xmin = seq_len(n_ct) - 1, xmax = seq_len(n_ct),
        fill = unname(pal[contrasts]))

    p <- ggplot(d, aes(x = theta, y = nlp)) +
        geom_rect(data = sectors, aes(xmin = xmin, xmax = xmax, ymin = 0,
            ymax = r_max, fill = fill), inherit.aes = FALSE,
            alpha = pc$sector_alpha, colour = "grey55", linewidth = 0.25) +
        scale_fill_identity()
    if (is.null(pc$colour)) {
        p <- p + geom_point(aes(colour = col), size = config$point_size,
            alpha = pc$point_alpha, stroke = 0) +
            scale_colour_identity()
    } else {
        stops <- palette_colors(pc$colour_palette)
        if (isTRUE(pc$colour_reverse)) stops <- rev(stops)
        lims <- if (is.null(pc$colour_limits)) range(d[[pc$colour]],
            na.rm = TRUE) else pc$colour_limits
        p <- p + geom_point(aes(colour = .data[[pc$colour]]),
            size = config$point_size, alpha = pc$point_alpha, stroke = 0) +
            scale_colour_gradientn(colours = stops, limits = lims,
                oob = scales::squish,
                name = if (is.null(pc$colour_title)) pc$colour else
                    pc$colour_title) +
            guides(colour = guide_colourbar(frame.colour = "grey35"))
    }
    p <- p +
        geom_point(data = ld, shape = 1, size = config$point_size + 2.2,
            colour = "grey15", stroke = 0.5) +
        ggrepel::geom_text_repel(
            data = ld, aes(label = text),
            size = config$label_size, colour = "grey15",
            fontface = if (isTRUE(config$label_italic)) "italic" else "plain",
            min.segment.length = 0, segment.size = 0.3,
            box.padding = 0.4, point.padding = 0.3, force = 2,
            bg.color = "white", bg.r = 0.08, max.overlaps = Inf, seed = 1
        ) +
        annotate("segment", x = to_theta(1, lo), xend = to_theta(1, hi),
            y = r_max * 1.2, yend = r_max * 1.2, linewidth = 0.35,
            colour = "grey25") +
        annotate("text", x = to_theta(1, c(lo, hi)), y = r_max * 1.28,
            label = format(c(lo, hi)), size = 3) +
        annotate("text", x = 0.5, y = r_max * 1.42,
            label = "Fold change\n(log2)", size = 3, lineheight = 0.95) +
        scale_x_continuous(limits = c(0, n_ct), breaks = seq_len(n_ct) - 0.5,
            labels = contrasts, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, r_max),
            breaks = pretty(c(0, r_max), n = 3), expand = c(0, 0),
            oob = scales::oob_keep) +
        coord_radial(theta = "x", expand = FALSE, r.axis.inside = FALSE,
            inner.radius = pc$inner_radius, clip = "off") +
        guides(theta = guide_axis_theta(angle = 0),
            r = guide_axis(cap = "both")) +
        labs(x = NULL, y = expression(-log[10](italic(P))), tag = pc$tag) +
        theme_classic(base_size = 11) +
        theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(size = 11, colour = "grey15"),
            axis.text.y = element_text(size = 10, colour = "grey15"),
            plot.tag = element_text(size = 12, face = "bold", hjust = 0,
                lineheight = 0.95),
            plot.tag.position = c(0.02, 0.06),
            plot.background = element_blank(),
            panel.background = element_blank(),
            plot.margin = margin(40, 20, 15, 30)
        )
    p
}

draw_strip <- function() {
    xmax <- max(df$nlp)
    lb <- lab
    lb$col <- unname(pal[as.character(lb$panel)])
    # label toward the open side of the panel
    right <- lb$nlp < 0.55 * xmax
    lb$nudge <- ifelse(right, 0.3, -0.3) * xmax
    lb$hjust <- ifelse(right, 0, 1)
    ggplot(df, aes(x = nlp, y = .data[[e_col]])) +
        geom_hline(yintercept = c(-1, 1) * config$fc_threshold,
            linetype = "dashed", linewidth = 0.4, colour = "grey55") +
        geom_point(colour = "grey78", size = config$point_size, stroke = 0) +
        geom_point(data = lb, aes(colour = col),
            size = config$point_size + 0.6, stroke = 0) +
        ggrepel::geom_text_repel(
            data = lb,
            aes(label = text, colour = col, hjust = hjust),
            nudge_x = lb$nudge,
            direction = "y",
            size = config$label_size,
            fontface = if (isTRUE(config$label_italic)) "italic" else "plain",
            segment.size = 0.35,
            min.segment.length = 0,
            box.padding = 0.15,
            force = 1,
            force_pull = 3,
            max.overlaps = Inf,
            seed = 1
        ) +
        scale_colour_identity() +
        scale_x_continuous(limits = c(0, xmax * 1.02),
            breaks = scales::breaks_pretty(3),
            expand = expansion(mult = c(0.03, 0.03))) +
        scale_y_continuous(breaks = function(l) seq(ceiling(l[[1]]),
            floor(l[[2]]))) +
        facet_wrap(~panel, nrow = 1) +
        labs(x = config$labels$x, y = config$labels$y) +
        theme_bw(base_size = 10) +
        theme(
            panel.grid = element_blank(),
            panel.spacing.x = grid::unit(3, "pt"),
            strip.background = element_rect(fill = "grey88", colour = "black",
                linewidth = 0.4),
            strip.text = element_text(size = 9, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.title = element_text(size = 10, face = "bold"),
            plot.background = element_blank()
        )
}

p <- if (identical(config$layout, "polar")) draw_polar() else
    if (identical(config$layout, "strip")) draw_strip() else
    ggplot(df, aes(x = x, y = .data[[e_col]])) +
    geom_blank(data = edge, aes(x = x, y = y)) +
    geom_vline(xintercept = -Inf, linewidth = 0.5, colour = "black") +
    geom_hline(yintercept = c(-1, 1) * config$fc_threshold,
        linetype = "dashed", linewidth = 0.5, colour = "grey20") +
    geom_point(aes(colour = col), size = config$point_size, stroke = 0) +
    geom_point(data = lab, aes(fill = col), shape = 21, colour = "grey15",
        size = config$point_size + 0.4, stroke = 0.3) +
    ggrepel::geom_text_repel(
        data = lab,
        aes(label = text, hjust = hjust),
        nudge_x = lab$nudge,
        direction = "y",
        size = config$label_size,
        fontface = if (isTRUE(config$label_italic)) "italic" else "plain",
        segment.size = 0.35,
        segment.colour = "grey15",
        min.segment.length = 0,
        box.padding = 0.25,
        point.padding = 0.15,
        force = 2,
        max.overlaps = Inf,
        seed = 1
    ) +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_x_continuous(breaks = sym_breaks,
        labels = function(b) format(round(b), trim = TRUE),
        expand = expansion(mult = 0.04)) +
    scale_y_continuous(breaks = function(l) seq(ceiling(l[[1]]), floor(l[[2]]))) +
    facet_wrap(~panel, nrow = 1, scales = "free_x", strip.position = "bottom") +
    labs(x = config$labels$x, y = config$labels$y) +
    theme_prism(base_size = 10) +
    theme(
        panel.spacing.x = grid::unit(10, "pt"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "plain",
            margin = margin(t = 4, b = 2)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10, face = "plain"),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width,
    height = config$size$height)
