#!/usr/bin/env Rscript

# Template-ID: scatter-volcano-barcode
#
# Purpose:
#   Volcano with gene-set context: shaded significant regions on each
#   side, optional dotted boxes around a supplied or ranked signature
#   (e.g. Top 200), gene-set members highlighted as filled points, a
#   barcode track per gene set under the x axis with supplied
#   enrichment statistics, and labelled key genes.
#
# Inputs:
#   de.tsv - one row per gene:
#     - gene, log2FC, padj
#     - label: optional display label (empty = no label)
#     - signature: optional up/down membership for the dotted boxes
#   genesets.tsv - one row per set member: set, gene (header-only = none)
#   Optional sidecar beside de.tsv: set_stats.tsv with set, ES, pval
#
# Output:
#   A PDF, PNG, or SVG volcano with gene-set barcodes.
#
# Dependencies:
#   ggplot2, readr, ggrepel
#
# Example:
#   Rscript plot.R example.tsv genesets.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, thresholds, which layers are on, region
#   and set labels, palettes, axis titles). Regions on/off, signature
#   boxes on/off and the barcode track on/off are CONFIG, not extra
#   ids. Do not run DE, GSEA or a gene-set test in this script; ES and
#   p-values come from set_stats.tsv.
#
# Scientific assumptions:
#   log2FC, padj, set membership and set statistics are supplied. The
#   shaded regions span genes passing |log2FC| > fc_cutoff and
#   padj < alpha. When no signature column is given, the boxes show the
#   top_n genes by log2FC per side among padj < signature_padj, which is
#   a display selection, not a new test.

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

io <- parse_io_args(n_input = 2L, input_names = c("de", "genesets"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        gene = "gene",
        effect = "log2FC",
        padj = "padj",
        label = "label",         # NULL when absent
        signature = NULL         # optional column with "up" / "down"
    ),
    set_columns = list(set = "set", gene = "gene"),
    stats_sidecar = "set_stats.tsv",
    fc_cutoff = 0.5,
    alpha = 0.05,
    show_regions = TRUE,
    show_signature = TRUE,
    signature_top_n = 200,
    signature_padj = 0.01,
    signature_label = "Top 200",
    show_barcode = TRUE,
    highlight_sets = TRUE,       # fill set members in the volcano
    sets = NULL,                 # NULL = first appearance
    set_palette = "Qualitative.Bold",
    set_index = c(4, 11),        # palette positions, one per set
    region_palette = "Diverging.RdBu",
    region_index = c(neg = 8, pos = 2),
    region_labels = c(
        neg = "PDAC KRAS UP\n(KRAS-dependent)",
        pos = "PDAC KRAS DOWN\n(KRAS-inhibited)"
    ),
    key_fill = "#88ccee",
    xlim = NULL,                 # NULL = symmetric, from the data
    labels = list(
        x = "KRAS vs NS siRNA (log2FC)",
        y = "Significance (-log10 adj. p-val.)"
    ),
    size = list(
        width = 7,
        height = 6
    )
)

load_packages(c("ggplot2", "readr", "ggrepel"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
de <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
cols <- config$columns
require_columns(de, cols[!vapply(cols, is.null, logical(1))])
g_col <- cols$gene
de[[g_col]] <- as.character(de[[g_col]])
de$fc <- as.numeric(de[[cols$effect]])
de$p <- as.numeric(de[[cols$padj]])
de <- de[is.finite(de$fc) & is.finite(de$p), , drop = FALSE]
if (any(de$p < 0 | de$p > 1)) stop("padj must be in [0, 1].", call. = FALSE)
de$nlp <- -log10(de$p)
fin <- is.finite(de$nlp)
de$nlp[!fin] <- max(de$nlp[fin])
ymax <- max(de$nlp)

gs <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
sc <- config$set_columns
if (nrow(gs)) {
    require_columns(gs, sc)
    gs <- data.frame(set = as.character(gs[[sc$set]]),
        gene = as.character(gs[[sc$gene]]), stringsAsFactors = FALSE)
    gs <- gs[gs$gene %in% de[[g_col]], , drop = FALSE]
} else {
    gs <- data.frame(set = character(0), gene = character(0))
}
sets <- if (is.null(config$sets)) unique(gs$set) else
    intersect(as.character(config$sets), gs$set)
set_pal <- palette_colors(config$set_palette)
set_col <- setNames(set_pal[config$set_index][seq_along(sets)], sets)
if (anyNA(set_col)) {
    set_col[is.na(set_col)] <- palette_colors("Qualitative.Safe",
        n = length(sets))[is.na(set_col)]
}

stats_path <- file.path(dirname(io$input[[1]]), config$stats_sidecar)
stats <- if (file.exists(stats_path)) {
    as.data.frame(read_table_auto(stats_path), stringsAsFactors = FALSE)
} else {
    NULL
}

# significant regions
sig_neg <- de[de$fc < -config$fc_cutoff & de$p < config$alpha, , drop = FALSE]
sig_pos <- de[de$fc > config$fc_cutoff & de$p < config$alpha, , drop = FALSE]
region_pal <- palette_colors(config$region_palette)
reg_col <- c(neg = region_pal[[config$region_index[["neg"]]]],
    pos = region_pal[[config$region_index[["pos"]]]])

# signature boxes
sig_box <- function(side) {
    if (!is.null(cols$signature)) {
        want <- if (side == "neg") "down" else "up"
        d <- de[tolower(as.character(de[[cols$signature]])) == want, ,
            drop = FALSE]
    } else {
        d <- de[de$p < config$signature_padj, , drop = FALSE]
        d <- if (side == "neg") d[order(d$fc), , drop = FALSE] else
            d[order(-d$fc), , drop = FALSE]
        d <- utils::head(d, config$signature_top_n)
    }
    if (!nrow(d)) return(NULL)
    data.frame(xmin = min(d$fc), xmax = max(d$fc), ymin = min(d$nlp),
        ymax = ymax, side = side)
}
boxes <- do.call(rbind, lapply(c("neg", "pos"), sig_box))

# point layers
member <- gs[gs$set %in% sets, , drop = FALSE]
member <- member[!duplicated(member$gene), , drop = FALSE]
de$set <- member$set[match(de[[g_col]], member$gene)]
keys <- if (!is.null(cols$label)) {
    lv <- as.character(de[[cols$label]])
    de[!is.na(lv) & nzchar(lv), , drop = FALSE]
} else {
    de[0, , drop = FALSE]
}
if (nrow(keys)) keys$text <- as.character(keys[[cols$label]])

xl <- if (is.null(config$xlim)) {
    m <- ceiling(max(abs(de$fc)) * 2) / 2
    c(-m, m)
} else {
    config$xlim
}

# barcode rows under y = 0
bar_on <- isTRUE(config$show_barcode) && length(sets) > 0
row_h <- 0.09 * ymax
bars <- if (bar_on) {
    do.call(rbind, lapply(seq_along(sets), function(i) {
        g <- gs$gene[gs$set == sets[[i]]]
        d <- de[de[[g_col]] %in% g, , drop = FALSE]
        top <- -row_h * (0.4 + (i - 1) * 1.15)
        data.frame(x = d$fc, ymin = top - row_h, ymax = top, set = sets[[i]])
    }))
} else {
    NULL
}
y_lo <- if (bar_on) -row_h * (0.4 + length(sets) * 1.15) else -0.03 * ymax
y_hi <- ymax * (if (isTRUE(config$show_regions)) 1.16 else 1.05)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(de, aes(x = fc, y = nlp))
if (isTRUE(config$show_regions)) {
    if (nrow(sig_neg)) {
        p <- p + annotate("rect", xmin = min(sig_neg$fc), xmax = max(sig_neg$fc),
            ymin = -log10(config$alpha), ymax = ymax, fill = reg_col[["neg"]],
            alpha = 0.18)
    }
    if (nrow(sig_pos)) {
        p <- p + annotate("rect", xmin = min(sig_pos$fc), xmax = max(sig_pos$fc),
            ymin = -log10(config$alpha), ymax = ymax, fill = reg_col[["pos"]],
            alpha = 0.18)
    }
    p <- p + annotate("text", x = xl * 0.62, y = ymax * 1.02,
        label = unname(config$region_labels[c("neg", "pos")]),
        colour = unname(reg_col[c("neg", "pos")]), size = 4.6,
        lineheight = 0.9, vjust = 0)
}
p <- p +
    geom_point(data = de[is.na(de$set) | !isTRUE(config$highlight_sets), ,
        drop = FALSE], shape = 21, colour = "black", fill = NA, alpha = 0.18,
        size = 1.3, stroke = 0.5)
if (isTRUE(config$highlight_sets) && length(sets)) {
    hl <- de[!is.na(de$set), , drop = FALSE]
    p <- p + geom_point(data = hl, aes(fill = set), shape = 21,
        colour = "black", size = 2, stroke = 0.6)
}
if (isTRUE(config$show_signature) && !is.null(boxes)) {
    p <- p +
        geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin,
            ymax = ymax), inherit.aes = FALSE, fill = NA, colour = "black",
            linetype = "dotted", linewidth = 0.6) +
        geom_text(data = boxes, aes(x = ifelse(side == "neg", xmin, xmax),
            y = ymax, label = config$signature_label,
            hjust = ifelse(side == "neg", -0.05, 1.05)), vjust = 1.3,
            inherit.aes = FALSE, size = 4.6)
}
if (nrow(keys)) {
    p <- p +
        geom_point(data = keys, shape = 21, colour = "black",
            fill = config$key_fill, size = 2.2, stroke = 0.7) +
        ggrepel::geom_label_repel(data = keys, aes(label = text), size = 4.5,
            fill = "white", alpha = 0.9, label.size = 0.3,
            box.padding = 0.4, point.padding = 0.4, min.segment.length = 0,
            segment.colour = "grey45", segment.size = 0.4, seed = 1)
}
p <- p +
    annotate("segment", x = 0, xend = 0, y = 0, yend = y_hi, colour = "grey60",
        linewidth = 0.5) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.5) +
    geom_hline(yintercept = -log10(config$alpha), colour = reg_col[["pos"]],
        linewidth = 0.6, linetype = "dotted") +
    annotate("text", x = xl[[2]], y = -log10(config$alpha),
        label = paste0("\u03b1 = ", config$alpha), colour = reg_col[["pos"]],
        hjust = 1, vjust = 1.4, size = 4)
if (bar_on) {
    p <- p + geom_linerange(data = bars, aes(x = x, ymin = ymin, ymax = ymax,
        colour = set), inherit.aes = FALSE, linewidth = 0.35, alpha = 0.8)
    lab_y <- -row_h * (0.4 + (seq_along(sets) - 1) * 1.15) - row_h / 2
    p <- p + annotate("text", x = xl[[1]], y = lab_y, label = sets,
        colour = unname(set_col), hjust = 0, size = 4)
    if (!is.null(stats)) {
        st <- stats[match(sets, as.character(stats$set)), , drop = FALSE]
        ok <- !is.na(st$set)
        if (any(ok)) {
            p <- p + annotate("text", x = xl[[2]], y = lab_y[ok],
                label = sprintf("ES: %.2f\npval: %s", st$ES[ok],
                    formatC(st$pval[ok], format = "g", digits = 2)),
                colour = unname(set_col[ok]), hjust = 1, size = 3.8,
                lineheight = 0.85)
        }
    }
}
p <- p +
    scale_fill_manual(values = set_col, guide = "none") +
    scale_colour_manual(values = set_col, guide = "none") +
    scale_x_continuous(breaks = seq(ceiling(xl[[1]]), floor(xl[[2]]), 1)) +
    scale_y_continuous(breaks = pretty(c(0, ymax))[pretty(c(0, ymax)) <= ymax]) +
    coord_cartesian(xlim = xl, ylim = c(y_lo, y_hi), expand = FALSE,
        clip = "off") +
    labs(x = config$labels$x, y = config$labels$y) +
    theme_classic(base_size = 13) +
    theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour = "black"),
        plot.background = element_blank(),
        plot.margin = margin(10, 14, 8, 8)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width,
    height = config$size$height)
