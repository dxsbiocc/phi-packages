#!/usr/bin/env Rscript

# Template-ID: sankey-enrichment
#
# Purpose:
#   Gene-to-term Sankey aligned with an enrichment panel. Left column:
#   genes; right column: enriched terms; ribbons: gene membership in a
#   term. A dot (or bar) panel on the right shares the term rows and
#   shows GeneRatio, Count, and -log10(p.adjust).
#
# Inputs:
#   One enrichment table (tsv/csv), e.g. clusterProfiler output:
#     - Description: term name
#     - GeneRatio: number or "k/n" string
#     - Count: genes in the term (full enrichment result)
#     - p.adjust: adjusted p-value
#     - geneID: genes separated by "/"
#
# Output:
#   A PDF, PNG, or SVG two-panel figure.
#
# Dependencies:
#   ggplot2, readr, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG: columns, top_n, max_genes_per_term, genes
#   (explicit gene list), right_panel ("dot" or "bar"), flow_colour
#   ("term" or "gene"), palettes, label widths.
#
# Scientific assumptions:
#   Enrichment statistics are supplied; the script does not run
#   enrichment or correct p-values. Terms are the top_n rows by
#   p.adjust. Genes are taken in geneID order (first
#   max_genes_per_term per term) or from CONFIG genes; never randomly
#   sampled. A term block's height counts the genes drawn in this
#   figure, not Count; Count is shown only in the right panel.

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
        gene_ratio = "GeneRatio",
        count = "Count",
        pvalue = "p.adjust",
        genes = "geneID"
    ),
    gene_sep = "/",
    top_n = 10L,
    max_genes_per_term = 6L,
    genes = NULL,               # character vector to restrict genes; NULL = all
    term_order = "gene_ratio",  # top-to-bottom: "gene_ratio" or "pvalue"
    right_panel = "dot",        # "dot" or "bar"
    flow_colour = "term",       # ribbons coloured by "term" or "gene"
    flow_alpha = 0.28,
    gap_fraction = 0.25,        # vertical gaps as a fraction of total flow height
    term_wrap = 30L,
    palettes = list(
        term = "Qualitative.Safe",
        gene = "Quantitative.BluGrn",
        stat = "Quantitative.Sunset"
    ),
    gene_block_colour = "grey35",
    labels = list(
        x = "GeneRatio",
        fill = "-log10(p.adjust)",
        size = "Count"
    ),
    size_range = c(2.5, 8),
    gene_label_size = 3.1,
    term_label_size = 3.4,
    widths = c(3.2, 1.25)
)

load_packages(c("ggplot2", "readr", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
cc <- config$columns
require_columns(df, cc)

parse_ratio <- function(x) {
    x <- as.character(x)
    out <- suppressWarnings(as.numeric(x))
    frac <- grepl("/", x, fixed = TRUE)
    if (any(frac)) {
        parts <- strsplit(x[frac], "/", fixed = TRUE)
        out[frac] <- vapply(parts, function(p) as.numeric(p[[1]]) / as.numeric(p[[2]]),
            numeric(1))
    }
    out
}
df$.term <- as.character(df[[cc$description]])
df$.ratio <- parse_ratio(df[[cc$gene_ratio]])
df$.count <- as.numeric(df[[cc$count]])
df$.p <- as.numeric(df[[cc$pvalue]])
df$.genes <- as.character(df[[cc$genes]])
ok <- nzchar(df$.term) & is.finite(df$.ratio) & is.finite(df$.count) &
    is.finite(df$.p) & df$.p > 0 & !is.na(df$.genes) & nzchar(df$.genes)
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with missing term, ratio, count, p, or genes.")
}
df <- df[ok, , drop = FALSE]
if (anyDuplicated(df$.term)) {
    stop("Description must be unique.", call. = FALSE)
}
df <- df[order(df$.p), , drop = FALSE]
df <- utils::head(df, as.integer(config$top_n))
if (!nrow(df)) {
    stop("No terms left to plot.", call. = FALSE)
}
df$.neglog <- -log10(df$.p)
term_ord <- if (identical(config$term_order, "pvalue")) {
    order(-df$.neglog, df$.term)
} else {
    order(-df$.ratio, df$.term)
}
df <- df[term_ord, , drop = FALSE]
terms <- df$.term

# Gene-term links: deterministic selection, never random.
links <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    g <- trimws(strsplit(df$.genes[[i]], config$gene_sep, fixed = TRUE)[[1]])
    g <- unique(g[nzchar(g)])
    if (!is.null(config$genes)) {
        g <- g[g %in% config$genes]
    } else {
        g <- utils::head(g, as.integer(config$max_genes_per_term))
    }
    if (!length(g)) return(NULL)
    data.frame(gene = g, term = df$.term[[i]], stringsAsFactors = FALSE)
}))
if (is.null(links) || !nrow(links)) {
    stop("No gene-term links after selection.", call. = FALSE)
}
# Within the selected genes, keep every membership among the shown terms.
sel_genes <- unique(links$gene)
links <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    g <- trimws(strsplit(df$.genes[[i]], config$gene_sep, fixed = TRUE)[[1]])
    g <- intersect(unique(g), sel_genes)
    if (!length(g)) return(NULL)
    data.frame(gene = g, term = df$.term[[i]], stringsAsFactors = FALSE)
}))
links$term_i <- match(links$term, terms)
terms <- terms[sort(unique(links$term_i))]
df <- df[df$.term %in% terms, , drop = FALSE]
links$term_i <- match(links$term, terms)

# Gene order: by the first (topmost) term each gene belongs to.
first_term <- tapply(links$term_i, links$gene, min)
genes <- names(first_term)[order(first_term, match(names(first_term), sel_genes))]
links$gene_i <- match(links$gene, genes)

# Vertical layout (top = largest y). Each link occupies one unit.
n_link <- nrow(links)
gap_total <- config$gap_fraction * n_link
H <- n_link + gap_total
stack <- function(sizes) {
    k <- length(sizes)
    gap <- if (k > 1) gap_total / (k - 1) else 0
    top <- H - c(0, cumsum(sizes + gap))[seq_len(k)]
    data.frame(ymax = top, ymin = top - sizes)
}
gene_n <- as.integer(table(factor(links$gene_i, levels = seq_along(genes))))
term_n <- as.integer(table(factor(links$term_i, levels = seq_along(terms))))
gene_blk <- cbind(stack(gene_n), name = genes, n = gene_n)
term_blk <- cbind(stack(term_n), name = terms, n = term_n)

# Slot each link inside its gene block (by term order) and term block
# (by gene order) so ribbons cross as little as possible.
links <- links[order(links$gene_i, links$term_i), , drop = FALSE]
links$g_slot <- stats::ave(links$term_i, links$gene_i, FUN = seq_along)
links <- links[order(links$term_i, links$gene_i), , drop = FALSE]
links$t_slot <- stats::ave(links$gene_i, links$term_i, FUN = seq_along)
links$y0_top <- gene_blk$ymax[links$gene_i] - (links$g_slot - 1)
links$y1_top <- term_blk$ymax[links$term_i] - (links$t_slot - 1)

x_gene <- c(0, 0.05)
x_term <- c(1, 1.05)
sigmoid <- function(a, b, t) a + (b - a) * (3 * t^2 - 2 * t^3)
tt <- seq(0, 1, length.out = 40)
xs <- x_gene[[2]] + (x_term[[1]] - x_gene[[2]]) * tt
ribbons <- do.call(rbind, lapply(seq_len(nrow(links)), function(k) {
    top <- sigmoid(links$y0_top[[k]], links$y1_top[[k]], tt)
    bot <- sigmoid(links$y0_top[[k]] - 1, links$y1_top[[k]] - 1, tt)
    data.frame(x = c(xs, rev(xs)), y = c(top, rev(bot)), id = k,
        term = links$term[[k]], gene = links$gene[[k]], stringsAsFactors = FALSE)
}))

term_cols <- palette_colors(config$palettes$term, n = length(terms))
names(term_cols) <- terms
gene_cols <- grDevices::colorRampPalette(palette_colors(config$palettes$gene))(length(genes))
names(gene_cols) <- genes
ribbons$col <- if (identical(config$flow_colour, "gene")) {
    unname(gene_cols[ribbons$gene])
} else {
    unname(term_cols[ribbons$term])
}
gene_blk$col <- if (identical(config$flow_colour, "gene")) {
    unname(gene_cols[gene_blk$name])
} else {
    config$gene_block_colour
}
term_blk$col <- unname(term_cols[term_blk$name])
term_blk$label <- vapply(term_blk$name, function(s) {
    paste(strwrap(s, width = config$term_wrap), collapse = "\n")
}, character(1))
term_blk$y <- (term_blk$ymin + term_blk$ymax) / 2
gene_blk$y <- (gene_blk$ymin + gene_blk$ymax) / 2

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p_sankey <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
        data = ribbons,
        ggplot2::aes(x = x, y = y, group = id, fill = col),
        alpha = config$flow_alpha, colour = NA
    ) +
    ggplot2::geom_rect(
        data = gene_blk,
        ggplot2::aes(xmin = x_gene[[1]], xmax = x_gene[[2]], ymin = ymin, ymax = ymax,
            fill = col)
    ) +
    ggplot2::geom_rect(
        data = term_blk,
        ggplot2::aes(xmin = x_term[[1]], xmax = x_term[[2]], ymin = ymin, ymax = ymax,
            fill = col)
    ) +
    ggplot2::geom_text(
        data = gene_blk,
        ggplot2::aes(x = x_gene[[1]] - 0.015, y = y, label = name),
        hjust = 1, size = config$gene_label_size
    ) +
    ggplot2::geom_text(
        data = term_blk,
        ggplot2::aes(x = x_term[[1]] - 0.015, y = y, label = label),
        hjust = 1, size = config$term_label_size, lineheight = 0.85
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(-0.32, x_term[[2]]), ylim = c(0, H),
        expand = FALSE, clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(6, 2, 6, 6)
    )

stat_cols <- palette_colors(config$palettes$stat)
right_y <- term_blk[, c("name", "y", "ymin", "ymax")]
right <- merge(right_y, df[, c(".term", ".ratio", ".count", ".neglog")],
    by.x = "name", by.y = ".term")

p_right <- ggplot2::ggplot(right)
if (identical(config$right_panel, "bar")) {
    p_right <- p_right +
        ggplot2::geom_rect(
            ggplot2::aes(xmin = 0, xmax = .ratio,
                ymin = ymin + 0.12 * (ymax - ymin), ymax = ymax - 0.12 * (ymax - ymin),
                fill = .neglog),
            colour = "grey25", linewidth = 0.25
        ) +
        ggplot2::geom_text(
            ggplot2::aes(x = .ratio, y = y, label = .count),
            hjust = -0.25, size = 3
        ) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18)))
} else {
    p_right <- p_right +
        ggplot2::geom_point(
            ggplot2::aes(x = .ratio, y = y, size = .count, fill = .neglog),
            shape = 21, colour = "grey20", stroke = 0.35
        ) +
        ggplot2::scale_size_continuous(name = config$labels$size,
            range = config$size_range) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.12)))
}
p_right <- p_right +
    ggplot2::scale_fill_gradientn(name = config$labels$fill, colours = stat_cols) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0, H), clip = "off") +
    ggplot2::labs(x = config$labels$x, y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 10),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(6, 6, 6, 2)
    )

p <- patchwork::wrap_plots(p_sankey, p_right, nrow = 1, widths = config$widths)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 11, height = max(5, 0.2 * length(genes) + 1.2))
