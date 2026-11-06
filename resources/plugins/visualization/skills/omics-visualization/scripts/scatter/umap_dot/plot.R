#!/usr/bin/env Rscript

# Template-ID: scatter-umap-dot
#
# Purpose:
#   Draw a numbered UMAP above a marker-gene dot plot. Each cluster has
#   a dashed density outline and a badge. The dot plot repeats those
#   numbers on a colour strip so the two panels share one identity key.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. embedding: cell, umap_1, umap_2, celltype
#     2. dots: celltype, gene, group, pct, avg
#        group is the facet / lineage bracket. pct is percent expressing.
#        avg is a supplied mean expression.
#
# Output:
#   A PDF, PNG, or SVG two-panel figure.
#
# Dependencies:
#   ggplot2, readr, patchwork, ggrepel, MASS
#
# Example:
#   Rscript plot.R umap.tsv dots.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG for column names and labels.
#   Pass the embedding table, then the dot table, on the CLI.
#   Edit PLOT only when the geometry must change.
#
# Scientific assumptions:
#   Coordinates and cell types are supplied. This script does not run
#   UMAP, PCA, clustering, or Seurat.
#   The dashed line is a display contour of the supplied points.
#   pct and avg are supplied. Colour is a within-gene scale of avg,
#   clamped to [-2, 2] for display.
#   Number order is the first-seen order of celltype in the embedding.

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

io <- parse_io_args(n_input = 2L, input_names = c("umap", "dots"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        cell = "cell",
        x = "umap_1",
        y = "umap_2",
        celltype = "celltype"
    ),
    dot_columns = list(
        celltype = "celltype",
        gene = "gene",
        group = "group",
        pct = "pct",
        avg = "avg"
    ),
    contour_level = 0.95,
    labels = list(
        pct = "% expressing",
        fill = "Scaled\nexpression"
    ),
    palettes = list(
        celltype = "Qualitative.Bold",
        expression = "Quantitative.BluGrn"
    )
)

load_packages(c("ggplot2", "readr", "patchwork", "ggrepel"))
if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Missing R package: MASS\nInstall it before running this script.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
emb <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
require_columns(emb, config$columns)
cell_col <- config$columns$cell
x_col <- config$columns$x
y_col <- config$columns$y
type_col <- config$columns$celltype
emb[[cell_col]] <- as.character(emb[[cell_col]])
emb[[type_col]] <- as.character(emb[[type_col]])
emb[[x_col]] <- as.numeric(emb[[x_col]])
emb[[y_col]] <- as.numeric(emb[[y_col]])
if (any(!is.finite(emb[[x_col]]) | !is.finite(emb[[y_col]]))) {
    stop("Embedding coordinates must be finite.", call. = FALSE)
}
type_levels <- unique(emb[[type_col]])
emb[[type_col]] <- factor(emb[[type_col]], levels = type_levels)
ids <- setNames(seq_along(type_levels), type_levels)
pal <- setNames(
    palette_colors(config$palettes$celltype, n = length(type_levels)),
    type_levels
)

set.seed(1)
emb <- emb[sample.int(nrow(emb)), , drop = FALSE]
lab <- aggregate(
    emb[, c(x_col, y_col)],
    by = list(celltype = emb[[type_col]]),
    FUN = median
)
lab$id <- ids[as.character(lab$celltype)]
counts <- table(emb[[type_col]])
legend_lab <- sprintf(
    "%d  %s  (%s)",
    ids[names(counts)],
    names(counts),
    format(as.integer(counts), big.mark = ",", trim = TRUE)
)
names(legend_lab) <- names(counts)

hull <- function(d, level = 0.95, grid = 80L) {
    if (nrow(d) < 8L) {
        return(NULL)
    }
    bw <- pmax(c(MASS::bandwidth.nrd(d[[x_col]]), MASS::bandwidth.nrd(d[[y_col]])), 0.35)
    pad <- 2 * max(bw)
    k <- MASS::kde2d(
        d[[x_col]], d[[y_col]],
        h = bw, n = grid,
        lims = c(range(d[[x_col]]) + c(-pad, pad), range(d[[y_col]]) + c(-pad, pad))
    )
    ix <- findInterval(d[[x_col]], k$x)
    iy <- findInterval(d[[y_col]], k$y)
    ix <- pmax(1L, pmin(ix, length(k$x)))
    iy <- pmax(1L, pmin(iy, length(k$y)))
    dens <- k$z[cbind(ix, iy)]
    dens <- dens[is.finite(dens)]
    if (!length(dens)) {
        return(NULL)
    }
    thr <- stats::quantile(dens, 1 - level)
    cl <- contourLines(k$x, k$y, k$z, levels = thr)
    if (!length(cl)) {
        return(NULL)
    }
    do.call(rbind, lapply(seq_along(cl), function(i) {
        data.frame(x = cl[[i]]$x, y = cl[[i]]$y, piece = i, stringsAsFactors = FALSE)
    }))
}
outline <- do.call(rbind, lapply(type_levels, function(k) {
    h <- hull(emb[emb[[type_col]] == k, , drop = FALSE], level = config$contour_level)
    if (is.null(h)) {
        return(NULL)
    }
    h$celltype <- k
    h$grp <- paste(k, h$piece)
    h
}))
if (!is.null(outline)) {
    outline$celltype <- factor(outline$celltype, levels = type_levels)
}

dots <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
require_columns(dots, config$dot_columns)
d_type <- config$dot_columns$celltype
d_gene <- config$dot_columns$gene
d_group <- config$dot_columns$group
d_pct <- config$dot_columns$pct
d_avg <- config$dot_columns$avg
dots[[d_type]] <- as.character(dots[[d_type]])
dots[[d_gene]] <- as.character(dots[[d_gene]])
dots[[d_group]] <- as.character(dots[[d_group]])
dots[[d_pct]] <- as.numeric(dots[[d_pct]])
dots[[d_avg]] <- as.numeric(dots[[d_avg]])
if (any(!is.finite(dots[[d_pct]]) | dots[[d_pct]] < 0 | dots[[d_pct]] > 100)) {
    stop("pct must be a percentage in [0, 100].", call. = FALSE)
}
if (any(!is.finite(dots[[d_avg]]))) {
    stop("avg must be finite.", call. = FALSE)
}
missing_types <- setdiff(unique(dots[[d_type]]), type_levels)
if (length(missing_types)) {
    stop(
        paste("Dot cell types missing from the embedding:", paste(missing_types, collapse = ", ")),
        call. = FALSE
    )
}
gene_levels <- unique(dots[[d_gene]])
group_levels <- unique(dots[[d_group]])
dots$z <- ave(dots[[d_avg]], dots[[d_gene]], FUN = function(v) {
    if (length(v) < 2L || stats::sd(v) < 1e-8) {
        return(rep(0, length(v)))
    }
    pmin(pmax(as.numeric(scale(v)), -2), 2)
})
dots[[d_gene]] <- factor(dots[[d_gene]], levels = gene_levels)
dots[[d_group]] <- factor(dots[[d_group]], levels = group_levels)
dots[[d_type]] <- factor(dots[[d_type]], levels = rev(type_levels))

xr <- range(emb[[x_col]])
yr <- range(emb[[y_col]])
pad_x <- 0.14 * diff(xr)
pad_y <- 0.14 * diff(yr)
x0 <- xr[1] - pad_x * 0.72
y0 <- yr[1] - pad_y * 0.55
arrow_len <- 0.12 * diff(xr)
expr_pal <- palette_colors(config$palettes$expression)
brackets <- do.call(rbind, lapply(split(dots, dots[[d_group]]), function(d) {
    genes <- unique(as.character(d[[d_gene]]))
    data.frame(
        group = d[[d_group]][1],
        gene = genes[1],
        gene_end = genes[length(genes)],
        stringsAsFactors = FALSE
    )
}))
brackets <- brackets[brackets$gene != brackets$gene_end, , drop = FALSE]
brackets$group <- factor(brackets$group, levels = group_levels)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p_umap <- ggplot(emb, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(colour = "grey20", size = 1.45, stroke = 0) +
    geom_point(aes(colour = .data[[type_col]]), size = 0.95, stroke = 0)
if (!is.null(outline)) {
    p_umap <- p_umap +
        geom_path(
            data = outline,
            aes(x = x, y = y, group = grp, colour = celltype),
            linetype = "22",
            linewidth = 0.5,
            show.legend = FALSE
        )
}
p_umap <- p_umap +
    geom_label_repel(
        data = lab,
        aes(x = .data[[x_col]], y = .data[[y_col]], label = id, fill = celltype),
        colour = "white",
        fontface = "bold",
        size = 3.2,
        label.size = 0,
        label.r = unit(0.3, "lines"),
        label.padding = unit(0.25, "lines"),
        min.segment.length = 0,
        segment.colour = "grey20",
        segment.size = 0.3,
        box.padding = 0.45,
        seed = 7,
        show.legend = FALSE
    ) +
    annotate(
        "segment",
        x = x0, y = y0, xend = x0 + arrow_len, yend = y0,
        linewidth = 0.4, colour = "black",
        arrow = arrow(length = unit(1.4, "mm"), type = "closed")
    ) +
    annotate(
        "segment",
        x = x0, y = y0, xend = x0, yend = y0 + arrow_len,
        linewidth = 0.4, colour = "black",
        arrow = arrow(length = unit(1.4, "mm"), type = "closed")
    ) +
    annotate("text", x = x0, y = y0, label = "UMAP 1", hjust = 0, vjust = 1.6, size = 2.6, colour = "grey35") +
    annotate("text", x = x0, y = y0, label = "UMAP 2", hjust = 0, vjust = -0.8, angle = 90, size = 2.6, colour = "grey35") +
    scale_colour_manual(values = pal, labels = legend_lab, name = NULL) +
    scale_fill_manual(values = pal, guide = "none") +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    coord_fixed(
        xlim = c(xr[1] - pad_x, xr[2] + pad_x * 0.08),
        ylim = c(yr[1] - pad_y, yr[2] + pad_y * 0.08),
        clip = "off"
    ) +
    theme_void(base_size = 10) +
    theme(
        legend.position = "right",
        legend.text = element_text(size = 8, colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(face = "bold", size = 14, colour = "black"),
        plot.margin = margin(6, 6, 8, 8)
    )

strip <- data.frame(
    celltype = factor(type_levels, levels = rev(type_levels)),
    stringsAsFactors = FALSE
)
p_strip <- ggplot(strip, aes(x = 1, y = celltype)) +
    geom_tile(aes(fill = celltype), width = 0.42, height = 0.72) +
    geom_text(
        aes(label = ids[as.character(celltype)]),
        colour = "white", fontface = "bold", size = 2.6
    ) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = c(0.55, 1.15))) +
    coord_fixed(xlim = c(0.72, 1.28), clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 8),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(18, 10, 8, 4)
    )

p_dot <- ggplot(dots, aes(x = .data[[d_gene]], y = .data[[d_type]])) +
    geom_point(aes(size = .data[[d_pct]], fill = z), shape = 21, colour = "grey35", stroke = 0.2) +
    geom_segment(
        data = brackets,
        aes(x = gene, xend = gene_end, y = length(type_levels) + 0.55, yend = length(type_levels) + 0.55),
        inherit.aes = FALSE,
        linewidth = 0.35,
        colour = "grey25"
    ) +
    geom_segment(
        data = brackets,
        aes(x = gene, xend = gene, y = length(type_levels) + 0.42, yend = length(type_levels) + 0.55),
        inherit.aes = FALSE,
        linewidth = 0.35,
        colour = "grey25"
    ) +
    geom_segment(
        data = brackets,
        aes(x = gene_end, xend = gene_end, y = length(type_levels) + 0.42, yend = length(type_levels) + 0.55),
        inherit.aes = FALSE,
        linewidth = 0.35,
        colour = "grey25"
    ) +
    facet_grid(cols = vars(.data[[d_group]]), scales = "free_x", space = "free_x") +
    scale_size_area(
        max_size = 6,
        limits = c(0, 100),
        breaks = c(25, 50, 75, 100),
        name = config$labels$pct
    ) +
    scale_fill_gradientn(
        colours = expr_pal,
        limits = c(-2, 2),
        name = config$labels$fill
    ) +
    guides(
        size = guide_legend(override.aes = list(fill = "grey70")),
        fill = guide_colourbar(barheight = unit(22, "mm"), barwidth = unit(3, "mm"))
    ) +
    scale_y_discrete(expand = expansion(add = c(0.55, 1.15))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
        panel.grid = element_blank(),
        panel.spacing.x = unit(3.2, "mm"),
        strip.clip = "off",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8, colour = "black", margin = margin(b = 2)),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic", colour = "black", size = 7.5),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 4, 4, 0)
    )

p_bottom <- patchwork::wrap_elements(
    p_strip + p_dot +
        patchwork::plot_layout(widths = c(0.24, 1), guides = "collect") +
        theme(legend.position = "right")
) +
    theme(plot.tag = element_text(face = "bold", size = 14, colour = "black"))
p <- p_umap / p_bottom +
    patchwork::plot_layout(heights = c(1.05, 1)) +
    patchwork::plot_annotation(tag_levels = "a")

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 12.2, height = 10.2)
