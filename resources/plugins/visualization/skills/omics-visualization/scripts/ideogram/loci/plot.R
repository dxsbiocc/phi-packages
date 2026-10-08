#!/usr/bin/env Rscript

# Template-ID: ideogram-loci
#
# Purpose:
#   Mark supplied gene loci on chromosome ideograms, coloured by a
#   category, with gene-name labels. Vertical (ggplot) and circular
#   (circlize) are config$layout, not two ids.
#
# Inputs:
#   One row per locus. Default example:
#     - chr, start, end: genomic interval
#     - gene: label
#     - category: fill / colour group
#   Sidecar beside the input (same directory):
#     - cytoband.tsv: Chrom, Start, End, Stain, Arm
#
# Output:
#   A PDF, PNG, or SVG loci-on-ideogram figure.
#
# Dependencies:
#   ggplot2, readr, ggrepel, ggideogram (GitHub: dxsbiocc/ggideogram);
#   circular layout also needs circlize and ComplexHeatmap.
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, layout, chromosomes).
#   layout = "circular" or "vertical" is CONFIG, not a second id.
#   Catalog preview is the circular figure. Vertical is the same
#   loci table with layout = "vertical".
#   hg19 vs hg38 is the cytoband sidecar plus matching gene
#   coordinates, not a new id. Do not add an id per gene set (m6A,
#   TF, surface marker). Transcript exon models are ideogram-gene.
#   G-banding without loci is ideogram-karyotype.
#
# Scientific assumptions:
#   Loci and cytobands are supplied. This script does not look up a
#   GTF or fetch a genome. Point position is the interval midpoint.

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
        chrom = "chr",
        start = "start",
        end = "end",
        gene = "gene",
        category = "category"
    ),
    sidecars = list(
        cytoband = "cytoband.tsv"
    ),
    cyto = list(
        chrom = "Chrom",
        start = "Start",
        end = "End",
        stain = "Stain",
        arm = "Arm"
    ),
    layout = "circular",
    chromosomes = "all",
    just = 1,
    radius_pt = 5,
    width = 0.3,
    labels = list(
        title = "",
        colour = NULL
    ),
    size = list(
        width = 8,
        height = 8
    )
)

# Catalog hex used in the source figure (Brand.Algolia / Brand.joomla).
category_colors <- c(
    Eraser = "#2dde98",
    Reader = "#f9a541",
    Writer = "#3369e7"
)

load_packages(c(
    "ggplot2", "readr", "ggrepel", "ggideogram", "circlize", "ComplexHeatmap"
))
patch_ggideogram_ggplot2()

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
genes <- read_table_auto(io$input)
require_columns(genes, config$columns)

cyto_path <- file.path(dirname(io$input), config$sidecars$cytoband)
if (!file.exists(cyto_path)) {
    stop(
        paste("Missing cytoband sidecar:", cyto_path),
        call. = FALSE
    )
}
cyto <- read_table_auto(cyto_path)
require_columns(cyto, config$cyto)

g_chr <- config$columns$chrom
g_st <- config$columns$start
g_en <- config$columns$end
g_nm <- config$columns$gene
g_cat <- config$columns$category
c_chr <- config$cyto$chrom
c_st <- config$cyto$start
c_en <- config$cyto$end
c_sn <- config$cyto$stain
c_ar <- config$cyto$arm

genes[[g_chr]] <- as.character(genes[[g_chr]])
genes[[g_st]] <- as.numeric(genes[[g_st]])
genes[[g_en]] <- as.numeric(genes[[g_en]])
genes[[g_nm]] <- as.character(genes[[g_nm]])
swap <- genes[[g_st]] > genes[[g_en]]
if (any(swap)) {
    tmp <- genes[[g_st]][swap]
    genes[[g_st]][swap] <- genes[[g_en]][swap]
    genes[[g_en]][swap] <- tmp
}
genes$mid <- (genes[[g_st]] + genes[[g_en]]) / 2

cat_levels <- names(category_colors)
extra_cat <- setdiff(unique(as.character(genes[[g_cat]])), cat_levels)
if (length(extra_cat)) {
    extra_cols <- palette_colors("Qualitative.Safe", n = length(extra_cat))
    names(extra_cols) <- extra_cat
    category_colors <- c(category_colors, extra_cols)
    cat_levels <- names(category_colors)
}
genes[[g_cat]] <- factor(genes[[g_cat]], levels = cat_levels)
genes <- genes[!is.na(genes[[g_cat]]), , drop = FALSE]

cyto[[c_chr]] <- as.character(cyto[[c_chr]])
cyto[[c_st]] <- as.numeric(cyto[[c_st]])
cyto[[c_en]] <- as.numeric(cyto[[c_en]])
cyto[[c_sn]] <- as.character(cyto[[c_sn]])
cyto[[c_ar]] <- as.character(cyto[[c_ar]])

chrom_order <- function(x) {
    x <- unique(as.character(x))
    key <- sub("^chr", "", x)
    rank <- match(key, c(as.character(seq_len(22)), "X", "Y", "M", "MT"))
    rank[is.na(rank)] <- 1000 + seq_len(sum(is.na(rank)))
    x[order(rank, x)]
}

if (identical(config$chromosomes, "loci")) {
    keep <- chrom_order(genes[[g_chr]])
} else {
    keep <- chrom_order(cyto[[c_chr]])
}
cyto <- cyto[cyto[[c_chr]] %in% keep, , drop = FALSE]
genes <- genes[genes[[g_chr]] %in% keep, , drop = FALSE]
if (!nrow(genes) || !nrow(cyto)) {
    stop("No overlapping chromosomes between loci and cytoband.", call. = FALSE)
}
cyto[[c_chr]] <- factor(cyto[[c_chr]], levels = keep)
genes[[g_chr]] <- factor(genes[[g_chr]], levels = keep)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_vertical <- function() {
    ggplot2::ggplot() +
        ggideogram::geom_ideogram(
            ggplot2::aes(
                x = .data[[c_chr]],
                ymin = .data[[c_st]],
                ymax = .data[[c_en]],
                chrom = .data[[c_chr]],
                fill = .data[[c_sn]],
                arm = .data[[c_ar]]
            ),
            data = cyto,
            width = config$width,
            radius = grid::unit(config$radius_pt, "pt"),
            chrom.col = "#888888",
            just = config$just,
            show.legend = FALSE
        ) +
        ggplot2::geom_point(
            ggplot2::aes(
                x = .data[[g_chr]],
                y = .data$mid,
                colour = .data[[g_cat]]
            ),
            data = genes,
            size = 2,
            alpha = 0.75
        ) +
        ggrepel::geom_text_repel(
            ggplot2::aes(
                x = .data[[g_chr]],
                y = .data$mid,
                colour = .data[[g_cat]],
                label = .data[[g_nm]]
            ),
            data = genes,
            nudge_x = 0.22,
            direction = "y",
            segment.curvature = -0.1,
            size = 3.2,
            angle = 90,
            max.overlaps = 100,
            force = 0.5,
            show.legend = FALSE,
            seed = 1
        ) +
        ggplot2::scale_fill_manual(
            values = ggideogram::cytoband_colors,
            na.value = "white"
        ) +
        ggplot2::scale_colour_manual(
            name = config$labels$colour,
            values = category_colors
        ) +
        ggplot2::labs(title = config$labels$title, x = NULL, y = NULL) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.background = ggplot2::element_blank(),
            legend.background = ggplot2::element_blank(),
            legend.position = "top",
            axis.text.x = ggplot2::element_text(colour = "black", size = 9),
            plot.margin = ggplot2::margin(8, 8, 8, 8)
        )
}

draw_circular <- function() {
    band <- data.frame(
        chr = as.character(cyto[[c_chr]]),
        start = cyto[[c_st]],
        end = cyto[[c_en]],
        value1 = cyto[[c_ar]],
        value2 = cyto[[c_sn]],
        stringsAsFactors = FALSE
    )
    loc <- data.frame(
        chr = as.character(genes[[g_chr]]),
        start = genes[[g_st]],
        end = genes[[g_en]],
        gene = genes[[g_nm]],
        category = genes[[g_cat]],
        stringsAsFactors = FALSE
    )
    n_chr <- length(keep)
    box_cols <- palette_colors("Qualitative.Safe", n = n_chr)
    names(box_cols) <- keep

    circlize::circos.par(
        start.degree = 90,
        gap.after = 2,
        cell.padding = c(0, 0, 0, 0),
        track.margin = c(0.005, 0.005)
    )
    circlize::circos.initializeWithIdeogram(
        cytoband = band,
        plotType = NULL,
        chromosome.index = keep
    )
    circlize::circos.track(
        ylim = c(0, 1),
        track.height = 0.06,
        panel.fun = function(x, y) {
            circlize::circos.text(
                circlize::CELL_META$xcenter,
                0.5,
                sub("^chr", "", circlize::CELL_META$sector.index),
                cex = 0.75,
                niceFacing = TRUE
            )
        },
        cell.padding = c(0, 0, 0, 0),
        bg.border = NA
    )
    circlize::circos.track(
        ylim = c(0, 1),
        panel.fun = function(x, y) {
            chr <- circlize::CELL_META$sector.index
            xlim <- circlize::CELL_META$xlim
            circlize::circos.rect(
                xlim[[1]], 0, xlim[[2]], 1,
                col = box_cols[[chr]],
                border = NA
            )
        },
        track.height = 0.035,
        cell.padding = c(0, 0, 0, 0),
        bg.border = NA
    )
    circlize::circos.genomicIdeogram(band, track.height = 0.1)
    circlize::circos.genomicTrackPlotRegion(
        loc[, c("chr", "start", "end", "category")],
        ylim = c(0, 1),
        track.height = 0.08,
        bg.border = NA,
        cell.padding = c(0, 0, 0, 0),
        panel.fun = function(region, value, ...) {
            mid <- (region[[1]] + region[[2]]) / 2
            circlize::circos.points(
                mid,
                rep(0.5, length(mid)),
                pch = 16,
                cex = 0.7,
                col = category_colors[as.character(value$category)]
            )
        }
    )
    circlize::circos.genomicLabels(
        loc,
        labels.column = 4,
        side = "inside",
        col = category_colors[as.character(loc$category)],
        line_col = category_colors[as.character(loc$category)],
        cex = 0.65
    )
    circlize::circos.clear()
    lgd <- ComplexHeatmap::Legend(
        at = names(category_colors),
        type = "points",
        legend_gp = grid::gpar(col = unname(category_colors)),
        title = NULL,
        nrow = 3,
        background = NA
    )
    ComplexHeatmap::draw(
        lgd,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0.5, "npc"),
        just = "center"
    )
}

if (identical(config$layout, "circular")) {
    out_path <- io$output
    fmt <- tolower(tools::file_ext(out_path))
    if (!fmt %in% c("pdf", "png", "svg")) {
        stop(sprintf("Unsupported output format: %s", fmt), call. = FALSE)
    }
    out_dir <- dirname(out_path)
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    circ_w <- max(config$size$width, config$size$height)
    if (identical(fmt, "png")) {
        grDevices::png(
            out_path,
            width = circ_w,
            height = circ_w,
            units = "in",
            res = 300,
            bg = "transparent"
        )
    } else if (identical(fmt, "pdf")) {
        grDevices::pdf(out_path, width = circ_w, height = circ_w)
    } else {
        grDevices::svg(
            out_path,
            width = circ_w,
            height = circ_w,
            bg = "transparent"
        )
    }
    graphics::par(mar = c(0, 0, 0, 0), bg = "transparent")
    tryCatch(
        draw_circular(),
        finally = grDevices::dev.off()
    )
    message("Saved: ", normalizePath(out_path, mustWork = TRUE))
} else {
    save_ggplot(
        draw_vertical(),
        io$output,
        width = config$size$width,
        height = config$size$height
    )
}
