#!/usr/bin/env Rscript

# Template-ID: ideogram-circos-heatmap
#
# Purpose:
#   Draw genomic intervals as a circular heatmap, with connector
#   lines from each locus to its tile stack. This is
#   circos.genomicHeatmap, not a group-split circular matrix.
#
# Inputs:
#   One row per interval. Default example:
#     - chr, start, end: genomic interval
#     - remaining numeric columns: heatmap layers (outside to inside)
#   Sidecar beside the input:
#     - cytoband.tsv: chr, start, end, band, stain
#
# Output:
#   A PDF, PNG, or SVG genomic Circos heatmap.
#
# Dependencies:
#   circlize, ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Extra numeric columns are extra heatmap layers, not extra ids.
#   side = "inside" vs "outside" is CONFIG. A Group-split gene
#   matrix is heatmap-circos-split. Histogram / scatter tracks
#   without tiles are ideogram-circos. Do not call peaks, z-score,
#   or fetch a genome in this script.
#
# Scientific assumptions:
#   Interval values are supplied. Connector lines are a display
#   of interval midpoints, not a statistical join. This script
#   does not window, impute, or rescale the matrix.

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
        end = "end"
    ),
    sidecars = list(
        cytoband = "cytoband.tsv"
    ),
    cyto = list(
        chrom = "chr",
        start = "start",
        end = "end",
        band = "band",
        stain = "stain"
    ),
    side = "inside",
    connection_height = 0.09,
    heatmap_height = 0.18,
    palettes = list(
        fill = "Diverging.TealRose",
        chromosome = "Qualitative.Bold"
    ),
    labels = list(
        fill = "score"
    ),
    size = list(
        width = 8,
        height = 8
    )
)

load_packages(c("circlize", "ComplexHeatmap", "grid", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
c_chr <- config$columns$chrom
c_st <- config$columns$start
c_en <- config$columns$end
df[[c_chr]] <- as.character(df[[c_chr]])
df[[c_st]] <- as.numeric(df[[c_st]])
df[[c_en]] <- as.numeric(df[[c_en]])
if (any(!nzchar(df[[c_chr]]))) {
    stop("chr must be non-empty.", call. = FALSE)
}
if (any(!is.finite(df[[c_st]]) | !is.finite(df[[c_en]]) | df[[c_en]] <= df[[c_st]])) {
    stop("start/end must be finite with end > start.", call. = FALSE)
}

value_cols <- setdiff(names(df), unlist(config$columns, use.names = FALSE))
if (length(value_cols) < 1L) {
    stop("Need at least one numeric heatmap column.", call. = FALSE)
}
for (nm in value_cols) {
    df[[nm]] <- as.numeric(df[[nm]])
    if (any(!is.finite(df[[nm]]))) {
        stop(paste("Non-finite value in column", nm), call. = FALSE)
    }
}

cyto_path <- file.path(dirname(io$input), config$sidecars$cytoband)
if (!file.exists(cyto_path)) {
    stop(paste("Missing cytoband sidecar:", cyto_path), call. = FALSE)
}
cyto <- as.data.frame(read_table_auto(cyto_path), stringsAsFactors = FALSE)
require_columns(cyto, config$cyto)
cyto[[config$cyto$chrom]] <- as.character(cyto[[config$cyto$chrom]])
cyto[[config$cyto$start]] <- as.numeric(cyto[[config$cyto$start]])
cyto[[config$cyto$end]] <- as.numeric(cyto[[config$cyto$end]])

keep <- unique(cyto[[config$cyto$chrom]])
n_drop <- sum(!df[[c_chr]] %in% keep)
df <- df[df[[c_chr]] %in% keep, , drop = FALSE]
if (n_drop) {
    message("Dropped ", n_drop, " interval row(s) whose chr is not in the cytoband.")
}
if (!nrow(df)) {
    stop("No interval rows match cytoband chromosomes.", call. = FALSE)
}
df[[c_chr]] <- factor(df[[c_chr]], levels = keep)

bed <- data.frame(
    chr = as.character(df[[c_chr]]),
    start = df[[c_st]],
    end = df[[c_en]],
    df[, value_cols, drop = FALSE],
    stringsAsFactors = FALSE
)

band <- data.frame(
    chr = cyto[[config$cyto$chrom]],
    start = cyto[[config$cyto$start]],
    end = cyto[[config$cyto$end]],
    value1 = cyto[[config$cyto$band]],
    value2 = cyto[[config$cyto$stain]],
    stringsAsFactors = FALSE
)

fill_pal <- palette_colors(config$palettes$fill)
mat <- as.matrix(bed[, value_cols, drop = FALSE])
lim <- max(abs(mat), na.rm = TRUE)
if (!is.finite(lim) || lim == 0) {
    lim <- 1
}
col_fun <- circlize::colorRamp2(
    seq(-lim, lim, length.out = length(fill_pal)),
    fill_pal
)
n_chr <- length(keep)
chr_cols <- palette_colors(config$palettes$chromosome, n = n_chr)
names(chr_cols) <- keep
outline <- palette_colors("Qualitative.Safe")[[12]]
side <- config$side
if (!side %in% c("inside", "outside")) {
    stop("side must be 'inside' or 'outside'.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_circos <- function() {
    circlize::circos.par(
        start.degree = 90,
        gap.after = 2,
        cell.padding = c(0, 0, 0, 0),
        track.margin = c(0.004, 0.004),
        circle.margin = 0.02,
        points.overflow.warning = FALSE
    )
    circlize::circos.initializeWithIdeogram(
        cytoband = band,
        plotType = NULL,
        chromosome.index = keep
    )
    if (identical(side, "inside")) {
        circlize::circos.track(
            ylim = c(0, 1),
            track.height = 0.05,
            bg.border = NA,
            panel.fun = function(x, y) {
                circlize::circos.text(
                    circlize::CELL_META$xcenter,
                    0.5,
                    sub("^chr", "", circlize::CELL_META$sector.index),
                    cex = 0.72,
                    facing = "downward"
                )
            }
        )
        circlize::circos.genomicIdeogram(band, track.height = 0.07)
        circlize::circos.genomicHeatmap(
            bed,
            col = col_fun,
            side = "inside",
            border = "white",
            connection_height = config$connection_height,
            heatmap_height = config$heatmap_height,
            line_col = unname(chr_cols[bed$chr]),
            line_lwd = 0.6
        )
    } else {
        circlize::circos.genomicHeatmap(
            bed,
            col = col_fun,
            side = "outside",
            border = "white",
            connection_height = config$connection_height,
            heatmap_height = config$heatmap_height,
            line_col = unname(chr_cols[bed$chr]),
            line_lwd = 0.6
        )
        circlize::circos.track(
            ylim = c(0, 1),
            track.height = 0.05,
            bg.border = NA,
            panel.fun = function(x, y) {
                circlize::circos.text(
                    circlize::CELL_META$xcenter,
                    0.5,
                    sub("^chr", "", circlize::CELL_META$sector.index),
                    cex = 0.72,
                    facing = "downward"
                )
            }
        )
        circlize::circos.genomicIdeogram(band, track.height = 0.07)
    }
    circlize::circos.clear()

    lgd <- ComplexHeatmap::Legend(
        title = config$labels$fill,
        col_fun = col_fun,
        title_gp = grid::gpar(fontsize = 9),
        labels_gp = grid::gpar(fontsize = 8),
        title_position = "topcenter",
        background = NA
    )
    ComplexHeatmap::draw(
        lgd,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0.5, "npc"),
        just = "center"
    )
}

out_path <- io$output
fmt <- tolower(tools::file_ext(out_path))
if (!fmt %in% c("pdf", "png", "svg")) {
    stop(sprintf("Unsupported output format: %s", fmt), call. = FALSE)
}
out_dir <- dirname(out_path)
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
width <- config$size$width
height <- config$size$height
if (identical(fmt, "png")) {
    grDevices::png(
        out_path, width = width, height = height, units = "in",
        res = 300, bg = "transparent"
    )
} else if (identical(fmt, "pdf")) {
    grDevices::pdf(out_path, width = width, height = height)
} else {
    grDevices::svg(out_path, width = width, height = height, bg = "transparent")
}
graphics::par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "transparent")
tryCatch(
    draw_circos(),
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
