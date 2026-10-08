#!/usr/bin/env Rscript

# Template-ID: ideogram-nested
#
# Purpose:
#   Draw a nested Circos: an outer genome (or sector run) and an
#   inner zoom of supplied windows, joined by correspondence
#   connectors. Inner points are supplied scores on those windows
#   (for example DMR methDiff). This is circos.nested.
#
# Inputs:
#   One row per zoom window on the outer scale. Default example:
#     - chr, start, end: locus on the outer genome
#     - name, zstart, zend: inner sector and its coordinates
#   Sidecars beside the input:
#     - cytoband.tsv: chr, start, end, band, stain
#     - windows.tsv: name, start, end, chr
#     - points.tsv: name, start, end, score
#
# Output:
#   A PDF, PNG, or SVG nested Circos figure.
#
# Dependencies:
#   circlize, ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Extra zoom windows are extra correspondence rows, not extra
#   ids. Inner point vs empty zoom is the points sidecar. Do not
#   add an id per assay (T-WGBS / DMR). A single Circos without
#   a nested ring is ideogram-circos. Do not call DMR, windows,
#   or fetch a genome in this script.
#
# Scientific assumptions:
#   Correspondence, window bounds, and point scores are supplied.
#   Connectors are a display of that table, not a liftOver.
#   Signed colour is a display of the supplied score, not a test.

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
        name = "name",
        zstart = "zstart",
        zend = "zend"
    ),
    sidecars = list(
        cytoband = "cytoband.tsv",
        windows = "windows.tsv",
        points = "points.tsv"
    ),
    cyto = list(
        chrom = "chr",
        start = "start",
        end = "end",
        band = "band",
        stain = "stain"
    ),
    window = list(
        name = "name",
        start = "start",
        end = "end",
        chrom = "chr"
    ),
    point = list(
        name = "name",
        start = "start",
        end = "end",
        score = "score"
    ),
    palettes = list(
        chromosome = "Qualitative.Bold",
        score = "Diverging.TealRose"
    ),
    labels = list(
        score = "score"
    ),
    show_axis = TRUE,
    size = list(
        width = 8,
        height = 8
    )
)

load_packages(c("circlize", "ComplexHeatmap", "grid", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
corr <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(corr, config$columns)
for (nm in c("chrom", "name")) {
    corr[[config$columns[[nm]]]] <- as.character(corr[[config$columns[[nm]]]])
}
for (nm in c("start", "end", "zstart", "zend")) {
    corr[[config$columns[[nm]]]] <- as.numeric(corr[[config$columns[[nm]]]])
}
if (any(!is.finite(corr[[config$columns$start]]) |
        !is.finite(corr[[config$columns$end]]) |
        corr[[config$columns$end]] <= corr[[config$columns$start]] |
        !is.finite(corr[[config$columns$zstart]]) |
        !is.finite(corr[[config$columns$zend]]) |
        corr[[config$columns$zend]] <= corr[[config$columns$zstart]])) {
    stop("correspondence intervals must be finite with end > start.", call. = FALSE)
}

in_dir <- dirname(io$input)
cyto_path <- file.path(in_dir, config$sidecars$cytoband)
win_path <- file.path(in_dir, config$sidecars$windows)
pt_path <- file.path(in_dir, config$sidecars$points)
if (!file.exists(cyto_path)) {
    stop(paste("Missing cytoband sidecar:", cyto_path), call. = FALSE)
}
if (!file.exists(win_path)) {
    stop(paste("Missing windows sidecar:", win_path), call. = FALSE)
}
if (!file.exists(pt_path)) {
    stop(paste("Missing points sidecar:", pt_path), call. = FALSE)
}

cyto <- as.data.frame(read_table_auto(cyto_path), stringsAsFactors = FALSE)
require_columns(cyto, config$cyto)
cyto[[config$cyto$chrom]] <- as.character(cyto[[config$cyto$chrom]])
cyto[[config$cyto$start]] <- as.numeric(cyto[[config$cyto$start]])
cyto[[config$cyto$end]] <- as.numeric(cyto[[config$cyto$end]])

windows <- as.data.frame(read_table_auto(win_path), stringsAsFactors = FALSE)
require_columns(windows, config$window)
windows[[config$window$name]] <- as.character(windows[[config$window$name]])
windows[[config$window$chrom]] <- as.character(windows[[config$window$chrom]])
windows[[config$window$start]] <- as.numeric(windows[[config$window$start]])
windows[[config$window$end]] <- as.numeric(windows[[config$window$end]])

points <- as.data.frame(read_table_auto(pt_path), stringsAsFactors = FALSE)
require_columns(points, config$point)
points[[config$point$name]] <- as.character(points[[config$point$name]])
points[[config$point$start]] <- as.numeric(points[[config$point$start]])
points[[config$point$end]] <- as.numeric(points[[config$point$end]])
points[[config$point$score]] <- as.numeric(points[[config$point$score]])
if (any(!is.finite(points[[config$point$score]]))) {
    stop("point score must be finite.", call. = FALSE)
}

keep <- unique(cyto[[config$cyto$chrom]])
n_drop <- sum(!corr[[config$columns$chrom]] %in% keep)
corr <- corr[corr[[config$columns$chrom]] %in% keep, , drop = FALSE]
if (n_drop) {
    message("Dropped ", n_drop, " correspondence row(s) whose chr is not in the cytoband.")
}
if (!nrow(corr)) {
    stop("No correspondence rows match cytoband chromosomes.", call. = FALSE)
}

win_keep <- windows[[config$window$name]] %in% corr[[config$columns$name]]
windows <- windows[win_keep, , drop = FALSE]
if (!nrow(windows)) {
    stop("No windows match correspondence names.", call. = FALSE)
}
points <- points[points[[config$point$name]] %in% windows[[config$window$name]], , drop = FALSE]
if (!nrow(points)) {
    stop("No points match window names.", call. = FALSE)
}

band <- data.frame(
    chr = cyto[[config$cyto$chrom]],
    start = cyto[[config$cyto$start]],
    end = cyto[[config$cyto$end]],
    value1 = cyto[[config$cyto$band]],
    value2 = cyto[[config$cyto$stain]],
    stringsAsFactors = FALSE
)
tagments <- data.frame(
    name = windows[[config$window$name]],
    start = windows[[config$window$start]],
    end = windows[[config$window$end]],
    chr = windows[[config$window$chrom]],
    stringsAsFactors = FALSE
)
dmr <- data.frame(
    chr = points[[config$point$name]],
    start = points[[config$point$start]],
    end = points[[config$point$end]],
    score = points[[config$point$score]],
    stringsAsFactors = FALSE
)
correspondance <- data.frame(
    chr = corr[[config$columns$chrom]],
    start = corr[[config$columns$start]],
    end = corr[[config$columns$end]],
    name = corr[[config$columns$name]],
    zstart = corr[[config$columns$zstart]],
    zend = corr[[config$columns$zend]],
    stringsAsFactors = FALSE
)

n_chr <- length(keep)
chr_cols <- palette_colors(config$palettes$chromosome, n = n_chr)
names(chr_cols) <- keep
score_pal <- palette_colors(config$palettes$score)
ylim_sc <- max(abs(dmr$score), na.rm = TRUE)
if (!is.finite(ylim_sc) || ylim_sc == 0) {
    ylim_sc <- 0.5
}
col_fun <- circlize::colorRamp2(
    seq(-ylim_sc, ylim_sc, length.out = length(score_pal)),
    score_pal
)
outline <- palette_colors("Qualitative.Safe")[[12]]
conn_col <- grDevices::adjustcolor(
    unname(chr_cols[correspondance$chr]),
    alpha.f = 0.55
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
f_outer <- function() {
    circlize::circos.par(
        start.degree = 90,
        gap.after = 2,
        cell.padding = c(0, 0, 0, 0),
        track.margin = c(0.002, 0.002),
        circle.margin = 0.12,
        points.overflow.warning = FALSE
    )
    circlize::circos.initializeWithIdeogram(
        cytoband = band,
        plotType = NULL,
        chromosome.index = keep
    )
    circlize::circos.track(
        ylim = c(0, 1),
        track.height = 0.14,
        bg.border = NA,
        panel.fun = function(x, y) {
            circlize::circos.text(
                circlize::CELL_META$xcenter,
                1 + circlize::mm_y(3),
                sub("^chr", "", circlize::CELL_META$sector.index),
                cex = 1.1,
                facing = "downward",
                adj = c(0.5, 0.5)
            )
        }
    )
    circlize::circos.genomicIdeogram(band, track.height = 0.08)
    if (isTRUE(config$show_axis)) {
        circlize::circos.track(
            track.index = circlize::get.current.track.index(),
            bg.border = NA,
            panel.fun = function(x, y) {
                circlize::circos.genomicAxis(
                    h = "top",
                    labels.cex = 0.6,
                    major.tick.length = circlize::mm_y(2),
                    labels.facing = "clockwise",
                    labels.niceFacing = TRUE,
                    labels.pos.adjust = TRUE,
                    col = outline,
                    labels.col = "black"
                )
            }
        )
    }
}

f_inner <- function() {
    n_win <- nrow(tagments)
    circlize::circos.par(
        cell.padding = c(0.002, 0, 0.002, 0),
        gap.after = c(rep(1, n_win - 1L), 8),
        track.margin = c(0.01, 0)
    )
    circlize::circos.genomicInitialize(tagments[, 1:3], plotType = NULL)
    circlize::circos.genomicTrack(
        dmr,
        ylim = c(-ylim_sc, ylim_sc),
        bg.col = grDevices::adjustcolor(
            unname(chr_cols[tagments$chr[match(
                circlize::get.all.sector.index(),
                tagments$name
            )]]),
            alpha.f = 0.22
        ),
        bg.border = outline,
        bg.lwd = 0.3,
        track.height = 0.28,
        panel.fun = function(region, value, ...) {
            circlize::circos.lines(
                circlize::CELL_META$xlim,
                c(0, 0),
                col = outline,
                lwd = 0.4,
                lty = 2
            )
            mid <- (region[[1]] + region[[2]]) / 2
            circlize::circos.points(
                mid,
                value[[1]],
                pch = 16,
                cex = 0.7,
                col = col_fun(value[[1]])
            )
        }
    )
}

draw_circos <- function() {
    circlize::circos.nested(
        f_outer,
        f_inner,
        correspondance,
        connection_col = conn_col
    )
    lgd <- ComplexHeatmap::Legend(
        title = config$labels$score,
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
graphics::par(mar = c(0.15, 0.15, 0.15, 0.15), bg = "transparent")
tryCatch(
    draw_circos(),
    finally = {
        circlize::circos.clear()
        grDevices::dev.off()
    }
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
