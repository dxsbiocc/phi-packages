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
#   Optional sidecars beside the input:
#     - tracks.tsv: chr, start, end, value, track, type
#       type is scatter / histogram / line; track order is file order.
#     - genes.tsv: chr, start, end, label
#   Coordinates are 0-based, half-open base-pair intervals on one assembly.
#
# Output:
#   A PDF, PNG, or SVG genomic Circos heatmap.
#
# Dependencies:
#   circlize, ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#   Rscript plot.R multitrack/example.tsv output.pdf
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
#   Optional track values and gene labels are supplied, never inferred.

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
        cytoband = "cytoband.tsv",
        tracks = "tracks.tsv", # optional; NULL disables this sidecar
        genes = "genes.tsv"    # optional; NULL disables this sidecar
    ),
    cyto = list(
        chrom = "chr",
        start = "start",
        end = "end",
        band = "band",
        stain = "stain"
    ),
    side = "inside",
    # auto: solid sectors with tracks, otherwise the original cytobands
    sector_style = "auto", # auto / solid / ideogram
    start_degree = 90,
    gap_degree = 2,
    connection_height = 0.06,
    heatmap_height = 0.14,
    tracks = list(
        columns = list(chrom = "chr", start = "start", end = "end",
                       value = "value", track = "track", type = "type"),
        height = 0.09,
        point_cex = 0.32
    ),
    genes = list(
        columns = list(chrom = "chr", start = "start", end = "end", label = "label"),
        cex = 0.48,
        connection_height = 0.025,
        height = 0.12
    ),
    palettes = list(
        fill = "Diverging.Tropic",
        chromosome = "Qualitative.Bold"
    ),
    labels = list(
        fill = "Supplied score",
        center = "Sample 1"
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
validate_intervals <- function(data, columns, table_name, check_bounds = TRUE) {
    require_columns(data, columns)
    if (!nrow(data)) stop(paste(table_name, "must not be empty."), call. = FALSE)
    chr <- as.character(data[[columns$chrom]])
    start <- suppressWarnings(as.numeric(data[[columns$start]]))
    end <- suppressWarnings(as.numeric(data[[columns$end]]))
    if (any(is.na(chr) | !nzchar(trimws(chr)))) {
        stop(paste(table_name, "chr must be non-empty."), call. = FALSE)
    }
    if (any(!is.finite(start) | !is.finite(end) | start < 0 | end <= start |
            start %% 1 != 0 | end %% 1 != 0)) {
        stop(paste(table_name, "requires integer start >= 0 and end > start."), call. = FALSE)
    }
    if (check_bounds) {
        if (any(!chr %in% keep)) {
            stop(paste(table_name, "contains unknown chromosome(s):",
                       paste(unique(chr[!chr %in% keep]), collapse = ", ")), call. = FALSE)
        }
        if (any(start < chr_start[chr] | end > chr_end[chr])) {
            stop(paste(table_name, "interval exceeds chromosome bounds."), call. = FALSE)
        }
    }
    data[[columns$chrom]] <- chr
    data[[columns$start]] <- start
    data[[columns$end]] <- end
    data
}

cyto_path <- file.path(dirname(io$input), config$sidecars$cytoband)
if (!file.exists(cyto_path)) {
    stop(paste("Missing cytoband sidecar:", cyto_path), call. = FALSE)
}
cyto <- as.data.frame(read_table_auto(cyto_path), stringsAsFactors = FALSE)
cyto <- validate_intervals(cyto, config$cyto, "cytoband", check_bounds = FALSE)
keep <- unique(cyto[[config$cyto$chrom]])
chr_start <- tapply(cyto[[config$cyto$start]], cyto[[config$cyto$chrom]], min)
chr_end <- tapply(cyto[[config$cyto$end]], cyto[[config$cyto$chrom]], max)

df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
df <- validate_intervals(df, config$columns, "heatmap")
c_chr <- config$columns$chrom
c_st <- config$columns$start
c_en <- config$columns$end

value_cols <- setdiff(names(df), unlist(config$columns, use.names = FALSE))
if (length(value_cols) < 1L) {
    stop("Need at least one numeric heatmap column.", call. = FALSE)
}
for (nm in value_cols) {
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    if (any(!is.finite(df[[nm]]))) {
        stop(paste("Non-finite value in column", nm), call. = FALSE)
    }
}

bed <- data.frame(
    chr = as.character(df[[c_chr]]),
    start = df[[c_st]],
    end = df[[c_en]],
    df[, value_cols, drop = FALSE],
    stringsAsFactors = FALSE
)
# circlize's stack=TRUE puts the first numeric column at the inner radius.
# Reverse display columns so input order and the key both read outside in.
heatmap_bed <- bed[, c("chr", "start", "end", rev(value_cols)), drop = FALSE]

band <- data.frame(
    chr = cyto[[config$cyto$chrom]],
    start = cyto[[config$cyto$start]],
    end = cyto[[config$cyto$end]],
    value1 = cyto[[config$cyto$band]],
    value2 = cyto[[config$cyto$stain]],
    stringsAsFactors = FALSE
)

read_optional <- function(filename) {
    if (is.null(filename)) return(NULL)
    path <- file.path(dirname(io$input), filename)
    if (!file.exists(path)) return(NULL)
    as.data.frame(read_table_auto(path), stringsAsFactors = FALSE)
}

tracks <- read_optional(config$sidecars$tracks)
track_order <- character()
track_limits <- list()
track_types <- character()
track_ranges <- character()
if (!is.null(tracks)) {
    tc <- config$tracks$columns
    tracks <- validate_intervals(tracks, tc, "tracks")
    tracks <- data.frame(
        chr = tracks[[tc$chrom]], start = tracks[[tc$start]], end = tracks[[tc$end]],
        value = suppressWarnings(as.numeric(tracks[[tc$value]])),
        track = as.character(tracks[[tc$track]]), type = as.character(tracks[[tc$type]])
    )
    if (any(!is.finite(tracks$value))) {
        stop("tracks value must be finite numeric data.", call. = FALSE)
    }
    if (any(is.na(tracks$track) | !nzchar(trimws(tracks$track)))) {
        stop("tracks track names must be non-empty.", call. = FALSE)
    }
    if (any(is.na(tracks$type) | !tracks$type %in% c("scatter", "histogram", "line"))) {
        stop("tracks type must be scatter, histogram, or line.", call. = FALSE)
    }
    track_order <- unique(tracks$track)
    for (tr in track_order) {
        rows <- tracks[tracks$track == tr, , drop = FALSE]
        types <- unique(rows$type)
        if (length(types) != 1L) {
            stop(paste("Each track must have one type:", tr), call. = FALSE)
        }
        track_types[tr] <- types
        limits <- range(c(0, rows$value))
        if (diff(limits) == 0) limits <- c(-1, 1)
        track_limits[[tr]] <- limits + diff(limits) * c(-0.08, 0.08)
        track_ranges[tr] <- paste(format(signif(range(rows$value), 2), trim = TRUE), collapse = " to ")
    }
}

genes <- read_optional(config$sidecars$genes)
if (!is.null(genes)) {
    gc <- config$genes$columns
    genes <- validate_intervals(genes, gc, "genes")
    genes <- data.frame(chr = genes[[gc$chrom]], start = genes[[gc$start]],
                        end = genes[[gc$end]], label = as.character(genes[[gc$label]]))
    if (any(is.na(genes$label) | !nzchar(trimws(genes$label)))) {
        stop("genes label must be non-empty.", call. = FALSE)
    }
}

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
ink <- palette_colors("Brand.Emma")[[1]]
negative <- fill_pal[[1]]
positive <- fill_pal[[length(fill_pal)]]
side <- config$side
if (!side %in% c("inside", "outside")) {
    stop("side must be 'inside' or 'outside'.", call. = FALSE)
}
sector_style <- config$sector_style
if (identical(sector_style, "auto")) {
    sector_style <- if (length(track_order)) "solid" else "ideogram"
}
if (!sector_style %in% c("solid", "ideogram")) {
    stop("sector_style must be auto, solid, or ideogram.", call. = FALSE)
}
# Reserve enough radial space for every supplied track; never drop a layer.
radial_height <- 0.14 + length(track_order) * (config$tracks$height + 0.008) +
    config$connection_height + config$heatmap_height + 0.016
if (!is.null(genes)) radial_height <- radial_height +
    config$genes$connection_height + config$genes$height + 0.016
if (radial_height > 0.90) {
    stop("Tracks exceed radial space; reduce CONFIG track / heatmap / label heights.", call. = FALSE)
}

draw_sectors <- function() {
    if (identical(sector_style, "solid")) {
        circlize::circos.track(
            ylim = c(0, 1), track.height = 0.065, bg.border = NA,
            panel.fun = function(x, y) {
                chr <- circlize::CELL_META$sector.index
                xr <- circlize::CELL_META$xlim
                fill <- chr_cols[[chr]]
                text_col <- if (.hex_luminance(fill) > 0.6) ink else "white"
                circlize::circos.rect(xr[[1]], 0, xr[[2]], 1, col = fill,
                                      border = outline, lwd = 0.35)
                circlize::circos.text(mean(xr), 0.5, chr, cex = 0.52,
                    col = text_col, facing = "bending.inside", niceFacing = TRUE)
                at <- seq(0, xr[[2]], by = 1e8)
                at <- at[at >= xr[[1]]]
                if (length(at)) {
                    circlize::circos.axis(h = "top", major.at = at,
                        labels = paste0(at / 1e6, "Mb"), labels.cex = 0.28,
                        minor.ticks = 4, major.tick.length = circlize::mm_y(0.8),
                        labels.niceFacing = TRUE, col = outline, labels.col = ink, lwd = 0.4)
                }
            }
        )
    } else {
        circlize::circos.track(
            ylim = c(0, 1), track.height = 0.05, bg.border = NA,
            panel.fun = function(x, y) {
                circlize::circos.text(circlize::CELL_META$xcenter, 0.5,
                    sub("^chr", "", circlize::CELL_META$sector.index),
                    cex = 0.72, facing = "downward")
            }
        )
        circlize::circos.genomicIdeogram(band, track.height = 0.07)
    }
}

draw_tracks <- function() {
    for (tr in track_order) {
        rows <- tracks[tracks$track == tr, , drop = FALSE]
        type <- track_types[[tr]]
        circlize::circos.track(ylim = track_limits[[tr]],
            track.height = config$tracks$height, bg.border = outline, bg.lwd = 0.35,
            panel.fun = function(x, y) {
                d <- rows[rows$chr == circlize::CELL_META$sector.index, , drop = FALSE]
                # Draw even for chromosomes without measurements in this track.
                circlize::circos.lines(circlize::CELL_META$xlim, c(0, 0),
                    col = outline, lty = 3, lwd = 0.45)
                if (!nrow(d)) return(invisible())
                mid <- (d$start + d$end) / 2
                if (identical(type, "histogram")) {
                    circlize::circos.rect(d$start, pmin(0, d$value), d$end, pmax(0, d$value),
                        col = ifelse(d$value < 0, negative, positive), border = NA)
                } else if (identical(type, "scatter")) {
                    circlize::circos.points(mid, d$value, col = ink,
                                           pch = 16, cex = config$tracks$point_cex)
                } else {
                    order <- order(mid)
                    if (length(mid) > 1L) {
                        circlize::circos.lines(mid[order], d$value[order], col = ink, lwd = 0.65)
                    } else {
                        circlize::circos.points(mid, d$value, col = ink, pch = 16,
                                               cex = config$tracks$point_cex)
                    }
                }
            }
        )
    }
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_circos <- function() {
    on.exit(circlize::circos.clear(), add = TRUE)
    circlize::circos.par(
        start.degree = config$start_degree,
        gap.after = config$gap_degree,
        cell.padding = c(0, 0, 0, 0),
        track.margin = c(0.004, 0.004),
        circle.margin = 0.04
    )
    circlize::circos.initializeWithIdeogram(
        cytoband = band,
        plotType = NULL,
        chromosome.index = keep
    )
    if (identical(side, "inside")) {
        draw_sectors()
        draw_tracks()
        circlize::circos.genomicHeatmap(
            heatmap_bed,
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
            heatmap_bed,
            col = col_fun,
            side = "outside",
            border = "white",
            connection_height = config$connection_height,
            heatmap_height = config$heatmap_height,
            line_col = unname(chr_cols[bed$chr]),
            line_lwd = 0.6
        )
        draw_sectors()
        draw_tracks()
    }
    if (!is.null(genes)) {
        circlize::circos.genomicLabels(
            genes, labels.column = 4, side = "inside",
            cex = config$genes$cex, col = ink, line_col = outline, line_lwd = 0.45,
            connection_height = config$genes$connection_height,
            labels_height = config$genes$height, padding = 0.25
        )
    }
    circlize::circos.clear()

    has_tracks <- length(track_order) > 0L
    if (has_tracks) {
        grid::grid.text(config$labels$center, x = 0.5, y = 0.57,
                        gp = grid::gpar(fontsize = 12, fontface = "bold", col = ink))
        track_key <- paste0(seq_along(track_order), ". ", track_order,
                            " [", track_ranges[track_order], "]")
        grid::grid.text(paste(c("Rings: outside to inside", track_key), collapse = "\n"),
                        x = 0.5, y = 0.525,
                        gp = grid::gpar(fontsize = 7, lineheight = 1.25, col = ink))
        layer_key <- paste("Tiles: outside to inside", paste(value_cols, collapse = " / "), sep = "\n")
        grid::grid.text(layer_key, x = 0.5, y = 0.47,
                        gp = grid::gpar(fontsize = 6.8, lineheight = 1.1, col = ink))
    }
    lgd <- ComplexHeatmap::Legend(
        title = config$labels$fill,
        col_fun = col_fun,
        at = c(-lim, 0, lim),
        labels = format(signif(c(-lim, 0, lim), 3), trim = TRUE),
        direction = if (has_tracks) "horizontal" else "vertical",
        legend_width = grid::unit(22, "mm"),
        title_gp = grid::gpar(fontsize = if (has_tracks) 7 else 9),
        labels_gp = grid::gpar(fontsize = if (has_tracks) 6.5 else 8),
        title_position = "topcenter",
        background = NA
    )
    ComplexHeatmap::draw(
        lgd,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(if (has_tracks) 0.425 else 0.5, "npc"),
        just = "center"
    )
    if (has_tracks && any(track_types == "histogram")) {
        grid::grid.text("Bars:", x = 0.435, y = 0.382,
                        gp = grid::gpar(fontsize = 6.5, col = ink))
        grid::grid.points(x = grid::unit(c(0.475, 0.535), "npc"),
                          y = grid::unit(c(0.382, 0.382), "npc"), pch = 15,
                          size = grid::unit(1.6, "mm"),
                          gp = grid::gpar(col = c(negative, positive)))
        grid::grid.text(c("< 0", "> 0"), x = c(0.487, 0.547), y = 0.382, just = "left",
                        gp = grid::gpar(fontsize = 6.5, col = ink))
    }
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
