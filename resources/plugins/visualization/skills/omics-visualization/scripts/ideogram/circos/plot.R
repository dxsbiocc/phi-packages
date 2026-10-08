#!/usr/bin/env Rscript

# Template-ID: ideogram-circos
#
# Purpose:
#   Draw a circular genome: chromosome ideogram, concentric
#   quantitative tracks, and optional rearrangement links. This is
#   the classic Circos genome figure (CNV / density / translocations).
#
# Inputs:
#   One row per genomic window. Default example:
#     - chr, start, end: interval
#     - value: track height or signed score
#     - track: track name (file order is draw order, outside in)
#     - type: histogram / scatter / line
#   Sidecars beside the input (same directory):
#     - cytoband.tsv: chr, start, end, band, stain
#     - links.tsv: chr1, start1, end1, chr2, start2, end2
#
# Output:
#   A PDF, PNG, or SVG multi-track genomic Circos figure.
#
# Dependencies:
#   circlize, ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Extra tracks are extra `track` levels, not extra ids. Histogram
#   vs scatter vs line is the type column. links on/off is CONFIG.
#   One chromosome (microbial / mito) is the cytoband table, not a
#   new id. Gene-loci points with names are ideogram-loci. Two
#   assemblies joined by homology blocks are ideogram-synteny.
#   Do not call CNV, SV, or density in this script.
#
# Scientific assumptions:
#   Window values and links are supplied. Signed scatter is a display
#   of the supplied value (for example log2 copy number), not a test.
#   Expanding windows or linking breakpoints is not performed here.

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
        value = "value",
        track = "track",
        type = "type"
    ),
    sidecars = list(
        cytoband = "cytoband.tsv",
        links = "links.tsv"
    ),
    cyto = list(
        chrom = "chr",
        start = "start",
        end = "end",
        band = "band",
        stain = "stain"
    ),
    link = list(
        chrom1 = "chr1",
        start1 = "start1",
        end1 = "end1",
        chrom2 = "chr2",
        start2 = "start2",
        end2 = "end2"
    ),
    show_links = TRUE,
    show_ideogram = TRUE,
    track_height = 0.16,
    palettes = list(
        histogram = "Quantitative.BluGrn",
        scatter = "Diverging.TealRose",
        chromosome = "Qualitative.Safe"
    ),
    labels = list(
        histogram = "density",
        scatter = "log2 CN"
    ),
    size = list(
        width = 8,
        height = 8.4
    )
)

load_packages(c("circlize", "ComplexHeatmap", "grid", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
df <- as.data.frame(df, stringsAsFactors = FALSE)

c_chr <- config$columns$chrom
c_st <- config$columns$start
c_en <- config$columns$end
c_val <- config$columns$value
c_tr <- config$columns$track
c_ty <- config$columns$type

df[[c_chr]] <- as.character(df[[c_chr]])
df[[c_st]] <- as.numeric(df[[c_st]])
df[[c_en]] <- as.numeric(df[[c_en]])
df[[c_val]] <- as.numeric(df[[c_val]])
df[[c_tr]] <- as.character(df[[c_tr]])
df[[c_ty]] <- tolower(as.character(df[[c_ty]]))
if (any(!nzchar(df[[c_chr]]) | !nzchar(df[[c_tr]]))) {
    stop("chr and track must be non-empty.", call. = FALSE)
}
if (any(!is.finite(df[[c_st]]) | !is.finite(df[[c_en]]) | df[[c_en]] <= df[[c_st]])) {
    stop("start/end must be finite with end > start.", call. = FALSE)
}
ok_type <- df[[c_ty]] %in% c("histogram", "scatter", "line")
if (any(!ok_type)) {
    stop("type must be histogram, scatter, or line.", call. = FALSE)
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
    message("Dropped ", n_drop, " track row(s) whose chr is not in the cytoband.")
}
if (!nrow(df)) {
    stop("No track rows match cytoband chromosomes.", call. = FALSE)
}
df[[c_chr]] <- factor(df[[c_chr]], levels = keep)
track_names <- unique(df[[c_tr]])
for (tn in track_names) {
    types <- unique(df[[c_ty]][df[[c_tr]] == tn])
    if (length(types) != 1L) {
        stop(
            paste0(
                "Track '", tn, "' has mixed type values: ",
                paste(types, collapse = ", "),
                ". One track, one type."
            ),
            call. = FALSE
        )
    }
}

link_path <- file.path(dirname(io$input), config$sidecars$links)
links <- NULL
if (isTRUE(config$show_links) && file.exists(link_path)) {
    links <- as.data.frame(read_table_auto(link_path), stringsAsFactors = FALSE)
    require_columns(links, config$link)
    for (nm in c("chrom1", "start1", "end1", "chrom2", "start2", "end2")) {
        col <- config$link[[nm]]
        if (startsWith(nm, "chrom")) {
            links[[col]] <- as.character(links[[col]])
        } else {
            links[[col]] <- as.numeric(links[[col]])
        }
    }
    keep_link <- links[[config$link$chrom1]] %in% keep &
        links[[config$link$chrom2]] %in% keep
    n_drop_link <- sum(!keep_link)
    links <- links[keep_link, , drop = FALSE]
    if (n_drop_link) {
        message(
            "Dropped ", n_drop_link,
            " link row(s) whose chr is not in the cytoband."
        )
    }
    if (!nrow(links)) {
        links <- NULL
    }
}

n_chr <- length(keep)
chr_cols <- palette_colors(config$palettes$chromosome, n = n_chr)
names(chr_cols) <- keep
outline <- palette_colors("Qualitative.Safe")[[12]]
hist_pal <- palette_colors(config$palettes$histogram)
scat_pal <- palette_colors(config$palettes$scatter)

ramp_from_palette <- function(pal, limits) {
    pal <- unname(as.character(pal))
    n <- length(pal)
    if (limits[[1]] >= limits[[2]]) {
        limits[[2]] <- limits[[1]] + 1e-6
    }
    circlize::colorRamp2(seq(limits[[1]], limits[[2]], length.out = n), pal)
}

band <- data.frame(
    chr = cyto[[config$cyto$chrom]],
    start = cyto[[config$cyto$start]],
    end = cyto[[config$cyto$end]],
    value1 = cyto[[config$cyto$band]],
    value2 = cyto[[config$cyto$stain]],
    stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_circos <- function() {
    circlize::circos.par(
        start.degree = 90,
        gap.after = 2,
        cell.padding = c(0.002, 0, 0.002, 0),
        track.margin = c(0.004, 0.004),
        circle.margin = c(0.02, 0.02, 0.28, 0.02),
        points.overflow.warning = FALSE
    )
    circlize::circos.initializeWithIdeogram(
        cytoband = band,
        plotType = NULL,
        chromosome.index = keep
    )
    if (isTRUE(config$show_ideogram)) {
        circlize::circos.genomicIdeogram(band, track.height = 0.07)
    }
    circlize::circos.track(
        ylim = c(0, 1),
        track.height = 0.05,
        bg.border = NA,
        cell.padding = c(0, 0, 0, 0),
        panel.fun = function(x, y) {
            circlize::circos.text(
                circlize::CELL_META$xcenter,
                0.5,
                sub("^chr", "", circlize::CELL_META$sector.index),
                cex = 0.78,
                facing = "downward"
            )
        }
    )

    legends <- list()
    for (tn in track_names) {
        sub <- df[df[[c_tr]] == tn, , drop = FALSE]
        typ <- sub[[c_ty]][[1]]
        bed <- data.frame(
            chr = as.character(sub[[c_chr]]),
            start = sub[[c_st]],
            end = sub[[c_en]],
            value = sub[[c_val]],
            stringsAsFactors = FALSE
        )
        if (identical(typ, "histogram")) {
            vmax <- max(bed$value, na.rm = TRUE)
            if (!is.finite(vmax) || vmax <= 0) {
                vmax <- 1
            }
            col_fun <- ramp_from_palette(hist_pal, c(0, vmax))
            circlize::circos.genomicTrack(
                bed,
                ylim = c(0, vmax),
                track.height = config$track_height,
                bg.border = outline,
                bg.lwd = 0.3,
                panel.fun = function(region, value, ...) {
                    circlize::circos.genomicRect(
                        region, value,
                        ytop.column = 1,
                        ybottom = 0,
                        col = col_fun(value[[1]]),
                        border = NA
                    )
                }
            )
            legends[[length(legends) + 1L]] <- ComplexHeatmap::Legend(
                title = if (nzchar(config$labels$histogram)) {
                    config$labels$histogram
                } else {
                    tn
                },
                col_fun = col_fun,
                title_gp = grid::gpar(fontsize = 8),
                labels_gp = grid::gpar(fontsize = 7),
                title_position = "leftcenter-rot",
                background = NA
            )
        } else if (identical(typ, "scatter")) {
            lim <- max(abs(bed$value), na.rm = TRUE)
            if (!is.finite(lim) || lim == 0) {
                lim <- 1
            }
            col_fun <- ramp_from_palette(scat_pal, c(-lim, lim))
            circlize::circos.genomicTrack(
                bed,
                ylim = c(-lim, lim),
                track.height = config$track_height,
                bg.border = outline,
                bg.lwd = 0.3,
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
                        cex = 0.38,
                        col = col_fun(value[[1]])
                    )
                }
            )
            legends[[length(legends) + 1L]] <- ComplexHeatmap::Legend(
                title = if (nzchar(config$labels$scatter)) {
                    config$labels$scatter
                } else {
                    tn
                },
                col_fun = col_fun,
                title_gp = grid::gpar(fontsize = 8),
                labels_gp = grid::gpar(fontsize = 7),
                title_position = "leftcenter-rot",
                background = NA
            )
        } else {
            yr <- range(bed$value, na.rm = TRUE)
            if (!all(is.finite(yr)) || yr[[1]] == yr[[2]]) {
                yr <- c(0, 1)
            }
            col_line <- palette_colors("Brand.Algolia")[[1]]
            circlize::circos.genomicTrack(
                bed,
                ylim = yr,
                track.height = config$track_height,
                bg.border = outline,
                bg.lwd = 0.3,
                panel.fun = function(region, value, ...) {
                    mid <- (region[[1]] + region[[2]]) / 2
                    ord <- order(mid)
                    circlize::circos.lines(
                        mid[ord],
                        value[[1]][ord],
                        col = col_line,
                        lwd = 0.9
                    )
                }
            )
        }
    }

    if (!is.null(links)) {
        bed1 <- data.frame(
            chr = links[[config$link$chrom1]],
            start = links[[config$link$start1]],
            end = links[[config$link$end1]],
            stringsAsFactors = FALSE
        )
        bed2 <- data.frame(
            chr = links[[config$link$chrom2]],
            start = links[[config$link$start2]],
            end = links[[config$link$end2]],
            stringsAsFactors = FALSE
        )
        link_col <- grDevices::adjustcolor(
            unname(chr_cols[bed1$chr]),
            alpha.f = 0.55
        )
        circlize::circos.genomicLink(
            bed1, bed2,
            col = link_col,
            border = NA,
            lwd = 1.1
        )
    }
    circlize::circos.clear()

    if (length(legends)) {
        packed <- do.call(
            ComplexHeatmap::packLegend,
            c(legends, list(
                direction = "horizontal",
                gap = grid::unit(8, "mm")
            ))
        )
        ComplexHeatmap::draw(
            packed,
            x = grid::unit(0.5, "npc"),
            y = grid::unit(0.02, "npc"),
            just = c("center", "bottom")
        )
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
graphics::par(mar = c(0.15, 0.15, 0.15, 0.15), bg = "transparent")
tryCatch(
    draw_circos(),
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
