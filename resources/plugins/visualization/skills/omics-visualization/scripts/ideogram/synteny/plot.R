#!/usr/bin/env Rscript

# Template-ID: ideogram-synteny
#
# Purpose:
#   Draw a two-genome Circos: each assembly occupies a run of
#   sectors, and homology blocks are ribbons between them.
#
# Inputs:
#   One row per aligned block. Default example:
#     - chr1, start1, end1: interval on genome 1
#     - chr2, start2, end2: homologous interval on genome 2
#     - group: ribbon colour (usually the query chromosome)
#   Sidecar beside the input:
#     - karyotype.tsv: genome, chr, start, end
#       chr names must be unique across genomes (Q1 / T1, not
#       chr1 on both assemblies). Prefix both the karyotype and
#       the block chr1/chr2/group columns.
#
# Output:
#   A PDF, PNG, or SVG comparative-synteny Circos figure.
#
# Dependencies:
#   circlize, ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Extra chromosomes are extra karyotype rows, not extra ids. A
#   third genome is another genome level in the sidecar. One genome
#   with quantitative tracks and SV links is ideogram-circos, not
#   this id. Do not run BLAST, minimap, or MCScanX here.
#
# Scientific assumptions:
#   Homology blocks and chromosome lengths are supplied. Ribbon
#   colour is the supplied group, not an inferred orthology class.
#   Sector order is karyotype file order.

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
        chrom1 = "chr1",
        start1 = "start1",
        end1 = "end1",
        chrom2 = "chr2",
        start2 = "start2",
        end2 = "end2",
        group = "group"
    ),
    sidecars = list(
        karyotype = "karyotype.tsv"
    ),
    karyotype = list(
        genome = "genome",
        chrom = "chr",
        start = "start",
        end = "end"
    ),
    gap_within = 1.2,
    gap_between = 14,
    palettes = list(
        ribbon = "Qualitative.Safe",
        genome = "Qualitative.Bold"
    ),
    labels = list(
        genome = "Genome"
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
blocks <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(blocks, config$columns)

kar_path <- file.path(dirname(io$input), config$sidecars$karyotype)
if (!file.exists(kar_path)) {
    stop(paste("Missing karyotype sidecar:", kar_path), call. = FALSE)
}
kar <- as.data.frame(read_table_auto(kar_path), stringsAsFactors = FALSE)
require_columns(kar, config$karyotype)

k_gen <- config$karyotype$genome
k_chr <- config$karyotype$chrom
k_st <- config$karyotype$start
k_en <- config$karyotype$end
kar[[k_gen]] <- as.character(kar[[k_gen]])
kar[[k_chr]] <- as.character(kar[[k_chr]])
kar[[k_st]] <- as.numeric(kar[[k_st]])
kar[[k_en]] <- as.numeric(kar[[k_en]])
if (any(duplicated(kar[[k_chr]]))) {
    stop("karyotype chr names must be unique across genomes.", call. = FALSE)
}
if (any(kar[[k_en]] <= kar[[k_st]])) {
    stop("karyotype end must be greater than start.", call. = FALSE)
}

for (nm in c("chrom1", "chrom2", "group")) {
    blocks[[config$columns[[nm]]]] <- as.character(blocks[[config$columns[[nm]]]])
}
for (nm in c("start1", "end1", "start2", "end2")) {
    blocks[[config$columns[[nm]]]] <- as.numeric(blocks[[config$columns[[nm]]]])
}
if (any(blocks[[config$columns$end1]] <= blocks[[config$columns$start1]] |
        blocks[[config$columns$end2]] <= blocks[[config$columns$start2]])) {
    stop("block end must be greater than start.", call. = FALSE)
}

sectors <- kar[[k_chr]]
genomes <- unique(kar[[k_gen]])
if (length(genomes) < 2L) {
    stop("karyotype needs at least two genome levels.", call. = FALSE)
}
chr_to_genome <- kar[[k_gen]]
names(chr_to_genome) <- kar[[k_chr]]

ok <- blocks[[config$columns$chrom1]] %in% sectors &
    blocks[[config$columns$chrom2]] %in% sectors
blocks <- blocks[ok, , drop = FALSE]
if (!nrow(blocks)) {
    stop("No synteny blocks match karyotype chromosomes.", call. = FALSE)
}

init <- data.frame(
    chr = kar[[k_chr]],
    start = kar[[k_st]],
    end = kar[[k_en]],
    stringsAsFactors = FALSE
)

n_per <- as.integer(table(factor(kar[[k_gen]], levels = genomes)))
gap_after <- unlist(lapply(n_per, function(n) {
    c(rep(config$gap_within, n - 1L), config$gap_between)
}), use.names = FALSE)

grp_levels <- unique(blocks[[config$columns$group]])
rib_cols <- palette_colors(config$palettes$ribbon, n = length(grp_levels))
names(rib_cols) <- grp_levels
gen_cols <- palette_colors(config$palettes$genome, n = length(genomes))
names(gen_cols) <- genomes
outline <- palette_colors("Qualitative.Safe")[[12]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_circos <- function() {
    circlize::circos.par(
        start.degree = 90,
        gap.after = gap_after,
        cell.padding = c(0, 0, 0, 0),
        track.margin = c(0.002, 0.008),
        circle.margin = c(0.04, 0.04, 0.14, 0.04),
        points.overflow.warning = FALSE
    )
    circlize::circos.genomicInitialize(
        init,
        plotType = NULL,
        tickLabelsStartFromZero = TRUE
    )
    circlize::circos.track(
        ylim = c(0, 1),
        track.height = 0.06,
        bg.border = NA,
        panel.fun = function(x, y) {
            circlize::circos.text(
                circlize::CELL_META$xcenter,
                0.45,
                circlize::CELL_META$sector.index,
                facing = "bending.inside",
                niceFacing = TRUE,
                cex = 0.72
            )
        }
    )
    circlize::circos.track(
        ylim = c(0, 1),
        track.height = 0.08,
        bg.border = NA,
        panel.fun = function(x, y) {
            sector <- circlize::CELL_META$sector.index
            g <- unname(chr_to_genome[[sector]])
            xlim <- circlize::CELL_META$xlim
            circlize::circos.rect(
                xlim[[1]], 0, xlim[[2]], 1,
                col = gen_cols[[g]],
                border = outline,
                lwd = 0.3
            )
        }
    )
    bed1 <- data.frame(
        chr = blocks[[config$columns$chrom1]],
        start = blocks[[config$columns$start1]],
        end = blocks[[config$columns$end1]],
        stringsAsFactors = FALSE
    )
    bed2 <- data.frame(
        chr = blocks[[config$columns$chrom2]],
        start = blocks[[config$columns$start2]],
        end = blocks[[config$columns$end2]],
        stringsAsFactors = FALSE
    )
    col <- grDevices::adjustcolor(
        unname(rib_cols[blocks[[config$columns$group]]]),
        alpha.f = 0.55
    )
    circlize::circos.genomicLink(
        bed1, bed2,
        col = col,
        border = NA
    )
    circlize::circos.clear()

    lgd_g <- ComplexHeatmap::Legend(
        title = config$labels$genome,
        at = genomes,
        legend_gp = grid::gpar(fill = unname(gen_cols[genomes])),
        title_gp = grid::gpar(fontsize = 8),
        labels_gp = grid::gpar(fontsize = 7),
        background = NA
    )
    ComplexHeatmap::draw(
        lgd_g,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0.02, "npc"),
        just = c("center", "bottom")
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
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
