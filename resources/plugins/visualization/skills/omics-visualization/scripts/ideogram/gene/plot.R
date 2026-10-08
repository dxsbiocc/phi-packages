#!/usr/bin/env Rscript

# Template-ID: ideogram-gene
#
# Purpose:
#   Draw transcript models as rounded ideogram tracks: one row per
#   exon / CDS / UTR, one track per transcript.
#
# Inputs:
#   One row per feature. Default example is TP53 from hg19 refGene:
#     - transcript: RefSeq accession (track)
#     - type: 5UTR, CDS, 3UTR
#     - start, end: genomic coordinates
#
# Output:
#   A PDF, PNG, or SVG gene-structure figure.
#
# Dependencies:
#   ggplot2, readr, ggideogram (GitHub: dxsbiocc/ggideogram)
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, fill, flip). Horizontal
#   vs vertical is config$flip, not a second id. Do not add an id per
#   gene, per GTF type filter, or per rounded-corner radius.
#   Karyotype G-banding is ideogram-karyotype. Density along a
#   chromosome is ideogram-density. Genes of interest as points on
#   chromosomes are ideogram-loci.
#
# Scientific assumptions:
#   Features are supplied intervals, not predicted by this script.
#   Overlapping exon and CDS rows will overplot; prefer CDS + UTR or
#   exon alone. The backbone is the transcript span, not an intron
#   table. Coordinates are already on one reference (genomic or local).

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
        transcript = "transcript",
        type = "type",
        start = "start",
        end = "end"
    ),
    flip = TRUE,
    radius_pt = 3,
    width = 0.35,
    labels = list(
        title = "",
        fill = "Feature"
    )
)

load_packages(c("ggplot2", "readr", "ggideogram"))
patch_ggideogram_ggplot2()

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

tx <- config$columns$transcript
ty <- config$columns$type
st <- config$columns$start
en <- config$columns$end

df[[tx]] <- factor(df[[tx]], levels = unique(df[[tx]]))
df[[ty]] <- factor(df[[ty]], levels = unique(df[[ty]]))
df[[st]] <- as.numeric(df[[st]])
df[[en]] <- as.numeric(df[[en]])
swap <- df[[st]] > df[[en]]
if (any(swap)) {
    tmp <- df[[st]][swap]
    df[[st]][swap] <- df[[en]][swap]
    df[[en]][swap] <- tmp
}

fill_cols <- palette_colors("Qualitative.Safe", n = nlevels(df[[ty]]))
names(fill_cols) <- levels(df[[ty]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot2::ggplot(df) +
    ggideogram::geom_ideogram(
        ggplot2::aes(
            x = .data[[tx]],
            ymin = .data[[st]],
            ymax = .data[[en]],
            chrom = .data[[tx]],
            fill = .data[[ty]]
        ),
        width = config$width,
        radius = grid::unit(config$radius_pt, "pt"),
        chrom.col = "black",
        linewidth = 0.3
    ) +
    ggplot2::scale_fill_manual(values = fill_cols, name = config$labels$fill) +
    ggplot2::labs(title = config$labels$title, x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(colour = "black", size = 10),
        plot.margin = ggplot2::margin(8, 12, 8, 8)
    )

if (isTRUE(config$flip)) {
    p <- p + ggplot2::coord_flip()
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 5.2)
