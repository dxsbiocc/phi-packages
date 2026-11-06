#!/usr/bin/env Rscript

# Template-ID: ideogram-synteny-linear
#
# Purpose:
#   Draw a horizontal multi-genome synteny (riparian) plot: one
#   chromosome bar track per genome, with ribbons between adjacent
#   genomes for supplied homology blocks. Large rearrangements can
#   be coloured separately via a supplied class column.
#
# Inputs:
#   1. Blocks table (one row per aligned interval pair):
#        genome1, chr1, start1, end1, genome2, chr2, start2, end2
#        optional class (synteny | rearrangement)
#   2. Sidecar karyotype.tsv beside the blocks file:
#        genome, chr, start, end
#
# Output:
#   A PDF, PNG, or SVG linear synteny figure.
#
# Dependencies:
#   ggplot2, readr, ggprism, dplyr (optional not required)
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for column names, genome order/colours, ribbon
#   colours, chromosome gap, and sizes. Extra genomes are extra
#   karyotype levels + blocks between adjacent pairs — not a new
#   id. Do not run BLAST, minimap, MCScanX, WGDI, or JCVI here.
#
# Scientific assumptions:
#   Chromosome lengths and homology blocks are supplied. Ribbon
#   class (synteny vs rearrangement) is supplied, not inferred from
#   alignment. Ribbons are drawn only between adjacent genome rows
#   in CONFIG genome_order. start2 > end2 is allowed and means an
#   inverted block on genome2 (drawn with swapped corners).

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
        genome1 = "genome1",
        chrom1 = "chr1",
        start1 = "start1",
        end1 = "end1",
        genome2 = "genome2",
        chrom2 = "chr2",
        start2 = "start2",
        end2 = "end2",
        class = "class"
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
    # Top → bottom. NULL = first-appearance order in karyotype.
    genome_order = c(
        "Potato (DM)",
        "Sam (SP1102)",
        "Eggplant (HQ-1315)"
    ),
    # Qualitative.Bold / Brand.Algolia hexes — do not invent colours.
    genome_colours = c(
        "Potato (DM)" = "#3969ac",
        "Sam (SP1102)" = "#e68310",
        "Eggplant (HQ-1315)" = "#11a579"
    ),
    ribbon_colours = c(
        synteny = "#888888",
        rearrangement = "#e68310"
    ),
    ribbon_alpha = c(
        synteny = 0.28,
        rearrangement = 0.55
    ),
    chr_gap = 0.035,
    track_height = 0.055,
    track_gap = 0.34,
    chr_label_size = 3.0,
    genome_label_size = 3.6,
    labels = list(
        title = ""
    ),
    size = list(width = 12.5, height = 5.2)
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
blocks <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
req_block <- config$columns[c(
    "genome1", "chrom1", "start1", "end1",
    "genome2", "chrom2", "start2", "end2"
)]
require_columns(blocks, req_block)

kar_path <- file.path(dirname(io$input), config$sidecars$karyotype)
if (!file.exists(kar_path)) {
    stop(paste("Missing karyotype sidecar:", kar_path), call. = FALSE)
}
kar <- as.data.frame(read_table_auto(kar_path), stringsAsFactors = FALSE)
require_columns(kar, config$karyotype)

g1c <- config$columns$genome1
c1c <- config$columns$chrom1
s1c <- config$columns$start1
e1c <- config$columns$end1
g2c <- config$columns$genome2
c2c <- config$columns$chrom2
s2c <- config$columns$start2
e2c <- config$columns$end2
cls_c <- config$columns$class

kg <- config$karyotype$genome
kc <- config$karyotype$chrom
ks <- config$karyotype$start
ke <- config$karyotype$end

kar[[kg]] <- as.character(kar[[kg]])
kar[[kc]] <- as.character(kar[[kc]])
kar[[ks]] <- as.numeric(kar[[ks]])
kar[[ke]] <- as.numeric(kar[[ke]])
ok_k <- is.finite(kar[[ks]]) & is.finite(kar[[ke]]) & kar[[ke]] > kar[[ks]]
if (any(!ok_k)) {
    message("Dropped ", sum(!ok_k), " karyotype row(s) with bad intervals.")
    kar <- kar[ok_k, , drop = FALSE]
}
if (!nrow(kar)) {
    stop("No karyotype rows left.", call. = FALSE)
}

blocks[[g1c]] <- as.character(blocks[[g1c]])
blocks[[c1c]] <- as.character(blocks[[c1c]])
blocks[[g2c]] <- as.character(blocks[[g2c]])
blocks[[c2c]] <- as.character(blocks[[c2c]])
s1_raw <- as.numeric(blocks[[s1c]])
e1_raw <- as.numeric(blocks[[e1c]])
s2_raw <- as.numeric(blocks[[s2c]])
e2_raw <- as.numeric(blocks[[e2c]])
# Invert flag: start2 > end2 on the second interval (JCVI / MCScanX style).
blocks$`.inv2` <- s2_raw > e2_raw
blocks[[s1c]] <- pmin(s1_raw, e1_raw)
blocks[[e1c]] <- pmax(s1_raw, e1_raw)
blocks[[s2c]] <- pmin(s2_raw, e2_raw)
blocks[[e2c]] <- pmax(s2_raw, e2_raw)

if (!is.null(cls_c) && cls_c %in% names(blocks)) {
    blocks[[cls_c]] <- as.character(blocks[[cls_c]])
    blocks[[cls_c]][is.na(blocks[[cls_c]]) | !nzchar(blocks[[cls_c]])] <- "synteny"
} else {
    cls_c <- ".class"
    blocks[[cls_c]] <- "synteny"
}

ok_b <- is.finite(blocks[[s1c]]) & is.finite(blocks[[e1c]]) &
    is.finite(blocks[[s2c]]) & is.finite(blocks[[e2c]]) &
    blocks[[e1c]] > blocks[[s1c]] & blocks[[e2c]] > blocks[[s2c]]
if (any(!ok_b)) {
    message("Dropped ", sum(!ok_b), " block(s) with bad intervals.")
    blocks <- blocks[ok_b, , drop = FALSE]
}
if (!nrow(blocks)) {
    stop("No blocks left to plot.", call. = FALSE)
}

g_order <- config$genome_order
if (is.null(g_order) || !length(g_order)) {
    g_order <- unique(kar[[kg]])
} else {
    g_order <- as.character(g_order)
    miss <- setdiff(unique(kar[[kg]]), g_order)
    g_order <- c(g_order[g_order %in% unique(kar[[kg]])], miss)
}
kar[[kg]] <- factor(kar[[kg]], levels = g_order)

# Chromosome order within each genome: numeric ChrN when possible.
chr_sort_key <- function(ch) {
    num <- suppressWarnings(as.numeric(sub(".*?([0-9]+)$", "\\1", ch)))
    ifelse(is.finite(num), num, match(ch, sort(unique(ch))))
}

# Layout chromosomes left→right per genome on a shared [0, 1] x scale.
chr_gap <- as.numeric(config$chr_gap)[[1]]
if (!is.finite(chr_gap) || chr_gap < 0) {
    chr_gap <- 0.03
}

layout_genome <- function(gdf) {
    gdf <- gdf[order(chr_sort_key(gdf[[kc]]), gdf[[kc]]), , drop = FALSE]
    lens <- gdf[[ke]] - gdf[[ks]]
    total <- sum(lens)
    n <- nrow(gdf)
    gap_total <- chr_gap * max(n - 1L, 0L)
    scale <- (1 - gap_total) / total
    x0 <- numeric(n)
    x1 <- numeric(n)
    cursor <- 0
    for (i in seq_len(n)) {
        w <- lens[[i]] * scale
        x0[[i]] <- cursor
        x1[[i]] <- cursor + w
        cursor <- x1[[i]] + chr_gap
    }
    data.frame(
        genome = as.character(gdf[[kg]]),
        chr = gdf[[kc]],
        chr_start = gdf[[ks]],
        chr_end = gdf[[ke]],
        x0 = x0,
        x1 = x1,
        stringsAsFactors = FALSE
    )
}

chr_layout <- do.call(rbind, lapply(g_order, function(g) {
    layout_genome(kar[as.character(kar[[kg]]) == g, , drop = FALSE])
}))

# Map a genomic coordinate to plot x within a chromosome bar.
map_x <- function(genome, chr, pos) {
    hit <- chr_layout$genome == genome & chr_layout$chr == chr
    if (!any(hit)) {
        return(NA_real_)
    }
    row <- chr_layout[which(hit)[[1]], , drop = FALSE]
    frac <- (pos - row$chr_start) / (row$chr_end - row$chr_start)
    frac <- pmin(pmax(frac, 0), 1)
    as.numeric(row$x0 + frac * (row$x1 - row$x0))
}

n_gen <- length(g_order)
track_h <- as.numeric(config$track_height)[[1]]
track_gap <- as.numeric(config$track_gap)[[1]]
if (!is.finite(track_h) || track_h <= 0) {
    track_h <- 0.05
}
if (!is.finite(track_gap) || track_gap <= 0) {
    track_gap <- 0.32
}

# y centres for genome tracks (top = first).
y_centre <- setNames(
    seq(1, 0, length.out = n_gen),
    g_order
)
chr_layout$y <- unname(y_centre[chr_layout$genome])
chr_layout$ymin <- chr_layout$y - track_h / 2
chr_layout$ymax <- chr_layout$y + track_h / 2

# Keep only blocks between adjacent genomes in genome_order.
adj_pairs <- cbind(g_order[-n_gen], g_order[-1L])
is_adj <- function(a, b) {
    any(adj_pairs[, 1] == a & adj_pairs[, 2] == b) ||
        any(adj_pairs[, 1] == b & adj_pairs[, 2] == a)
}
keep_adj <- vapply(seq_len(nrow(blocks)), function(i) {
    is_adj(blocks[[g1c]][[i]], blocks[[g2c]][[i]])
}, logical(1))
if (any(!keep_adj)) {
    message("Dropped ", sum(!keep_adj), " block(s) not spanning adjacent genomes.")
    blocks <- blocks[keep_adj, , drop = FALSE]
}
if (!nrow(blocks)) {
    stop("No adjacent-genome blocks left to plot.", call. = FALSE)
}

# Orient each block so genome A is the upper track.
orient_blocks <- blocks
swap <- y_centre[orient_blocks[[g1c]]] < y_centre[orient_blocks[[g2c]]]
if (any(swap)) {
    tmp <- orient_blocks[swap, , drop = FALSE]
    orient_blocks[swap, g1c] <- tmp[[g2c]]
    orient_blocks[swap, c1c] <- tmp[[c2c]]
    orient_blocks[swap, s1c] <- tmp[[s2c]]
    orient_blocks[swap, e1c] <- tmp[[e2c]]
    orient_blocks[swap, g2c] <- tmp[[g1c]]
    orient_blocks[swap, c2c] <- tmp[[c1c]]
    orient_blocks[swap, s2c] <- tmp[[s1c]]
    orient_blocks[swap, e2c] <- tmp[[e1c]]
    # Inversion flag follows the original genome2 endpoint order; after
    # swap, invert the upper interval presentation when needed.
    orient_blocks$`.inv2`[swap] <- FALSE
}

# Smooth sigmoid ribbon (left and right edges interpolated).
smooth_ribbon <- function(x1a, x1b, ya, x2a, x2b, yb, n = 24L) {
    t <- seq(0, 1, length.out = n)
    s <- t * t * (3 - 2 * t)
    x_left <- (1 - s) * x1a + s * x2a
    x_right <- (1 - s) * x1b + s * x2b
    y <- (1 - s) * ya + s * yb
    data.frame(
        x = c(x_left, rev(x_right)),
        y = c(y, rev(y))
    )
}

# Build ribbon polygons (invert flips bottom edge).
ribbon_list <- vector("list", nrow(orient_blocks))
for (i in seq_len(nrow(orient_blocks))) {
    ga <- orient_blocks[[g1c]][[i]]
    gb <- orient_blocks[[g2c]][[i]]
    ya <- y_centre[[ga]] - track_h / 2
    yb <- y_centre[[gb]] + track_h / 2
    x1a <- map_x(ga, orient_blocks[[c1c]][[i]], orient_blocks[[s1c]][[i]])
    x1b <- map_x(ga, orient_blocks[[c1c]][[i]], orient_blocks[[e1c]][[i]])
    x2a <- map_x(gb, orient_blocks[[c2c]][[i]], orient_blocks[[s2c]][[i]])
    x2b <- map_x(gb, orient_blocks[[c2c]][[i]], orient_blocks[[e2c]][[i]])
    if (any(!is.finite(c(x1a, x1b, x2a, x2b)))) {
        next
    }
    if (isTRUE(orient_blocks$`.inv2`[[i]])) {
        tmp <- x2a
        x2a <- x2b
        x2b <- tmp
    }
    cls <- orient_blocks[[cls_c]][[i]]
    poly <- smooth_ribbon(x1a, x1b, ya, x2a, x2b, yb)
    poly$id <- paste0("r", i)
    poly$class <- cls
    ribbon_list[[i]] <- poly
}
ribbons <- do.call(rbind, ribbon_list)
if (is.null(ribbons) || !nrow(ribbons)) {
    stop("Could not place any ribbons (chr ids mismatch karyotype?).", call. = FALSE)
}

# Genome fill colours.
gcols <- config$genome_colours
if (is.null(names(gcols))) {
    gcols <- setNames(
        palette_colors("Qualitative.Bold", n = max(n_gen, 3L))[seq_len(n_gen)],
        g_order
    )
} else {
    miss <- setdiff(g_order, names(gcols))
    if (length(miss)) {
        extra <- palette_colors("Qualitative.Bold", n = max(length(miss), 3L))
        gcols <- c(gcols, setNames(extra[seq_along(miss)], miss))
    }
    gcols <- gcols[g_order]
}
chr_layout$fill <- unname(gcols[chr_layout$genome])

rcols <- config$ribbon_colours
ralpha <- config$ribbon_alpha
cls_lv <- unique(ribbons$class)
for (cl in cls_lv) {
    if (is.null(rcols[[cl]]) && is.null(rcols[cl])) {
        rcols[[cl]] <- if (identical(cl, "rearrangement")) "#e68310" else "#888888"
    }
    if (is.null(ralpha[[cl]]) && is.null(ralpha[cl])) {
        ralpha[[cl]] <- if (identical(cl, "rearrangement")) 0.55 else 0.28
    }
}
ribbons$colour <- unname(unlist(rcols[ribbons$class]))
ribbons$alpha <- as.numeric(unlist(ralpha[ribbons$class]))

# Draw rearrangements on top of synteny.
ribbons$class <- factor(
    ribbons$class,
    levels = c(
        setdiff(unique(as.character(ribbons$class)), "rearrangement"),
        "rearrangement"
    )
)
ribbons <- ribbons[order(ribbons$class, ribbons$id), , drop = FALSE]

# Chromosome labels from the top genome.
top_g <- g_order[[1]]
chr_labs <- chr_layout[chr_layout$genome == top_g, , drop = FALSE]
chr_labs$xmid <- (chr_labs$x0 + chr_labs$x1) / 2
chr_labs$label <- sub("^Chr\\.?", "Chr. ", chr_labs$chr)
chr_labs$label <- sub("^Chr\\. Chr", "Chr.", chr_labs$label)
# Normalize Chr1 → Chr. 1
chr_labs$label <- vapply(chr_labs$chr, function(ch) {
    num <- sub(".*?([0-9]+)$", "\\1", ch)
    if (grepl("^[0-9]+$", num)) paste0("Chr. ", num) else ch
}, character(1))

genome_labs <- data.frame(
    genome = g_order,
    x = -0.02,
    y = unname(y_centre[g_order]),
    stringsAsFactors = FALSE
)

ink <- palette_colors("Brand.Emma")[[1]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_polygon(
        data = ribbons,
        aes(x = x, y = y, group = id),
        fill = ribbons$colour,
        alpha = ribbons$alpha,
        colour = NA
    ) +
    geom_rect(
        data = chr_layout,
        aes(xmin = x0, xmax = x1, ymin = ymin, ymax = ymax),
        fill = chr_layout$fill,
        colour = ink,
        linewidth = 0.25
    ) +
    geom_text(
        data = chr_labs,
        aes(x = xmid, y = ymax + 0.045, label = label),
        size = as.numeric(config$chr_label_size),
        colour = ink,
        vjust = 0
    ) +
    geom_text(
        data = genome_labs,
        aes(x = x, y = y, label = genome),
        hjust = 1,
        size = as.numeric(config$genome_label_size),
        colour = ink,
        fontface = "bold"
    ) +
    coord_cartesian(
        xlim = c(-0.28, 1.02),
        ylim = c(min(chr_layout$ymin) - 0.06, max(chr_layout$ymax) + 0.10),
        clip = "off"
    ) +
    labs(title = config$labels$title) +
    theme_void(base_size = 11) +
    theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(12, 12, 12, 12)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = config$size$width,
    height = config$size$height
)
