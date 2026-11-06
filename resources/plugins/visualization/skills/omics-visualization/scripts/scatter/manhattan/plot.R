#!/usr/bin/env Rscript

# Template-ID: scatter-manhattan
#
# Purpose:
#   Genome-wide association Manhattan plot. Chromosomes are laid end to
#   end on x, y = -log10(p). Genome-wide and suggestive threshold lines,
#   highlighted SNPs (past the threshold, from a list, or inside
#   regions), and gene / SNP labels with leader lines.
#
# Inputs:
#   A SNP table (tsv/csv), one row per variant (per trait). Default example:
#     - snp: variant id
#     - chr: chromosome (1..22, X, Y, MT; a "chr" prefix is stripped)
#     - bp: physical position
#     - p: supplied association p-value (or a score with logp = FALSE)
#     - trait: optional; one panel per trait (example_traits.tsv)
#     - gene: optional label column; non-empty rows are labelled
#   Optional sidecars next to the input (CONFIG gene_file, regions_file):
#     - genes.tsv (snp, gene): SNPs to highlight and label
#     - regions.tsv (chr, start, end): candidate regions boxed across
#       every panel and/or used to highlight the SNPs inside
#
# Layouts and options:
#   layout = "auto" draws one panel, or one stacked panel per trait.
#   layout = "miami" mirrors two traits: the first up, the second down.
#   layout = "overlay" draws every trait in one panel, coloured by trait,
#   with significant SNPs enlarged in their trait colour. Optional
#   chr_bands and threshold_label support multi-tissue haQTL overlays
#   (preview_haqtl.png). Optional
#   chr_bands and threshold_label support multi-tissue haQTL overlays
#   (preview_haqtl.png).
#   layout = "circular" bends the genome into a ring, one concentric
#   track per trait (gene / SNP labels are not drawn in this layout).
#   signal_tiers = TRUE colours SNPs between the suggestive and the
#   genome-wide line with a second colour.
#   qq_inset = TRUE adds an observed vs expected -log10(p) QQ inset with
#   a 95% band and lambda (single panel or overlay layouts).
#   color_mode = "value" fades points from grey to orange by y.
#   density_track = TRUE adds a SNP-density strip above the plot.
#   label_mode = "top_per_chr" labels the top SNP of each chromosome
#   past the suggestive line (else the genome-wide line).
#
# Output:
#   A PDF, PNG, or SVG Manhattan plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggrepel, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#   Rscript plot.R example_traits.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, thresholds, layout, colour mode,
#   highlight and label modes, sidecars). Use scatter-qq-pvalue for the
#   matching QQ plot.
#
# Scientific assumptions:
#   p-values come from GAPIT, rMVP, GEMMA, GCTA, PLINK, REGENIE, ...
#   This script runs no association test. "bonferroni" only draws a
#   line at 0.05 / number of SNPs. With logp = FALSE the value column is
#   plotted as is and thresholds are read on that scale. thin_p_above
#   drops a deterministic share of non-significant points for display;
#   points below it are always kept. SNP density is a count per window.

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
        snp = "snp",
        chr = "chr",
        bp = "bp",
        p = "p",
        trait = "trait",
        gene = "gene"
    ),
    # Sidecars relative to the input; NULL or a missing file = none
    gene_file = "genes.tsv",
    regions_file = NULL,
    labels = list(
        title = NULL,
        x = "Chromosome",
        y = NULL                # NULL: -log10(p); e.g. "-log10(FDR)"
    ),
    # "auto" (one panel per trait), "miami" (exactly two traits),
    # "overlay" (all traits in one panel), or "circular" (ring per trait)
    layout = "auto",
    trait_palette = "Qualitative.Safe",
    traits = NULL,              # subset / order of traits; NULL keeps all
    # TRUE: y = -log10(p). FALSE: plot the value column as is
    logp = TRUE,
    # A p-value, "bonferroni" (0.05 / N), or NULL for no line
    genomewide = 1e-8,
    suggestive = NULL,
    # Optional text next to the genome-wide line (e.g. "emp P 0.005")
    threshold_label = NULL,
    # Alternating chromosome background bands (cartesian layouts)
    chr_bands = FALSE,
    chr_band_colors = c("#EFEFEF", "#FFFFFF"),
    # overlay only: p >= overlay_sig_p drawn grey; p < coloured by trait.
    # NULL = colour every point by trait.
    overlay_sig_p = NULL,
    overlay_ns_color = "#8A8A8A",
    ymax = NULL,                # fixed y upper limit; NULL = from data
    # "per_chr": one palette colour per chromosome, cycled
    # "alternate": two colours alternating between chromosomes
    # "value": grey (low) to orange (high) by y
    color_mode = "per_chr",
    palette = "Qualitative.Prism",
    alternate_colors = c("#1f78b4", "#a6cee3"),
    value_colors = c("#b3b3b3", "#e68310"),
    # "both": past genomewide + SNPs with a gene / in highlight_snps
    # "threshold", "list", or "none"
    highlight_mode = "both",
    highlight_snps = character(),
    highlight_color = "#e41a1c",
    # TRUE: SNPs past suggestive but not genomewide get tier_color
    signal_tiers = FALSE,
    tier_color = "#33a02c",
    # Regions as "chr3:1-30000000" strings, added to regions_file rows
    highlight_regions = character(),
    region_box = TRUE,          # outline regions across every panel
    region_points = FALSE,      # highlight SNPs inside regions
    region_color = "#e41a1c",
    # "genes", "top_per_chr", "both", or "none"
    label_mode = "genes",
    genomewide_color = "#e41a1c",
    suggestive_color = "#1f78b4",
    density_track = FALSE,
    density_window = 10e6,
    density_palette = "Gradient.YlOrRd",
    density_position = "top",   # "top" or "bottom" (cartesian layouts)
    outward = TRUE,             # circular: larger -log10(p) points outward
    circular_hole = 2.4,        # circular: inner hole, in ring heights
    chr_gap = 0.012,            # gap between chromosomes, share of genome length
    point_size = 1.3,
    point_shape = "solid",      # "solid" or "hollow"
    qq_inset = FALSE,
    qq_inset_box = c(0.76, 0.52, 0.995, 0.99),  # left, bottom, right, top
    highlight_size = 2.4,
    thin_p_above = NULL,        # e.g. 0.01 thins points with p above it
    thin_keep = 0.2,            # share of those points kept
    width = 12,
    height = NULL               # NULL: sized from panel rows and tracks
)

load_packages(c("ggplot2", "readr", "ggprism", "ggrepel", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
cols <- config$columns
has_gene <- !is.null(cols$gene) && cols$gene %in% names(df)
has_trait <- !is.null(cols$trait) && cols$trait %in% names(df)
require_columns(df, cols[c("snp", "chr", "bp", "p")])

d <- data.frame(
    snp = as.character(df[[cols$snp]]),
    chr = sub("^chr", "", as.character(df[[cols$chr]]), ignore.case = TRUE),
    bp = as.numeric(df[[cols$bp]]),
    p = as.numeric(df[[cols$p]]),
    trait = if (has_trait) as.character(df[[cols$trait]]) else "all",
    gene = if (has_gene) as.character(df[[cols$gene]]) else "",
    stringsAsFactors = FALSE
)
d$gene[is.na(d$gene)] <- ""
d <- d[is.finite(d$p) & is.finite(d$bp), , drop = FALSE]
if (isTRUE(config$logp)) {
    d <- d[d$p > 0, , drop = FALSE]
    d$y <- -log10(d$p)
} else {
    d$y <- d$p
}

sidecar_path <- function(name) {
    if (is.null(name)) return(NULL)
    path <- if (file.exists(name)) name else file.path(dirname(io$input), name)
    if (file.exists(path)) path else NULL
}
gene_path <- sidecar_path(config$gene_file)
if (!is.null(gene_path)) {
    genes <- read_table_auto(gene_path)
    require_columns(genes, list(snp = "snp", gene = "gene"))
    hit <- match(d$snp, as.character(genes$snp))
    d$gene[!is.na(hit)] <- as.character(genes$gene[hit[!is.na(hit)]])
}

trait_levels <- if (is.null(config$traits)) unique(d$trait) else config$traits
d <- d[d$trait %in% trait_levels, , drop = FALSE]
d$trait <- factor(d$trait, levels = trait_levels)
is_miami <- identical(config$layout, "miami")
is_circular <- identical(config$layout, "circular")
is_overlay <- identical(config$layout, "overlay")
if (is_miami && length(trait_levels) != 2L) {
    stop("layout = \"miami\" needs exactly two traits.", call. = FALSE)
}

n_snp <- length(unique(d$snp))
to_y <- function(value) {
    if (is.null(value)) return(NULL)
    if (identical(value, "bonferroni")) value <- 0.05 / n_snp
    if (isTRUE(config$logp)) -log10(as.numeric(value)) else as.numeric(value)
}
gw_y <- to_y(config$genomewide)
sg_y <- to_y(config$suggestive)

chr_rank <- suppressWarnings(as.numeric(d$chr))
special <- c(X = 1001, Y = 1002, XY = 1003, MT = 1004, M = 1004)
chr_rank[is.na(chr_rank)] <- special[toupper(d$chr[is.na(chr_rank)])]
chr_rank[is.na(chr_rank)] <- 2000
chr_levels <- unique(d$chr[order(chr_rank)])
d$chr <- factor(d$chr, levels = chr_levels)

chr_len <- tapply(d$bp, d$chr, max)
gap <- config$chr_gap * sum(chr_len)
chr_offset <- c(0, cumsum(as.numeric(chr_len) + gap))[seq_along(chr_len)]
names(chr_offset) <- names(chr_len)
d$pos <- d$bp + chr_offset[as.character(d$chr)]
axis_df <- data.frame(
    chr = names(chr_len),
    mid = chr_offset + as.numeric(chr_len) / 2
)
x_limits <- c(0, max(chr_offset + as.numeric(chr_len)))

region_rows <- list()
reg_path <- sidecar_path(config$regions_file)
if (!is.null(reg_path)) {
    reg <- read_table_auto(reg_path)
    require_columns(reg, list(chr = "chr", start = "start", end = "end"))
    region_rows[[1]] <- data.frame(chr = as.character(reg$chr),
        start = as.numeric(reg$start), end = as.numeric(reg$end))
}
if (length(config$highlight_regions)) {
    parts <- regmatches(config$highlight_regions,
        regexec("^(?:chr)?([^:]+):([0-9.e+]+)-([0-9.e+]+)$",
            config$highlight_regions, perl = TRUE))
    region_rows[[2]] <- do.call(rbind, lapply(parts, function(m) {
        data.frame(chr = m[[2]], start = as.numeric(m[[3]]), end = as.numeric(m[[4]]))
    }))
}
regions <- do.call(rbind, region_rows)
in_region <- rep(FALSE, nrow(d))
if (!is.null(regions) && nrow(regions)) {
    regions$chr <- sub("^chr", "", regions$chr, ignore.case = TRUE)
    regions <- regions[regions$chr %in% chr_levels, , drop = FALSE]
    regions$xmin <- regions$start + chr_offset[regions$chr]
    regions$xmax <- regions$end + chr_offset[regions$chr]
    for (i in seq_len(nrow(regions))) {
        in_region <- in_region | (as.character(d$chr) == regions$chr[[i]] &
            d$bp >= regions$start[[i]] & d$bp <= regions$end[[i]])
    }
}

is_threshold <- if (is.null(gw_y)) rep(FALSE, nrow(d)) else d$y >= gw_y
is_tier2 <- rep(FALSE, nrow(d))
if (isTRUE(config$signal_tiers) && !is.null(sg_y)) {
    is_tier2 <- d$y >= sg_y & !is_threshold
}
is_listed <- nzchar(d$gene) | d$snp %in% config$highlight_snps
d$highlight <- switch(config$highlight_mode,
    both = is_threshold | is_tier2 | is_listed,
    threshold = is_threshold | is_tier2,
    list = is_listed,
    rep(FALSE, nrow(d))
) | (isTRUE(config$region_points) & in_region)
tier2_only <- is_tier2 & !is_listed & !(isTRUE(config$region_points) & in_region)
d$hl_fill <- ifelse(tier2_only, config$tier_color, config$highlight_color)
d$hl_border <- ifelse(tier2_only, "#0b4d0b", "#5a0000")

qq_source <- d[, c("trait", "p")]
if (!is.null(config$thin_p_above)) {
    thin_cut <- to_y(config$thin_p_above)
    thin_pool <- which(d$y < thin_cut & !d$highlight)
    step <- max(1L, round(1 / config$thin_keep))
    drop <- thin_pool[seq_along(thin_pool) %% step != 0]
    if (length(drop)) d <- d[-drop, , drop = FALSE]
}

label_line <- if (!is.null(sg_y)) sg_y else gw_y
gene_rows <- nzchar(d$gene) &
    (d$highlight | identical(config$highlight_mode, "none"))
top_rows <- rep(FALSE, nrow(d))
if (config$label_mode %in% c("top_per_chr", "both") && !is.null(label_line)) {
    key <- paste(d$trait, d$chr)
    best <- tapply(seq_len(nrow(d)), key, function(i) i[which.max(d$y[i])])
    best <- unlist(best)
    top_rows[best[d$y[best] >= label_line]] <- TRUE
}
label_rows <- switch(config$label_mode,
    genes = gene_rows,
    top_per_chr = top_rows,
    both = gene_rows | top_rows,
    rep(FALSE, nrow(d))
)
labels_df <- d[label_rows, , drop = FALSE]
labels_df$label <- ifelse(nzchar(labels_df$gene), labels_df$gene, labels_df$snp)
labels_df$italic <- nzchar(labels_df$gene)

y_top <- if (!is.null(config$ymax)) config$ymax else
    max(c(d$y, gw_y, sg_y), na.rm = TRUE)

sign_of <- function(trait) if (is_miami) ifelse(trait == trait_levels[[2]], -1, 1) else 1
d$y_plot <- d$y * sign_of(d$trait)
labels_df$y_plot <- labels_df$y * sign_of(labels_df$trait)

chr_cols <- switch(config$color_mode,
    alternate = rep_len(config$alternate_colors, length(chr_levels)),
    rep_len(palette_colors(config$palette), length(chr_levels))
)
names(chr_cols) <- chr_levels
trait_cols <- setNames(palette_colors(config$trait_palette, length(trait_levels)),
    trait_levels)
if (is_overlay) {
    over_hl <- d$hl_fill == config$highlight_color
    d$hl_fill[over_hl] <- trait_cols[as.character(d$trait[over_hl])]
    d$hl_border[over_hl] <- "#222222"
}
pt_shape <- if (identical(config$point_shape, "hollow")) 1 else 19
pt_stroke <- if (identical(config$point_shape, "hollow")) 0.45 else 0
y_title <- if (!is.null(config$labels$y)) config$labels$y else
    if (isTRUE(config$logp)) expression(-log[10](italic(p))) else "Value"

dens <- NULL
if (isTRUE(config$density_track)) {
    snp_pos <- unique(d[, c("snp", "chr", "bp")])
    snp_pos$win <- floor(snp_pos$bp / config$density_window)
    dens <- aggregate(snp ~ chr + win, data = snp_pos, FUN = length)
    names(dens)[names(dens) == "snp"] <- "count"
    dens$xmin <- dens$win * config$density_window + chr_offset[as.character(dens$chr)]
    dens$xmax <- pmin((dens$win + 1) * config$density_window,
        as.numeric(chr_len[as.character(dens$chr)])) +
        chr_offset[as.character(dens$chr)]
}
density_scale <- function() {
    scale_fill_gradientn(colours = palette_colors(config$density_palette),
        name = "SNPs /\nwindow")
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
n_rows <- 1L
if (is_circular) {
    n_tr <- length(trait_levels)
    ring_h <- 1
    ring_gap <- 0.3
    wedge <- 0.05 * x_limits[[2]]
    base_r <- config$circular_hole + (seq_len(n_tr) - 1) * (ring_h + ring_gap)
    names(base_r) <- trait_levels
    to_r <- function(y, trait) {
        frac <- pmin(y, y_top) / y_top * 0.95
        b <- base_r[as.character(trait)]
        if (isTRUE(config$outward)) b + frac * ring_h else b + ring_h - frac * ring_h
    }
    d$r <- to_r(d$y, d$trait)
    outer_r <- max(base_r) + ring_h
    band_lo <- outer_r + 0.15
    band_hi <- band_lo + 0.22
    dens_lo <- band_hi + 0.08
    dens_hi <- dens_lo + 0.3
    label_r <- (if (is.null(dens)) band_hi else dens_hi) + 0.35

    arc <- function(r, trait, kind) {
        xs <- seq(x_limits[[1]], x_limits[[2]], length.out = 360)
        data.frame(x = xs, r = unname(r), trait = trait, kind = kind)
    }
    guide_arcs <- do.call(rbind, lapply(trait_levels, function(tr) {
        rbind(
            arc(base_r[[tr]] + (if (isTRUE(config$outward)) 0 else ring_h), tr, "base"),
            if (!is.null(sg_y)) arc(to_r(sg_y, tr), tr, "suggestive"),
            if (!is.null(gw_y)) arc(to_r(gw_y, tr), tr, "genomewide")
        )
    }))
    chr_band <- data.frame(chr = factor(names(chr_len), levels = chr_levels),
        xmin = chr_offset, xmax = chr_offset + as.numeric(chr_len))

    p <- ggplot() +
        geom_rect(data = chr_band,
            aes(xmin = xmin, xmax = xmax, ymin = band_lo, ymax = band_hi),
            fill = unname(chr_cols[as.character(chr_band$chr)]), colour = NA) +
        geom_path(data = guide_arcs[guide_arcs$kind == "base", ],
            aes(x = x, y = r, group = trait), colour = "#c8c8c8", linewidth = 0.3)
    base_pts <- d[!d$highlight, , drop = FALSE]
    if (identical(config$color_mode, "value")) {
        p <- p +
            geom_point(data = base_pts, aes(x = pos, y = r, colour = y, alpha = y),
                size = config$point_size * 0.55, stroke = 0) +
            scale_colour_gradient(low = config$value_colors[[1]],
                high = config$value_colors[[2]], guide = "none") +
            scale_alpha(range = c(0.25, 1), guide = "none")
    } else {
        p <- p +
            geom_point(data = base_pts, aes(x = pos, y = r, colour = chr),
                size = config$point_size * 0.55, shape = pt_shape,
                stroke = pt_stroke * 0.6) +
            scale_colour_manual(values = chr_cols, guide = "none")
    }
    for (kind in c("suggestive", "genomewide")) {
        arcs <- guide_arcs[guide_arcs$kind == kind, , drop = FALSE]
        if (!nrow(arcs)) next
        p <- p + geom_path(data = arcs, aes(x = x, y = r, group = trait),
            linetype = "dashed", linewidth = 0.4,
            colour = if (kind == "genomewide") config$genomewide_color else
                config$suggestive_color)
    }
    hl_pts <- d[d$highlight, , drop = FALSE]
    for (fill_hex in unique(hl_pts$hl_fill)) {
        grp <- hl_pts[hl_pts$hl_fill == fill_hex, , drop = FALSE]
        p <- p + geom_point(data = grp, aes(x = pos, y = r), shape = 21,
            size = config$highlight_size * 0.75, stroke = 0.3,
            fill = fill_hex, colour = grp$hl_border[[1]])
    }
    p <- p +
        geom_text(data = axis_df, aes(x = mid, y = label_r, label = chr),
            size = 3.6, colour = "#222222") +
        geom_text(
            data = data.frame(label = trait_levels,
                r = unname(base_r) + ring_h / 2),
            aes(x = x_limits[[1]] - wedge / 2, y = r, label = label),
            size = 3.4, fontface = "bold", colour = "#222222"
        )
    if (!is.null(dens)) {
        p <- p +
            geom_rect(data = dens,
                aes(xmin = xmin, xmax = xmax, ymin = dens_lo, ymax = dens_hi,
                    fill = count)) +
            density_scale()
    }
    p <- p +
        scale_x_continuous(limits = c(x_limits[[1]] - wedge, x_limits[[2]]),
            expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, label_r + 0.3), expand = c(0, 0)) +
        coord_polar(theta = "x", start = -pi * wedge / (x_limits[[2]] + wedge),
            clip = "off") +
        labs(title = config$labels$title) +
        theme_void() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "inside",
            legend.position.inside = c(0.5, 0.5),
            legend.direction = "horizontal",
            legend.title.position = "top",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            legend.key.height = unit(2.5, "mm"),
            legend.key.width = unit(6, "mm")
        )
} else {
    step_y <- max(1, ceiling(y_top / 6 / 2) * 2)
    y_breaks <- seq(0, ceiling(y_top), by = step_y)
    if (is_miami) y_breaks <- sort(unique(c(-y_breaks, y_breaks)))
    y_range <- if (is_miami) c(-y_top, y_top) else c(0, y_top)

    p <- ggplot(d, aes(x = pos, y = y_plot))
    if (isTRUE(config$chr_bands)) {
        band_cols <- rep_len(config$chr_band_colors, length(chr_levels))
        band_df <- data.frame(
            xmin = unname(chr_offset),
            xmax = unname(chr_offset) + as.numeric(chr_len),
            fill = band_cols,
            stringsAsFactors = FALSE
        )
        for (i in seq_len(nrow(band_df))) {
            p <- p + annotate(
                "rect",
                xmin = band_df$xmin[[i]],
                xmax = band_df$xmax[[i]],
                ymin = -Inf,
                ymax = Inf,
                fill = band_df$fill[[i]],
                colour = NA
            )
        }
    }
    if (!is.null(regions) && nrow(regions) && isTRUE(config$region_box)) {
        p <- p + geom_rect(
            data = regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = config$region_color, alpha = 0.06,
            colour = config$region_color, linewidth = 0.6
        )
    }
    base_pts <- d[!d$highlight, , drop = FALSE]
    if (is_overlay) {
        base_pts <- base_pts[order(base_pts$y), , drop = FALSE]
        if (!is.null(config$overlay_sig_p)) {
            sig_cut <- to_y(config$overlay_sig_p)
            ns_pts <- base_pts[base_pts$y < sig_cut, , drop = FALSE]
            sig_pts <- base_pts[base_pts$y >= sig_cut, , drop = FALSE]
            if (nrow(ns_pts)) {
                p <- p + geom_point(
                    data = ns_pts,
                    colour = config$overlay_ns_color,
                    size = config$point_size,
                    shape = pt_shape,
                    stroke = pt_stroke,
                    show.legend = FALSE
                )
            }
            colour_pts <- sig_pts
        } else {
            colour_pts <- base_pts
        }
        if (nrow(colour_pts)) {
            p <- p +
                geom_point(
                    data = colour_pts, aes(colour = trait),
                    size = config$point_size, shape = pt_shape, stroke = pt_stroke
                )
        } else {
            # Keep a colour scale for legend when all sig points are highlights.
            p <- p + geom_point(
                data = d[0, , drop = FALSE], aes(colour = trait),
                size = config$point_size
            )
        }
        p <- p + scale_colour_manual(
            values = trait_cols, name = NULL, drop = FALSE,
            guide = guide_legend(override.aes = list(size = 3.5, shape = 15))
        )
    } else if (identical(config$color_mode, "value")) {
        p <- p +
            geom_point(data = base_pts, aes(colour = y, alpha = y),
                size = config$point_size, shape = pt_shape, stroke = pt_stroke) +
            scale_colour_gradient(low = config$value_colors[[1]],
                high = config$value_colors[[2]], guide = "none") +
            scale_alpha(range = c(0.25, 1), guide = "none")
    } else {
        p <- p +
            geom_point(data = base_pts, aes(colour = chr),
                size = config$point_size, shape = pt_shape, stroke = pt_stroke) +
            scale_colour_manual(values = chr_cols, guide = "none")
    }
    line_signs <- if (is_miami) c(1, -1) else 1
    if (!is.null(sg_y)) {
        p <- p + geom_hline(yintercept = sg_y * line_signs, linetype = "dashed",
            linewidth = 0.45, colour = config$suggestive_color)
    }
    if (!is.null(gw_y)) {
        p <- p + geom_hline(yintercept = gw_y * line_signs, linetype = "dashed",
            linewidth = 0.45, colour = config$genomewide_color)
        if (!is.null(config$threshold_label) && nzchar(config$threshold_label)) {
            p <- p + annotate(
                "text",
                x = x_limits[[1]],
                y = gw_y,
                label = config$threshold_label,
                hjust = 0,
                vjust = -0.6,
                size = 3.2,
                colour = "grey25"
            )
        }
    }
    if (is_miami) {
        p <- p + geom_hline(yintercept = 0, linewidth = 0.4, colour = "#222222")
    }
    hl_pts <- d[d$highlight, , drop = FALSE]
    for (fill_hex in unique(hl_pts$hl_fill)) {
        grp <- hl_pts[hl_pts$hl_fill == fill_hex, , drop = FALSE]
        p <- p + geom_point(
            data = grp, shape = 21, size = config$highlight_size, stroke = 0.35,
            fill = fill_hex, colour = grp$hl_border[[1]]
        )
    }
    for (direction in c(1, -1)) {
        sub <- labels_df[sign(labels_df$y_plot) == direction |
            (direction == 1 & labels_df$y_plot == 0), , drop = FALSE]
        if (!nrow(sub)) next
        p <- p + geom_text_repel(
            data = sub, aes(label = label,
                fontface = ifelse(italic, "italic", "plain")),
            size = 3.6, colour = "#222222",
            nudge_y = direction * 0.06 * y_top, nudge_x = 0.015 * x_limits[[2]],
            segment.colour = "#444444", segment.size = 0.35,
            min.segment.length = 0, box.padding = 0.4, seed = 1
        )
    }
    p <- p +
        scale_x_continuous(breaks = axis_df$mid, labels = axis_df$chr,
            limits = x_limits, expand = expansion(mult = 0.01)) +
        scale_y_continuous(breaks = y_breaks, labels = abs(y_breaks),
            expand = expansion(mult = c(if (is_miami) 0.08 else 0, 0.1))) +
        coord_cartesian(ylim = y_range, clip = "off") +
        labs(title = config$labels$title, x = config$labels$x,
            y = y_title) +
        theme_prism(base_size = 12) +
        theme(
            axis.title.x = element_text(margin = margin(t = 6)),
            axis.text.x = element_text(size = 10),
            legend.position = "inside",
            legend.position.inside = c(0.995, 0.99),
            legend.justification = c(1, 1),
            legend.background = element_rect(fill = "white", colour = "#bdbdbd",
                linewidth = 0.3),
            legend.text = element_text(size = 10, face = "plain"),
            legend.key.spacing.y = unit(0.5, "mm")
        )

    if (is_miami) {
        trait_tags <- data.frame(label = trait_levels, y = c(y_top, -y_top),
            vjust = c(1, 0))
        p <- p + geom_text(data = trait_tags,
            aes(x = x_limits[[1]], y = y, label = label, vjust = vjust),
            inherit.aes = FALSE, hjust = 0, size = 4.4, fontface = "bold")
    } else if (length(trait_levels) > 1L && !is_overlay) {
        n_rows <- length(trait_levels)
        p <- p + facet_grid(trait ~ ., switch = "y") +
            theme(
                strip.placement = "outside",
                strip.background = element_blank(),
                strip.text.y.left = element_text(angle = 90, face = "bold", size = 12),
                panel.spacing.y = unit(4, "mm")
            )
    }

    if (isTRUE(config$qq_inset) && isTRUE(config$logp) && n_rows == 1L &&
        !is_miami) {
        qq_df <- do.call(rbind, lapply(trait_levels, function(tr) {
            pv <- sort(qq_source$p[qq_source$trait == tr & qq_source$p > 0])
            n <- length(pv)
            i <- seq_len(n)
            out <- data.frame(trait = tr, expected = -log10(ppoints(n)),
                observed = -log10(pv),
                lo = -log10(stats::qbeta(0.975, i, n - i + 1)),
                hi = -log10(stats::qbeta(0.025, i, n - i + 1)))
            step <- max(1L, n %/% 2000L)
            out[out$expected > 2 | i %% step == 0, , drop = FALSE]
        }))
        qq_df$trait <- factor(qq_df$trait, levels = trait_levels)
        qq_xmax <- max(qq_df$expected)
        qq_ymax <- max(qq_df$observed)
        p_qq <- ggplot(qq_df, aes(x = expected, y = observed)) +
            geom_ribbon(data = qq_df[qq_df$trait == trait_levels[[1]], ],
                aes(ymin = lo, ymax = hi), fill = "#d9d9d9") +
            geom_abline(slope = 1, intercept = 0, colour = config$genomewide_color,
                linewidth = 0.4) +
            geom_point(aes(colour = trait), size = 0.9, shape = pt_shape,
                stroke = if (pt_shape == 1) 0.35 else 0) +
            scale_colour_manual(values = if (is_overlay) trait_cols else
                setNames(rep(config$alternate_colors[[1]], length(trait_levels)),
                    trait_levels), guide = "none") +
            coord_cartesian(xlim = c(0, qq_xmax), ylim = c(0, qq_ymax)) +
            labs(x = expression(Expected ~ -log[10](italic(p))),
                y = expression(Observed ~ -log[10](italic(p)))) +
            theme_prism(base_size = 8, border = TRUE) +
            theme(plot.background = element_rect(fill = "white", colour = NA),
                axis.title = element_text(size = 7.5))
        if (!is_overlay) {
            pv_all <- qq_source$p[qq_source$p > 0]
            lambda <- stats::median(stats::qchisq(1 - pv_all, 1)) /
                stats::qchisq(0.5, 1)
            p_qq <- p_qq + annotate("text", x = 0, y = qq_ymax, hjust = 0,
                vjust = 1, size = 2.8,
                label = sprintf("lambda == %.3f", lambda), parse = TRUE)
        }
        box <- config$qq_inset_box
        p <- p + inset_element(p_qq, left = box[[1]], bottom = box[[2]],
            right = box[[3]], top = box[[4]])
    }

    if (isTRUE(config$density_track)) {
        p_dens <- ggplot(dens) +
            geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = count)) +
            density_scale() +
            scale_x_continuous(limits = x_limits, expand = expansion(mult = 0.01)) +
            scale_y_continuous(expand = c(0, 0)) +
            labs(y = "Density") +
            theme_void() +
            theme(
                axis.title.y = element_text(angle = 90, size = 11,
                    margin = margin(r = 6)),
                legend.title = element_text(size = 9),
                legend.text = element_text(size = 8),
                legend.key.height = unit(3, "mm"),
                legend.key.width = unit(10, "mm"),
                legend.direction = "horizontal"
            )
        p <- if (identical(config$density_position, "bottom")) {
            p / p_dens + plot_layout(heights = c(n_rows, 0.08 * n_rows),
                guides = "collect")
        } else {
            p_dens / p + plot_layout(heights = c(0.08 * n_rows, n_rows),
                guides = "collect")
        }
        p <- p + plot_annotation(theme = theme(legend.position = "bottom"))
    }
}

plot_width <- if (is_circular && is.null(config$height)) 9 else config$width
plot_height <- if (!is.null(config$height)) config$height else
    if (is_circular) 9 else
    (if (n_rows > 1L) 1.2 + 2.5 * n_rows else 4.8) *
    (if (is_miami) 1.4 else 1) +
    (if (isTRUE(config$density_track)) 0.9 else 0)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
# patchwork 1.3 flags ggplot2 4 themes in plot_annotation() as invalid but
# still applies them; only that warning is muffled.
withCallingHandlers(
    save_ggplot(p, io$output, width = plot_width, height = plot_height),
    warning = function(w) {
        if (grepl("annotation$theme is not a valid theme", conditionMessage(w),
            fixed = TRUE)) invokeRestart("muffleWarning")
    }
)
