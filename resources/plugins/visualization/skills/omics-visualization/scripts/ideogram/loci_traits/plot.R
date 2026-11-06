#!/usr/bin/env Rscript

# Template-ID: ideogram-loci-traits
#
# Purpose:
#   Multi-trait loci on a circular ideogram: chromosome arcs and
#   cytobands outside, one band label per locus with its association
#   count, branching leaders to a column of trait markers per locus,
#   with one concentric marker ring per organ (shape = organ, fill =
#   trait) and an organ-grouped legend in the centre.
#
# Inputs:
#   One row per locus-by-trait association. Default example:
#     - chr, start, end: locus interval
#     - locus: locus id (rows sharing it are one locus)
#     - band: display label (e.g. cytoband "17q21.31")
#     - organ: marker ring / shape group
#     - trait: marker fill
#   Sidecar beside the input (same directory):
#     - cytoband.tsv: Chrom, Start, End, Stain, Arm
#
# Output:
#   A PDF, PNG, or SVG multi-trait loci circos.
#
# Dependencies:
#   circlize, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, organ shapes / palettes / ring order,
#   spacing, radii, chromosome set, label counts). Extra loci, organs
#   or traits are extra rows, not extra ids. hg19 vs hg38 is the
#   cytoband sidecar. Single-category gene labels are ideogram-loci.
#   Do not run GWAS, LD clumping, colocalisation, or locus definition
#   in this script.
#
# Scientific assumptions:
#   Loci, organs and traits are supplied. The count after each label
#   is the number of supplied association rows for that locus. Marker
#   and label angles are spread for legibility (order preserved,
#   isotonic least-squares shift); the leader's anchor stays at the
#   true genomic midpoint.

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
        locus = "locus",
        label = "band",
        organ = "organ",
        trait = "trait"
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
    chromosomes = "loci",        # "loci" (only chromosomes with loci) or "all"
    # Ring order = list order (outermost first). NULL = first appearance.
    organs = list(
        Abdominal = list(pch = 21, palette = "Quantitative.Sunset"),
        Cardiac = list(pch = 22, palette = "Quantitative.BrwnYl"),
        Brain = list(pch = 24, palette = "Gradient.YlGnBu")
    ),
    chrom_palette = "Qualitative.Bold",
    legend_suffix = " traits",
    show_count = TRUE,
    spacing = list(
        marker_within = 2.9,     # degrees between columns of one locus
        marker_between = 3.5,    # degrees between loci
        label = 2.4,             # degrees between labels
        margin = 2               # keep clear of the 12 o'clock seam
    ),
    radius = list(
        anchor = 0.925,
        label_start = 0.918,
        label_end = 0.790,
        hub = 0.770,
        marker_top = 0.745,
        ring_step = 0.058
    ),
    label_cex = 0.46,
    marker_cex = 1.45,
    chrom_cex = 1.0,
    legend_cex = 0.55,
    size = list(
        width = 10,
        height = 10
    )
)

load_packages(c("circlize", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
cc <- config$columns
for (k in c("chrom", "locus", "label", "organ", "trait")) {
    df[[cc[[k]]]] <- as.character(df[[cc[[k]]]])
}
df[[cc$start]] <- as.numeric(df[[cc$start]])
df[[cc$end]] <- as.numeric(df[[cc$end]])
if (any(!is.finite(df[[cc$start]]) | !is.finite(df[[cc$end]]))) {
    stop("start and end must be numeric.", call. = FALSE)
}
df <- unique(df)

cyto_path <- file.path(dirname(io$input), config$sidecars$cytoband)
if (!file.exists(cyto_path)) {
    stop(paste("Missing cytoband sidecar:", cyto_path), call. = FALSE)
}
cyto <- as.data.frame(read_table_auto(cyto_path), stringsAsFactors = FALSE)
require_columns(cyto, config$cyto)
cy <- config$cyto
band <- data.frame(
    chr = as.character(cyto[[cy$chrom]]),
    start = as.numeric(cyto[[cy$start]]),
    end = as.numeric(cyto[[cy$end]]),
    name = as.character(cyto[[cy$arm]]),
    gieStain = as.character(cyto[[cy$stain]]),
    stringsAsFactors = FALSE
)

chrom_order <- function(x) {
    x <- unique(as.character(x))
    key <- sub("^chr", "", x)
    rank <- match(key, c(as.character(seq_len(22)), "X", "Y", "M", "MT"))
    rank[is.na(rank)] <- 1000 + seq_len(sum(is.na(rank)))
    x[order(rank, x)]
}
keep <- if (identical(config$chromosomes, "loci")) {
    chrom_order(intersect(df[[cc$chrom]], band$chr))
} else {
    chrom_order(band$chr)
}
miss <- setdiff(unique(df[[cc$chrom]]), band$chr)
if (length(miss)) {
    stop(sprintf("Chromosome(s) not in cytoband: %s",
        paste(miss, collapse = ", ")), call. = FALSE)
}
band <- band[band$chr %in% keep, , drop = FALSE]
df <- df[df[[cc$chrom]] %in% keep, , drop = FALSE]
if (!nrow(df)) stop("No loci on the selected chromosomes.", call. = FALSE)
chr_len <- tapply(band$end, band$chr, max)[keep]

# organs, shapes and trait colours
organ_names <- if (is.null(config$organs)) unique(df[[cc$organ]]) else
    names(config$organs)
extra <- setdiff(unique(df[[cc$organ]]), organ_names)
if (length(extra)) {
    stop(sprintf("Organ(s) missing from CONFIG$organs: %s",
        paste(extra, collapse = ", ")), call. = FALSE)
}
organ_names <- organ_names[organ_names %in% df[[cc$organ]]]
default_pch <- c(21, 22, 24, 23, 25)
default_pal <- c("Quantitative.Sunset", "Quantitative.BrwnYl",
    "Quantitative.TealGrn", "Quantitative.Purp", "Quantitative.Burg")
organ_style <- lapply(seq_along(organ_names), function(i) {
    o <- organ_names[[i]]
    s <- if (!is.null(config$organs)) config$organs[[o]] else list()
    if (is.null(s$pch)) s$pch <- default_pch[[(i - 1) %% 5 + 1]]
    if (is.null(s$palette)) s$palette <- default_pal[[(i - 1) %% 5 + 1]]
    traits <- unique(df[[cc$trait]][df[[cc$organ]] == o])
    stops <- palette_colors(s$palette)
    # drop the palest stops so markers stay visible on white
    stops <- stops[-seq_len(max(1L, floor(length(stops) / 4)))]
    s$colors <- setNames(grDevices::colorRampPalette(stops)(length(traits)),
        traits)
    s$radius <- config$radius$marker_top - (i - 1) * config$radius$ring_step
    s
})
names(organ_style) <- organ_names

# one row per locus
loci <- df[!duplicated(df[[cc$locus]]),
    c(cc$locus, cc$chrom, cc$start, cc$end, cc$label), drop = FALSE]
names(loci) <- c("locus", "chr", "start", "end", "label")
loci$n <- as.integer(table(df[[cc$locus]])[loci$locus])
loci$mid <- (loci$start + loci$end) / 2
if (isTRUE(config$show_count)) {
    loci$label <- paste0(loci$label, " (", loci$n, ")")
}

# marker columns: column j of a locus holds the j-th trait of each organ
markers <- do.call(rbind, lapply(loci$locus, function(l) {
    sub <- df[df[[cc$locus]] == l, , drop = FALSE]
    per <- lapply(organ_names, function(o) unique(sub[[cc$trait]][sub[[cc$organ]] == o]))
    k <- max(lengths(per))
    out <- data.frame(locus = rep(l, k), col = seq_len(k), stringsAsFactors = FALSE)
    for (i in seq_along(organ_names)) {
        v <- per[[i]]
        out[[organ_names[[i]]]] <- c(v, rep(NA_character_, k - length(v)))
    }
    out
}))

spread <- function(target, step, margin) {
    off <- cumsum(c(0, step))
    span <- 360 - 2 * margin
    if (off[[length(off)]] > span * 0.98) {
        off <- off * span * 0.98 / off[[length(off)]]
    }
    fit <- stats::isoreg(seq_along(target), target - off)$yf
    fit <- pmax(margin, pmin(360 - margin - off[[length(off)]], fit))
    fit + off
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_plot <- function() {
    chr_cols <- setNames(palette_colors(config$chrom_palette, n = length(keep)),
        keep)
    circlize::circos.par(
        start.degree = 90, gap.degree = 1.2,
        cell.padding = c(0, 0, 0, 0),
        canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1),
        points.overflow.warning = FALSE
    )
    circlize::circos.initialize(keep, xlim = cbind(0, as.numeric(chr_len)))
    circlize::circos.track(
        ylim = c(0, 1), track.height = 0.007, track.margin = c(0.016, 0),
        bg.border = NA, bg.col = unname(chr_cols),
        panel.fun = function(x, y) {
            circlize::circos.text(
                circlize::CELL_META$xcenter,
                circlize::convert_y(4.5, "mm"),
                sub("^chr", "", circlize::CELL_META$sector.index),
                facing = "clockwise", niceFacing = TRUE,
                cex = config$chrom_cex, font = 2
            )
        }
    )
    circlize::circos.genomicIdeogram(band, track.height = 0.047,
        track.margin = c(0, 0))

    loci$angle <- vapply(seq_len(nrow(loci)), function(i) {
        th <- circlize::circlize(loci$mid[[i]], 0,
            sector.index = loci$chr[[i]], track.index = 1)[1, "theta"]
        (90 - th) %% 360
    }, numeric(1))
    circlize::circos.clear()

    loci <- loci[order(loci$angle), , drop = FALSE]
    mk <- markers[order(match(markers$locus, loci$locus), markers$col), ,
        drop = FALSE]
    sp <- config$spacing
    step <- ifelse(diff(match(mk$locus, loci$locus)) == 0,
        sp$marker_within, sp$marker_between)
    mk$angle <- spread(loci$angle[match(mk$locus, loci$locus)], step,
        sp$margin)
    loci$hub <- as.numeric(tapply(mk$angle, mk$locus, mean)[loci$locus])
    loci$label_angle <- spread(loci$angle, rep(sp$label, nrow(loci) - 1L),
        sp$margin)

    graphics::par(new = TRUE)
    circlize::circos.par(
        start.degree = 90, gap.degree = 0, cell.padding = c(0, 0, 0, 0),
        canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1),
        points.overflow.warning = FALSE
    )
    circlize::circos.initialize("display", xlim = c(0, 360))
    r <- config$radius
    r_low <- min(vapply(organ_style, `[[`, numeric(1), "radius")) - 0.03
    circlize::circos.track(
        ylim = c(r_low, 0.93), track.height = 0.93 - r_low,
        track.margin = c(0, 0.07), bg.border = NA
    )
    lc <- "grey20"
    circlize::circos.segments(loci$angle, r$anchor, loci$label_angle,
        r$label_start, straight = TRUE, lwd = 0.6, col = lc)
    circlize::circos.text(loci$label_angle, r$label_start - 0.004, loci$label,
        facing = "clockwise", niceFacing = TRUE, adj = c(1, 0.5),
        cex = config$label_cex, font = 2)
    text_len <- graphics::strwidth(loci$label, cex = config$label_cex,
        font = 2)
    tail_r <- pmax(r$hub + 0.02, r$label_start - 0.004 - text_len - 0.008)
    circlize::circos.segments(loci$label_angle, tail_r,
        loci$label_angle, tail_r - 0.012, straight = TRUE,
        lwd = 0.6, col = lc)
    circlize::circos.segments(loci$label_angle, tail_r - 0.012,
        loci$hub, r$hub, straight = TRUE, lwd = 0.6, col = lc)
    circlize::circos.segments(loci$hub[match(mk$locus, loci$locus)], r$hub,
        mk$angle, r$marker_top + 0.012, straight = TRUE, lwd = 0.6, col = lc)

    for (o in organ_names) {
        s <- organ_style[[o]]
        tr <- mk[[o]]
        ok <- !is.na(tr)
        if (!any(ok)) next
        circlize::circos.points(mk$angle[ok], rep(s$radius, sum(ok)),
            pch = s$pch, cex = config$marker_cex, col = NA,
            bg = unname(s$colors[tr[ok]]))
    }
    circlize::circos.clear()

    # centred legend columns sized from measured text width
    lg_text <- lapply(organ_names, function(o) c(
        paste0(o, config$legend_suffix), names(organ_style[[o]]$colors)))
    col_w <- vapply(lg_text, function(t) {
        max(graphics::strwidth(t, cex = config$legend_cex)) +
            graphics::strwidth("MMM", cex = config$legend_cex)
    }, numeric(1))
    lx <- -sum(col_w) / 2 + c(0, cumsum(col_w))[seq_along(col_w)]
    n_rows <- max(lengths(lg_text))
    row_h <- graphics::strheight("Xg", cex = config$legend_cex) * 1.35 * 1.5
    ly <- n_rows * row_h / 2
    for (i in seq_along(organ_names)) {
        s <- organ_style[[i]]
        graphics::legend(
            lx[[i]], ly,
            legend = lg_text[[i]],
            pch = s$pch, pt.bg = c("black", unname(s$colors)), col = NA,
            bty = "n", cex = config$legend_cex, pt.cex = 1.2,
            y.intersp = 1.35, x.intersp = 0.75
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
w <- max(config$size$width, config$size$height)
if (identical(fmt, "png")) {
    grDevices::png(out_path, width = w, height = w, units = "in", res = 300,
        bg = "white")
} else if (identical(fmt, "pdf")) {
    grDevices::pdf(out_path, width = w, height = w)
} else {
    grDevices::svg(out_path, width = w, height = w, bg = "white")
}
graphics::par(mar = c(0.5, 0.5, 0.5, 0.5), bg = "white", family = "sans")
tryCatch(
    draw_plot(),
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
