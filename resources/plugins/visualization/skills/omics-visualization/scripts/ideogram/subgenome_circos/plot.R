#!/usr/bin/env Rscript

# Template-ID: ideogram-subgenome-circos
#
# Purpose:
#   Multi-species subgenome Circos (polyploid comparative genomics):
#   the reference species' subgenomes are sectors with an inner
#   chromosome ring, stacked heatmap tracks and central homeolog
#   links; outside, each sector carries a radial stack of species
#   rows (one arc per species x subgenome) joined by supplied
#   synteny blocks between adjacent rows ("circular riparian").
#   Subgenomes absent from the reference (e.g. H, D) are rows
#   anchored on a reference sector; a dashed arc marks missing
#   positions at a chosen level.
#
# Inputs:
#   example.tsv - karyotype: species, subgenome, chr, length
#                 (chr unique within a species; order = draw order)
#   Sidecars beside the input (filenames in CONFIG):
#     rows.tsv   - species, subgenome, anchor, level, group, ploidy
#                  anchor = reference subgenome sector it is drawn on;
#                  level = 1 (innermost) ...; group / ploidy feed the
#                  ploidy legend (empty group = not in legend)
#     blocks.tsv - species1, chr1, start1, end1, species2, chr2,
#                  start2, end2 (start2 > end2 = inverted block)
#     tracks.tsv - chr, start, end, track, value on reference chr
#     links.tsv  - chr1, start1, end1, chr2, start2, end2 on
#                  reference chr (central homeolog ribbons)
#
# Output:
#   A PDF, PNG, or SVG multi-subgenome synteny Circos.
#
# Dependencies:
#   circlize, ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (reference species, sector order and gaps,
#   subgenome colours, heatmap tracks and palettes, stack height,
#   dashed level, which subgenomes get outer chr labels). Extra
#   species are extra rows.tsv / karyotype rows; extra heatmap
#   tracks are extra CONFIG entries. Set a sidecar to NULL to drop
#   it. Do not run MCScanX, WGDI, JCVI, BLAST, or compute
#   expression / methylation windows in this script.
#
# Scientific assumptions:
#   Synteny blocks, homeolog links, chromosome lengths and window
#   values are supplied. Rows are scaled to their anchor sector, so
#   arc length is relative within a row, not absolute bp across
#   species. Heatmap colour spans each track's supplied range.

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
    reference = "PWB sp1",
    sectors = c("C", "B", "A"),     # reference subgenomes, clockwise from top
    gap_after = c(C = 6, B = 10, A = 26),  # last gap = species-label gap at top
    start_degree = 87,
    sidecars = list(
        rows = "rows.tsv",
        blocks = "blocks.tsv",
        tracks = "tracks.tsv",
        links = "links.tsv"
    ),
    subgenome_palette = "Qualitative.Bold",
    subgenome_index = c(H = 1, C = 3, B = 7, A = 6, D = 2),
    heat = list(
        list(track = "Expression", palette = "Gradient.Blues"),
        list(track = "mCG", palette = "Diverging.Spectral", reverse = TRUE),
        list(track = "Gene", palette = "Gradient.YlOrRd"),
        list(track = "Synteny", palette = "Gradient.Greys")
    ),
    chr_gap = 0.012,                # gap between chromosomes, fraction of sector
    stack_height = 0.40,
    row_line = 0.10,                # row line thickness, level units
    block_alpha = 1,
    link_alpha = 0.55,
    dashed_level = 9,               # NULL = no dashed missing arc
    outer_label_subgenomes = c("H", "D"),
    species_italic = TRUE,
    cex = list(species = 0.55, chr = 0.5, ref_chr = 0.55),
    legend = list(x = 0.02, y = 0.98),   # npc, top-left anchor
    size = list(width = 10, height = 10)
)

load_packages(c("circlize", "ComplexHeatmap", "grid", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
read_sidecar <- function(key, cols) {
    fn <- config$sidecars[[key]]
    if (is.null(fn)) return(NULL)
    path <- file.path(dirname(io$input), fn)
    if (!file.exists(path)) {
        stop(paste("Missing sidecar:", path), call. = FALSE)
    }
    d <- as.data.frame(read_table_auto(path), stringsAsFactors = FALSE)
    require_columns(d, as.list(stats::setNames(cols, cols)))
    d
}

kar <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(kar, list(species = "species", subgenome = "subgenome",
    chr = "chr", length = "length"))
kar$species <- as.character(kar$species)
kar$subgenome <- as.character(kar$subgenome)
kar$chr <- as.character(kar$chr)
kar$length <- as.numeric(kar$length)
if (any(!is.finite(kar$length) | kar$length <= 0)) {
    stop("karyotype length must be positive.", call. = FALSE)
}
kar$key <- paste(kar$species, kar$chr, sep = "\r")
if (anyDuplicated(kar$key)) {
    stop("chr must be unique within a species.", call. = FALSE)
}

rows <- read_sidecar("rows", c("species", "subgenome", "anchor", "level",
    "group", "ploidy"))
for (nm in c("species", "subgenome", "anchor", "group", "ploidy")) {
    rows[[nm]] <- as.character(rows[[nm]])
    rows[[nm]][is.na(rows[[nm]])] <- ""
}
rows$level <- as.integer(rows$level)
if (!all(rows$anchor %in% config$sectors)) {
    stop("rows.tsv anchor must be one of CONFIG sectors.", call. = FALSE)
}
if (anyDuplicated(paste(rows$anchor, rows$level))) {
    stop("rows.tsv level must be unique within an anchor.", call. = FALSE)
}
rows$rkey <- paste(rows$species, rows$subgenome, sep = "\r")
if (anyDuplicated(rows$rkey)) {
    stop("one row per species x subgenome.", call. = FALSE)
}
miss <- setdiff(rows$rkey, paste(kar$species, kar$subgenome, sep = "\r"))
if (length(miss)) {
    stop(sprintf("rows without karyotype: %s",
        gsub("\r", " / ", miss[[1]])), call. = FALSE)
}

ref_kar <- kar[kar$species == config$reference, , drop = FALSE]
if (!all(config$sectors %in% ref_kar$subgenome)) {
    stop("reference lacks a CONFIG sector subgenome.", call. = FALSE)
}

sg_all <- unique(c(config$sectors, rows$subgenome))
pal <- palette_colors(config$subgenome_palette)
sg_col <- stats::setNames(pal[config$subgenome_index[sg_all]], sg_all)
if (anyNA(sg_col)) {
    sg_col[is.na(sg_col)] <- palette_colors(config$subgenome_palette,
        n = length(sg_all))[is.na(sg_col)]
}

# Sector length = reference subgenome bp + chromosome gaps
sector_len <- vapply(config$sectors, function(s) {
    L <- ref_kar$length[ref_kar$subgenome == s]
    sum(L) / (1 - config$chr_gap * (length(L) - 1))
}, numeric(1))

# Per (species, subgenome) layout: chr offset and scale on its anchor sector
layout_row <- function(species, subgenome, anchor) {
    k <- kar[kar$species == species & kar$subgenome == subgenome, ]
    slen <- sector_len[[anchor]]
    g <- config$chr_gap * slen
    sc <- (slen - g * (nrow(k) - 1)) / sum(k$length)
    off <- c(0, cumsum(k$length * sc + g))[seq_len(nrow(k))]
    data.frame(chr = k$chr, start = off, end = off + k$length * sc,
        scale = sc, stringsAsFactors = FALSE)
}
row_lay <- lapply(seq_len(nrow(rows)), function(i) {
    layout_row(rows$species[[i]], rows$subgenome[[i]], rows$anchor[[i]])
})
names(row_lay) <- rows$rkey
ref_lay <- lapply(config$sectors, function(s) {
    layout_row(config$reference, s, s)
})
names(ref_lay) <- config$sectors
ref_sector <- stats::setNames(ref_kar$subgenome, ref_kar$chr)

ref_x <- function(chr, pos) {
    s <- unname(ref_sector[chr])
    out <- numeric(length(chr))
    for (i in seq_along(chr)) {
        lay <- ref_lay[[s[[i]]]]
        j <- match(chr[[i]], lay$chr)
        out[[i]] <- lay$start[[j]] + pos[[i]] * lay$scale[[j]]
    }
    out
}

check_ref_chr <- function(chr, what) {
    bad <- setdiff(unique(chr), names(ref_sector))
    if (length(bad)) {
        stop(sprintf("%s chr not on reference: %s", what, bad[[1]]),
            call. = FALSE)
    }
}

sg_of <- stats::setNames(kar$subgenome, kar$key)
blocks <- read_sidecar("blocks", c("species1", "chr1", "start1", "end1",
    "species2", "chr2", "start2", "end2"))
blk <- NULL
if (!is.null(blocks) && nrow(blocks)) {
    k1 <- paste(blocks$species1, blocks$chr1, sep = "\r")
    k2 <- paste(blocks$species2, blocks$chr2, sep = "\r")
    if (!all(c(k1, k2) %in% kar$key)) {
        stop("blocks.tsv species/chr not in karyotype.", call. = FALSE)
    }
    r1 <- match(paste(blocks$species1, sg_of[k1], sep = "\r"), rows$rkey)
    r2 <- match(paste(blocks$species2, sg_of[k2], sep = "\r"), rows$rkey)
    keep <- !is.na(r1) & !is.na(r2) & rows$anchor[r1] == rows$anchor[r2]
    if (any(!keep)) {
        stop("each block must join two drawn rows on the same anchor.",
            call. = FALSE)
    }
    to_x <- function(r, chr, pos) {
        vapply(seq_along(r), function(i) {
            lay <- row_lay[[r[[i]]]]
            j <- match(chr[[i]], lay$chr)
            lay$start[[j]] + as.numeric(pos[[i]]) * lay$scale[[j]]
        }, numeric(1))
    }
    blk <- data.frame(
        anchor = rows$anchor[r1],
        y1 = rows$level[r1], y2 = rows$level[r2],
        x1a = to_x(r1, blocks$chr1, blocks$start1),
        x1b = to_x(r1, blocks$chr1, blocks$end1),
        x2a = to_x(r2, blocks$chr2, blocks$start2),
        x2b = to_x(r2, blocks$chr2, blocks$end2),
        col = unname(sg_col[rows$subgenome[r1]]),
        stringsAsFactors = FALSE
    )
}

tracks <- read_sidecar("tracks", c("chr", "start", "end", "track", "value"))
heat <- list()
if (!is.null(tracks)) {
    tracks$chr <- as.character(tracks$chr)
    check_ref_chr(tracks$chr, "tracks.tsv")
    for (h in config$heat) {
        d <- tracks[tracks$track == h$track, , drop = FALSE]
        if (!nrow(d)) next
        v <- as.numeric(d$value)
        stops <- palette_colors(h$palette)
        if (isTRUE(h$reverse)) stops <- rev(stops)
        lims <- range(v, na.rm = TRUE)
        heat[[length(heat) + 1L]] <- list(
            title = h$track,
            sector = unname(ref_sector[d$chr]),
            xs = ref_x(d$chr, as.numeric(d$start)),
            xe = ref_x(d$chr, as.numeric(d$end)),
            value = v,
            col_fun = circlize::colorRamp2(
                seq(lims[[1]], lims[[2]], length.out = length(stops)), stops),
            lims = lims
        )
    }
}

links <- read_sidecar("links", c("chr1", "start1", "end1", "chr2",
    "start2", "end2"))
if (!is.null(links) && nrow(links)) {
    links$chr1 <- as.character(links$chr1)
    links$chr2 <- as.character(links$chr2)
    check_ref_chr(c(links$chr1, links$chr2), "links.tsv")
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
ribbon_poly <- function(x1a, x1b, y1, x2a, x2b, y2, n = 24) {
    t <- seq(0, 1, length.out = n)
    list(
        x = c(x1a + (x1b - x1a) * t, x1b + (x2b - x1b) * t,
            x2b + (x2a - x2b) * t, x2a + (x1a - x2a) * t),
        y = c(rep(y1, n), y1 + (y2 - y1) * t, rep(y2, n), y2 + (y1 - y2) * t)
    )
}

draw_group_legend <- function() {
    grp <- rows[nzchar(rows$group), , drop = FALSE]
    if (!nrow(grp)) return(invisible())
    lv <- unique(grp$group)
    cell <- grid::unit(4.2, "mm")
    head_h <- grid::unit(3.6, "mm")
    y <- grid::unit(1, "npc")
    w_mm <- max(vapply(lv, function(g) {
        pl <- grp$ploidy[grp$group == g][[1]]
        tw <- grid::convertWidth(grid::stringWidth(paste0(g, "  ", pl)),
            "mm", valueOnly = TRUE) * 7 / 12
        max(tw + 3, 4.2 * length(unique(grp$subgenome[grp$group == g])))
    }, numeric(1)))
    w <- grid::unit(w_mm, "mm")
    for (g in lv) {
        sgs <- sg_all[sg_all %in% grp$subgenome[grp$group == g]]
        pl <- grp$ploidy[grp$group == g][[1]]
        grid::grid.rect(x = 0, y = y, width = w, height = head_h,
            just = c("left", "top"),
            gp = grid::gpar(fill = "white", col = "grey30", lwd = 0.6))
        grid::grid.text(paste0(g, "  ", pl), x = grid::unit(1, "mm"),
            y = y - head_h * 0.5, just = "left",
            gp = grid::gpar(fontsize = 7))
        cw <- w * (1 / length(sgs))
        for (i in seq_along(sgs)) {
            grid::grid.rect(x = cw * (i - 1), y = y - head_h,
                width = cw, height = cell, just = c("left", "top"),
                gp = grid::gpar(fill = sg_col[[sgs[[i]]]], col = "grey30",
                    lwd = 0.6))
            grid::grid.text(sgs[[i]], x = cw * (i - 0.5),
                y = y - head_h - cell * 0.5,
                gp = grid::gpar(fontsize = 7, col = "white",
                    fontface = "bold"))
        }
        y <- y - head_h - cell - grid::unit(2.5, "mm")
    }
}

draw_figure <- function() {
    circlize::circos.par(
        start.degree = config$start_degree,
        gap.after = unname(config$gap_after[config$sectors]),
        cell.padding = c(0, 0, 0, 0),
        track.margin = c(0, 0),
        circle.margin = 0.04,
        points.overflow.warning = FALSE
    )
    circlize::circos.initialize(
        sectors = factor(config$sectors, levels = config$sectors),
        xlim = cbind(0, unname(sector_len))
    )

    max_level <- max(rows$level)
    stack_ylim <- c(-0.4, max_level + 1.6)
    circlize::circos.track(
        ylim = stack_ylim, track.height = config$stack_height,
        bg.border = NA,
        panel.fun = function(x, y) {
            s <- circlize::CELL_META$sector.index
            if (!is.null(blk)) {
                b <- blk[blk$anchor == s, , drop = FALSE]
                for (i in seq_len(nrow(b))) {
                    pp <- ribbon_poly(b$x1a[[i]], b$x1b[[i]], b$y1[[i]],
                        b$x2a[[i]], b$x2b[[i]], b$y2[[i]])
                    circlize::circos.polygon(pp$x, pp$y,
                        col = grDevices::adjustcolor(b$col[[i]],
                            alpha.f = config$block_alpha),
                        border = "white", lwd = 0.35)
                }
            }
            rs <- rows[rows$anchor == s, , drop = FALSE]
            for (i in seq_len(nrow(rs))) {
                lay <- row_lay[[rs$rkey[[i]]]]
                lv <- rs$level[[i]]
                circlize::circos.rect(lay$start, lv - config$row_line,
                    lay$end, lv + config$row_line,
                    col = "grey10", border = NA)
            }
            for (sg in config$outer_label_subgenomes) {
                rr <- rs[rs$subgenome == sg, , drop = FALSE]
                if (!nrow(rr)) next
                top <- rr[which.max(rr$level), ]
                lay <- row_lay[[top$rkey]]
                circlize::circos.text((lay$start + lay$end) / 2,
                    top$level + 0.9, lay$chr,
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = config$cex$chr)
            }
        }
    )

    circlize::circos.track(
        ylim = c(0, 1), track.height = 0.045, bg.border = NA,
        panel.fun = function(x, y) {
            lay <- ref_lay[[circlize::CELL_META$sector.index]]
            circlize::circos.text((lay$start + lay$end) / 2, 0.4, lay$chr,
                facing = "clockwise", niceFacing = TRUE,
                cex = config$cex$ref_chr)
        }
    )
    circlize::circos.track(
        ylim = c(0, 1), track.height = 0.022, bg.border = NA,
        panel.fun = function(x, y) {
            s <- circlize::CELL_META$sector.index
            lay <- ref_lay[[s]]
            circlize::circos.rect(lay$start, 0, lay$end, 1,
                col = sg_col[[s]], border = NA)
        }
    )
    for (h in heat) {
        circlize::circos.track(
            ylim = c(0, 1), track.height = 0.032, bg.border = NA,
            track.margin = c(0.003, 0.003),
            panel.fun = function(x, y) {
                s <- circlize::CELL_META$sector.index
                i <- which(h$sector == s)
                if (!length(i)) return()
                circlize::circos.rect(h$xs[i], 0, h$xe[i], 1,
                    col = h$col_fun(h$value[i]), border = NA)
            }
        )
    }

    if (!is.null(links) && nrow(links)) {
        s1 <- unname(ref_sector[links$chr1])
        s2 <- unname(ref_sector[links$chr2])
        a1 <- ref_x(links$chr1, as.numeric(links$start1))
        b1 <- ref_x(links$chr1, as.numeric(links$end1))
        a2 <- ref_x(links$chr2, as.numeric(links$start2))
        b2 <- ref_x(links$chr2, as.numeric(links$end2))
        for (i in seq_len(nrow(links))) {
            circlize::circos.link(s1[[i]], c(a1[[i]], b1[[i]]),
                s2[[i]], c(a2[[i]], b2[[i]]),
                col = grDevices::adjustcolor(sg_col[[s1[[i]]]],
                    alpha.f = config$link_alpha),
                border = NA)
        }
    }

    # species labels at the start of the first sector
    first <- config$sectors[[1]]
    rs <- rows[rows$anchor == first, , drop = FALSE]
    for (i in seq_len(nrow(rs))) {
        pt <- circlize::circlize(0, rs$level[[i]], sector.index = first,
            track.index = 1)
        th <- pt[1, "theta"] / 180 * pi
        graphics::text(pt[1, "rou"] * cos(th) - 0.012,
            pt[1, "rou"] * sin(th), rs$species[[i]], adj = c(1, 0.5),
            cex = config$cex$species,
            font = if (isTRUE(config$species_italic)) 3 else 1)
    }

    # dashed arc at a level wherever no row is drawn (not the label gap)
    if (!is.null(config$dashed_level)) {
        lv <- config$dashed_level
        n <- length(config$sectors)
        deg <- t(vapply(config$sectors, function(s) {
            c(circlize::get.cell.meta.data("cell.start.degree", s, 1),
                circlize::get.cell.meta.data("cell.end.degree", s, 1))
        }, numeric(2)))
        rou <- circlize::circlize(0, lv, sector.index = first,
            track.index = 1)[1, "rou"]
        arc <- function(d0, d1) {
            while (d1 > d0) d1 <- d1 - 360
            th <- seq(d0, d1, length.out = 120) / 180 * pi
            graphics::lines(rou * cos(th), rou * sin(th), lty = 2,
                lwd = 0.9, col = "grey10")
        }
        for (i in seq_len(n)) {
            s <- config$sectors[[i]]
            if (!any(rows$anchor == s & rows$level == lv)) {
                arc(deg[i, 1], deg[i, 2])
            }
            if (i < n) arc(deg[i, 2], deg[i + 1L, 1])
        }
    }
    circlize::circos.clear()

    lg <- config$legend
    grid::pushViewport(grid::viewport(x = lg$x, y = lg$y,
        width = 0.2, height = 0.3, just = c("left", "top")))
    draw_group_legend()
    grid::popViewport()

    if (length(heat)) {
        hl <- lapply(heat, function(h) {
            ComplexHeatmap::Legend(col_fun = h$col_fun, at = h$lims,
                labels = c("Low", "High"), title = h$title,
                title_position = "topleft", direction = "horizontal",
                legend_width = grid::unit(16, "mm"),
                grid_height = grid::unit(2.5, "mm"),
                tick_length = grid::unit(0, "mm"),
                border = "grey40",
                title_gp = grid::gpar(fontsize = 7),
                labels_gp = grid::gpar(fontsize = 6))
        })
        pk <- do.call(ComplexHeatmap::packLegend,
            c(hl, list(direction = "vertical", gap = grid::unit(1, "mm"))))
        ComplexHeatmap::draw(pk, x = grid::unit(lg$x, "npc") +
            grid::unit(26, "mm"), y = grid::unit(lg$y, "npc"),
            just = c("left", "top"))
    }
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
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
    grDevices::png(out_path, width = width, height = height, units = "in",
        res = 300, bg = "transparent")
} else if (identical(fmt, "pdf")) {
    grDevices::pdf(out_path, width = width, height = height)
} else {
    grDevices::svg(out_path, width = width, height = height,
        bg = "transparent")
}
graphics::par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "transparent")
tryCatch(draw_figure(), finally = grDevices::dev.off())
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
