#!/usr/bin/env Rscript

# Template-ID: heatmap-circos-split
#
# Purpose:
#   Draw a circular heatmap split into sectors by a supplied group,
#   with an optional q-value diamond track and an optional Euler
#   diagram of gene-set overlap in the hole.
#
# Inputs:
#   One row per gene within a group. Default example:
#     - Gene: row label (may repeat across groups)
#     - Group: sector split
#     - qval: supplied adjusted p (diamond fill)
#     - remaining numeric columns: heatmap values (sites / samples)
#
# Output:
#   A PDF, PNG, or SVG split circular heatmap.
#
# Dependencies:
#   circlize, ComplexHeatmap, grid, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, inner, show_qval,
#   start.degree, gaps, palette, value range, borders, row names,
#   sector / column label style, legend position). inner = "euler"
#   vs "none", the q-value track, outside sector labels, repelled
#   column labels and a centre legend are CONFIG, not extra ids.
#   layout = "locus" (example_locus.tsv) adds outer gene labels
#   coloured by a class column, a band ring, a merged locus-number
#   ring and several binary / continuous track blocks with names in
#   the top gap; it is a variant, not a new id.
#   Number of groups or sites is the table, not a new id. A rectangular heatmap with row_split is
#   not a second template. Do not add an id per gene set (Acceleration
#   / Curl / Divergence).
#
# Scientific assumptions:
#   Heatmap cells and qval are supplied. This script does not cluster,
#   rescale, or run a test. Within-group row order is a display sort
#   (qval or file order). The Euler diagram counts gene-name overlap
#   across Group; it is not a statistical test. Euler inset supports
#   2 or 3 groups.

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
        gene = "Gene",
        group = "Group",
        qval = "qval"
    ),
    inner = "euler",
    show_qval = TRUE,            # FALSE also makes the qval column optional
    sort = "qval",
    start.degree = 60,
    gap.within = 5,
    gap.last = 30,
    heat_palette = "Quantitative.BluGrn",
    heat_ramp = "ends",          # "ends" (two-colour) or "full" (every stop)
    heat_reverse = FALSE,
    value_range = NULL,          # NULL = data range
    cell_border = "white",
    sector_border = "black",
    track_height = 0.28,
    show_rownames = TRUE,
    sector_label = "band",       # "band" (inner tinted band) or "outside"
    column_label = "ticks",      # "ticks" or "repel" (beside each ring, de-overlapped)
    legend_position = "corner",  # "corner" or "center"
    legend_labels = NULL,        # e.g. c("Low", "High") at the range ends
    # "split" = one heatmap block (default); "locus" = outer coloured
    # gene labels + band + locus-number ring + several track blocks
    # (binary / continuous), track names in the top gap, sector names
    # inside. Colours are "Family.Name:index" catalog references.
    layout = "split",
    locus = list(
        band = "Band",
        locus = "Locus",
        label_class = "Class",
        band_colors = c("All-cause HF" = "Qualitative.Safe:12",
            DCM = "Qualitative.Vivid:1", HFrEF = "Qualitative.Bold:2",
            HFpEF = "Qualitative.Bold:8"),
        label_colors = c(Novel = "Diverging.RdBu:1"),  # other classes black
        blocks = list(
            list(title = "Locus Classifier",
                columns = c("Coloc", "Mendel", "ABC"),
                type = "binary", palette = "Gradient.Blues:7"),
            list(title = "Predictor",
                columns = c("Overall", "V2G", "TWAS", "PoPS"),
                type = "continuous", palette = "Gradient.Blues")
        ),
        ring_height = 0.034,
        label_track = 0.13,
        label_cex = 0.52
    ),
    labels = list(
        fill = "Value",
        qval = "qval",
        inner = "Morphogenic\ngenes"
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
df <- as.data.frame(read_table_auto(io$input), check.names = FALSE,
    stringsAsFactors = FALSE)
g_gene <- config$columns$gene
g_grp <- config$columns$group
g_q <- config$columns$qval
has_q <- g_q %in% names(df)
if (isTRUE(config$show_qval) || has_q) {
    require_columns(df, config$columns)
} else {
    require_columns(df, config$columns[c("gene", "group")])
}

df[[g_gene]] <- as.character(df[[g_gene]])
df[[g_grp]] <- as.character(df[[g_grp]])
if (any(!nzchar(df[[g_gene]]) | !nzchar(df[[g_grp]]))) {
    stop("Gene and Group must be non-empty.", call. = FALSE)
}
if (has_q) {
    df[[g_q]] <- as.numeric(df[[g_q]])
    if (any(!is.finite(df[[g_q]]) | df[[g_q]] < 0 | df[[g_q]] > 1)) {
        stop("qval must be numeric in [0, 1].", call. = FALSE)
    }
}

is_locus <- identical(config$layout, "locus")
meta_cols <- if (is_locus) {
    unlist(config$locus[c("band", "locus", "label_class")], use.names = FALSE)
} else {
    character()
}
if (is_locus) {
    require_columns(df, as.list(stats::setNames(meta_cols, meta_cols)))
    blk_cols <- unlist(lapply(config$locus$blocks, `[[`, "columns"))
    require_columns(df, as.list(stats::setNames(blk_cols, blk_cols)))
}
value_cols <- setdiff(names(df),
    c(unlist(config$columns, use.names = FALSE), meta_cols))
if (length(value_cols) < 2L) {
    stop("Need at least two numeric site / sample columns.", call. = FALSE)
}
for (nm in value_cols) {
    df[[nm]] <- as.numeric(df[[nm]])
    if (any(!is.finite(df[[nm]]))) {
        stop(paste("Non-finite value in column", nm), call. = FALSE)
    }
}

grp_levels <- unique(df[[g_grp]])
df[[g_grp]] <- factor(df[[g_grp]], levels = grp_levels)
if (identical(config$sort, "qval") && has_q) {
    df <- df[order(df[[g_grp]], df[[g_q]], df[[g_gene]]), , drop = FALSE]
} else {
    df <- df[order(df[[g_grp]]), , drop = FALSE]
}

mat <- as.matrix(df[, value_cols, drop = FALSE])
storage.mode(mat) <- "double"
rownames(mat) <- df[[g_gene]]

heat_pal <- palette_colors(config$heat_palette)
if (isTRUE(config$heat_reverse)) heat_pal <- rev(heat_pal)
if (!identical(config$heat_ramp, "full")) {
    heat_pal <- heat_pal[c(1L, length(heat_pal))]
}
q_pal <- palette_colors("Quantitative.RedOr")
q_light <- q_pal[[1]]
q_dark <- q_pal[[length(q_pal)]]

v_rng <- if (is.null(config$value_range)) range(mat, na.rm = TRUE) else
    as.numeric(config$value_range)
if (v_rng[[1]] >= v_rng[[2]]) {
    v_rng[[2]] <- v_rng[[1]] + 1e-6
}
heat_col <- circlize::colorRamp2(
    seq(v_rng[[1]], v_rng[[2]], length.out = length(heat_pal)),
    heat_pal
)
# Small qval is darker (more significant).
q_col <- circlize::colorRamp2(c(0, 1), c(q_dark, q_light))

n_grp <- length(grp_levels)
group_cols <- palette_colors("Qualitative.Safe", n = n_grp)
names(group_cols) <- grp_levels

draw_euler <- function(sets, fills, title) {
    n <- length(sets)
    if (n < 2L || n > 3L) {
        grid::grid.text(
            "Euler inset needs 2 or 3 groups",
            gp = grid::gpar(fontsize = 8, col = "#333333")
        )
        return(invisible())
    }
    fills <- unname(fills[names(sets)])
    fills <- grDevices::adjustcolor(fills, alpha.f = 0.35)
    A <- unique(sets[[1]])
    B <- unique(sets[[2]])
    if (n == 2L) {
        circ <- list(
            list(x = 0.38, y = 0.52, r = 0.30),
            list(x = 0.62, y = 0.52, r = 0.30)
        )
        only1 <- setdiff(A, B)
        only2 <- setdiff(B, A)
        both <- intersect(A, B)
        qty <- list(
            list(x = 0.28, y = 0.52, n = length(only1)),
            list(x = 0.72, y = 0.52, n = length(only2)),
            list(x = 0.50, y = 0.52, n = length(both))
        )
    } else {
        C <- unique(sets[[3]])
        circ <- list(
            list(x = 0.38, y = 0.58, r = 0.28),
            list(x = 0.62, y = 0.58, r = 0.28),
            list(x = 0.50, y = 0.36, r = 0.28)
        )
        abc <- Reduce(intersect, list(A, B, C))
        ab <- setdiff(intersect(A, B), C)
        ac <- setdiff(intersect(A, C), B)
        bc <- setdiff(intersect(B, C), A)
        only1 <- setdiff(setdiff(A, B), C)
        only2 <- setdiff(setdiff(B, A), C)
        only3 <- setdiff(setdiff(C, A), B)
        qty <- list(
            list(x = 0.30, y = 0.64, n = length(only1)),
            list(x = 0.70, y = 0.64, n = length(only2)),
            list(x = 0.50, y = 0.26, n = length(only3)),
            list(x = 0.50, y = 0.64, n = length(ab)),
            list(x = 0.36, y = 0.42, n = length(ac)),
            list(x = 0.64, y = 0.42, n = length(bc)),
            list(x = 0.50, y = 0.48, n = length(abc))
        )
    }
    for (i in seq_along(circ)) {
        grid::grid.circle(
            x = circ[[i]]$x,
            y = circ[[i]]$y,
            r = circ[[i]]$r,
            default.units = "npc",
            gp = grid::gpar(
                fill = fills[[i]],
                col = "black",
                lty = 2,
                lwd = 1.2
            )
        )
    }
    for (q in qty) {
        if (q$n <= 0) {
            next
        }
        grid::grid.text(
            as.character(q$n),
            x = q$x,
            y = q$y,
            default.units = "npc",
            gp = grid::gpar(fontsize = 8, col = "#222222")
        )
    }
    if (nzchar(title)) {
        grid::grid.text(
            title,
            x = 0.5,
            y = -0.04,
            default.units = "npc",
            gp = grid::gpar(fontsize = 9, lineheight = 0.9, col = "#222222")
        )
    }
}

pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

draw_locus <- function() {
    lc <- config$locus
    bands <- as.character(df[[lc$band]])
    loci <- as.character(df[[lc$locus]])
    cls <- as.character(df[[lc$label_class]])
    band_lv <- unique(bands)
    band_col <- vapply(band_lv, function(b) {
        if (b %in% names(lc$band_colors)) pal_ref(lc$band_colors[[b]]) else
            NA_character_
    }, character(1))
    if (anyNA(band_col)) {
        band_col[is.na(band_col)] <- palette_colors("Qualitative.Safe",
            n = length(band_lv))[is.na(band_col)]
    }
    lab_col <- ifelse(cls %in% names(lc$label_colors),
        vapply(cls, function(k) {
            if (k %in% names(lc$label_colors)) pal_ref(lc$label_colors[[k]])
            else "black"
        }, character(1)), "black")
    blocks <- lapply(lc$blocks, function(b) {
        m <- as.matrix(df[, b$columns, drop = FALSE])
        if (identical(b$type, "binary")) {
            if (any(!m %in% c(0, 1))) {
                stop(sprintf("Block '%s' must be 0/1.", b$title),
                    call. = FALSE)
            }
            b$fill <- function(v) ifelse(v == 1, pal_ref(b$palette), "white")
        } else {
            stops <- pal_ref(b$palette)
            rng <- range(m)
            if (rng[[1]] >= rng[[2]]) rng[[2]] <- rng[[1]] + 1e-6
            b$col_fun <- circlize::colorRamp2(
                seq(rng[[1]], rng[[2]], length.out = length(stops)), stops)
            b$fill <- b$col_fun
            b$rng <- rng
        }
        b$mat <- m
        b
    })
    grp <- df[[g_grp]]
    idx_in <- split(seq_len(nrow(df)), grp)
    xlim <- cbind(0, vapply(grp_levels, function(g) length(idx_in[[g]]),
        numeric(1)))

    circlize::circos.par(
        start.degree = 90 - config$gap.last / 2,
        gap.after = c(rep(config$gap.within, n_grp - 1L), config$gap.last),
        track.margin = c(0, 0),
        cell.padding = c(0, 0, 0, 0),
        circle.margin = 0.02,
        points.overflow.warning = FALSE
    )
    circlize::circos.initialize(sectors = factor(grp_levels,
        levels = grp_levels), xlim = xlim)
    rows_of <- function() idx_in[[circlize::CELL_META$sector.index]]

    circlize::circos.track(ylim = c(0, 1), track.height = lc$label_track,
        bg.border = NA, panel.fun = function(x, y) {
            i <- rows_of()
            circlize::circos.text(seq_along(i) - 0.5, rep(0.02, length(i)),
                df[[g_gene]][i], facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = lc$label_cex, col = lab_col[i])
        })
    circlize::circos.track(ylim = c(0, 1), track.height = 0.035,
        bg.border = NA, panel.fun = function(x, y) {
            i <- rows_of()
            circlize::circos.rect(seq_along(i) - 1, 0, seq_along(i), 1,
                col = band_col[bands[i]], border = "grey25", lwd = 0.4)
        })
    circlize::circos.track(ylim = c(0, 1), track.height = 0.04,
        bg.border = NA, track.margin = c(0.004, 0.004),
        panel.fun = function(x, y) {
            i <- rows_of()
            r <- rle(loci[i])
            end <- cumsum(r$lengths)
            start <- end - r$lengths
            circlize::circos.text((start + end) / 2, rep(0.5, length(end)),
                r$values, facing = "bending.inside", niceFacing = TRUE,
                cex = lc$label_cex * 0.9)
        })
    track_y <- list()
    for (bi in seq_along(blocks)) {
        b <- blocks[[bi]]
        for (ci in seq_along(b$columns)) {
            top_gap <- if (ci == 1L) 0.03 else 0
            circlize::circos.track(ylim = c(0, 1),
                track.height = lc$ring_height, bg.border = NA,
                track.margin = c(0, top_gap),
                panel.fun = function(x, y) {
                    i <- rows_of()
                    circlize::circos.rect(seq_along(i) - 1, 0, seq_along(i),
                        1, col = b$fill(b$mat[i, ci]), border = "grey15",
                        lwd = 0.45)
                })
            track_y[[length(track_y) + 1L]] <- list(
                track = circlize::get.current.track.index(),
                label = b$columns[[ci]], title = if (ci == 1L) b$title)
        }
    }
    circlize::circos.track(ylim = c(0, 1), track.height = 0.06,
        bg.border = NA, track.margin = c(0, 0.01),
        panel.fun = function(x, y) {
            circlize::circos.text(circlize::CELL_META$xcenter, 0.4,
                circlize::CELL_META$sector.index, facing = "bending.inside",
                niceFacing = TRUE, cex = 0.95)
        })

    first <- grp_levels[[1]]
    for (t in track_y) {
        rou <- circlize::circlize(0, 0.5, sector.index = first,
            track.index = t$track)[1, "rou"]
        graphics::text(0, rou, t$label, cex = 0.62)
        if (!is.null(t$title)) {
            rt <- circlize::circlize(0, 1, sector.index = first,
                track.index = t$track)[1, "rou"]
            graphics::text(0, rt + 0.018, t$title, cex = 0.72, font = 2)
        }
    }
    circlize::circos.clear()

    lg <- list(
        ComplexHeatmap::Legend(labels = band_lv,
            legend_gp = grid::gpar(fill = unname(band_col[band_lv])),
            title = lc$band, title_gp = grid::gpar(fontsize = 8,
                fontface = "bold"), labels_gp = grid::gpar(fontsize = 7)),
        ComplexHeatmap::Legend(labels = unique(cls), type = "points",
            pch = 15, background = NA,
            legend_gp = grid::gpar(col = unname(lab_col[match(unique(cls),
                cls)])),
            title = lc$label_class, title_gp = grid::gpar(fontsize = 8,
                fontface = "bold"), labels_gp = grid::gpar(fontsize = 7))
    )
    for (b in blocks) {
        lg[[length(lg) + 1L]] <- if (identical(b$type, "binary")) {
            ComplexHeatmap::Legend(labels = c("Yes", "No"),
                legend_gp = grid::gpar(fill = c(pal_ref(b$palette),
                    "white")), border = "grey15", title = b$title,
                title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                labels_gp = grid::gpar(fontsize = 7))
        } else {
            ComplexHeatmap::Legend(col_fun = b$col_fun, at = b$rng,
                labels = c("Low", "High"), title = b$title,
                direction = "horizontal", tick_length = grid::unit(0, "mm"),
                legend_width = grid::unit(18, "mm"), border = "grey15",
                title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                labels_gp = grid::gpar(fontsize = 7))
        }
    }
    pk <- do.call(ComplexHeatmap::packLegend,
        c(lg, list(direction = "vertical", gap = grid::unit(2.5, "mm"))))
    ComplexHeatmap::draw(pk, x = grid::unit(0.99, "npc"),
        y = grid::unit(0.01, "npc"), just = c("right", "bottom"))
}

draw_circos <- function() {
    if (is_locus) return(draw_locus())
    circlize::circos.par(
        start.degree = config$start.degree,
        gap.after = c(rep(config$gap.within, n_grp - 1L), config$gap.last),
        track.margin = c(0, 0.01),
        cell.padding = c(0, 0, 0, 0),
        circle.margin = 0.02
    )
    outside_labels <- identical(config$sector_label, "outside")
    if (outside_labels) {
        circlize::circos.par(circle.margin = 0.12,
            points.overflow.warning = FALSE)
    }
    circlize::circos.heatmap(
        mat,
        split = df[[g_grp]],
        cluster = FALSE,
        bg.border = config$sector_border,
        bg.lwd = 1,
        cell.border = config$cell_border,
        cell.lwd = 0.4,
        rownames.side = if (isTRUE(config$show_rownames)) "outside" else "none",
        rownames.cex = 0.55,
        col = heat_col,
        track.height = config$track_height
    )
    heat_track <- circlize::get.current.track.index()
    if (outside_labels) {
        circlize::circos.track(
            track.index = heat_track,
            bg.border = NA,
            panel.fun = function(x, y) {
                circlize::circos.text(
                    circlize::CELL_META$xcenter,
                    circlize::CELL_META$cell.ylim[[2]] +
                        circlize::convert_y(4, "mm"),
                    circlize::CELL_META$sector.index,
                    facing = "bending.inside",
                    niceFacing = TRUE,
                    cex = 0.95,
                    adj = c(0.5, 0),
                    col = "black"
                )
            }
        )
    }
    if (identical(config$column_label, "repel")) {
        cn <- colnames(mat)
        n <- length(cn)
        first <- grp_levels[[1]]
        ylim <- circlize::get.cell.meta.data("cell.ylim",
            sector.index = first, track.index = heat_track)
        x0 <- circlize::get.cell.meta.data("cell.xlim",
            sector.index = first, track.index = heat_track)[[1]]
        cell_h <- diff(ylim) / n
        y_mid <- seq(ylim[[1]] + cell_h / 2, ylim[[2]] - cell_h / 2,
            length.out = n)
        pt <- circlize::circlize(rep(x0, n), y_mid,
            sector.index = first, track.index = heat_track)
        th <- pt[, "theta"] * pi / 180
        ex <- pt[, "rou"] * cos(th)
        ey <- pt[, "rou"] * sin(th)
        cex <- 0.62
        # 1-D repel along y (ggrepel direction = "y"): merge overlapping
        # runs and centre each run on the mean of its anchors.
        gap <- graphics::strheight("Xg", cex = cex) * 1.25
        ord <- order(ey)
        blocks <- lapply(ord, function(i) list(idx = i, target = ey[[i]]))
        place <- function(b) {
            k <- length(b$idx)
            b$target - (k - 1) * gap / 2 + (seq_len(k) - 1) * gap
        }
        repeat {
            merged <- FALSE
            j <- 1L
            while (j < length(blocks)) {
                a <- blocks[[j]]
                b <- blocks[[j + 1L]]
                if (min(place(b)) - max(place(a)) < gap - 1e-9) {
                    idx <- c(a$idx, b$idx)
                    blocks[[j]] <- list(idx = idx, target = mean(ey[idx]))
                    blocks[[j + 1L]] <- NULL
                    merged <- TRUE
                } else {
                    j <- j + 1L
                }
            }
            if (!merged) break
        }
        ty <- numeric(n)
        for (b in blocks) ty[b$idx] <- place(b)
        pad <- graphics::strwidth("M", cex = cex) * 1.4
        tx <- ex - pad
        moved <- abs(ty - ey) > gap * 0.15
        graphics::segments(ex, ey, tx + pad * 0.25, ifelse(moved, ty, ey),
            col = "grey45", lwd = 0.5)
        graphics::text(tx, ty, cn, adj = c(1, 0.5), cex = cex,
            col = "black")
    }
    circlize::circos.track(
        track.index = heat_track,
        bg.border = NA,
        panel.fun = function(x, y) {
            if (!identical(config$column_label, "ticks") ||
                circlize::CELL_META$sector.numeric.index != n_grp) {
                return(invisible())
            }
            cn <- colnames(mat)
            n <- length(cn)
            cell_h <- diff(circlize::CELL_META$cell.ylim) / n
            y_coords <- seq(
                circlize::CELL_META$cell.ylim[[1]] + cell_h / 2,
                circlize::CELL_META$cell.ylim[[2]] - cell_h / 2,
                length.out = n
            )
            x0 <- circlize::CELL_META$cell.xlim[[2]]
            x1 <- x0 + circlize::convert_x(1, "mm")
            for (i in seq_len(n)) {
                circlize::circos.lines(
                    c(x0, x1),
                    c(y_coords[[i]], y_coords[[i]]),
                    col = "black",
                    lwd = 1
                )
            }
            circlize::circos.text(
                rep(x0, n) + circlize::convert_x(1.6, "mm"),
                y_coords,
                cn,
                cex = 0.55,
                adj = c(0, 0.5),
                facing = "inside",
                col = "black"
            )
        }
    )

    if (isTRUE(config$show_qval)) {
        circlize::circos.track(
            ylim = c(0, 1),
            track.height = 0.06,
            bg.border = NA,
            panel.fun = function(x, y) {
                sector <- as.character(circlize::CELL_META$sector.index)
                sub <- df[df[[g_grp]] == sector, , drop = FALSE]
                n <- nrow(sub)
                if (!n) {
                    return(invisible())
                }
                x_at <- circlize::CELL_META$xlim[[1]] +
                    (circlize::CELL_META$xlim[[2]] - circlize::CELL_META$xlim[[1]]) *
                    (seq_len(n) - 0.5) / n
                circlize::circos.points(
                    x_at,
                    rep(0.5, n),
                    pch = 18,
                    cex = 0.9,
                    col = q_col(sub[[g_q]])
                )
                if (circlize::CELL_META$sector.numeric.index == n_grp) {
                    x0 <- circlize::CELL_META$cell.xlim[[2]]
                    circlize::circos.lines(
                        c(x0, x0 + circlize::convert_x(1, "mm")),
                        c(0.5, 0.5),
                        col = "black",
                        lwd = 1
                    )
                    circlize::circos.text(
                        x0 + circlize::convert_x(1.6, "mm"),
                        0.5,
                        config$labels$qval,
                        cex = 0.55,
                        adj = c(0, 0.5),
                        facing = "inside",
                        col = "black"
                    )
                }
            }
        )
    }

    if (!outside_labels) {
        circlize::circos.track(
            ylim = c(0, 1),
            track.height = 0.07,
            bg.col = grDevices::adjustcolor(group_cols[grp_levels], alpha.f = 0.35),
            bg.border = NA,
            panel.fun = function(x, y) {
                circlize::circos.text(
                    circlize::CELL_META$xcenter,
                    0.45,
                    circlize::CELL_META$sector.index,
                    facing = "bending.inside",
                    niceFacing = TRUE,
                    cex = 0.7,
                    adj = c(0.5, 0.5),
                    col = "black"
                )
            }
        )
    }
    circlize::circos.clear()

    centre_lgd <- identical(config$legend_position, "center")
    lgd_args <- list(
        title = config$labels$fill,
        col_fun = heat_col,
        title_position = "leftcenter-rot",
        title_gp = grid::gpar(fontsize = 9),
        labels_gp = grid::gpar(fontsize = 8),
        background = NA
    )
    if (is.null(config$legend_labels)) {
        lgd_args$at <- pretty(v_rng, n = 4)
    } else {
        lgd_args$at <- v_rng
        lgd_args$labels <- config$legend_labels
    }
    if (centre_lgd) {
        lgd_args$legend_height <- grid::unit(32, "mm")
        lgd_args$grid_width <- grid::unit(4, "mm")
    }
    heat_lgd <- do.call(ComplexHeatmap::Legend, lgd_args)
    if (centre_lgd) {
        ComplexHeatmap::draw(heat_lgd, x = grid::unit(0.5, "npc"),
            y = grid::unit(0.5, "npc"), just = c("center", "center"))
    } else {
        ComplexHeatmap::draw(
            heat_lgd,
            x = grid::unit(0.98, "npc") - grid::unit(2, "mm"),
            y = grid::unit(0.92, "npc"),
            just = c("right", "top")
        )
    }
    if (isTRUE(config$show_qval)) {
        q_lgd <- ComplexHeatmap::Legend(
            title = config$labels$qval,
            col_fun = q_col,
            at = c(0, 0.25, 0.5, 0.75, 1),
            title_position = "leftcenter-rot",
            title_gp = grid::gpar(fontsize = 9),
            labels_gp = grid::gpar(fontsize = 8),
            background = NA
        )
        ComplexHeatmap::draw(
            q_lgd,
            x = grid::unit(0.98, "npc") - grid::unit(2, "mm"),
            y = grid::unit(0.22, "npc"),
            just = c("right", "bottom")
        )
    }

    if (identical(config$inner, "euler")) {
        sets <- lapply(grp_levels, function(g) {
            unique(df[[g_gene]][df[[g_grp]] == g])
        })
        names(sets) <- grp_levels
        grid::pushViewport(grid::viewport(
            x = 0.5,
            y = 0.5,
            width = 0.26,
            height = 0.26
        ))
        draw_euler(sets, group_cols, config$labels$inner)
        grid::popViewport()
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
dpi <- 300
if (identical(fmt, "png")) {
    grDevices::png(
        out_path, width = width, height = height, units = "in",
        res = dpi, bg = "transparent"
    )
} else if (identical(fmt, "pdf")) {
    grDevices::pdf(out_path, width = width, height = height)
} else {
    grDevices::svg(out_path, width = width, height = height, bg = "transparent")
}
old_par <- graphics::par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "transparent")
tryCatch(
    draw_circos(),
    finally = {
        graphics::par(old_par)
        grDevices::dev.off()
    }
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
