#!/usr/bin/env Rscript

# Template-ID: ideogram-coverage
#
# Purpose:
#   Draw a genomic locus browser: stacked coverage tracks, optional
#   chromatin-loop arcs, and gene intervals on one shared window.
#
# Inputs:
#   One row per coverage bin. Default example:
#     - chrom, start, end: genomic bin
#     - score: supplied signal
#     - track: lane name (cell type or sample)
#     - group: assay / class strip (scATAC vs HiChIP)
#   Optional sidecars beside the input:
#     - loops.tsv: start, end, padj (or a supplied score)
#     - genes.tsv: start, end, gene
#
# Output:
#   A PDF, PNG, or SVG multi-track locus figure.
#
# Dependencies:
#   ggplot2, readr, ggprism, patchwork, ggh4x, gtable, grid
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, region, highlight).
#   Extra tracks are extra rows. Loops and genes are sidecars, not
#   a second id. Group-to-group gap is config$group_gap. Do not add
#   an id per assay (ATAC / HiChIP) or for highlight on/off.
#   Whole-chromosome window fill is ideogram-density.
#   layout = "linked" (linked/example.tsv) is the peak-to-gene
#   variant: gene models, line tracks with gene-linked tile strips,
#   a coefficient track, SNP markers, scale bar and right-hand
#   expression violins, all from sidecars in linked/. Do not
#   compute links, z scores, coefficients, or colocalisation here.
#
# Scientific assumptions:
#   Coverage, loops, and gene intervals are supplied. This script
#   does not call peaks, loops, MACS, or Hi-C. -log10(p) is a
#   display transform of a supplied p / padj. Coordinates are
#   already on one reference.

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
        chrom = "chrom",
        start = "start",
        end = "end",
        score = "score",
        track = "track",
        group = "group"
    ),
    sidecars = list(
        loops = "loops.tsv",
        genes = "genes.tsv"
    ),
    loop = list(
        start = "start",
        end = "end",
        p = "padj"
    ),
    gene = list(
        start = "start",
        end = "end",
        name = "gene"
    ),
    highlight = c(127400000, 127418000),
    group_gap = 0.4,
    # "stack" = grouped filled tracks + loops (default). "linked" =
    # gene models on top, line tracks each with a gene-linked tile
    # strip (supplied z, optional flag ring), a coefficient track,
    # SNP / tile markers with a scale bar, and per-track expression
    # violins on the right. Sidecars sit beside the input; colours
    # are "Family.Name:index" catalog references.
    layout = "stack",
    linked = list(
        sidecars = list(
            models = "gene_models.tsv",     # gene, feature (gene|exon), start, end, lane
            tiles = "tiles.tsv",            # track, start, end, z, optional flag
            markers = "markers.tsv",        # pos, label, shape (triangle|tick), colour, lane
            expression = "expression.tsv"   # track, value (one row per cell)
        ),
        coef_track = "Coefficients",
        focus_gene = "BCL2",
        gene_colors = c(focus = "Diverging.RdBu:3", other = "Diverging.RdBu:7"),
        track_palette = "Qualitative.Bold",
        coef_color = "black",
        z_palette = "Gradient.Blues",
        z_limits = c(0, 2),
        flag_color = "Diverging.RdBu:1",
        scale_bar = 100000,
        title_gene = "BCL2",
        title_rest = "in pituitary (\u03c1 = 0.35)",
        y_left = "Accessibility",
        y_right = "Gene expression",
        heights = c(gene = 1.5, track = 1.5, tiles = 0.32, coef = 1.9,
            markers = 0.55),
        widths = c(3.2, 1),
        size = c(width = 9, height = 10.5)
    ),
    labels = list(
        title = NULL,
        x = "Mb",
        loop = "-\u202flog10(P adj)"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "patchwork", "ggh4x", "gtable"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
need <- config$columns[c("start", "end", "score", "track")]
require_columns(df, need)

st <- config$columns$start
en <- config$columns$end
sc <- config$columns$score
tr <- config$columns$track
gr <- config$columns$group
chr_col <- config$columns$chrom

df[[st]] <- as.numeric(df[[st]])
df[[en]] <- as.numeric(df[[en]])
df[[sc]] <- as.numeric(df[[sc]])
ok <- is.finite(df[[st]]) & is.finite(df[[en]]) & is.finite(df[[sc]]) &
    df[[sc]] >= 0
if (any(!ok)) {
    message("Dropped ", sum(!ok), " coverage row(s) with invalid bins.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No coverage rows left to plot.", call. = FALSE)
}
swap <- df[[st]] > df[[en]]
if (any(swap)) {
    tmp <- df[[st]][swap]
    df[[st]][swap] <- df[[en]][swap]
    df[[en]][swap] <- tmp
}
df$mid <- (df[[st]] + df[[en]]) / 2
df[[tr]] <- factor(df[[tr]], levels = unique(as.character(df[[tr]])))
if (!is.null(gr) && gr %in% names(df)) {
    df[[gr]] <- factor(df[[gr]], levels = unique(as.character(df[[gr]])))
} else {
    gr <- NULL
}

x_min <- min(df[[st]], na.rm = TRUE)
x_max <- max(df[[en]], na.rm = TRUE)

if (identical(config$layout, "linked")) {
    lk <- config$linked
    pref <- function(ref) {
        if (!grepl(".", ref, fixed = TRUE)) return(ref)
        parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
        cols <- palette_colors(parts[[1]])
        if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
    }
    read_side <- function(key, cols) {
        fn <- lk$sidecars[[key]]
        if (is.null(fn)) return(NULL)
        path <- file.path(dirname(io$input), fn)
        if (!file.exists(path)) return(NULL)
        d <- as.data.frame(read_table_auto(path), stringsAsFactors = FALSE)
        require_columns(d, as.list(stats::setNames(cols, cols)))
        d
    }
    models <- read_side("models", c("gene", "feature", "start", "end", "lane"))
    tiles <- read_side("tiles", c("track", "start", "end", "z"))
    markers <- read_side("markers", c("pos", "label", "shape", "colour",
        "lane"))
    expr <- read_side("expression", c("track", "value"))

    trk <- levels(df[[tr]])
    has_coef <- lk$coef_track %in% trk
    cell_tracks <- setdiff(trk, lk$coef_track)
    tile_row <- function(t) paste0(t, "\rtiles")
    rows <- c(
        if (!is.null(models)) "gene",
        unlist(lapply(cell_tracks, function(t) {
            c(t, if (!is.null(tiles)) tile_row(t))
        })),
        if (has_coef) lk$coef_track,
        if (!is.null(markers)) "markers"
    )
    kind <- ifelse(rows == "gene", "gene", ifelse(rows == "markers", "markers",
        ifelse(grepl("\rtiles$", rows), "tiles",
            ifelse(rows == lk$coef_track, "coef", "track"))))
    heights <- unname(lk$heights[kind])
    row_lab <- ifelse(kind %in% c("track", "coef"), rows, "")
    names(row_lab) <- rows
    tcol <- stats::setNames(palette_colors(lk$track_palette,
        n = length(cell_tracks)), cell_tracks)
    tcol[[lk$coef_track]] <- pref(lk$coef_color)
    rowf <- function(x) factor(x, levels = rows)

    cov <- data.frame(row = rowf(as.character(df[[tr]])), x = df$mid,
        y = df[[sc]], col = unname(tcol[as.character(df[[tr]])]))
    ax <- stats::aggregate(y ~ row, cov, max)
    left <- ggplot() +
        geom_segment(data = ax, aes(x = x_min, xend = x_min, y = 0,
            yend = y), linewidth = 0.35) +
        geom_line(data = cov, aes(x = x, y = y, colour = col, group = row),
            linewidth = 0.35) +
        scale_colour_identity()

    if (!is.null(tiles)) {
        tiles$z <- as.numeric(tiles$z)
        tiles$flag <- if ("flag" %in% names(tiles)) as.integer(tiles$flag) else 0L
        zl <- as.numeric(lk$z_limits)
        td <- data.frame(row = rowf(tile_row(as.character(tiles$track))),
            x = (as.numeric(tiles$start) + as.numeric(tiles$end)) / 2,
            z = pmin(pmax(tiles$z, zl[[1]]), zl[[2]]), flag = tiles$flag)
        td <- td[!is.na(td$row), , drop = FALSE]
        td <- td[order(td$z), , drop = FALSE]
        zpal <- pref(lk$z_palette)
        bg <- data.frame(row = rowf(tile_row(cell_tracks)), y0 = 0.1,
            y1 = 0.9)
        left <- left +
            geom_rect(data = bg, aes(xmin = -Inf, xmax = Inf, ymin = y0,
                ymax = y1), fill = zpal[[2]], alpha = 0.45) +
            geom_point(data = td, aes(x = x, y = 0.5, fill = z), shape = 21,
                colour = "transparent", size = 2.1) +
            geom_point(data = td[td$flag == 1L, , drop = FALSE],
                aes(x = x, y = 0.5), shape = 21, fill = NA,
                colour = pref(lk$flag_color), stroke = 0.8, size = 2.7) +
            scale_fill_gradientn(colours = zpal, limits = zl, guide = "none")
    }

    if (!is.null(models)) {
        models$start <- as.numeric(models$start)
        models$end <- as.numeric(models$end)
        models$lane <- as.integer(models$lane)
        models$ly <- ifelse(models$lane == 1L, 0.74, 0.3)
        models$col <- ifelse(models$gene == lk$focus_gene,
            pref(lk$gene_colors[["focus"]]), pref(lk$gene_colors[["other"]]))
        models$row <- rowf("gene")
        gl <- models[models$feature == "gene", , drop = FALSE]
        ex <- models[models$feature == "exon", , drop = FALSE]
        gl$ty <- ifelse(gl$lane == 1L, gl$ly + 0.17, gl$ly - 0.17)
        gl$vj <- ifelse(gl$lane == 1L, 0, 1)
        left <- left +
            geom_rect(data = gl, aes(xmin = start, xmax = end,
                ymin = ly - 0.035, ymax = ly + 0.035), fill = gl$col) +
            geom_rect(data = ex, aes(xmin = start, xmax = end,
                ymin = ly - 0.11, ymax = ly + 0.11), fill = ex$col) +
            geom_text(data = gl, aes(x = start, y = ty, label = gene,
                vjust = vj), hjust = 0, size = 3, fontface = "italic")
    }

    if (!is.null(markers)) {
        markers$pos <- as.numeric(markers$pos)
        markers$col <- vapply(markers$colour, pref, character(1))
        markers$y <- ifelse(as.integer(markers$lane) == 1L, 0.72, 0.3)
        markers$pch <- ifelse(markers$shape == "tick", 124L, 17L)
        markers$row <- rowf("markers")
        markers <- markers[rev(seq_len(nrow(markers))), , drop = FALSE]
        sb <- data.frame(row = rowf("markers"), x = x_min + (x_max - x_min) *
            0.01, xe = x_min + (x_max - x_min) * 0.01 + lk$scale_bar)
        left <- left +
            geom_point(data = markers, aes(x = pos, y = y),
                shape = markers$pch, colour = markers$col, size = 2.4) +
            geom_segment(data = sb, aes(x = x, xend = xe, y = 0.5, yend = 0.5),
                linewidth = 2, colour = "grey30") +
            geom_text(data = sb, aes(x = xe, y = 0.5,
                label = paste0(" ", lk$scale_bar / 1000, " kb")),
                hjust = 0, size = 3)
    }

    brk <- function(lim) {
        b <- pretty(c(0, lim[[2]]), 2)
        b[b >= 0 & b <= lim[[2]]]
    }
    yscales <- lapply(kind, function(k) {
        if (k %in% c("track", "coef")) {
            scale_y_continuous(breaks = brk, limits = c(0, NA),
                expand = expansion(mult = c(0.02, 0.08)))
        } else {
            scale_y_continuous(limits = c(0, 1), breaks = NULL,
                expand = c(0, 0))
        }
    })
    chr_lab <- if (!is.null(chr_col) && chr_col %in% names(df)) {
        as.character(df[[chr_col]][[1]])
    } else {
        "chr"
    }
    xb <- pretty(c(x_min, x_max), n = 6)
    xb <- xb[xb >= x_min & xb <= x_max]
    lab_fun <- function(x) unname(row_lab[as.character(x)])
    th <- theme_classic(base_size = 9) +
        theme(
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9),
            axis.line = element_blank(),
            panel.spacing.y = unit(1.5, "pt"),
            plot.background = element_blank(),
            panel.background = element_blank()
        )
    left <- left +
        facet_grid(row ~ ., scales = "free_y", switch = "y",
            labeller = labeller(row = lab_fun)) +
        ggh4x::force_panelsizes(rows = heights) +
        ggh4x::facetted_pos_scales(y = yscales) +
        scale_x_continuous(breaks = xb, labels = sprintf("%.2f", xb / 1e6),
            expand = c(0, 0)) +
        coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
        labs(x = sprintf("%s:%s\u2013%s (Mb)", chr_lab,
            format(round(x_min), big.mark = ","),
            format(round(x_max), big.mark = ",")), y = NULL) +
        th +
        theme(axis.line.x = element_line(linewidth = 0.35))

    right <- ggplot()
    xr <- c(0, 1)
    if (!is.null(expr)) {
        expr$value <- as.numeric(expr$value)
        ed <- data.frame(row = rowf(as.character(expr$track)),
            v = expr$value, col = unname(tcol[as.character(expr$track)]))
        ed <- ed[!is.na(ed$row), , drop = FALSE]
        xr <- range(c(0, ed$v))
        med <- stats::aggregate(v ~ row, ed, stats::median)
        right <- right +
            geom_violin(data = ed, aes(x = v, y = 0.5, group = row,
                fill = col), orientation = "y", width = 0.9, scale = "width",
                colour = "grey20", linewidth = 0.3, trim = TRUE) +
            geom_boxplot(data = ed, aes(x = v, y = 0.5, group = row),
                orientation = "y", width = 0.1, outlier.shape = NA,
                fill = "grey20", colour = "grey20", linewidth = 0.3) +
            geom_point(data = med, aes(x = v, y = 0.5), shape = 21,
                fill = "white", colour = "transparent", size = 1.3) +
            scale_fill_identity()
    }
    span <- diff(xr)
    if (!is.null(tiles) && "gene" %in% rows) {
        zpal <- pref(lk$z_palette)
        nz <- 40
        gx <- seq(xr[[1]], xr[[2]], length.out = nz + 1)
        zb <- data.frame(row = rowf("gene"), xmin = gx[-(nz + 1)],
            xmax = gx[-1], f = seq(0, 1, length.out = nz))
        zl <- as.numeric(lk$z_limits)
        right <- right +
            geom_rect(data = zb, aes(xmin = xmin, xmax = xmax, ymin = 0.62,
                ymax = 0.8), fill = grDevices::colorRampPalette(zpal)(nz),
                colour = NA) +
            geom_rect(data = data.frame(row = rowf("gene"), a = xr[[1]],
                b = xr[[2]]), aes(xmin = a, xmax = b, ymin = 0.62,
                ymax = 0.8), fill = NA, colour = "grey20", linewidth = 0.3) +
            geom_text(data = data.frame(row = rowf("gene"),
                x = c(xr[[1]], xr[[2]], mean(xr)),
                lab = c(sprintf("\u2264%g", zl[[1]]),
                    sprintf("\u2265%.1f", zl[[2]]), "Z score"),
                hj = c(0, 1, 0.5)),
                aes(x = x, y = 0.48, label = lab, hjust = hj), vjust = 1,
                size = 2.8)
    }
    if (!is.null(markers) && has_coef) {
        ml <- markers[rev(seq_len(nrow(markers))), , drop = FALSE]
        ml <- ml[!duplicated(ml$label), , drop = FALSE]
        n <- nrow(ml)
        ml$row <- rowf(lk$coef_track)
        ml$ly <- seq(0.8, 0.2, length.out = max(n, 2))[seq_len(n)]
        right <- right +
            geom_point(data = ml, aes(x = xr[[1]] - span * 0.02, y = ly),
                shape = ml$pch, colour = ml$col, size = 2) +
            geom_text(data = ml, aes(x = xr[[1]] + span * 0.04, y = ly,
                label = label), hjust = 0, size = 2.7)
    }
    right <- right +
        facet_grid(row ~ ., drop = FALSE) +
        ggh4x::force_panelsizes(rows = heights) +
        scale_x_continuous(breaks = pretty(xr, 3)) +
        scale_y_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0),
            name = lk$y_left, sec.axis = dup_axis(name = lk$y_right)) +
        coord_cartesian(xlim = xr, clip = "off") +
        labs(x = NULL) +
        th +
        theme(strip.text = element_blank(),
            axis.line.x = element_line(linewidth = 0.35),
            axis.title.y.left = element_text(size = 10, margin = margin(r = 4)),
            axis.title.y.right = element_text(size = 10, angle = 90,
                margin = margin(l = 6)))

    g <- lk$title_gene
    rest <- lk$title_rest
    p <- patchwork::wrap_plots(left, right, nrow = 1, widths = lk$widths) +
        patchwork::plot_annotation(
            title = bquote(italic(.(g)) ~ .(rest)),
            theme = theme(plot.title = element_text(hjust = 0.5, size = 13)))
    save_ggplot(p, io$output, width = lk$size[["width"]],
        height = lk$size[["height"]])
    quit(save = "no", status = 0)
}

hl <- as.numeric(config$highlight)
hl <- hl[is.finite(hl)]
if (length(hl) != 2L) {
    hl <- NULL
}

side_path <- function(name) {
    if (is.null(name) || !nzchar(name)) {
        return(NA_character_)
    }
    file.path(dirname(io$input), name)
}

loops <- NULL
loop_path <- side_path(config$sidecars$loops)
if (!is.na(loop_path) && file.exists(loop_path)) {
    loops <- read_table_auto(loop_path)
    require_columns(loops, config$loop)
    ls <- config$loop$start
    le <- config$loop$end
    lp <- config$loop$p
    loops[[ls]] <- as.numeric(loops[[ls]])
    loops[[le]] <- as.numeric(loops[[le]])
    loops[[lp]] <- as.numeric(loops[[lp]])
    ok_l <- is.finite(loops[[ls]]) & is.finite(loops[[le]]) &
        is.finite(loops[[lp]]) & loops[[lp]] > 0
    if (any(!ok_l)) {
        message("Dropped ", sum(!ok_l), " loop row(s).")
        loops <- loops[ok_l, , drop = FALSE]
    }
    if (!nrow(loops)) {
        loops <- NULL
    } else {
        swap_l <- loops[[ls]] > loops[[le]]
        if (any(swap_l)) {
            tmp <- loops[[ls]][swap_l]
            loops[[ls]][swap_l] <- loops[[le]][swap_l]
            loops[[le]][swap_l] <- tmp
        }
        loops$mlogp <- -log10(pmax(loops[[lp]], .Machine$double.xmin))
    }
}

genes <- NULL
gene_path <- side_path(config$sidecars$genes)
if (!is.na(gene_path) && file.exists(gene_path)) {
    genes <- read_table_auto(gene_path)
    require_columns(genes, config$gene)
    gs <- config$gene$start
    ge <- config$gene$end
    gn <- config$gene$name
    genes[[gs]] <- as.numeric(genes[[gs]])
    genes[[ge]] <- as.numeric(genes[[ge]])
    genes[[gn]] <- as.character(genes[[gn]])
    ok_g <- is.finite(genes[[gs]]) & is.finite(genes[[ge]])
    if (any(!ok_g)) {
        genes <- genes[ok_g, , drop = FALSE]
    }
    if (!nrow(genes)) {
        genes <- NULL
    }
}

chr_lab <- if (!is.null(chr_col) && chr_col %in% names(df)) {
    as.character(df[[chr_col]][[1]])
} else {
    "chr"
}
title <- config$labels$title
if (is.null(title) || !nzchar(title)) {
    title <- sprintf(
        "%s: %s\u2013%s",
        chr_lab,
        format(round(x_min), big.mark = ",", scientific = FALSE),
        format(round(x_max), big.mark = ",", scientific = FALSE)
    )
}

alg <- palette_colors("Brand.Algolia")
safe <- palette_colors("Qualitative.Safe")
blu <- palette_colors("Quantitative.BluGrn")
df$fill_key <- as.character(df[[tr]])
if (!is.null(gr)) {
    hichip <- grepl("hichip|h3k27", as.character(df[[gr]]), ignore.case = TRUE)
    df$fill_key[hichip] <- "HiChIP"
}
fill_lv <- unique(df$fill_key)
fill_vals <- palette_colors("Brand.Algolia", n = length(fill_lv))
names(fill_vals) <- fill_lv
if ("HiChIP" %in% names(fill_vals)) {
    fill_vals[["HiChIP"]] <- blu[[5]]
}
hl_col <- alg[[5]]
ink <- palette_colors("Brand.Emma")[[1]]
outline <- safe[[12]]

x_breaks <- pretty(c(x_min, x_max), n = 4)
x_labs <- sprintf("%.1f", x_breaks / 1e6)
base_th <- theme_prism() +
    theme(
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = outline, fill = NA, linewidth = 0.4),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        legend.key = element_blank(),
        panel.spacing.y = unit(0, "pt")
    )

hl_df <- NULL
if (!is.null(hl)) {
    hl_df <- data.frame(xmin = min(hl), xmax = max(hl))
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
add_highlight <- function(p) {
    if (is.null(hl_df)) {
        return(p)
    }
    p + geom_rect(
        data = hl_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = hl_col,
        alpha = 0.18,
        colour = NA
    )
}

tighten_left <- function(gt) {
    hit <- which(gt$layout$name == "guide-box-left")
    if (length(hit)) {
        gt$grobs[[hit[[1]]]] <- ggplot2::zeroGrob()
        col <- unique(gt$layout$l[gt$layout$name == "guide-box-left"])
        gt$widths[col] <- grid::unit(0, "pt")
    }
    lay <- gt$layout
    inner_i <- which(grepl("^strip-l-", lay$name) & lay$t == lay$b)
    outer_i <- which(grepl("^strip-l-", lay$name) & lay$t != lay$b)
    if (!length(inner_i) || !length(outer_i)) {
        return(gt)
    }
    strip_col <- unique(lay$l[c(inner_i, outer_i)])[[1]]
    for (i in outer_i) {
        sg <- gt$grobs[[i]]
        if (inherits(sg, "gtable") && length(sg$widths) >= 2L) {
            nest <- which(sg$layout$name == "nester")
            if (length(nest)) {
                sg$layout$l[nest] <- ncol(sg)
                sg$layout$r[nest] <- ncol(sg)
            }
            sg$widths[[length(sg$widths)]] <- grid::unit(7, "mm")
            gt$grobs[[i]] <- sg
        }
    }
    gt$widths[strip_col] <- grid::grobWidth(gt$grobs[[outer_i[[1]]]])
    for (i in outer_i) {
        sg <- gt$grobs[[i]]
        gt$grobs[[i]] <- grid::grobTree(
            sg,
            vp = grid::viewport(
                x = 1,
                just = "right",
                width = grid::grobWidth(sg),
                height = grid::unit(1, "npc")
            )
        )
    }
    track_names <- as.character(levels(df[[tr]]))
    inner_w <- max(grid::stringWidth(track_names)) + grid::unit(8, "pt")
    gt <- gtable::gtable_add_cols(gt, inner_w, pos = strip_col)
    gt$layout$l[inner_i] <- strip_col + 1L
    gt$layout$r[inner_i] <- strip_col + 1L
    for (i in inner_i) {
        sg <- gt$grobs[[i]]
        gt$grobs[[i]] <- grid::grobTree(
            sg,
            vp = grid::viewport(
                x = 1,
                just = "right",
                width = grid::grobWidth(sg),
                height = grid::unit(1, "npc")
            )
        )
    }
    gt
}

place_loop_legend <- function(gt, legend) {
    if (is.null(legend)) {
        return(gt)
    }
    pan <- gt$layout[grepl("^panel-", gt$layout$name), , drop = FALSE]
    pan <- pan[order(pan$t, pan$l), , drop = FALSE]
    if (!nrow(pan)) {
        return(gt)
    }
    outer <- which(
        grepl("^strip-l-", gt$layout$name) & gt$layout$t != gt$layout$b
    )
    if (!length(outer)) {
        return(gt)
    }
    strip_col <- unique(gt$layout$l[outer])[[1]]
    loops_t <- pan$t[nrow(pan)]
    xlab <- gt$layout$t[gt$layout$name == "xlab-b"]
    loops_b <- if (length(xlab)) max(pan$b[nrow(pan)], xlab[[1]]) else pan$b[nrow(pan)]
    wrapped <- grid::grobTree(
        legend,
        vp = grid::viewport(
            x = 0.5,
            y = 0.08,
            just = c("center", "bottom"),
            width = grid::grobWidth(legend),
            height = grid::grobHeight(legend)
        )
    )
    gtable::gtable_add_grob(
        gt,
        wrapped,
        t = loops_t,
        b = loops_b,
        l = strip_col,
        r = strip_col,
        name = "loop-legend",
        clip = "off",
        z = Inf
    )
}

insert_group_gaps <- function(gt, n_per_group, gap) {
    if (length(n_per_group) < 2L || gap <= 0) {
        return(gt)
    }
    pan <- gt$layout[grepl("^panel-", gt$layout$name), , drop = FALSE]
    pan <- pan[order(pan$t, pan$l), , drop = FALSE]
    if (nrow(pan) < sum(n_per_group)) {
        return(gt)
    }
    ends <- cumsum(n_per_group)
    ends <- ends[-length(ends)]
    for (i in rev(ends)) {
        gt <- gtable::gtable_add_rows(
            gt,
            grid::unit(gap, "null"),
            pos = pan$t[[i]]
        )
    }
    gt
}

tracks_per_group <- function() {
    if (is.null(gr)) {
        return(length(levels(df[[tr]])))
    }
    vapply(levels(df[[gr]]), function(g) {
        length(unique(as.character(df[[tr]][df[[gr]] == g])))
    }, integer(1))
}

loop_group <- "\u00a0"
n_per_group <- tracks_per_group()
group_lv <- if (is.null(gr)) {
    NULL
} else if (is.null(loops)) {
    levels(df[[gr]])
} else {
    c(as.character(levels(df[[gr]])), loop_group)
}
track_lv <- if (is.null(loops)) {
    levels(df[[tr]])
} else {
    c(as.character(levels(df[[tr]])), "Loops")
}

if (!is.null(gr) && !is.null(group_lv)) {
    df[[gr]] <- factor(as.character(df[[gr]]), levels = group_lv)
}
df[[tr]] <- factor(as.character(df[[tr]]), levels = track_lv)

blank_group <- function(x) {
    x[x == loop_group] <- ""
    x
}
facet_lab <- if (!is.null(gr)) {
    args <- list()
    args[[gr]] <- blank_group
    args[[tr]] <- identity
    do.call(ggplot2::labeller, args)
} else {
    "label_value"
}

p <- ggplot(df, aes(x = mid, y = .data[[sc]]))
p <- add_highlight(p)
p <- p +
    geom_ribbon(
        aes(ymin = 0, ymax = .data[[sc]], fill = fill_key),
        colour = NA
    ) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_x_continuous(
        limits = c(x_min, x_max),
        breaks = x_breaks,
        labels = x_labs,
        expand = c(0, 0)
    ) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    labs(title = title, x = config$labels$x, y = NULL) +
    base_th +
    theme(
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.text.x = element_text(size = 8),
        axis.ticks.x = element_line(colour = ink),
        axis.title.x = element_text(size = 9),
        legend.position = "left",
        legend.justification = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 8),
        plot.margin = margin(t = 4, r = 8, b = 4, l = 4)
    )

drew_genes <- FALSE
if (!is.null(loops)) {
    ls <- config$loop$start
    le <- config$loop$end
    span <- loops[[le]] - loops[[ls]]
    max_span <- max(span, na.rm = TRUE)
    if (!is.finite(max_span) || max_span <= 0) {
        max_span <- 1
    }
    y_floor <- if (is.null(genes)) 0 else 0.22
    arcs <- do.call(rbind, lapply(seq_len(nrow(loops)), function(i) {
        th <- seq(0, pi, length.out = 64)
        h <- (1 - y_floor - 0.04) * (0.45 + 0.55 * (span[[i]] / max_span))
        out <- data.frame(
            id = i,
            x = loops[[ls]][[i]] + span[[i]] * (1 - cos(th)) / 2,
            y = 1 - h * sin(th),
            mlogp = loops$mlogp[[i]]
        )
        out[[tr]] <- "Loops"
        if (!is.null(gr)) {
            out[[gr]] <- loop_group
        }
        out
    }))
    arcs[[tr]] <- factor(arcs[[tr]], levels = track_lv)
    if (!is.null(gr)) {
        arcs[[gr]] <- factor(arcs[[gr]], levels = group_lv)
        n_per_group <- c(n_per_group, 1L)
    }
    p <- p +
        geom_line(
            data = arcs,
            aes(x = x, y = y, group = id, colour = mlogp),
            inherit.aes = FALSE,
            linewidth = 0.55,
            lineend = "round"
        ) +
        scale_colour_gradientn(
            name = config$labels$loop,
            colours = blu,
            guide = guide_colorbar(
                title.position = "top",
                title.hjust = 0.5,
                direction = "horizontal",
                barwidth = unit(1.7, "cm"),
                barheight = unit(0.28, "cm"),
                order = 1
            )
        )
    if (!is.null(genes)) {
        gs <- config$gene$start
        ge <- config$gene$end
        gn <- config$gene$name
        gene_df <- genes
        gene_df[[tr]] <- factor("Loops", levels = track_lv)
        if (!is.null(gr)) {
            gene_df[[gr]] <- factor(loop_group, levels = group_lv)
        }
        p <- p +
            geom_rect(
                data = gene_df,
                aes(
                    xmin = .data[[gs]],
                    xmax = .data[[ge]],
                    ymin = 0.04,
                    ymax = 0.16
                ),
                inherit.aes = FALSE,
                fill = alg[[4]],
                colour = NA
            ) +
            geom_text(
                data = gene_df,
                aes(
                    x = .data[[gs]],
                    y = 0.10,
                    label = .data[[gn]]
                ),
                inherit.aes = FALSE,
                colour = alg[[4]],
                size = 3.2,
                fontface = "bold.italic",
                hjust = 1.15,
                vjust = 0.5
            )
        drew_genes <- TRUE
    }
}

if (!is.null(gr)) {
    p <- p + ggh4x::facet_nested(
        stats::as.formula(paste(gr, "+", tr, "~ .")),
        scales = "free_y",
        switch = "y",
        axes = "margins",
        labeller = facet_lab,
        nest_line = element_line(colour = outline, linewidth = 0.35),
        solo_line = FALSE,
        strip = ggh4x::strip_nested(size = "variable")
    )
} else {
    p <- p + facet_grid(
        stats::as.formula(paste(tr, "~ .")),
        scales = "free_y",
        switch = "y"
    )
}

gt <- ggplotGrob(p)
loop_legend <- NULL
leg_hit <- which(gt$layout$name == "guide-box-left")
if (length(leg_hit)) {
    loop_legend <- gt$grobs[[leg_hit[[1]]]]
}
gap <- as.numeric(config$group_gap)
if (!is.finite(gap) || gap < 0) {
    gap <- 0
}
if (!is.null(gr)) {
    n_gap <- n_per_group
    if (!is.null(loops) && length(n_gap) > 1L) {
        n_gap <- n_gap[-length(n_gap)]
    }
    gt <- insert_group_gaps(gt, n_gap, gap)
}
gt <- tighten_left(gt)
if (!is.null(loops)) {
    gt <- place_loop_legend(gt, loop_legend)
}

if (!is.null(genes) && !drew_genes) {
    gs <- config$gene$start
    ge <- config$gene$end
    gn <- config$gene$name
    p_gene <- ggplot(genes)
    p_gene <- add_highlight(p_gene)
    p_gene <- p_gene +
        geom_rect(
            aes(xmin = .data[[gs]], xmax = .data[[ge]], ymin = 0.28, ymax = 0.72),
            fill = alg[[4]],
            colour = NA
        ) +
        geom_text(
            aes(
                x = (.data[[gs]] + .data[[ge]]) / 2,
                y = 0.12,
                label = .data[[gn]]
            ),
            colour = alg[[4]],
            size = 3.2,
            fontface = "bold.italic",
            hjust = 0.5
        ) +
        scale_x_continuous(
            limits = c(x_min, x_max),
            breaks = x_breaks,
            labels = x_labs,
            expand = c(0, 0)
        ) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
        labs(x = config$labels$x, y = NULL) +
        base_th +
        theme(
            axis.text.x = element_text(size = 8),
            axis.ticks.x = element_line(colour = ink),
            axis.title.x = element_text(size = 9),
            plot.margin = margin(t = 0, r = 4, b = 4, l = 4)
        )
    p <- patchwork::wrap_plots(gt, p_gene, ncol = 1, heights = c(11, 0.7))
} else {
    p <- gt
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8.8, height = 7.8)
