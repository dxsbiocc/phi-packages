#!/usr/bin/env Rscript

# Template-ID: bar-circos-links
#
# Purpose:
#   Circos sector bar chart: one sector per category (e.g. disease), an
#   outer header band with the sector name, a ring of item labels (e.g.
#   genes) coloured by a supplied direction, an inner bar track of a
#   supplied score, and optional central links between items in
#   different sectors (e.g. shared genes).
#
# Inputs:
#   genes.tsv — one row per item-in-sector:
#     - disease: sector name
#     - gene: item label
#     - score: non-negative bar height
#     - direction: two-level label (default "positive" / "negative")
#   links.tsv — one row per link between two items:
#     - from_disease, from_gene, to_disease, to_gene
#     Each endpoint must match a row of genes.tsv. An empty table
#     (header only) or draw_links = FALSE draws no links.
#
# Output:
#   A PDF, PNG, or SVG circular sector bar chart.
#
# Dependencies:
#   circlize, readr
#
# Example:
#   Rscript plot.R genes.tsv links.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, sector order, within-sector sort,
#   direction colours, track heights, score axis, links on/off, title). Extra
#   diseases or genes are extra rows, not extra ids. Do not run GWAS,
#   Mendelian randomization, or an overlap test in this script.
#
# Scientific assumptions:
#   score and direction are supplied. Bars share one radial scale
#   (0 to score_max). Links are drawn only between supplied endpoints;
#   this script does not infer shared genes or link weights.

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

io <- parse_io_args(n_input = 2L, input_names = c("genes", "links"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        sector = "disease",
        item = "gene",
        score = "score",
        direction = "direction"
    ),
    link_columns = list(
        from_sector = "from_disease",
        from_item = "from_gene",
        to_sector = "to_disease",
        to_item = "to_gene"
    ),
    # NULL keeps first-appearance order
    sectors = NULL,
    sort_within = "score",       # "score" (descending) or "input"
    directions = c(positive = "positive", negative = "negative"),
    direction_labels = c(positive = "Positive", negative = "Negative"),
    palettes = list(
        direction = "Diverging.RdBu",
        link = "Qualitative.Safe"
    ),
    # palette indices: RdBu 8 = blue, 2 = red-orange; Safe 11 = steel blue
    direction_index = c(positive = 8L, negative = 2L),
    link_index = 11L,
    link_alpha = 0.35,
    draw_links = TRUE,
    score_max = NULL,            # NULL = max(score)
    score_breaks = NULL,         # NULL = pretty()
    legend_title = "Direction",
    track_height = list(
        header = 0.065,
        label = 0.13,
        bar = 0.24
    ),
    header_fill = "#e5e5e5",
    header_border = "#b0b0b0",
    header_cex = 0.62,
    label_cex = 0.55,
    gap_degree = 2,
    axis_gap_degree = 9,
    title = NULL,
    size = list(
        width = 9,
        height = 9
    )
)

load_packages(c("circlize", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
genes <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
require_columns(genes, config$columns)
s_col <- config$columns$sector
i_col <- config$columns$item
v_col <- config$columns$score
d_col <- config$columns$direction
genes[[s_col]] <- as.character(genes[[s_col]])
genes[[i_col]] <- as.character(genes[[i_col]])
genes[[v_col]] <- as.numeric(genes[[v_col]])
genes[[d_col]] <- as.character(genes[[d_col]])
if (any(!nzchar(genes[[s_col]]) | !nzchar(genes[[i_col]]))) {
    stop("Sector and item columns must be non-empty.", call. = FALSE)
}
if (any(!is.finite(genes[[v_col]]) | genes[[v_col]] < 0)) {
    stop("score must be finite and non-negative.", call. = FALSE)
}
bad_dir <- setdiff(unique(genes[[d_col]]), unname(config$directions))
if (length(bad_dir)) {
    stop(sprintf("Unknown direction value(s): %s",
        paste(bad_dir, collapse = ", ")), call. = FALSE)
}
if (anyDuplicated(paste(genes[[s_col]], genes[[i_col]], sep = "\t"))) {
    stop("Each item may appear only once per sector.", call. = FALSE)
}

sectors <- if (is.null(config$sectors)) unique(genes[[s_col]]) else
    as.character(config$sectors)
genes <- genes[genes[[s_col]] %in% sectors, , drop = FALSE]
if (!nrow(genes)) stop("No rows left after sector filter.", call. = FALSE)
sectors <- sectors[sectors %in% genes[[s_col]]]

ord <- if (identical(config$sort_within, "score")) {
    order(match(genes[[s_col]], sectors), -genes[[v_col]])
} else {
    order(match(genes[[s_col]], sectors), seq_len(nrow(genes)))
}
genes <- genes[ord, , drop = FALSE]
genes$x <- stats::ave(seq_len(nrow(genes)), genes[[s_col]], FUN = seq_along)
n_per <- table(factor(genes[[s_col]], levels = sectors))

dir_pal <- palette_colors(config$palettes$direction)
dir_cols <- c(
    dir_pal[[config$direction_index[["positive"]]]],
    dir_pal[[config$direction_index[["negative"]]]]
)
names(dir_cols) <- unname(config$directions[c("positive", "negative")])
genes$col <- unname(dir_cols[genes[[d_col]]])
link_col <- grDevices::adjustcolor(
    palette_colors(config$palettes$link)[[config$link_index]],
    alpha.f = config$link_alpha
)

links <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
if (isTRUE(config$draw_links) && nrow(links)) {
    require_columns(links, config$link_columns)
    lc <- config$link_columns
    key <- paste(genes[[s_col]], genes[[i_col]], sep = "\t")
    fk <- paste(links[[lc$from_sector]], links[[lc$from_item]], sep = "\t")
    tk <- paste(links[[lc$to_sector]], links[[lc$to_item]], sep = "\t")
    miss <- unique(c(fk[!fk %in% key], tk[!tk %in% key]))
    if (length(miss)) {
        stop(sprintf("Link endpoint(s) not in genes table: %s",
            paste(sub("\t", ":", miss), collapse = ", ")), call. = FALSE)
    }
    links$from_x <- genes$x[match(fk, key)]
    links$to_x <- genes$x[match(tk, key)]
    links$from_s <- links[[lc$from_sector]]
    links$to_s <- links[[lc$to_sector]]
    links <- links[links$from_s %in% sectors & links$to_s %in% sectors, ,
        drop = FALSE]
} else {
    links <- links[0, , drop = FALSE]
}

score_max <- if (is.null(config$score_max)) max(genes[[v_col]]) else
    config$score_max
breaks <- if (is.null(config$score_breaks)) {
    b <- pretty(c(0, score_max), n = 5)
    b[b <= score_max * 1.0001]
} else {
    config$score_breaks
}
score_top <- max(score_max, max(breaks)) * 1.05
grid_col <- "#d9d9d9"
axis_col <- "#4d4d4d"

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_plot <- function() {
    circlize::circos.par(
        start.degree = 90,
        gap.after = c(rep(config$gap_degree, length(sectors) - 1L),
            config$axis_gap_degree),
        track.margin = c(0.008, 0.008),
        cell.padding = c(0, 0, 0, 0),
        circle.margin = 0.02,
        points.overflow.warning = FALSE
    )
    sec_f <- factor(genes[[s_col]], levels = sectors)
    circlize::circos.initialize(
        sectors = sec_f,
        xlim = cbind(rep(0.5, length(sectors)), as.numeric(n_per) + 0.5)
    )

    # Track 1: sector header band
    circlize::circos.track(
        sectors = sec_f, ylim = c(0, 1),
        track.height = config$track_height$header,
        bg.col = config$header_fill, bg.border = config$header_border,
        panel.fun = function(x, y) {
            circlize::circos.text(
                circlize::CELL_META$xcenter, 0.5,
                circlize::CELL_META$sector.index,
                facing = "bending.inside", niceFacing = TRUE,
                cex = config$header_cex, font = 2, col = "grey15"
            )
        }
    )

    # Track 2: item labels coloured by direction
    circlize::circos.track(
        sectors = sec_f, ylim = c(0, 1),
        track.height = config$track_height$label,
        bg.col = NA, bg.border = NA,
        panel.fun = function(x, y) {
            sub <- genes[genes[[s_col]] == circlize::CELL_META$sector.index, ,
                drop = FALSE]
            circlize::circos.text(
                sub$x, rep(0, nrow(sub)), sub[[i_col]],
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = config$label_cex, col = sub$col
            )
        }
    )

    # Track 3: score bars
    circlize::circos.track(
        sectors = sec_f, ylim = c(0, score_top),
        track.height = config$track_height$bar,
        bg.col = NA, bg.border = grid_col,
        panel.fun = function(x, y) {
            xl <- circlize::CELL_META$xlim
            for (b in breaks[breaks > 0]) {
                circlize::circos.lines(xl, c(b, b), col = grid_col, lwd = 0.5)
            }
            sub <- genes[genes[[s_col]] == circlize::CELL_META$sector.index, ,
                drop = FALSE]
            circlize::circos.rect(
                sub$x - 0.38, 0, sub$x + 0.38, sub[[v_col]],
                col = sub$col, border = NA
            )
            if (identical(circlize::CELL_META$sector.index, sectors[[1]])) {
                circlize::circos.yaxis(
                    side = "left", at = breaks,
                    labels = format(breaks, trim = TRUE),
                    tick.length = circlize::convert_x(0.6, "mm"),
                    labels.cex = 0.48, col = axis_col, labels.col = axis_col
                )
            }
        }
    )

    if (nrow(links)) {
        for (i in seq_len(nrow(links))) {
            circlize::circos.link(
                links$from_s[[i]], links$from_x[[i]],
                links$to_s[[i]], links$to_x[[i]],
                col = link_col, lwd = 1, h.ratio = 0.6
            )
        }
    }
    circlize::circos.clear()

    graphics::legend(
        "bottomright",
        legend = unname(config$direction_labels[c("positive", "negative")]),
        fill = unname(dir_cols), border = NA, bty = "n",
        cex = 0.7, inset = c(0.01, 0.01), title = config$legend_title
    )
    if (!is.null(config$title)) {
        graphics::title(config$title, cex.main = 1.1, line = -1)
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
    grDevices::png(out_path, width = width, height = height, units = "in",
        res = 300, bg = "white")
} else if (identical(fmt, "pdf")) {
    grDevices::pdf(out_path, width = width, height = height)
} else {
    grDevices::svg(out_path, width = width, height = height, bg = "white")
}
graphics::par(mar = c(0.3, 0.3, 0.3, 0.3), bg = "white")
tryCatch(
    draw_plot(),
    finally = grDevices::dev.off()
)
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
