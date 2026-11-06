#!/usr/bin/env Rscript

# Template-ID: bar-radial-donut
#
# Purpose:
#   Radial stacked bars around an independent summary doughnut: one
#   equal-angle bar per item (e.g. cell type), grouped and coloured by
#   a supplied group (e.g. lineage) with a darker tip for the second
#   stack part, rotated item labels outside, a radial scale axis at
#   12 o'clock, and an inner doughnut of a second table (e.g. eGenes by
#   number of sharing cell types) with its total in the centre.
#
# Inputs:
#   bars.tsv - one row per item:
#     - cell_type: item label; "^{...}" marks a superscript
#     - lineage: group (contiguous colour run)
#     - shared, specific: non-negative stack parts (inner -> outer)
#   ring.tsv - one row per doughnut category, in display order:
#     - n_cell_types: category
#     - n_eGenes: non-negative value
#
# Output:
#   A PDF, PNG, or SVG radial bar + doughnut figure.
#
# Dependencies:
#   ggplot2, readr, grid
#
# Example:
#   Rscript plot.R bars.tsv ring.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, stack parts, group order, palettes, gaps,
#   radii, axis tick, ring labels, centre text, title). Extra items or
#   groups are extra rows, not extra ids. Do not run eQTL mapping,
#   colocalisation, or sharing tests in this script.
#
# Scientific assumptions:
#   Bar parts and doughnut values are supplied. Bar length is linear in
#   the stacked total; doughnut angle is proportional to value. The
#   centre number is the sum of the doughnut values.

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

io <- parse_io_args(n_input = 2L, input_names = c("bars", "ring"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    bar_columns = list(
        item = "cell_type",
        group = "lineage"
    ),
    # stack parts, inner -> outer; names are legend labels
    parts = c(
        "Cell Type-Shared" = "shared",
        "Cell Type-Specific" = "specific"
    ),
    ring_columns = list(
        category = "n_cell_types",
        value = "n_eGenes"
    ),
    groups = NULL,               # NULL = first appearance
    group_palette = "Qualitative.Safe",
    tip_darken = 0.35,           # darker shade for the outer stack part
    ring_palette = "Gradient.PuBu",
    ring_labels = c("1", "5", "10+"),
    ring_label_map = c("10+" = "\u226510"),
    ring_legend_title = "Number of shared Cell Types",
    centre_title = "eGene",
    centre_arrow = TRUE,
    title = 'italic("cis")*"-eGene Counts"',   # plotmath; NULL = none
    gap_slots = 1,               # empty slots between groups
    axis_gap_slots = 3,          # empty slots at 12 o'clock for the axis
    axis_tick = NULL,            # NULL = largest pretty value <= max
    radius = list(
        ring_in = 0.30,
        ring_out = 0.56,
        base = 0.62,
        bar_span = 0.66
    ),
    bar_width = 0.82,            # fraction of a slot
    label_size = 2.6,
    size = list(
        width = 9,
        height = 9
    )
)

load_packages(c("ggplot2", "readr", "grid"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
bars <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE,
    check.names = FALSE)
ring <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE,
    check.names = FALSE)
require_columns(bars, c(config$bar_columns, as.list(unname(config$parts))))
require_columns(ring, config$ring_columns)

it_col <- config$bar_columns$item
g_col <- config$bar_columns$group
parts <- config$parts
bars[[it_col]] <- as.character(bars[[it_col]])
bars[[g_col]] <- as.character(bars[[g_col]])
for (pc in parts) {
    bars[[pc]] <- as.numeric(bars[[pc]])
    if (any(!is.finite(bars[[pc]]) | bars[[pc]] < 0)) {
        stop(sprintf("Column %s must be finite and non-negative.", pc),
            call. = FALSE)
    }
}
groups <- if (is.null(config$groups)) unique(bars[[g_col]]) else
    as.character(config$groups)
bars <- bars[bars[[g_col]] %in% groups, , drop = FALSE]
bars <- bars[order(match(bars[[g_col]], groups), seq_len(nrow(bars))), ,
    drop = FALSE]
groups <- groups[groups %in% bars[[g_col]]]
bars$total <- rowSums(bars[, unname(parts), drop = FALSE])

rc <- config$ring_columns
ring[[rc$category]] <- as.character(ring[[rc$category]])
ring[[rc$value]] <- as.numeric(ring[[rc$value]])
if (any(!is.finite(ring[[rc$value]]) | ring[[rc$value]] < 0)) {
    stop("Ring values must be finite and non-negative.", call. = FALSE)
}

# slots: axis gap split around 12 o'clock, gap_slots between groups
n_grp_gap <- (length(groups) - 1L) * config$gap_slots
n_slots <- nrow(bars) + n_grp_gap + config$axis_gap_slots
w <- 360 / n_slots
grp_idx <- match(bars[[g_col]], groups)
slot <- config$axis_gap_slots / 2 + (seq_len(nrow(bars)) - 1) +
    (grp_idx - 1) * config$gap_slots
bars$angle <- (slot + 0.5) * w

r <- config$radius
max_total <- max(bars$total)
to_r <- function(v) r$base + v / max_total * r$bar_span
tick <- if (is.null(config$axis_tick)) {
    pt <- pretty(c(0, max_total))
    max(pt[pt <= max_total])
} else {
    config$axis_tick
}

grp_col <- setNames(palette_colors(config$group_palette, n = length(groups)),
    groups)
darken <- function(col, p) {
    m <- grDevices::col2rgb(col) / 255 * (1 - p)
    grDevices::rgb(m[1, ], m[2, ], m[3, ])
}
part_shade <- function(base, k) {
    if (k == 1L) base else darken(base, config$tip_darken * (k - 1) /
        max(1, length(parts) - 1))
}
stack <- do.call(rbind, lapply(seq_along(parts), function(k) {
    lo <- if (k == 1L) 0 else rowSums(bars[, unname(parts)[seq_len(k - 1)],
        drop = FALSE])
    hi <- lo + bars[[parts[[k]]]]
    data.frame(
        xmin = bars$angle - w * config$bar_width / 2,
        xmax = bars$angle + w * config$bar_width / 2,
        ymin = to_r(lo), ymax = to_r(hi),
        fill = vapply(bars[[g_col]], function(g) part_shade(grp_col[[g]], k),
            character(1)),
        stringsAsFactors = FALSE
    )
}))
stack <- stack[stack$ymax > stack$ymin, , drop = FALSE]

to_plotmath <- function(s) {
    q <- function(x) paste0('"', gsub('(["\\\\])', "\\\\\\1", x), '"')
    m <- gregexpr("\\^\\{[^}]*\\}", s)[[1]]
    if (m[[1]] == -1) return(q(s))
    sups <- regmatches(s, list(m))[[1]]
    plain <- regmatches(s, list(m), invert = TRUE)[[1]]
    out <- q(plain[[1]])
    for (i in seq_along(sups)) {
        out <- paste0(out, "^", q(substr(sups[[i]], 3, nchar(sups[[i]]) - 1)))
        if (nzchar(plain[[i + 1]])) out <- paste0(out, "*", q(plain[[i + 1]]))
    }
    out
}
bars$lab <- vapply(bars[[it_col]], to_plotmath, character(1))
right <- bars$angle <= 180
bars$lab_angle <- ifelse(right, 90 - bars$angle, 270 - bars$angle)
bars$hjust <- ifelse(right, 0, 1)
bars$lab_r <- to_r(bars$total) + 0.02
bars$lab_col <- unname(grp_col[bars[[g_col]]])

ring_stops <- palette_colors(config$ring_palette)
ring_stops <- rev(ring_stops[-seq_len(max(1L, floor(length(ring_stops) / 4)))])
ring$fill <- grDevices::colorRampPalette(ring_stops)(nrow(ring))
tot <- sum(ring[[rc$value]])
ring$width <- 360 * ring[[rc$value]] / tot
ring$end <- cumsum(ring$width)
ring$start <- ring$end - ring$width
ring$start <- ring$start + 0.6
ring$end <- ring$end - 0.6
ring$mid <- (ring$start + ring$end) / 2
ring$tag <- ifelse(ring[[rc$category]] %in% names(config$ring_label_map),
    unname(config$ring_label_map[ring[[rc$category]]]), ring[[rc$category]])
ring_lab <- ring[ring[[rc$category]] %in% config$ring_labels, , drop = FALSE]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
y_top <- r$base + r$bar_span + 0.03
p <- ggplot() +
    geom_rect(data = ring, aes(xmin = start, xmax = end, ymin = r$ring_in,
        ymax = r$ring_out, fill = fill), colour = "white", linewidth = 0.45) +
    geom_text(data = ring_lab, aes(x = mid, y = (r$ring_in + r$ring_out) / 2,
        label = tag), colour = "white", size = 5) +
    geom_rect(data = stack, aes(xmin = xmin, xmax = xmax, ymin = ymin,
        ymax = ymax, fill = fill), colour = "white", linewidth = 0.3) +
    annotate("segment", x = 0, xend = 0, y = r$base, yend = to_r(tick),
        colour = "grey25", linewidth = 0.45) +
    annotate("segment", x = 0, xend = w * 0.6, y = to_r(tick),
        yend = to_r(tick), colour = "grey25", linewidth = 0.45) +
    annotate("text", x = 0, y = to_r(tick), label = format(tick, big.mark = ""),
        hjust = 1.15, vjust = 0.3, size = 4, colour = "grey10") +
    geom_text(data = bars, aes(x = angle, y = lab_r, label = lab,
        angle = lab_angle, hjust = hjust, colour = lab_col),
        parse = TRUE, size = config$label_size) +
    annotate("text", x = 0, y = 0, label = paste0(config$centre_title, "\n",
        format(tot, big.mark = "")), size = 6, colour = "grey15",
        lineheight = 0.95) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0, 360), expand = c(0, 0), breaks = NULL) +
    scale_y_continuous(limits = c(0, y_top), expand = c(0, 0), breaks = NULL) +
    coord_polar(theta = "x", start = 0, direction = 1, clip = "off") +
    theme_void() +
    theme(
        plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(hjust = 0.5, size = 16,
            margin = margin(b = 4)),
        plot.margin = margin(14, 14, 14, 14, unit = "pt")
    )
if (isTRUE(config$centre_arrow)) {
    arc <- data.frame(x = seq(120, 230, length.out = 50),
        y = r$ring_in * 0.82)
    p <- p + geom_path(data = arc, aes(x = x, y = y), colour = "grey60",
        linewidth = 0.5,
        arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "open"))
}
if (!is.null(config$title)) {
    p <- p + labs(title = parse(text = config$title)[[1]])
}

draw_legend <- function() {
    grid::pushViewport(grid::viewport(x = 0.985, y = 0.02, width = 0.25,
        height = 0.17, just = c("right", "bottom")))
    n <- nrow(ring)
    gp_t <- grid::gpar(fontsize = 9, col = "grey10")
    grid::grid.text(config$ring_legend_title, x = 0, y = 0.95, just = c(0, 1),
        gp = gp_t)
    bw <- 1 / n
    grid::grid.rect(x = (seq_len(n) - 0.5) * bw, y = 0.72, width = bw * 0.94,
        height = 0.11, gp = grid::gpar(fill = ring$fill, col = NA))
    grid::grid.text(ring$tag[[1]], x = 0.01, y = 0.55, just = c(0, 0.5),
        gp = gp_t)
    grid::grid.text(ring$tag[[n]], x = 0.99, y = 0.55, just = c(1, 0.5),
        gp = gp_t)
    grid::grid.lines(x = c(0.08, 0.84), y = 0.55, gp = grid::gpar(col = "grey10"),
        arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "open"))
    k <- length(parts)
    shades <- vapply(seq_len(k), function(i) part_shade("#9a9a9a", i),
        character(1))
    for (i in seq_len(k)) {
        y <- 0.32 - (i - 1) * 0.16
        grid::grid.rect(x = 0.02, y = y, width = 0.08, height = 0.1,
            just = c(0, 0.5), gp = grid::gpar(fill = shades[[i]], col = NA))
        grid::grid.text(names(parts)[[i]], x = 0.13, y = y, just = c(0, 0.5),
            gp = grid::gpar(fontsize = 8.5, col = "grey35"))
    }
    grid::grid.text("(inner \u2192 outer)", x = 0.13,
        y = 0.32 - k * 0.16 + 0.04, just = c(0, 0.5),
        gp = grid::gpar(fontsize = 7.5, col = "grey55"))
    grid::popViewport()
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
wd <- config$size$width
ht <- config$size$height
if (identical(fmt, "png")) {
    grDevices::png(out_path, width = wd, height = ht, units = "in", res = 300,
        bg = "white")
} else if (identical(fmt, "pdf")) {
    grDevices::cairo_pdf(out_path, width = wd, height = ht)
} else {
    grDevices::svg(out_path, width = wd, height = ht, bg = "white")
}
tryCatch({
    print(p)
    draw_legend()
}, finally = grDevices::dev.off())
message("Saved: ", normalizePath(out_path, mustWork = TRUE))
