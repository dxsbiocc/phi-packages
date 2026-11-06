#!/usr/bin/env Rscript

# Template-ID: sankey-alluvial-map
#
# Purpose:
#   Draw a three-axis alluvial / bipartite association map: thin
#   coloured sigmoid links from proteins (left) through clinical
#   traits (middle) to disease categories (right). Line colour is
#   the disease category. Link width does not encode flow weight.
#
# Inputs:
#   One row per protein–trait association with its disease group.
#   Default example:
#     - protein, trait, category
#
# Output:
#   A PDF, PNG, or SVG alluvial association map.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggforce
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (column names, category order, colours, sizes).
#   Edit DATA PREPARATION to change protein / trait ordering.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one supplied association (e.g. colocalised + MR
#   support). This script does not run coloc, MR, GWAS, or enrichment.
#   Curves are equal-width presence links, not conserved flow mass.
#   Category colours and axis orders are CONFIG / supplied.

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
        source = "protein",
        middle = "trait",
        target = "category",
        fill = "category"
    ),
    # Top → bottom on the right axis. NULL = first-appearance order.
    target_order = c("ASCVD", "MASH", "Diabetes", "Obesity", "CKD"),
    # Optional explicit orders (character vectors). NULL = auto.
    middle_order = NULL,
    source_order = NULL,
    # Brand.Algolia / Qualitative.Safe hexes (do not invent colours).
    fill_colours = c(
        ASCVD = "#3369e7",
        MASH = "#cc6677",
        Diabetes = "#ffc168",
        Obesity = "#b84592",
        CKD = "#1cc7d0"
    ),
    palette = "Brand.Algolia",
    curve_linewidth = 0.42,
    curve_alpha = 0.78,
    # Horizontal pull of bezier control points (0–1 of axis gap).
    strength = 0.55,
    tick_length = 0.035,
    label_size = 2.55,
    label_size_target = 3.4,
    axis_x = c(0, 1, 2),
    labels = list(
        title = "",
        x = NULL,
        y = NULL
    ),
    size = list(width = 9.2, height = 10.5)
)

load_packages(c("ggplot2", "readr", "ggprism", "ggforce"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
src_col <- config$columns$source
mid_col <- config$columns$middle
tgt_col <- config$columns$target
fill_col <- config$columns$fill
require_columns(df, c(src_col, mid_col, tgt_col, fill_col))

df[[src_col]] <- as.character(df[[src_col]])
df[[mid_col]] <- as.character(df[[mid_col]])
df[[tgt_col]] <- as.character(df[[tgt_col]])
df[[fill_col]] <- as.character(df[[fill_col]])
ok <- nzchar(df[[src_col]]) & nzchar(df[[mid_col]]) & nzchar(df[[tgt_col]])
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with empty labels.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

# Drop exact duplicate associations.
df <- unique(df[c(src_col, mid_col, tgt_col, fill_col)])

tgt_order <- config$target_order
if (is.null(tgt_order) || !length(tgt_order)) {
    tgt_order <- unique(df[[tgt_col]])
} else {
    tgt_order <- as.character(tgt_order)
    miss <- setdiff(unique(df[[tgt_col]]), tgt_order)
    tgt_order <- c(tgt_order[tgt_order %in% unique(df[[tgt_col]])], miss)
}
df[[tgt_col]] <- factor(df[[tgt_col]], levels = tgt_order)

# Traits: group by their primary (modal) category, preserve file order
# within each block unless middle_order is supplied.
if (!is.null(config$middle_order) && length(config$middle_order)) {
    mid_order <- as.character(config$middle_order)
    mid_order <- c(mid_order, setdiff(unique(df[[mid_col]]), mid_order))
} else {
    mid_primary <- vapply(unique(df[[mid_col]]), function(tr) {
        tab <- table(df[[tgt_col]][df[[mid_col]] == tr])
        names(tab)[which.max(tab)]
    }, character(1))
    mid_df <- data.frame(
        trait = unique(df[[mid_col]]),
        cat = unname(mid_primary[unique(df[[mid_col]])]),
        idx = seq_along(unique(df[[mid_col]])),
        stringsAsFactors = FALSE
    )
    mid_df$cat <- factor(mid_df$cat, levels = tgt_order)
    mid_df <- mid_df[order(mid_df$cat, mid_df$idx), , drop = FALSE]
    mid_order <- mid_df$trait
}
df[[mid_col]] <- factor(df[[mid_col]], levels = mid_order)

# Proteins: barycentre of connected trait ranks (reduces crossings),
# unless source_order is supplied.
mid_rank <- setNames(seq_along(mid_order), mid_order)
if (!is.null(config$source_order) && length(config$source_order)) {
    src_order <- as.character(config$source_order)
    src_order <- c(src_order, setdiff(unique(df[[src_col]]), src_order))
} else {
    src_lv <- unique(df[[src_col]])
    src_score <- vapply(src_lv, function(p) {
        median(mid_rank[as.character(df[[mid_col]][df[[src_col]] == p])])
    }, numeric(1))
    src_order <- src_lv[order(src_score, src_lv)]
}
df[[src_col]] <- factor(df[[src_col]], levels = src_order)

fill_lv <- levels(df[[tgt_col]])
n_fill <- length(fill_lv)
if (!is.null(config$fill_colours) && length(config$fill_colours)) {
    fill_vals <- config$fill_colours
    if (is.null(names(fill_vals))) {
        if (length(fill_vals) < n_fill) {
            stop("fill_colours is shorter than the number of categories.", call. = FALSE)
        }
        fill_vals <- setNames(fill_vals[seq_len(n_fill)], fill_lv)
    } else {
        miss <- setdiff(fill_lv, names(fill_vals))
        if (length(miss)) {
            pal <- palette_colors(config$palette, n = max(length(miss), 3L))
            extra <- setNames(pal[seq_along(miss)], miss)
            fill_vals <- c(fill_vals, extra)
        }
        fill_vals <- fill_vals[fill_lv]
    }
} else {
    fill_vals <- setNames(
        palette_colors(config$palette, n = max(n_fill, 3L))[seq_len(n_fill)],
        fill_lv
    )
}

axis_x <- as.numeric(config$axis_x)
if (length(axis_x) != 3L || any(!is.finite(axis_x))) {
    axis_x <- c(0, 1, 2)
}
names(axis_x) <- c("source", "middle", "target")

# Equal vertical spacing within each axis (top = high y).
place_axis <- function(levels_chr) {
    n <- length(levels_chr)
    if (!n) {
        return(numeric(0))
    }
    y <- seq(1, 0, length.out = n)
    setNames(y, levels_chr)
}

y_src <- place_axis(levels(df[[src_col]]))
y_mid <- place_axis(levels(df[[mid_col]]))
y_tgt <- place_axis(levels(df[[tgt_col]]))

# Sigmoid bezier control: four points per link.
make_bezier <- function(x0, y0, x1, y1, group, colour_key, strength) {
    dx <- (x1 - x0) * strength
    data.frame(
        x = c(x0, x0 + dx, x1 - dx, x1),
        y = c(y0, y0, y1, y1),
        group = group,
        fill_key = colour_key,
        stringsAsFactors = FALSE
    )
}

strength <- as.numeric(config$strength)[[1]]
if (!is.finite(strength) || strength <= 0 || strength > 1) {
    strength <- 0.55
}

link_rows <- vector("list", nrow(df) * 2L)
for (i in seq_len(nrow(df))) {
    p <- as.character(df[[src_col]][[i]])
    t <- as.character(df[[mid_col]][[i]])
    c <- as.character(df[[tgt_col]][[i]])
    fk <- as.character(df[[fill_col]][[i]])
    g1 <- paste0("L", i, "a")
    g2 <- paste0("L", i, "b")
    link_rows[[2L * i - 1L]] <- make_bezier(
        axis_x[["source"]], y_src[[p]],
        axis_x[["middle"]], y_mid[[t]],
        g1, fk, strength
    )
    link_rows[[2L * i]] <- make_bezier(
        axis_x[["middle"]], y_mid[[t]],
        axis_x[["target"]], y_tgt[[c]],
        g2, fk, strength
    )
}
curves <- do.call(rbind, link_rows)
curves$fill_key <- factor(curves$fill_key, levels = fill_lv)

tick_len <- as.numeric(config$tick_length)[[1]]
if (!is.finite(tick_len) || tick_len <= 0) {
    tick_len <- 0.035
}

nodes <- rbind(
    data.frame(
        label = names(y_src),
        x = axis_x[["source"]],
        y = unname(y_src),
        side = "source",
        stringsAsFactors = FALSE
    ),
    data.frame(
        label = names(y_mid),
        x = axis_x[["middle"]],
        y = unname(y_mid),
        side = "middle",
        stringsAsFactors = FALSE
    ),
    data.frame(
        label = names(y_tgt),
        x = axis_x[["target"]],
        y = unname(y_tgt),
        side = "target",
        stringsAsFactors = FALSE
    )
)
nodes$x0 <- nodes$x - tick_len
nodes$x1 <- nodes$x + tick_len

lab_src <- nodes[nodes$side == "source", , drop = FALSE]
lab_mid <- nodes[nodes$side == "middle", , drop = FALSE]
lab_tgt <- nodes[nodes$side == "target", , drop = FALSE]

ink <- palette_colors("Brand.Emma")[[1]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_bezier(
        data = curves,
        aes(
            x = x,
            y = y,
            group = group,
            colour = fill_key
        ),
        linewidth = as.numeric(config$curve_linewidth),
        alpha = as.numeric(config$curve_alpha),
        show.legend = FALSE
    ) +
    geom_segment(
        data = nodes,
        aes(x = x0, xend = x1, y = y, yend = y),
        colour = ink,
        linewidth = 0.45,
        lineend = "butt"
    ) +
    geom_text(
        data = lab_src,
        aes(x = x0, y = y, label = label),
        hjust = 1,
        nudge_x = -0.02,
        size = as.numeric(config$label_size),
        colour = ink
    ) +
    geom_text(
        data = lab_mid,
        aes(x = x0, y = y, label = label),
        hjust = 1,
        nudge_x = -0.02,
        size = as.numeric(config$label_size),
        colour = ink
    ) +
    geom_text(
        data = lab_tgt,
        aes(x = x1, y = y, label = label),
        hjust = 0,
        nudge_x = 0.03,
        size = as.numeric(config$label_size_target),
        colour = ink,
        fontface = "bold"
    ) +
    scale_colour_manual(values = fill_vals, guide = "none") +
    coord_cartesian(
        xlim = c(axis_x[["source"]] - 0.55, axis_x[["target"]] + 0.55),
        ylim = c(-0.02, 1.02),
        clip = "off"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_void(base_size = 11) +
    theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 28, 10, 10)
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
