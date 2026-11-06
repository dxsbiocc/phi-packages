#!/usr/bin/env Rscript

# Template-ID: sunburst-petal
#
# Purpose:
#   Draw a rounded "petal" sunburst of a two- or three-level taxonomy.
#   Each leaf is one tapered rounded capsule (equal angle) on the outer
#   ring; each category is a rounded wedge spanning its leaves; an
#   optional super-group is a thick arc with a label in the centre.
#
# Inputs:
#   A long table (tsv/csv), one row per leaf. Default example:
#     - group: optional super-group drawn as a centre arc (e.g. Gen / Edit)
#     - category: inner wedge
#     - category_label: optional wedge text; "\n" breaks lines
#     - colour: optional hex per category
#     - subcategory: outer petal text
#
# Output:
#   A PDF, PNG, or SVG petal sunburst.
#
# Dependencies:
#   ggplot2, readr, ggforce
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, radii, gaps, corner radius, text sizes,
#   group arc colours). Row order is the clockwise drawing order.
#
# Scientific assumptions:
#   Every leaf has equal angle, so wedge angle encodes the number of
#   leaves in a category, not a sample count. Counts such as "(13k)"
#   belong in category_label as text. The script does not count,
#   classify, or reorder the taxonomy.

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
        group = "group",
        category = "category",
        category_label = "category_label",
        colour = "colour",
        leaf = "subcategory"
    ),
    palettes = list(
        category = "Qualitative.Bold"
    ),
    group_colours = c(
        Gen = "#C03A3A",
        Edit = "#4A81BA"
    ),
    radius = list(
        group_arc = 0.20,
        group_label = 0.135,
        wedge_inner = 0.27,
        wedge_outer = 0.50,
        petal_inner = 0.54,
        petal_outer = 1.00
    ),
    start_degree = 0,          # 0 = 12 o'clock, clockwise
    wedge_gap_degree = 1.4,
    petal_gap = 0.16,          # fraction of each petal's angle left empty
    group_gap_degree = 7,
    wedge_corner = ggplot2::unit(3.2, "mm"),
    petal_fill_mix = 0.72,     # mix category colour toward white
    petal_edge_mix = 0.38,
    group_arc_width = 3.6,
    text = list(
        leaf_size = 2.55,
        category_size = 3.0,
        group_size = 5.2,
        dark = "#1F1F1F",
        light = "#FFFFFF"
    )
)

load_packages(c("ggplot2", "readr", "ggforce"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns[c("category", "leaf")])
cat_col <- config$columns$category
leaf_col <- config$columns$leaf
grp_col <- config$columns$group
lab_col <- config$columns$category_label
col_col <- config$columns$colour

df[[cat_col]] <- as.character(df[[cat_col]])
df[[leaf_col]] <- as.character(df[[leaf_col]])
ok <- !is.na(df[[cat_col]]) & nzchar(df[[cat_col]]) &
    !is.na(df[[leaf_col]]) & nzchar(df[[leaf_col]])
df <- df[ok, , drop = FALSE]
if (!nrow(df)) {
    stop("No category/leaf rows to plot.", call. = FALSE)
}
cat_levels <- unique(df[[cat_col]])
if (any(duplicated(rle(df[[cat_col]])$values))) {
    stop("Rows of one category must be contiguous (row order = drawing order).",
        call. = FALSE)
}

has_group <- !is.null(grp_col) && grp_col %in% names(df)
if (has_group) {
    df[[grp_col]] <- as.character(df[[grp_col]])
    df[[grp_col]][is.na(df[[grp_col]])] <- ""
}

# Category display labels and colours.
cat_label <- setNames(cat_levels, cat_levels)
if (!is.null(lab_col) && lab_col %in% names(df)) {
    first <- df[!duplicated(df[[cat_col]]), , drop = FALSE]
    lab <- as.character(first[[lab_col]])
    lab <- gsub("\\\\n", "\n", lab)
    use <- !is.na(lab) & nzchar(lab)
    cat_label[first[[cat_col]][use]] <- lab[use]
}
cat_fill <- palette_colors(config$palettes$category, n = length(cat_levels))
names(cat_fill) <- cat_levels
if (!is.null(col_col) && col_col %in% names(df)) {
    first <- df[!duplicated(df[[cat_col]]), , drop = FALSE]
    hex <- as.character(first[[col_col]])
    use <- !is.na(hex) & grepl("^#[0-9A-Fa-f]{6}", hex)
    cat_fill[first[[cat_col]][use]] <- hex[use]
}

mix_white <- function(hex, w) {
    rgb <- grDevices::col2rgb(hex) / 255
    out <- rgb * (1 - w) + w
    grDevices::rgb(out[1, ], out[2, ], out[3, ])
}
contrast_text <- function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
    ifelse(lum > 0.6, config$text$dark, config$text$light)
}

# Angles: theta measured clockwise from 12 o'clock, in radians.
n_leaf <- nrow(df)
step <- 2 * pi / n_leaf
start <- config$start_degree * pi / 180
df$theta0 <- start + (seq_len(n_leaf) - 1) * step
df$theta1 <- df$theta0 + step
df$theta <- (df$theta0 + df$theta1) / 2
to_xy <- function(r, th) cbind(x = r * sin(th), y = r * cos(th))

# Readable rotation for text whose baseline is radial / tangential.
radial_angle <- function(th) {
    a <- 90 - th * 180 / pi
    a <- ((a + 180) %% 360) - 180
    ifelse(a > 90 | a < -90, a + 180, a)
}
tangent_angle <- function(th) {
    a <- -th * 180 / pi
    a <- ((a + 180) %% 360) - 180
    ifelse(a > 90 | a < -90, a + 180, a)
}

rad <- config$radius

# Outer petals: tapered capsule = convex hull of an inner and outer circle.
petal_poly <- function(th, id) {
    half <- step / 2 * (1 - config$petal_gap)
    wi <- rad$petal_inner * sin(half)
    wo <- rad$petal_outer * sin(half)
    ci <- rad$petal_inner + wi
    co <- rad$petal_outer - wo
    a <- seq(0, 2 * pi, length.out = 49)[-49]
    u <- c(ci + wi * cos(a), co + wo * cos(a))
    v <- c(wi * sin(a), wo * sin(a))
    h <- grDevices::chull(u, v)
    u <- u[h]
    v <- v[h]
    data.frame(
        x = u * sin(th) + v * cos(th),
        y = u * cos(th) - v * sin(th),
        id = id,
        stringsAsFactors = FALSE
    )
}
petals <- do.call(rbind, lapply(seq_len(n_leaf), function(i) {
    petal_poly(df$theta[[i]], i)
}))
petals$category <- df[[cat_col]][petals$id]

# Inner wedges: annular sectors, corners rounded by ggforce::geom_shape.
gap_w <- config$wedge_gap_degree * pi / 180 / 2
wedge_tbl <- do.call(rbind, lapply(cat_levels, function(cg) {
    idx <- which(df[[cat_col]] == cg)
    data.frame(
        category = cg,
        theta0 = min(df$theta0[idx]) + gap_w,
        theta1 = max(df$theta1[idx]) - gap_w,
        stringsAsFactors = FALSE
    )
}))
wedges <- do.call(rbind, lapply(seq_len(nrow(wedge_tbl)), function(i) {
    th <- seq(wedge_tbl$theta0[[i]], wedge_tbl$theta1[[i]], length.out = 60)
    outer <- to_xy(rad$wedge_outer, th)
    inner <- to_xy(rad$wedge_inner, rev(th))
    data.frame(
        x = c(outer[, "x"], inner[, "x"]),
        y = c(outer[, "y"], inner[, "y"]),
        category = wedge_tbl$category[[i]],
        stringsAsFactors = FALSE
    )
}))
wedge_tbl$theta <- (wedge_tbl$theta0 + wedge_tbl$theta1) / 2
wedge_xy <- to_xy((rad$wedge_inner + rad$wedge_outer) / 2, wedge_tbl$theta)
wedge_tbl$x <- wedge_xy[, "x"]
wedge_tbl$y <- wedge_xy[, "y"]
wedge_tbl$angle <- tangent_angle(wedge_tbl$theta)
wedge_tbl$label <- unname(cat_label[wedge_tbl$category])
wedge_tbl$text_col <- contrast_text(cat_fill[wedge_tbl$category])

leaf_xy <- to_xy((rad$petal_inner + rad$petal_outer) / 2, df$theta)
leaf_tbl <- data.frame(
    x = leaf_xy[, "x"],
    y = leaf_xy[, "y"],
    angle = radial_angle(df$theta),
    label = df[[leaf_col]],
    stringsAsFactors = FALSE
)

# Centre super-group arcs.
arcs <- NULL
arc_lab <- NULL
if (has_group) {
    grp_levels <- unique(df[[grp_col]][nzchar(df[[grp_col]])])
    gap_g <- config$group_gap_degree * pi / 180 / 2
    gcols <- config$group_colours
    miss <- setdiff(grp_levels, names(gcols))
    if (length(miss)) {
        extra <- palette_colors("Qualitative.Safe", n = length(miss))
        names(extra) <- miss
        gcols <- c(gcols, extra)
    }
    arcs <- do.call(rbind, lapply(grp_levels, function(g) {
        idx <- which(df[[grp_col]] == g)
        th <- seq(min(df$theta0[idx]) + gap_g, max(df$theta1[idx]) - gap_g,
            length.out = 120)
        xy <- to_xy(rad$group_arc, th)
        data.frame(x = xy[, "x"], y = xy[, "y"], group = g, stringsAsFactors = FALSE)
    }))
    arc_lab <- do.call(rbind, lapply(grp_levels, function(g) {
        idx <- which(df[[grp_col]] == g)
        th <- (min(df$theta0[idx]) + max(df$theta1[idx])) / 2
        xy <- to_xy(rad$group_label, th)
        data.frame(
            x = xy[, "x"], y = xy[, "y"], group = g,
            angle = tangent_angle(th) + ifelse(cos(th) >= 0, 90, -90),
            stringsAsFactors = FALSE
        )
    }))
    arc_lab$angle <- ((arc_lab$angle + 180) %% 360) - 180
    arc_lab$angle <- ifelse(arc_lab$angle > 90 | arc_lab$angle < -90,
        arc_lab$angle + 180, arc_lab$angle)
}

petals$fill <- mix_white(cat_fill[petals$category], config$petal_fill_mix)
petals$edge <- mix_white(cat_fill[petals$category], config$petal_edge_mix)
wedges$fill <- unname(cat_fill[wedges$category])
if (!is.null(arcs)) {
    arcs$col <- unname(gcols[arcs$group])
    arc_lab$col <- unname(gcols[arc_lab$group])
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
        data = petals,
        ggplot2::aes(x = x, y = y, group = id, fill = fill, colour = edge),
        linewidth = 0.9,
        linejoin = "round"
    ) +
    ggforce::geom_shape(
        data = wedges,
        ggplot2::aes(x = x, y = y, group = category, fill = fill),
        colour = "grey25",
        linewidth = 0.35,
        radius = config$wedge_corner
    ) +
    ggplot2::geom_text(
        data = leaf_tbl,
        ggplot2::aes(x = x, y = y, label = label, angle = angle),
        size = config$text$leaf_size,
        fontface = "bold",
        colour = config$text$dark
    ) +
    ggplot2::geom_text(
        data = wedge_tbl,
        ggplot2::aes(x = x, y = y, label = label, angle = angle, colour = text_col),
        size = config$text$category_size,
        fontface = "bold",
        lineheight = 0.85
    )

if (!is.null(arcs)) {
    p <- p +
        ggplot2::geom_path(
            data = arcs,
            ggplot2::aes(x = x, y = y, group = group, colour = col),
            linewidth = config$group_arc_width,
            lineend = "round"
        ) +
        ggplot2::geom_text(
            data = arc_lab,
            ggplot2::aes(x = x, y = y, label = group, angle = angle, colour = col),
            size = config$text$group_size,
            fontface = "bold"
        )
}

p <- p +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::coord_fixed(
        xlim = c(-1.04, 1.04), ylim = c(-1.04, 1.04), expand = FALSE, clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 11, height = 11)
