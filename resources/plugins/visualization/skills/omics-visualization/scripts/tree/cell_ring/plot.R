#!/usr/bin/env Rscript

# Template-ID: tree-cell-ring
#
# Purpose:
#   Draw an equal-angle circular cell tree: lineage branches from the
#   centre, L3 radial capsules, numbered tips, and concentric diverging
#   tracks (for example Sex/Age by ATAC/RNA). Leaf angle is row order,
#   not descendant size.
#
# Inputs:
#   One row per cell, already in circular order:
#     id, l1, l2, l3, cell_type, percent
#     plus one score column per track, and optional star columns.
#   Track definitions live in CONFIG$tracks, inner ring first.
#
# Output:
#   A PDF, PNG, or SVG circular annotation tree.
#
# Dependencies:
#   ggplot2, readr, ggforce, cowplot, scales
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG for column names, track list, and radii.
#   Edit PLOT only when the geometry must change.
#
# Scientific assumptions:
#   Hierarchy, proportions, scores, and stars are supplied.
#   This script does not cluster cells or estimate proportions.
#   Score colours use Diverging.RdBu on [-1, 1]. Values outside that
#   range are clamped for display. ATAC-style tracks use a lower alpha.
#   Lineage colours follow first-seen l1 order in Qualitative.Bold.

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
        id = "id",
        l1 = "l1",
        l2 = "l2",
        l3 = "l3",
        cell = "cell_type",
        percent = "percent"
    ),
    # Inner ring first. alpha < 1 marks the lighter assay (ATAC).
    tracks = list(
        list(score = "sex_score", stars = "sex_atac_stars", alpha = 0.4, label = "Sex: ATAC / RNA"),
        list(score = "sex_score", stars = "sex_rna_stars", alpha = 1, label = ""),
        list(score = "age_score", stars = "age_atac_stars", alpha = 0.4, label = "Age: ATAC / RNA"),
        list(score = "age_score", stars = "age_rna_stars", alpha = 1, label = "")
    ),
    layout = list(
        l1 = 0.20,
        l2 = 0.42,
        l3 = 0.78,
        tip = 1.00,
        tip_width = 0.045,
        ring_width = 0.055,
        ring_gap = 0.012,
        ring_start = 1.07,
        label_radius = 1.38
    ),
    labels = list(
        legend = "Cell annotation"
    ),
    palettes = list(
        lineage = "Qualitative.Bold",
        score = "Diverging.RdBu"
    )
)

load_packages(c("ggplot2", "readr", "ggforce", "cowplot", "scales"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
id_col <- config$columns$id
l1_col <- config$columns$l1
l2_col <- config$columns$l2
l3_col <- config$columns$l3
cell_col <- config$columns$cell
pct_col <- config$columns$percent

df[[id_col]] <- as.integer(df[[id_col]])
df[[l1_col]] <- as.character(df[[l1_col]])
df[[l2_col]] <- as.character(df[[l2_col]])
df[[l3_col]] <- as.character(df[[l3_col]])
df[[cell_col]] <- as.character(df[[cell_col]])
df[[pct_col]] <- as.numeric(df[[pct_col]])
if (anyNA(df[[id_col]]) || anyDuplicated(df[[id_col]])) {
    stop("id must be unique integers.", call. = FALSE)
}
if (any(!is.finite(df[[pct_col]]) | df[[pct_col]] < 0)) {
    stop("percent must be finite and non-negative.", call. = FALSE)
}
df <- df[order(df[[id_col]]), , drop = FALSE]
n_cell <- nrow(df)
if (n_cell < 2L) {
    stop("Need at least two cells.", call. = FALSE)
}

df$theta <- (seq_len(n_cell) - 0.5) * 2 * pi / n_cell
half <- pi / n_cell
pad <- 0.004
df$start <- df$theta - half + pad
df$end <- df$theta + half - pad

l1_levels <- unique(df[[l1_col]])
lineage_pal <- palette_colors(config$palettes$lineage, n = length(l1_levels))
names(lineage_pal) <- l1_levels
df$lineage_colour <- unname(lineage_pal[df[[l1_col]]])

score_pal <- rev(palette_colors(config$palettes$score))
score_ramp <- grDevices::colorRampPalette(score_pal)(101)

run_id <- function(keys) {
    cumsum(c(1L, as.integer(keys[-1] != keys[-length(keys)])))
}
df$run_l1 <- run_id(df[[l1_col]])
df$run_l2 <- run_id(paste(df[[l1_col]], df[[l2_col]], sep = "\r"))
df$run_l3 <- run_id(paste(df[[l1_col]], df[[l2_col]], df[[l3_col]], sep = "\r"))

lay <- config$layout
radii <- c(l1 = lay$l1, l2 = lay$l2, l3 = lay$l3, cell = lay$tip)

node_at <- function(run, radius) {
    split_df <- split(df, run)
    do.call(rbind, lapply(split_df, function(sub) {
        data.frame(
            run = sub[[names(run)[1]]][[1]],
            theta = mean(sub$theta),
            radius = radius,
            label = sub[[ncol(sub)]][[1]],
            colour = sub$lineage_colour[[1]],
            percent = sum(sub[[pct_col]]),
            stringsAsFactors = FALSE
        )
    }))
}

# Rebuild nodes from explicit columns so labels are not accidental.
summarise_run <- function(run_values, label_values, radius) {
    pieces <- split(seq_len(n_cell), run_values)
    do.call(rbind, lapply(pieces, function(idx) {
        data.frame(
            theta = mean(df$theta[idx]),
            radius = radius,
            label = label_values[idx][[1]],
            colour = df$lineage_colour[idx][[1]],
            percent = sum(df[[pct_col]][idx]),
            stringsAsFactors = FALSE
        )
    }))
}
nodes_l1 <- summarise_run(df$run_l1, df[[l1_col]], radii[["l1"]])
nodes_l2 <- summarise_run(df$run_l2, df[[l2_col]], radii[["l2"]])
nodes_l3 <- summarise_run(df$run_l3, df[[l3_col]], radii[["l3"]])
nodes_l1$x <- nodes_l1$radius * sin(nodes_l1$theta)
nodes_l1$y <- nodes_l1$radius * cos(nodes_l1$theta)
nodes_l2$x <- nodes_l2$radius * sin(nodes_l2$theta)
nodes_l2$y <- nodes_l2$radius * cos(nodes_l2$theta)
nodes_l3$x <- nodes_l3$radius * sin(nodes_l3$theta)
nodes_l3$y <- nodes_l3$radius * cos(nodes_l3$theta)

# Map each leaf to its ancestor node row.
parent_row <- function(run_values, nodes) {
    match(run_values, as.integer(names(split(seq_len(n_cell), run_values))))
}
# names(split) is the run id as character in appearance order, which matches
# the row order of summarise_run. Match run value to that order.
run_index <- function(run_values) {
    keys <- unique(run_values)
    match(run_values, keys)
}
df$i_l1 <- run_index(df$run_l1)
df$i_l2 <- run_index(df$run_l2)
df$i_l3 <- run_index(df$run_l3)

polar_edge <- function(th0, r0, th1, r1, id, colour, n = 28L) {
    t <- seq(0, 1, length.out = n)
    th <- th0 + (th1 - th0) * t
    r <- r0 + (r1 - r0) * t
    data.frame(
        edge = id,
        x = r * sin(th),
        y = r * cos(th),
        colour = colour,
        stringsAsFactors = FALSE
    )
}
edge_rows <- list()
k <- 1L
for (i in seq_len(nrow(nodes_l1))) {
    edge_rows[[k]] <- polar_edge(0, 0, nodes_l1$theta[i], nodes_l1$radius[i], k, nodes_l1$colour[i])
    k <- k + 1L
}
for (i in seq_len(nrow(nodes_l2))) {
    parent <- nodes_l1[df$i_l1[df$i_l2 == i][[1]], ]
    edge_rows[[k]] <- polar_edge(parent$theta, parent$radius, nodes_l2$theta[i], nodes_l2$radius[i], k, nodes_l2$colour[i])
    k <- k + 1L
}
edges <- do.call(rbind, edge_rows)

upright <- function(angle) {
    ifelse(angle > 90, angle - 180, ifelse(angle < -90, angle + 180, angle))
}
label_on <- function(th0, r0, th1, r1, text, colour) {
    th <- th0 + (th1 - th0) * 0.55
    r <- r0 + (r1 - r0) * 0.55
    ang <- upright(atan2(cos(th1) * r1 - cos(th0) * r0, sin(th1) * r1 - sin(th0) * r0) * 180 / pi)
    data.frame(
        x = r * sin(th), y = r * cos(th), label = text, angle = ang,
        stringsAsFactors = FALSE
    )
}
branch_labels <- rbind(
    do.call(rbind, lapply(seq_len(nrow(nodes_l1)), function(i) {
        label_on(0, 0, nodes_l1$theta[i], nodes_l1$radius[i], nodes_l1$label[i], nodes_l1$colour[i])
    })),
    do.call(rbind, lapply(seq_len(nrow(nodes_l2)), function(i) {
        parent <- nodes_l1[df$i_l1[df$i_l2 == i][[1]], ]
        text <- nodes_l2$label[i]
        if (identical(text, nodes_l3$label[df$i_l3[df$i_l2 == i][[1]]])) {
            text <- ""
        }
        label_on(parent$theta, parent$radius, nodes_l2$theta[i], nodes_l2$radius[i], text, nodes_l2$colour[i])
    }))
)
branch_labels <- branch_labels[nzchar(branch_labels$label), , drop = FALSE]

capsules <- do.call(rbind, lapply(seq_len(nrow(nodes_l3)), function(i) {
    parent <- nodes_l2[df$i_l2[df$i_l3 == i][[1]], ]
    x0 <- parent$radius * sin(parent$theta)
    y0 <- parent$radius * cos(parent$theta)
    x1 <- nodes_l3$radius[i] * sin(nodes_l3$theta[i])
    y1 <- nodes_l3$radius[i] * cos(nodes_l3$theta[i])
    ang <- upright(atan2(y1 - y0, x1 - x0) * 180 / pi)
    data.frame(
        x = x0, y = y0, xend = x1, yend = y1,
        colour = nodes_l3$colour[i],
        lx = x0 + 0.55 * (x1 - x0),
        ly = y0 + 0.55 * (y1 - y0),
        angle = ang,
        label = sprintf("%s  %.1f%%", nodes_l3$label[i], nodes_l3$percent[i]),
        stringsAsFactors = FALSE
    )
}))

sectors <- do.call(rbind, lapply(split(df, df$run_l1), function(sub) {
    data.frame(
        start = min(sub$start) - 0.01,
        end = max(sub$end) + 0.01,
        r0 = 0.03,
        r = lay$l2,
        fill = sub$lineage_colour[[1]],
        alpha = 0.16,
        stringsAsFactors = FALSE
    )
}))
l3_bands <- do.call(rbind, lapply(split(df, df$run_l3), function(sub) {
    data.frame(
        start = min(sub$theta) - half - 0.006,
        end = max(sub$theta) + half + 0.006,
        r0 = lay$l2,
        r = lay$l3 + 0.03,
        fill = sub$lineage_colour[[1]],
        alpha = 0.22,
        stringsAsFactors = FALSE
    )
}))
bands <- rbind(sectors, l3_bands)

df$tip_x <- sin(df$theta)
df$tip_y <- cos(df$theta)
df$label_angle <- 90 - df$theta * 180 / pi
df$hjust <- ifelse(df$tip_x < 0, 1, 0)
df$label_angle <- ifelse(df$tip_x < 0, df$label_angle + 180, df$label_angle)

track_defs <- config$tracks
if (!length(track_defs)) {
    stop("CONFIG$tracks needs at least one track.", call. = FALSE)
}
ring_rows <- list()
for (i in seq_along(track_defs)) {
    spec <- track_defs[[i]]
    if (!spec$score %in% names(df)) {
        stop(paste("Missing score column:", spec$score), call. = FALSE)
    }
    score <- as.numeric(df[[spec$score]])
    if (any(!is.finite(score))) {
        stop(paste("Score column has non-finite values:", spec$score), call. = FALSE)
    }
    score <- pmax(pmin(score, 1), -1)
    idx <- pmax(1L, pmin(101L, as.integer(round((score + 1) / 2 * 100)) + 1L))
    inner <- lay$ring_start + (i - 1) * (lay$ring_width + lay$ring_gap)
    stars <- rep(NA_character_, n_cell)
    if (!is.null(spec$stars) && nzchar(spec$stars) && spec$stars %in% names(df)) {
        stars <- as.character(df[[spec$stars]])
        stars[!nzchar(stars) | stars == "NA"] <- NA_character_
    }
    ring_rows[[i]] <- data.frame(
        start = df$start,
        end = df$end,
        r0 = inner,
        r = inner + lay$ring_width,
        fill = score_ramp[idx],
        alpha = spec$alpha,
        theta = df$theta,
        stars = stars,
        label = if (!is.null(spec$label)) spec$label else "",
        stringsAsFactors = FALSE
    )
}
rings <- do.call(rbind, ring_rows)
rings$rx <- ((rings$r0 + rings$r) / 2) * sin(rings$theta)
rings$ry <- ((rings$r0 + rings$r) / 2) * cos(rings$theta)
rings$star_angle <- -rings$theta * 180 / pi
rings$star_colour <- ifelse(abs(match(rings$fill, score_ramp) - 51) > 25, "white", "#202020")

track_index <- rep(seq_along(track_defs), each = n_cell)
guide_rows <- list()
for (i in seq_along(track_defs)) {
    lab <- track_defs[[i]]$label
    if (is.null(lab) || !nzchar(lab)) {
        next
    }
    take <- i
    if (i < length(track_defs)) {
        next_lab <- track_defs[[i + 1]]$label
        if (is.null(next_lab) || !nzchar(next_lab)) {
            take <- c(i, i + 1)
        }
    }
    rows <- rings[track_index %in% take, , drop = FALSE]
    guide_rows[[length(guide_rows) + 1L]] <- data.frame(
        x = 0,
        y = mean(c(rows$r0, rows$r)),
        label = lab,
        stringsAsFactors = FALSE
    )
}
guide_labels <- do.call(rbind, guide_rows)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
base_layers <- function(show_labels) {
    ggplot() +
        ggforce::geom_arc_bar(
            data = bands,
            aes(x0 = 0, y0 = 0, r0 = r0, r = r, start = start, end = end, fill = fill, alpha = alpha),
            colour = NA
        ) +
        geom_path(
            data = edges,
            aes(x = x, y = y, group = edge, colour = colour),
            linewidth = 1.15,
            lineend = "round"
        ) +
        geom_segment(
            data = capsules,
            aes(x = x, y = y, xend = xend, yend = yend, colour = colour),
            linewidth = 3.4,
            lineend = "round"
        ) +
        geom_text(
            data = capsules,
            aes(x = lx, y = ly, label = label, angle = angle),
            colour = "white",
            size = 2.4
        ) +
        geom_text(
            data = branch_labels,
            aes(x = x, y = y, label = label, angle = angle),
            colour = "white",
            size = 2.6,
            fontface = "bold"
        ) +
        ggforce::geom_arc_bar(
            data = df,
            aes(x0 = 0, y0 = 0, r0 = lay$tip - lay$tip_width / 2, r = lay$tip + lay$tip_width / 2,
                start = start, end = end, fill = lineage_colour),
            colour = NA, alpha = 1
        ) +
        geom_text(
            data = df,
            aes(x = tip_x * lay$tip, y = tip_y * lay$tip, label = .data[[id_col]]),
            colour = "white",
            size = 2.2
        ) +
        ggforce::geom_arc_bar(
            data = rings,
            aes(x0 = 0, y0 = 0, r0 = r0, r = r, start = start, end = end, fill = fill, alpha = alpha),
            colour = NA
        ) +
        geom_text(
            data = rings[!is.na(rings$stars), , drop = FALSE],
            aes(x = rx, y = ry, label = stars, angle = star_angle, colour = star_colour),
            size = 2.1,
            fontface = "bold"
        ) +
        scale_fill_identity() +
        scale_colour_identity() +
        scale_alpha_identity() +
        coord_fixed(xlim = c(-2.15, 2.15), ylim = c(-2.15, 2.15), clip = "off") +
        theme_void(base_size = 11) +
        theme(
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.margin = margin(10, 10, 10, 10)
        )
}

p <- base_layers(TRUE) +
    geom_text(
        data = df,
        aes(
            x = tip_x * lay$label_radius,
            y = tip_y * lay$label_radius,
            label = sprintf("%s\n(%.1f%%)", .data[[cell_col]], .data[[pct_col]]),
            colour = lineage_colour,
            angle = label_angle,
            hjust = hjust
        ),
        size = 2.5,
        lineheight = 0.9
    ) +
    geom_label(
        data = guide_labels,
        aes(x = x, y = y, label = label),
        fill = "white",
        colour = "black",
        linewidth = 0,
        size = 2.7,
        label.padding = unit(0.12, "lines")
    ) +
    annotate(
        "segment",
        x = -1.82, xend = -1.82,
        y = -0.55, yend = 0.55,
        linewidth = 0.3, linetype = "22", colour = "#888888"
    ) +
    annotate(
        "text",
        x = -1.82,
        y = seq(0.55, -0.55, length.out = 4),
        label = paste0("L", 1:4),
        size = 2.5,
        colour = "#333333"
    )

p_mini <- ggplot() +
    geom_path(
        data = edges,
        aes(x = x, y = y, group = edge, colour = colour),
        linewidth = 0.4,
        lineend = "round"
    ) +
    ggforce::geom_arc_bar(
        data = df,
        aes(
            x0 = 0, y0 = 0,
            r0 = lay$tip - lay$tip_width / 2,
            r = lay$tip + lay$tip_width / 2,
            start = start, end = end, fill = lineage_colour
        ),
        colour = NA
    ) +
    ggforce::geom_arc_bar(
        data = rings,
        aes(x0 = 0, y0 = 0, r0 = r0, r = r, start = start, end = end, fill = fill, alpha = alpha),
        colour = NA
    ) +
    annotate(
        "label", x = 0, y = 0, label = config$labels$legend,
        fill = "white", size = 2.4
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_alpha_identity() +
    coord_fixed(xlim = c(-1.7, 1.7), ylim = c(-1.7, 1.7), clip = "on") +
    theme_void() +
    theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = NA, colour = NA)
    )

final_plot <- cowplot::ggdraw(p) +
    cowplot::draw_plot(p_mini, x = 0.78, y = 0.02, width = 0.20, height = 0.20)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(final_plot, io$output, width = 10, height = 10)
