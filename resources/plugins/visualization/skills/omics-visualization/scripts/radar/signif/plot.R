#!/usr/bin/env Rscript

# Template-ID: radar-signif
#
# Purpose:
#   Radial group profiles across many features: one closed line per
#   group with an optional uncertainty ribbon, plus outer rings of dots
#   marking significant pairwise comparisons per feature. Optional
#   side-by-side panels (e.g. tissue vs plasma).
#
# Inputs:
#   Two tables (tsv/csv):
#     1. profile: feature, group, mean; optional lower, upper, panel
#     2. signif:  feature, comparison, padj; optional panel
#
# Output:
#   A PDF, PNG, or SVG radial profile with significance rings.
#
# Dependencies:
#   ggplot2, readr, patchwork
#
# Example:
#   Rscript plot.R profile.tsv signif.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, colours, alpha, comparison order, hole).
#   Pass the profile table then the signif table on the CLI.
#   Row order of features in the profile table is the clockwise order.
#
# Scientific assumptions:
#   mean, lower/upper, and padj are supplied. The script does not
#   average cells, run differential tests, or correct p-values. A dot
#   means padj < alpha for that comparison; no dot means not
#   significant or not tested. Radius is shared within a panel, so
#   features are comparable within a panel, not across panels.

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

io <- parse_io_args(n_input = 2L, input_names = c("profile", "signif"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        panel = "panel",
        feature = "feature",
        group = "group",
        mean = "mean",
        lower = "lower",
        upper = "upper"
    ),
    sig_columns = list(
        panel = "panel",
        feature = "feature",
        comparison = "comparison",
        padj = "padj"
    ),
    labels = list(
        group = "Group",
        comparison = "Padj < 0.05"
    ),
    alpha = 0.05,
    group_colours = c(
        Control = "#2E8B57",
        HFpEF = "#7B68B5",
        HFrEF = "#F08A4B"
    ),
    comparison_colours = c(
        "HFpEF vs Control" = "#7B68B5",
        "HFrEF vs Control" = "#F08A4B",
        "HFpEF vs HFrEF" = "#1F8AC0"
    ),
    palettes = list(
        group = "Qualitative.Safe",
        comparison = "Qualitative.Bold"
    ),
    hole = 0.12,           # inner radius fraction kept empty
    n_rings = 5L,
    ribbon_alpha = 0.22,
    line_width = 0.7,
    dot_size = 2.6,
    dot_ring_start = 1.08,
    dot_ring_step = 0.075,
    label_size = 2.7,
    ncol = 2L
)

load_packages(c("ggplot2", "readr", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
prof <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
sig <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
cc <- config$columns
sc <- config$sig_columns
require_columns(prof, cc[c("feature", "group", "mean")])
require_columns(sig, sc[c("feature", "comparison", "padj")])

has_panel <- cc$panel %in% names(prof)
if (!has_panel) {
    prof[[cc$panel]] <- ""
}
if (!sc$panel %in% names(sig)) {
    sig[[sc$panel]] <- if (has_panel) NA_character_ else ""
}
has_ci <- cc$lower %in% names(prof) && cc$upper %in% names(prof)

prof$.panel <- as.character(prof[[cc$panel]])
prof$.feature <- as.character(prof[[cc$feature]])
prof$.group <- as.character(prof[[cc$group]])
prof$.mean <- as.numeric(prof[[cc$mean]])
if (has_ci) {
    prof$.lower <- as.numeric(prof[[cc$lower]])
    prof$.upper <- as.numeric(prof[[cc$upper]])
}
if (any(!is.finite(prof$.mean))) {
    stop("mean must be finite for every row.", call. = FALSE)
}
key <- paste(prof$.panel, prof$.feature, prof$.group, sep = "\r")
if (anyDuplicated(key)) {
    stop("Each panel x feature x group must appear once.", call. = FALSE)
}

sig$.panel <- as.character(sig[[sc$panel]])
sig$.feature <- as.character(sig[[sc$feature]])
sig$.comparison <- as.character(sig[[sc$comparison]])
sig$.padj <- as.numeric(sig[[sc$padj]])
sig <- sig[is.finite(sig$.padj) & sig$.padj < config$alpha, , drop = FALSE]

panel_levels <- unique(prof$.panel)
group_levels <- unique(prof$.group)
comp_levels <- unique(c(
    intersect(names(config$comparison_colours), unique(sig$.comparison)),
    unique(sig$.comparison)
))

fill_colours <- function(fixed, levels, palette) {
    out <- fixed[names(fixed) %in% levels]
    miss <- setdiff(levels, names(out))
    if (length(miss)) {
        extra <- palette_colors(palette, n = length(miss))
        names(extra) <- miss
        out <- c(out, extra)
    }
    out[levels]
}
group_cols <- fill_colours(config$group_colours, group_levels, config$palettes$group)
comp_cols <- if (length(comp_levels)) {
    fill_colours(config$comparison_colours, comp_levels, config$palettes$comparison)
} else {
    character(0)
}

n_comp <- length(comp_levels)
dot_r <- config$dot_ring_start + (seq_len(max(n_comp, 1L)) - 1) * config$dot_ring_step
label_r <- max(dot_r) + 0.06
radial_angle <- function(theta) {
    a <- theta * 180 / pi
    ifelse(cos(theta) < 0, a + 180, a)
}

# -----------------------------------------------------------------------------
# PLOT (one panel per panel level)
# -----------------------------------------------------------------------------
make_panel <- function(pn) {
    d <- prof[prof$.panel == pn, , drop = FALSE]
    feats <- unique(d$.feature)
    n <- length(feats)
    if (n < 3L) {
        stop("Panel '", pn, "' needs at least three features.", call. = FALSE)
    }
    theta <- pi / 2 - 2 * pi * (seq_len(n) - 1) / n
    names(theta) <- feats
    vmax <- max(c(d$.mean, if (has_ci) d$.upper), na.rm = TRUE)
    vmax <- max(pretty(c(0, vmax)))
    to_r <- function(v) config$hole + (1 - config$hole) * pmax(v, 0) / vmax
    xy <- function(r, th) data.frame(x = r * cos(th), y = r * sin(th))

    loop <- function(v, th) {
        p <- xy(to_r(v), th)
        rbind(p, p[1, , drop = FALSE])
    }

    rings <- do.call(rbind, lapply(seq_len(config$n_rings), function(i) {
        a <- seq(0, 2 * pi, length.out = 181)
        r <- config$hole + (1 - config$hole) * i / config$n_rings
        data.frame(x = r * cos(a), y = r * sin(a), ring = i)
    }))
    spokes <- data.frame(
        x = config$hole * cos(theta), y = config$hole * sin(theta),
        xend = cos(theta), yend = sin(theta)
    )

    lines <- do.call(rbind, lapply(group_levels, function(g) {
        dg <- d[d$.group == g, , drop = FALSE]
        dg <- dg[match(feats, dg$.feature), , drop = FALSE]
        if (anyNA(dg$.mean)) {
            stop("Group '", g, "' lacks some features in panel '", pn, "'.", call. = FALSE)
        }
        out <- loop(dg$.mean, theta[feats])
        out$group <- g
        out
    }))
    lines$group <- factor(lines$group, levels = group_levels)

    p <- ggplot2::ggplot() +
        ggplot2::geom_path(
            data = rings, ggplot2::aes(x = x, y = y, group = ring),
            colour = "grey88", linewidth = 0.35
        ) +
        ggplot2::geom_segment(
            data = spokes, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
            colour = "grey90", linewidth = 0.3
        )

    if (has_ci) {
        for (g in group_levels) {
            dg <- d[d$.group == g, , drop = FALSE]
            dg <- dg[match(feats, dg$.feature), , drop = FALSE]
            up <- loop(dg$.upper, theta[feats])
            lo <- loop(dg$.lower, theta[feats])
            # Annulus as one simple polygon: outer loop, then inner loop reversed.
            band <- rbind(up, lo[rev(seq_len(nrow(lo))), , drop = FALSE])
            p <- p + ggplot2::geom_polygon(
                data = band,
                ggplot2::aes(x = x, y = y),
                fill = group_cols[[g]],
                alpha = config$ribbon_alpha,
                colour = NA
            )
        }
    }

    p <- p + ggplot2::geom_path(
        data = lines,
        ggplot2::aes(x = x, y = y, group = group, colour = group),
        linewidth = config$line_width,
        linejoin = "round"
    )

    s <- sig[(sig$.panel == pn | is.na(sig$.panel)) & sig$.feature %in% feats, , drop = FALSE]
    if (nrow(s)) {
        k <- match(s$.comparison, comp_levels)
        th <- theta[s$.feature]
        s$x <- dot_r[k] * cos(th)
        s$y <- dot_r[k] * sin(th)
        s$comparison <- factor(s$.comparison, levels = comp_levels)
        p <- p + ggplot2::geom_point(
            data = s,
            ggplot2::aes(x = x, y = y, fill = comparison),
            shape = 21, colour = "white", stroke = 0.3,
            size = config$dot_size
        )
    }

    labs_df <- data.frame(
        x = label_r * cos(theta), y = label_r * sin(theta),
        label = feats,
        angle = radial_angle(theta),
        hjust = ifelse(cos(theta) < 0, 1, 0),
        stringsAsFactors = FALSE
    )
    lim <- label_r + 0.32

    p +
        ggplot2::geom_text(
            data = labs_df,
            ggplot2::aes(x = x, y = y, label = label, angle = angle, hjust = hjust),
            size = config$label_size, colour = "grey15"
        ) +
        ggplot2::scale_colour_manual(
            name = config$labels$group, values = group_cols, drop = FALSE,
            guide = ggplot2::guide_legend(order = 1, override.aes = list(linewidth = 1.2))
        ) +
        ggplot2::scale_fill_manual(
            name = config$labels$comparison, values = comp_cols, drop = FALSE,
            limits = comp_levels,
            guide = ggplot2::guide_legend(order = 2, override.aes = list(size = 3.6))
        ) +
        ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE, clip = "off") +
        ggplot2::labs(title = if (nzchar(pn)) pn else NULL) +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0,
                margin = ggplot2::margin(b = 6)),
            legend.title = ggplot2::element_text(face = "bold", size = 10),
            legend.text = ggplot2::element_text(size = 9),
            legend.key = ggplot2::element_rect(fill = "white", colour = NA),
            plot.background = ggplot2::element_rect(fill = "white", colour = NA),
            plot.margin = ggplot2::margin(8, 8, 8, 8)
        )
}

panels <- lapply(panel_levels, make_panel)
ncol_p <- max(1L, min(as.integer(config$ncol), length(panels)))
p <- patchwork::wrap_plots(panels, ncol = ncol_p) +
    patchwork::plot_layout(guides = "collect")

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
nrow_p <- ceiling(length(panels) / ncol_p)
save_ggplot(p, io$output, width = 6.2 * ncol_p + 1.8, height = 6.4 * nrow_p)
