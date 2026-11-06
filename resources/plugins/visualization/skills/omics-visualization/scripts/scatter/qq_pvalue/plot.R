#!/usr/bin/env Rscript

# Template-ID: scatter-qq-pvalue
#
# Purpose:
#   P-value Q-Q plot for GWAS / eQTL / DE results: expected vs observed
#   -log10(p) against the y = x null, with a pointwise null band and
#   the genomic inflation factor (lambda) per group.
#
# Inputs:
#   A long table (tsv/csv). Default example:
#     - p: p-values
#     - model: optional group (one coloured series each)
#
# Output:
#   A PDF, PNG, or SVG p-value Q-Q plot.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG: columns, conf, show_lambda, thin, sig_line
#   (e.g. 5e-8), palette, labels.
#
# Scientific assumptions:
#   p-values are supplied and assumed uniform under the null; the
#   script does not run tests or correct p-values. lambda =
#   median(qchisq(1 - p, 1)) / qchisq(0.5, 1), meaningful for 1-df
#   tests on many mostly-null features. The null band uses
#   Beta(i, n - i + 1) order-statistic quantiles for the largest
#   group. Thinning only drops overlapping low points; it never
#   changes lambda (computed on all p-values).

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
        p = "p",
        group = "model"
    ),
    conf = 0.95,
    show_lambda = TRUE,
    thin = TRUE,                # drop overlapping points with expected -log10(p) < 2
    sig_line = NULL,            # e.g. 5e-8 for genome-wide significance
    palettes = list(group = "Qualitative.Prism"),
    # Palette positions in use order, so overlaid groups contrast.
    palette_order = c(2, 7, 4, 8, 1, 6, 3, 9, 5, 10, 11),
    point_size = 1.6,
    point_alpha = 0.85,
    band_fill = "grey85",
    labels = list(
        title = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
cc <- config$columns
require_columns(df, cc["p"])
has_group <- !is.null(cc$group) && cc$group %in% names(df)
df$.p <- as.numeric(df[[cc$p]])
df$.g <- if (has_group) as.character(df[[cc$group]]) else "All tests"
bad <- !is.finite(df$.p) | df$.p < 0 | df$.p > 1 | is.na(df$.g)
if (any(bad)) {
    message("Dropped ", sum(bad), " row(s) with missing or out-of-range p.")
    df <- df[!bad, , drop = FALSE]
}
zero <- df$.p == 0
if (any(zero)) {
    message(sum(zero), " p-value(s) equal 0; drawn at the smallest double.")
    df$.p[zero] <- .Machine$double.xmin
}
g_levels <- unique(df$.g)

lambda_of <- function(p) stats::median(stats::qchisq(1 - p, 1)) / stats::qchisq(0.5, 1)
lam <- vapply(g_levels, function(g) lambda_of(df$.p[df$.g == g]), numeric(1))

pts <- do.call(rbind, lapply(g_levels, function(g) {
    p <- sort(df$.p[df$.g == g])
    n <- length(p)
    d <- data.frame(
        expected = -log10(stats::ppoints(n)),
        observed = -log10(p),
        group = g
    )
    if (isTRUE(config$thin) && n > 2000) {
        low <- d$expected < 2
        key <- paste(round(d$expected[low], 2), round(d$observed[low], 2))
        keep_low <- which(low)[!duplicated(key)]
        d <- d[sort(c(keep_low, which(!low))), , drop = FALSE]
    }
    d
}))
# plotmath, so lambda renders on the PDF device too.
lab <- if (isTRUE(config$show_lambda)) {
    lapply(g_levels, function(g) {
        bquote(.(g) ~~ (lambda == .(sprintf("%.3f", lam[[g]]))))
    })
} else {
    as.list(g_levels)
}
pts$group <- factor(pts$group, levels = g_levels)

n_max <- max(table(df$.g))
i <- seq_len(n_max)
alpha <- (1 - config$conf) / 2
band <- data.frame(
    expected = -log10(stats::ppoints(n_max)),
    lo = -log10(stats::qbeta(1 - alpha, i, n_max - i + 1)),
    hi = -log10(stats::qbeta(alpha, i, n_max - i + 1))
)
if (isTRUE(config$thin) && n_max > 2000) {
    band <- band[!duplicated(round(band$expected, 3)), , drop = FALSE]
}

pal <- palette_colors(config$palettes$group)
ord <- config$palette_order[config$palette_order <= length(pal)]
pal <- c(pal[ord], pal[setdiff(seq_along(pal), ord)])
gcols <- expand_palette(pal, length(g_levels))
names(gcols) <- g_levels
lim <- max(c(pts$observed, pts$expected, band$hi)) * 1.03

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
        data = band, ggplot2::aes(x = expected, ymin = lo, ymax = hi),
        fill = config$band_fill, alpha = 0.8
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey20",
        linetype = "dashed", linewidth = 0.6)
if (!is.null(config$sig_line)) {
    p <- p + ggplot2::geom_hline(yintercept = -log10(config$sig_line),
        colour = "#cc503e", linetype = "dotted", linewidth = 0.6)
}
p <- p +
    ggplot2::geom_point(
        data = pts, ggplot2::aes(x = expected, y = observed, colour = group),
        size = config$point_size, alpha = config$point_alpha, stroke = 0
    ) +
    ggplot2::scale_colour_manual(values = gcols, labels = lab, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, max(pts$expected) * 1.03),
        expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::scale_y_continuous(limits = c(0, lim),
        expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::labs(
        x = expression(Expected ~ -log[10](italic(P))),
        y = expression(Observed ~ -log[10](italic(P))),
        title = config$labels$title
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
    ggprism::theme_prism(base_size = 12, border = TRUE) +
    ggplot2::theme(
        legend.position = "inside",
        legend.position.inside = c(0.03, 0.97),
        legend.justification = c(0, 1),
        legend.text = ggplot2::element_text(size = 10, face = "plain"),
        legend.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(10, 14, 8, 8)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 5.6, height = 5.4)
