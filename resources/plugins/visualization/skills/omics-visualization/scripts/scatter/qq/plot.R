#!/usr/bin/env Rscript

# Template-ID: scatter-qq
#
# Purpose:
#   Distribution check: Q-Q plot (sample vs theoretical quantiles) or
#   P-P plot (empirical vs theoretical cumulative probability), with a
#   reference line and a pointwise confidence band, one panel per
#   group. detrend = TRUE subtracts the reference line.
#
# Inputs:
#   A long table (tsv/csv). Default example:
#     - value: observations
#     - group: optional; one panel per group
#
# Output:
#   A PDF, PNG, or SVG Q-Q / P-P figure.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG: columns, type ("qq" or "pp"), distribution and
#   dparams, band / conf, detrend, facet, palette, labels.
#
# Scientific assumptions:
#   The theoretical distribution and its parameters are supplied in
#   CONFIG; the script does not fit them. With dparams = NULL a normal
#   Q-Q line passes through the quartiles (qqline), so no parameters
#   are needed; a normal P-P with dparams = NULL standardises each group
#   by its own mean and sd. The band is pointwise from order-statistic
#   Beta(i, n - i + 1) quantiles; it is a visual guide, not a test.

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
        value = "value",
        group = "group"
    ),
    type = "qq",                # "qq" or "pp"
    distribution = "norm",      # any R d/p/q family: norm, t, gamma, beta, exp, lnorm, unif ...
    dparams = NULL,             # e.g. list(df = 5) for t, list(shape = 1, rate = 0.2) for gamma
    band = TRUE,
    conf = 0.95,
    detrend = FALSE,
    facet = TRUE,               # FALSE = overlay groups in one panel
    ncol = 3L,
    free_scales = TRUE,
    palettes = list(group = "Qualitative.Prism"),
    # Palette positions in use order, so overlaid groups contrast.
    palette_order = c(2, 7, 4, 8, 1, 6, 3, 9, 5, 10, 11),
    point_size = 1.7,
    point_alpha = 0.85,
    band_alpha = 0.18,
    line_colour = "grey15",
    labels = list(
        x = NULL,               # NULL = automatic from type / distribution
        y = NULL,
        title = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
cc <- config$columns
require_columns(df, cc["value"])
has_group <- !is.null(cc$group) && cc$group %in% names(df)
df$.v <- as.numeric(df[[cc$value]])
df$.g <- if (has_group) as.character(df[[cc$group]]) else "Sample"
df <- df[is.finite(df$.v) & !is.na(df$.g), , drop = FALSE]
g_levels <- unique(df$.g)

type <- match.arg(config$type, c("qq", "pp"))
dist <- config$distribution
qfun <- get(paste0("q", dist), mode = "function")
pfun <- get(paste0("p", dist), mode = "function")
alpha <- (1 - config$conf) / 2
call_d <- function(f, x, params) do.call(f, c(list(x), params))

per_group <- lapply(g_levels, function(g) {
    v <- sort(df$.v[df$.g == g])
    n <- length(v)
    if (n < 3) {
        stop("Group '", g, "' has fewer than three values.", call. = FALSE)
    }
    i <- seq_len(n)
    pp <- stats::ppoints(n)
    lo_p <- stats::qbeta(alpha, i, n - i + 1)
    hi_p <- stats::qbeta(1 - alpha, i, n - i + 1)
    params <- config$dparams
    if (identical(type, "qq")) {
        if (is.null(params)) params <- list()
        theo <- call_d(qfun, pp, params)
        # Reference line through the quartiles (as stats::qqline).
        qs <- stats::quantile(v, c(0.25, 0.75), names = FALSE)
        qt <- call_d(qfun, c(0.25, 0.75), params)
        slope <- diff(qs) / diff(qt)
        int <- qs[[1]] - slope * qt[[1]]
        ref <- int + slope * theo
        lo <- int + slope * call_d(qfun, lo_p, params)
        hi <- int + slope * call_d(qfun, hi_p, params)
        d <- data.frame(x = theo, y = v, ref = ref, lo = lo, hi = hi)
    } else {
        if (is.null(params) && identical(dist, "norm")) {
            params <- list(mean = mean(v), sd = stats::sd(v))
        } else if (is.null(params)) {
            params <- list()
        }
        emp <- call_d(pfun, v, params)
        d <- data.frame(x = pp, y = emp, ref = pp, lo = lo_p, hi = hi_p)
    }
    d$out <- d$y < d$lo | d$y > d$hi
    if (isTRUE(config$detrend)) {
        d$y <- d$y - d$ref
        d$lo <- d$lo - d$ref
        d$hi <- d$hi - d$ref
        d$ref <- 0
    }
    d$group <- g
    d
})
plot_df <- do.call(rbind, per_group)
plot_df$group <- factor(plot_df$group, levels = g_levels)
frac_out <- tapply(plot_df$out, plot_df$group, mean)

pal <- palette_colors(config$palettes$group)
ord <- config$palette_order[config$palette_order <= length(pal)]
pal <- c(pal[ord], pal[setdiff(seq_along(pal), ord)])
gcols <- expand_palette(pal, length(g_levels))
names(gcols) <- g_levels

dist_name <- c(norm = "Normal", t = "t", gamma = "Gamma", beta = "Beta",
    exp = "Exponential", lnorm = "Log-normal", unif = "Uniform",
    chisq = "Chi-squared", logis = "Logistic", weibull = "Weibull")
dname <- if (dist %in% names(dist_name)) dist_name[[dist]] else dist
x_lab <- config$labels$x
y_lab <- config$labels$y
if (is.null(x_lab)) {
    x_lab <- if (identical(type, "qq")) paste(dname, "theoretical quantiles") else
        "Theoretical cumulative probability"
}
if (is.null(y_lab)) {
    y_lab <- if (identical(type, "qq")) "Sample quantiles" else "Empirical cumulative probability"
    if (isTRUE(config$detrend)) y_lab <- paste(y_lab, "- reference")
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x))
if (isTRUE(config$band)) {
    p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lo, ymax = hi, fill = group, group = group),
        alpha = config$band_alpha, colour = NA
    )
}
p <- p +
    ggplot2::geom_line(
        ggplot2::aes(y = ref, group = group),
        colour = config$line_colour, linewidth = 0.6, linetype = "dashed"
    ) +
    ggplot2::geom_point(
        ggplot2::aes(y = y, colour = group),
        size = config$point_size, alpha = config$point_alpha, stroke = 0
    ) +
    ggplot2::scale_colour_manual(values = gcols, name = NULL) +
    ggplot2::scale_fill_manual(values = gcols, name = NULL) +
    ggplot2::labs(x = x_lab, y = y_lab, title = config$labels$title)

if (isTRUE(config$facet) && length(g_levels) > 1L) {
    strip_lab <- setNames(
        sprintf("%s\n%.0f%% outside band", g_levels, 100 * frac_out[g_levels]),
        g_levels
    )
    if (!isTRUE(config$band)) strip_lab <- setNames(g_levels, g_levels)
    p <- p +
        ggplot2::facet_wrap(
            ggplot2::vars(group), ncol = as.integer(config$ncol),
            scales = if (isTRUE(config$free_scales)) "free" else "fixed",
            labeller = ggplot2::as_labeller(strip_lab)
        )
}

p <- p +
    ggprism::theme_prism(base_size = 12, border = TRUE) +
    ggplot2::theme(
        legend.position = if (!isTRUE(config$facet) && length(g_levels) > 1L) "right" else "none",
        strip.text = ggplot2::element_text(size = 11, lineheight = 0.95),
        panel.spacing = ggplot2::unit(6, "mm"),
        plot.margin = ggplot2::margin(10, 12, 8, 8)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
n_panel <- if (isTRUE(config$facet)) length(g_levels) else 1L
n_col <- min(as.integer(config$ncol), n_panel)
n_row <- ceiling(n_panel / n_col)
save_ggplot(p, io$output, width = 3.6 * n_col + 1, height = 3.4 * n_row + 0.8)
