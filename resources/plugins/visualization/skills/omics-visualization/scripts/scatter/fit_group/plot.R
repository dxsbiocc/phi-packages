#!/usr/bin/env Rscript

# Template-ID: scatter-fit-group
#
# Purpose:
#   Draw a grouped scatter with one overall linear fit and one fit per
#   group, each with a confidence band, plus equation / R2_adj / P labels.
#
# Inputs:
#   A table with one row per observation. Default example:
#     - flipper_mm: x
#     - bill_mm: y
#     - species: group
#
# Output:
#   A PDF, PNG, or SVG multi-category fit scatter.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG for column names, labels, and whether the overall
#   fit is drawn. Edit PLOT only when the geometry must change.
#
# Scientific assumptions:
#   Each row is one observation of two numeric variables and a group.
#   Fits are ordinary least-squares (y ~ x). R2_adj and P come from
#   those display fits; they are not causal claims. Do not run
#   multivariate models or residual diagnostics in the script.

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
        x = "flipper_mm",
        y = "bill_mm",
        group = "species"
    ),
    labels = list(
        x = "Flipper length (mm)",
        y = "Bill length (mm)",
        group = NULL
    ),
    palettes = list(
        group = "Qualitative.Bold"
    ),
    show_overall = TRUE,
    overall_colour = "grey15",
    ribbon_fill = "grey70",
    ribbon_alpha = 0.28,
    point_alpha = 0.55,
    point_size = 2.2
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
g_col <- config$columns$group
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df[[g_col]] <- as.character(df[[g_col]])
ok <- is.finite(df[[x_col]]) & is.finite(df[[y_col]]) & !is.na(df[[g_col]]) & nzchar(df[[g_col]])
df <- df[ok, , drop = FALSE]
if (nrow(df) < 6L) {
    stop("Need at least 6 finite rows for grouped fits.", call. = FALSE)
}
group_levels <- unique(df[[g_col]])
if (length(group_levels) < 2L) {
    stop("Need at least two groups. Use scatter-correlation for one series.", call. = FALSE)
}
df[[g_col]] <- factor(df[[g_col]], levels = group_levels)
counts <- table(df[[g_col]])
if (any(counts < 3L)) {
    stop("Each group needs at least 3 points for a linear fit.", call. = FALSE)
}

fit_stats <- function(x, y, label) {
    fit <- stats::lm(y ~ x)
    sm <- summary(fit)
    coefs <- stats::coef(fit)
    a <- unname(coefs[[1]])
    b <- unname(coefs[[2]])
    r2 <- unname(sm$adj.r.squared)
    p <- stats::anova(fit)[["Pr(>F)"]][[1]]
    list(
        label = label,
        a = a,
        b = b,
        r2 = r2,
        p = p,
        n = length(x)
    )
}

format_eq <- function(st) {
    a_txt <- sprintf("%.3g", st$a)
    b_txt <- sprintf("%.3g", st$b)
    sign_b <- if (st$b >= 0) paste0("+ ", b_txt) else paste0("- ", sprintf("%.3g", abs(st$b)))
    r2_txt <- sprintf("%.3f", st$r2)
    if (!is.finite(st$p)) {
        p_txt <- "NA"
    } else if (st$p < 0.001) {
        p_txt <- "< 0.001"
    } else {
        p_txt <- sprintf("= %.3g", st$p)
    }
    sprintf(
        "y = %s %sx; R\u00b2_adj = %s; P %s",
        a_txt, sign_b, r2_txt, p_txt
    )
}

overall <- fit_stats(df[[x_col]], df[[y_col]], "Overall")
group_stats <- lapply(group_levels, function(g) {
    sub <- df[df[[g_col]] == g, , drop = FALSE]
    fit_stats(sub[[x_col]], sub[[y_col]], g)
})
names(group_stats) <- group_levels

pal <- setNames(
    palette_colors(config$palettes$group, n = length(group_levels)),
    group_levels
)

# Equation anchors: overall near top; groups stacked at lower-right.
xr <- range(df[[x_col]])
yr <- range(df[[y_col]])
pad_x <- 0.02 * diff(xr)
pad_y <- 0.04 * diff(yr)
eq_overall <- data.frame(
    x = mean(xr),
    y = yr[2] - pad_y,
    text = format_eq(overall),
    colour = config$overall_colour,
    stringsAsFactors = FALSE
)
eq_group <- do.call(rbind, lapply(seq_along(group_levels), function(i) {
    g <- group_levels[[i]]
    data.frame(
        x = xr[2] - pad_x,
        y = yr[1] + pad_y * (length(group_levels) - i + 1.2),
        text = format_eq(group_stats[[g]]),
        group = g,
        colour = pal[[g]],
        stringsAsFactors = FALSE
    )
}))

# Direct labels near the high-x end of each group cloud.
label_df <- do.call(rbind, lapply(group_levels, function(g) {
    sub <- df[df[[g_col]] == g, , drop = FALSE]
    # Prefer a high-x tip of the group, nudged slightly above the local mean y.
    ord <- order(sub[[x_col]], decreasing = TRUE)
    tip <- sub[ord[[1]], , drop = FALSE]
    data.frame(
        x = tip[[x_col]],
        y = tip[[y_col]] + 0.035 * diff(yr),
        label = g,
        group = g,
        stringsAsFactors = FALSE
    )
}))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(
        aes(colour = .data[[g_col]]),
        size = config$point_size,
        alpha = config$point_alpha,
        stroke = 0
    )

if (isTRUE(config$show_overall)) {
    p <- p +
        geom_smooth(
            method = "lm",
            formula = y ~ x,
            se = TRUE,
            colour = config$overall_colour,
            fill = config$ribbon_fill,
            alpha = config$ribbon_alpha,
            linetype = "22",
            linewidth = 0.7,
            inherit.aes = TRUE
        )
}

p <- p +
    geom_smooth(
        aes(colour = .data[[g_col]], group = .data[[g_col]]),
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        fill = config$ribbon_fill,
        alpha = config$ribbon_alpha,
        linetype = "22",
        linewidth = 0.65,
        show.legend = FALSE
    ) +
    geom_text(
        data = label_df,
        aes(x = x, y = y, label = label, colour = group),
        size = 3.2,
        fontface = "bold",
        show.legend = FALSE
    ) +
    geom_text(
        data = eq_overall,
        aes(x = x, y = y, label = text),
        colour = config$overall_colour,
        size = 2.7,
        hjust = 0.5,
        vjust = 1,
        inherit.aes = FALSE
    ) +
    geom_text(
        data = eq_group,
        aes(x = x, y = y, label = text, colour = group),
        size = 2.55,
        hjust = 1,
        vjust = 0.5,
        show.legend = FALSE,
        inherit.aes = FALSE
    ) +
    scale_colour_manual(values = pal, name = config$labels$group) +
    guides(colour = "none") +
    labs(x = config$labels$x, y = config$labels$y) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey55", fill = NA, linewidth = 0.5),
        axis.title = element_text(colour = "black", size = 10),
        axis.text = element_text(colour = "black", size = 9),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 12, 8, 8)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 6.4, height = 5.2)
