#!/usr/bin/env Rscript

# Template-ID: line-roc
#
# Purpose:
#   Draw multi-group ROC curves with optional confidence ribbons, a
#   diagonal chance line, and supplied AUC text in the legend labels.
#
# Inputs:
#   One row per curve point (example.tsv):
#     - group: short series / contrast label
#     - fpr: 1 - specificity (0–1)
#     - tpr: sensitivity (0–1)
#     - ci_lo, ci_hi: optional ribbon bounds on the tpr scale
#   Optional auc.tsv beside the input: group, auc_label
#     (short AUC text appended to legend keys)
#
# Output:
#   A PDF, PNG, or SVG ROC figure.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, palette, step vs path, diagonal, and
#   legend. Do not compute ROC points, AUC, or CIs in this script;
#   those values are supplied. Survival curves stay line-survival.
#
# Scientific assumptions:
#   fpr / tpr / ribbon bounds and AUC labels are supplied display
#   coordinates. The diagonal is a chance reference (AUC = 0.5), not
#   a fitted model.

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
        fpr = "fpr",
        tpr = "tpr",
        ci_lo = "ci_lo",
        ci_hi = "ci_hi"
    ),
    # Optional sidecar: group, auc_label  (merged into legend text).
    auc_file = "auc.tsv",
    group_order = c(
        "A vs NT +MCHAT", "A vs NT", "A vs Non-A", "A+LD vs NT", "LD vs NT"
    ),
    labels = list(
        title = NULL,
        x = "1 – specificity",
        y = "Sensitivity",
        caption = NULL
    ),
    # "step" = empirical stairs; "path" = connected polyline.
    line_type = "step",
    linewidth = 1.05,
    ribbon_alpha = 0.22,
    show_ribbon = TRUE,
    show_diagonal = TRUE,
    diagonal_linetype = "dashed",
    diagonal_linewidth = 0.55,
    # Soft chance line (catalog coral) rather than harsh primary red.
    palettes = list(
        group = "Qualitative.Bold",
        diagonal = "Qualitative.Bold:5"
    ),
    legend_position = "right",
    aspect_ratio = 1,
    base_size = 13,
    size = list(width = 9.2, height = 6.2)
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols[[1]]
}

order_levels <- function(vals, preferred) {
    lv <- unique(as.character(vals))
    if (is.null(preferred)) return(lv)
    c(intersect(as.character(preferred), lv), setdiff(lv, as.character(preferred)))
}

df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
g_col <- config$columns$group
x_col <- config$columns$fpr
y_col <- config$columns$tpr
lo_col <- config$columns$ci_lo
hi_col <- config$columns$ci_hi
require_columns(df, list(g_col, x_col, y_col))

df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
has_ci <- lo_col %in% names(df) && hi_col %in% names(df)
if (has_ci) {
    df[[lo_col]] <- as.numeric(df[[lo_col]])
    df[[hi_col]] <- as.numeric(df[[hi_col]])
}
df <- df[is.finite(df[[x_col]]) & is.finite(df[[y_col]]), , drop = FALSE]
if (!nrow(df)) stop("No finite ROC points.", call. = FALSE)

df[[g_col]] <- as.character(df[[g_col]])

# Optional AUC labels
auc_path <- config$auc_file
if (!is.null(auc_path) && nzchar(as.character(auc_path))) {
    if (!grepl("^(/|[A-Za-z]:)", auc_path) && !file.exists(auc_path)) {
        cand <- file.path(dirname(io$input), auc_path)
        if (file.exists(cand)) auc_path <- cand
    }
    if (file.exists(auc_path)) {
        auc <- as.data.frame(read_table_auto(auc_path), stringsAsFactors = FALSE)
        require_columns(auc, list("group", "auc_label"))
        map <- stats::setNames(
            as.character(auc$auc_label),
            as.character(auc$group)
        )
        df$legend <- vapply(df[[g_col]], function(g) {
            if (!is.null(map[[g]]) && nzchar(map[[g]])) {
                paste0(g, " | ", map[[g]])
            } else {
                g
            }
        }, character(1))
    } else {
        df$legend <- df[[g_col]]
    }
} else {
    df$legend <- df[[g_col]]
}

leg_lv <- order_levels(df$legend, config$group_order)
# If group_order matches raw group names, map through legend.
if (!is.null(config$group_order)) {
    raw <- order_levels(df[[g_col]], config$group_order)
    mapped <- unique(df$legend[match(raw, df[[g_col]])])
    mapped <- mapped[!is.na(mapped)]
    leg_lv <- c(mapped, setdiff(leg_lv, mapped))
}
df$legend <- factor(df$legend, levels = leg_lv)

# Stable within-group order along fpr
df <- df[order(df$legend, df[[x_col]], df[[y_col]]), , drop = FALSE]

cols <- expand_palette(
    palette_colors(config$palettes$group), nlevels(df$legend)
)
names(cols) <- levels(df$legend)
diag_col <- pal_ref(config$palettes$diagonal)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    colour = .data$legend,
    fill = .data$legend,
    group = .data$legend
))

if (isTRUE(config$show_ribbon) && has_ci) {
    p <- p + geom_ribbon(
        aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]),
        alpha = config$ribbon_alpha,
        colour = NA,
        linewidth = 0
    )
}

if (identical(config$line_type, "path")) {
    p <- p + geom_path(linewidth = config$linewidth, lineend = "round")
} else {
    p <- p + geom_step(
        linewidth = config$linewidth,
        direction = "hv",
        lineend = "round"
    )
}

if (isTRUE(config$show_diagonal)) {
    p <- p + geom_abline(
        intercept = 0, slope = 1,
        linetype = config$diagonal_linetype,
        colour = diag_col,
        linewidth = config$diagonal_linewidth,
        alpha = 0.85
    )
}

p <- p +
    scale_colour_manual(values = cols, name = NULL) +
    scale_fill_manual(values = cols, name = NULL) +
    scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = c(0.01, 0.01)
    ) +
    coord_fixed(ratio = config$aspect_ratio) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        caption = config$labels$caption
    ) +
    guides(
        colour = guide_legend(ncol = 1, override.aes = list(
            linewidth = 1.4, alpha = 1, fill = NA
        )),
        fill = "none"
    ) +
    theme_prism(base_size = config$base_size) +
    theme(
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = config$legend_position,
        legend.justification = if (identical(config$legend_position, "right")) {
            "center"
        } else {
            "left"
        },
        legend.text = element_text(size = rel(0.8)),
        legend.key.width = grid::unit(1.1, "lines"),
        legend.spacing.y = grid::unit(3, "pt"),
        legend.background = element_blank(),
        plot.caption = element_text(
            hjust = 0, size = rel(0.72), colour = "grey40",
            margin = margin(t = 10)
        ),
        plot.margin = margin(10, 10, 8, 10)
    )

if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
    p <- p + theme(plot.title = element_blank())
}
if (is.null(config$labels$caption) || !nzchar(as.character(config$labels$caption))) {
    p <- p + theme(plot.caption = element_blank())
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width, height = config$size$height)
