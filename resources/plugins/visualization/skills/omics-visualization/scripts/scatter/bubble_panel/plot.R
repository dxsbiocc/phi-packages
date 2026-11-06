#!/usr/bin/env Rscript

# Template-ID: scatter-bubble-panel
#
# Purpose:
#   Dual-panel bubble scatter with a log10 x-axis, bubble size from a
#   supplied weight, colour by region/group, horizontal and vertical
#   threshold dashed lines, and optional bordered country/id labels with
#   leader segments. Points outside the display window are clipped.
#
# Inputs:
#   A long table (tsv/csv). Default example:
#     - code: point id / country code
#     - region: colour group
#     - gdp_per_capita: x (log10)
#     - y: loss share or other magnitude
#     - bubble_weight: size mapping input
#     - panel: panel title / facet key
#     - show_label: 1 to label
#     - label_x, label_y: optional manual label anchors
#
# Output:
#   A PDF, PNG, or SVG dual-panel bubble figure.
#
# Dependencies:
#   ggplot2, readr, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, thresholds, colours, size range, limits).
#
# Scientific assumptions:
#   All plotted values are supplied. This script does not compute GDP,
#   loss shares, attribution, or risk scores. sqrt + rescale of
#   bubble_weight is a display transform only.

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
        id = "code",
        region = "region",
        x = "gdp_per_capita",
        y = "y",
        size = "bubble_weight",
        panel = "panel",
        show_label = "show_label",
        label_x = "label_x",
        label_y = "label_y"
    ),
    labels = list(
        x = "GDP per capita (USD)",
        y = "Loss share of GDP (%)",
        region = "Region",
        size = "Weight"
    ),
    region_colours = c(
        Africa = "#EA5D86",
        Americas = "#F1CC5C",
        Asia = "#2FCEB0",
        Europe = "#4AA3C8",
        Oceania = "#3E6670"
    ),
    palettes = list(
        region = "Qualitative.Safe"
    ),
    # Shared vertical threshold (same x on all panels).
    vline = 9574,
    # Horizontal threshold per panel title (exact match to panel column).
    hline = c(
        "(a) Attributable direct losses" = 0.075,
        "(b) Attributable direct and indirect losses" = 0.16
    ),
    # Default hline if panel name not in hline map.
    hline_default = 0.10,
    size_range = c(2.2, 16),
    point_alpha = 0.78,
    x_breaks = c(658, 1453, 3533, 9574, 29298, 102771),
    xlim = c(500, 120000),
    ylim = c(0, 0.95),
    label_size = 2.4,
    label_fill = "white",
    seed = 1L
)

load_packages(c("ggplot2", "readr", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns[c("id", "region", "x", "y", "size", "panel")])

id_col <- config$columns$id
reg_col <- config$columns$region
x_col <- config$columns$x
y_col <- config$columns$y
sz_col <- config$columns$size
panel_col <- config$columns$panel
lab_flag <- config$columns$show_label
lx_col <- config$columns$label_x
ly_col <- config$columns$label_y

df[[id_col]] <- as.character(df[[id_col]])
df[[reg_col]] <- as.character(df[[reg_col]])
df[[panel_col]] <- as.character(df[[panel_col]])
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df[[sz_col]] <- as.numeric(df[[sz_col]])

ok <- is.finite(df[[x_col]]) & df[[x_col]] > 0 &
    is.finite(df[[y_col]]) & is.finite(df[[sz_col]]) & df[[sz_col]] >= 0 &
    !is.na(df[[id_col]]) & nzchar(df[[id_col]]) &
    !is.na(df[[reg_col]]) & nzchar(df[[reg_col]]) &
    !is.na(df[[panel_col]]) & nzchar(df[[panel_col]])
df <- df[ok, , drop = FALSE]
if (!nrow(df)) {
    stop("No finite bubble rows to plot.", call. = FALSE)
}

panel_levels <- unique(df[[panel_col]])
df[[panel_col]] <- factor(df[[panel_col]], levels = panel_levels)
df[[reg_col]] <- factor(df[[reg_col]], levels = unique(df[[reg_col]]))

# Display size: sqrt then linear rescale into CONFIG range.
sz_raw <- sqrt(pmax(df[[sz_col]], 0))
rng <- range(sz_raw, finite = TRUE)
if (!all(is.finite(rng)) || diff(rng) < 1e-12) {
    df$`.bubble_size` <- mean(as.numeric(config$size_range))
} else {
    to <- as.numeric(config$size_range)
    df$`.bubble_size` <- to[[1]] + (sz_raw - rng[[1]]) / diff(rng) * diff(to)
}

# Region colours.
reg_levels <- levels(df[[reg_col]])
reg_cols <- config$region_colours
reg_cols <- reg_cols[names(reg_cols) %in% reg_levels]
miss <- setdiff(reg_levels, names(reg_cols))
if (length(miss)) {
    extra <- palette_colors(config$palettes$region, n = length(miss))
    names(extra) <- miss
    reg_cols <- c(reg_cols, extra)
}
reg_cols <- reg_cols[reg_levels]

# Labels.
df$`.label` <- ""
if (!is.null(lab_flag) && lab_flag %in% names(df)) {
    show <- as.integer(df[[lab_flag]])
    show[is.na(show)] <- 0L
    df$`.label`[show > 0L] <- df[[id_col]][show > 0L]
}
has_lx <- !is.null(lx_col) && lx_col %in% names(df)
has_ly <- !is.null(ly_col) && ly_col %in% names(df)
if (has_lx) {
    df[[lx_col]] <- as.numeric(df[[lx_col]])
}
if (has_ly) {
    df[[ly_col]] <- as.numeric(df[[ly_col]])
}

hline_map <- config$hline
if (is.null(hline_map)) {
    hline_map <- setNames(
        rep(as.numeric(config$hline_default), length(panel_levels)),
        panel_levels
    )
}

# -----------------------------------------------------------------------------
# PLOT (one ggplot per panel, then patchwork)
# -----------------------------------------------------------------------------
make_panel <- function(panel_name) {
    d <- df[df[[panel_col]] == panel_name, , drop = FALSE]
    h0 <- hline_map[[panel_name]]
    if (is.null(h0) || !is.finite(h0)) {
        h0 <- as.numeric(config$hline_default)
    }

    lab <- d[nzchar(d$`.label`), , drop = FALSE]
    if (nrow(lab)) {
        if (has_lx && has_ly) {
            use_manual <- is.finite(lab[[lx_col]]) & is.finite(lab[[ly_col]])
            lab$`.lx` <- ifelse(use_manual, lab[[lx_col]], lab[[x_col]] * 1.35)
            lab$`.ly` <- ifelse(use_manual, lab[[ly_col]], lab[[y_col]] + 0.06)
        } else {
            lab$`.lx` <- lab[[x_col]] * 1.35
            lab$`.ly` <- lab[[y_col]] + 0.06
        }
    }

    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
        ggplot2::geom_vline(
            xintercept = as.numeric(config$vline),
            linetype = "dashed",
            colour = "grey35",
            linewidth = 0.45
        ) +
        ggplot2::geom_hline(
            yintercept = h0,
            linetype = "dashed",
            colour = "grey35",
            linewidth = 0.45
        ) +
        ggplot2::geom_point(
            ggplot2::aes(size = .data$`.bubble_size`, colour = .data[[reg_col]]),
            alpha = as.numeric(config$point_alpha),
            shape = 16
        )

    if (nrow(lab)) {
        p <- p +
            ggplot2::geom_segment(
                data = lab,
                ggplot2::aes(
                    x = .data[[x_col]], y = .data[[y_col]],
                    xend = .data$`.lx`, yend = .data$`.ly`,
                    colour = .data[[reg_col]]
                ),
                linewidth = 0.35,
                alpha = 0.85,
                show.legend = FALSE
            ) +
            ggplot2::geom_label(
                data = lab,
                ggplot2::aes(
                    x = .data$`.lx`, y = .data$`.ly`,
                    label = .data$`.label`,
                    colour = .data[[reg_col]]
                ),
                fill = config$label_fill,
                size = as.numeric(config$label_size),
                linewidth = 0.35,
                label.padding = ggplot2::unit(0.12, "lines"),
                label.r = ggplot2::unit(0.08, "lines"),
                show.legend = FALSE
            )
    }

    p +
        ggplot2::scale_x_log10(
            breaks = as.numeric(config$x_breaks),
            labels = scales::label_comma(accuracy = 1)
        ) +
        ggplot2::scale_colour_manual(
            name = config$labels$region,
            values = reg_cols
        ) +
        ggplot2::scale_size_identity(guide = "none") +
        ggplot2::coord_cartesian(
            xlim = as.numeric(config$xlim),
            ylim = as.numeric(config$ylim),
            clip = "on"
        ) +
        ggplot2::labs(
            title = panel_name,
            x = config$labels$x,
            y = config$labels$y
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = 11, face = "bold", hjust = 0, colour = "#222222"
            ),
            panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.title = ggplot2::element_text(size = 9, face = "bold"),
            legend.text = ggplot2::element_text(size = 8.5),
            legend.key = ggplot2::element_rect(fill = "white", colour = NA),
            axis.title = ggplot2::element_text(size = 10),
            axis.text = ggplot2::element_text(size = 8),
            plot.background = ggplot2::element_rect(fill = "white", colour = NA),
            plot.margin = ggplot2::margin(8, 10, 6, 8)
        )
}

# scales::label_comma — load scales lightly via ggplot2 imports; ensure present.
if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required for axis labels.", call. = FALSE)
}

panels <- lapply(seq_along(panel_levels), function(i) {
    p0 <- make_panel(panel_levels[[i]])
    # Keep legend only on the last panel; patchwork collects it.
    if (i < length(panel_levels)) {
        p0 <- p0 + ggplot2::theme(legend.position = "none")
    }
    p0
})
p <- patchwork::wrap_plots(panels, nrow = 1L) +
    patchwork::plot_layout(guides = "collect")

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = max(10, 5.2 * length(panel_levels)),
    height = 5.6
)
