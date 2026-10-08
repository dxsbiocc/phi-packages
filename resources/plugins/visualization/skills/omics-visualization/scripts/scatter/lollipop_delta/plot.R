#!/usr/bin/env Rscript

# Template-ID: scatter-lollipop-delta
#
# Purpose:
#   Draw a signed horizontal lollipop of a supplied delta (for example
#   post − pre ssGSEA), with stem to zero, point fill by signature set,
#   point size by -log10(p), and stroke by a p-value bin.
#
# Inputs:
#   One row per item. Default example:
#     - celltype: immune signature / cell-type label
#     - delta: signed difference (already computed)
#     - pvalue: supplied p-value
#     - group: signature set (Cibersort, xCell, …)
#
# Output:
#   A PDF, PNG, or SVG signed lollipop plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, labels, label_strip).
#   The left coloured-name strip is CONFIG, not a second template.
#   Faceting by a sample-set column is CONFIG. Do not add an id per
#   immune deconvolution method. Do not run Wilcoxon or ssGSEA here.
#
# Scientific assumptions:
#   delta and pvalue are supplied results. This script does not test
#   groups or score signatures. -log10(pvalue) is a display transform.
#   Point size encodes significance, not effect size. Stem length is
#   the supplied delta. A Cleveland plot without p-size or a signed
#   zero line is scatter-cleveland.

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
        y = "celltype",
        x = "delta",
        p = "pvalue",
        group = "group"
    ),
    label_strip = TRUE,
    labels = list(
        title = "",
        x = "delta (post − pre)",
        y = "",
        fill = "Set",
        size = "-log10(p)",
        stroke = "p"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

y_col <- config$columns$y
x_col <- config$columns$x
p_col <- config$columns$p
g_col <- config$columns$group

df[[y_col]] <- as.character(df[[y_col]])
df[[g_col]] <- as.character(df[[g_col]])
df[[x_col]] <- as.numeric(df[[x_col]])
df[[p_col]] <- as.numeric(df[[p_col]])
ok <- is.finite(df[[x_col]]) & is.finite(df[[p_col]]) & df[[p_col]] > 0
if (any(!ok)) {
    message(
        "Dropped ", sum(!ok),
        " row(s) with non-finite delta or pvalue <= 0."
    )
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

# File order is top-to-bottom.
df[[y_col]] <- factor(df[[y_col]], levels = rev(unique(df[[y_col]])))
df[[g_col]] <- factor(df[[g_col]], levels = unique(df[[g_col]]))
df$neg_log10_p <- -log10(df[[p_col]])
df$p_bin <- ifelse(
    df[[p_col]] < 0.05,
    "<0.05",
    ifelse(df[[p_col]] < 0.1, "<0.1", ">=0.1")
)
df$p_bin <- factor(df$p_bin, levels = c("<0.05", "<0.1", ">=0.1"))

fill_cols <- palette_colors("Qualitative.Safe", n = nlevels(df[[g_col]]))
names(fill_cols) <- levels(df[[g_col]])
stroke_cols <- c(
    "<0.05" = "black",
    "<0.1" = "#888888",
    ">=0.1" = "#888888"
)
y_text_cols <- fill_cols[as.character(df[[g_col]][order(df[[y_col]])])]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p_main <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])
) +
    ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        colour = "black",
        linewidth = 0.4
    ) +
    ggplot2::geom_col(
        ggplot2::aes(fill = .data[[g_col]]),
        width = 0.12,
        show.legend = TRUE
    ) +
    ggplot2::geom_point(
        ggplot2::aes(
            size = neg_log10_p,
            fill = .data[[g_col]],
            colour = p_bin
        ),
        shape = 21,
        stroke = 1.4
    ) +
    ggplot2::scale_fill_manual(
        values = fill_cols,
        name = config$labels$fill
    ) +
    ggplot2::scale_colour_manual(
        values = stroke_cols,
        name = config$labels$stroke,
        drop = FALSE
    ) +
    ggplot2::scale_size_continuous(
        range = c(2, 6),
        name = config$labels$size
    ) +
    ggplot2::labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    ggplot2::guides(
        fill = ggplot2::guide_legend(override.aes = list(size = 4, colour = NA)),
        colour = ggplot2::guide_legend(override.aes = list(size = 4)),
        size = ggplot2::guide_legend()
    ) +
    ggprism::theme_prism() +
    ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank()
    )

if (isTRUE(config$label_strip)) {
    p_main <- p_main +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank()
        )
    p_lab <- ggplot2::ggplot(
        df,
        ggplot2::aes(y = .data[[y_col]], colour = .data[[g_col]])
    ) +
        ggplot2::geom_text(
            ggplot2::aes(x = 0, label = .data[[y_col]]),
            hjust = 1,
            fontface = "bold",
            size = 3.2,
            show.legend = FALSE
        ) +
        ggplot2::scale_colour_manual(values = fill_cols) +
        ggplot2::coord_cartesian(xlim = c(-0.55, -0.02), clip = "off") +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.background = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(5.5, 0, 5.5, 5.5)
        )
    p <- p_lab + p_main + patchwork::plot_layout(widths = c(1.05, 2.35))
} else {
    p_main <- p_main +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(colour = y_text_cols)
        )
    p <- p_main
}

n_row <- nlevels(df[[y_col]])
fig_h <- min(14, max(6, 0.32 * n_row + 2.2))

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 9, height = fig_h)
