#!/usr/bin/env Rscript

# Template-ID: scatter-taxon-manhattan
#
# Purpose:
#   Taxon-ordered "Manhattan" bubble plot of differential OTUs / ASVs.
#   One stacked panel per comparison; each point is one OTU at a fixed
#   x position grouped by Order (blocks) and Family (within a block).
#   y = supplied log2 fold change, size = relative abundance, colour =
#   Order, filled = significant, hollow = not significant. A coloured
#   Order label strip sits under the last panel.
#
# Inputs:
#   A long table (tsv/csv), one row per OTU-in-comparison. Default example:
#     - otu: feature id (same id keeps the same x in every panel)
#     - order: colour / block taxon (e.g. Order)
#     - family: optional within-block sort key (e.g. Family)
#     - comparison: panel id (e.g. "B_AZM vs B")
#     - log2fc: supplied log2 fold change
#     - abundance: relative abundance used for point size
#     - significant: yes/no, TRUE/FALSE, or 1/0
#
# Output:
#   A PDF, PNG, or SVG multi-panel taxon Manhattan plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, top_n_orders, others, palette, y limits,
#   size breaks, labels). Panel count follows the comparison column.
#
# Scientific assumptions:
#   log2fc, abundance and significance are supplied (edgeR, DESeq2,
#   ANCOM-BC, ...). This script runs no test. Keeping the top N Orders
#   by summed abundance is a display filter. Points are never jittered
#   in y, so the plotted fold change is the supplied one; x positions
#   are deterministic, so the same OTU sits at the same x in each panel.

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
        otu = "otu",
        order = "order",
        family = "family",
        comparison = "comparison",
        log2fc = "log2fc",
        abundance = "abundance",
        significant = "significant"
    ),
    labels = list(
        title = NULL,
        y = NULL,               # NULL draws Log[2] FC with a subscript
        size = "Relative abundance",
        status = NULL
    ),
    # Keep the N Orders with the largest summed abundance; NULL keeps all
    top_n_orders = 25,
    # Orders outside the top N: "drop" or "merge" (into others_label)
    others = "drop",
    others_label = "Others",
    # Catalog palette ids, concatenated in order, then extended if short
    palette = c("Qualitative.Paired", "Qualitative.Pastel", "Qualitative.Vivid"),
    ylim = c(-4, 4),
    zero_line = TRUE,
    order_bands = TRUE,         # alternating light background per Order
    band_fill = "#f3f4f6",
    block_gap = 3,              # empty x slots between Order blocks
    size_range = c(0.4, 9),
    size_breaks = c(0.001, 0.01, 0.03, 0.05, 0.1),
    stroke = 0.5,
    comparison_label_size = 4.6,
    order_label_size = 3.4,
    order_label_angle = 60,
    strip_height = 0.12,        # Order strip height relative to one panel
    width = 13,
    height = 10
)

load_packages(c("ggplot2", "readr", "ggprism", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
cols <- config$columns
if (!cols$family %in% names(df)) {
    cols$family <- NULL
}
require_columns(df, cols)

d <- data.frame(
    otu = as.character(df[[cols$otu]]),
    order = as.character(df[[cols$order]]),
    family = if (is.null(cols$family)) "" else as.character(df[[cols$family]]),
    comparison = as.character(df[[cols$comparison]]),
    log2fc = as.numeric(df[[cols$log2fc]]),
    abundance = as.numeric(df[[cols$abundance]]),
    stringsAsFactors = FALSE
)
sig_raw <- tolower(as.character(df[[cols$significant]]))
d$significant <- sig_raw %in% c("yes", "true", "1", "sig", "significant")

order_total <- sort(tapply(d$abundance, d$order, sum, na.rm = TRUE),
    decreasing = TRUE)
order_levels <- names(order_total)
if (!is.null(config$top_n_orders) && length(order_levels) > config$top_n_orders) {
    keep <- order_levels[seq_len(config$top_n_orders)]
    if (identical(config$others, "merge")) {
        d$order[!d$order %in% keep] <- config$others_label
        order_levels <- c(keep, config$others_label)
    } else {
        d <- d[d$order %in% keep, , drop = FALSE]
        order_levels <- keep
    }
}
d$order <- factor(d$order, levels = order_levels)
d$comparison <- factor(d$comparison, levels = unique(d$comparison))

otu_tab <- unique(d[, c("otu", "order", "family")])
otu_tab <- otu_tab[!duplicated(otu_tab$otu), , drop = FALSE]
otu_abund <- tapply(d$abundance, d$otu, sum, na.rm = TRUE)
otu_tab$total <- as.numeric(otu_abund[otu_tab$otu])
family_total <- tapply(otu_tab$total, otu_tab$family, sum)
otu_tab$family_total <- as.numeric(family_total[otu_tab$family])
otu_tab <- otu_tab[order(otu_tab$order, -otu_tab$family_total,
    otu_tab$family, otu_tab$otu), , drop = FALSE]
block_index <- as.integer(otu_tab$order) - 1L
otu_tab$x <- seq_len(nrow(otu_tab)) + block_index * config$block_gap
d$x <- otu_tab$x[match(d$otu, otu_tab$otu)]

blocks <- do.call(rbind, lapply(split(otu_tab, otu_tab$order, drop = TRUE),
    function(b) data.frame(order = b$order[[1]],
        xmin = min(b$x) - 0.5, xmax = max(b$x) + 0.5)))
blocks$order <- factor(blocks$order, levels = order_levels)
blocks <- blocks[order(blocks$order), , drop = FALSE]
blocks$xmid <- (blocks$xmin + blocks$xmax) / 2
blocks$band <- seq_len(nrow(blocks)) %% 2 == 0
x_limits <- c(min(blocks$xmin) - config$block_gap, max(blocks$xmax) + config$block_gap)

palette_ids <- config$palette
order_cols <- unique(unlist(lapply(palette_ids, palette_colors)))
order_cols <- setNames(expand_palette(order_cols, length(order_levels)),
    order_levels)

d$fill <- ifelse(d$significant, order_cols[as.character(d$order)], NA)
d$status <- factor(ifelse(d$significant, "Significant", "Not significant"),
    levels = c("Significant", "Not significant"))
d <- d[order(d$significant, -d$abundance), , drop = FALSE]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p_main <- ggplot(d, aes(x = x, y = log2fc))
if (isTRUE(config$order_bands)) {
    p_main <- p_main + geom_rect(
        data = blocks[blocks$band, , drop = FALSE],
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE, fill = config$band_fill, colour = NA
    )
}
if (isTRUE(config$zero_line)) {
    p_main <- p_main + geom_hline(yintercept = 0, linewidth = 0.35,
        linetype = "dashed", colour = "#9a9a9a")
}
p_main <- p_main +
    geom_point(
        aes(size = abundance, colour = order, fill = fill),
        shape = 21, stroke = config$stroke
    ) +
    geom_point(
        aes(shape = status), x = x_limits[[1]], y = 0,
        alpha = 0, size = 3, colour = "#555555", fill = "#555555"
    ) +
    facet_wrap(~comparison, ncol = 1) +
    scale_fill_identity() +
    scale_colour_manual(values = order_cols, guide = "none") +
    scale_shape_manual(
        values = c("Significant" = 21, "Not significant" = 21),
        name = config$labels$status,
        guide = guide_legend(order = 2, override.aes = list(
            alpha = 1, size = 4, stroke = 0.8,
            fill = c("#555555", NA)
        ))
    ) +
    scale_size_continuous(
        range = config$size_range, breaks = config$size_breaks,
        labels = format(config$size_breaks, scientific = FALSE,
            drop0trailing = TRUE),
        name = config$labels$size,
        guide = guide_legend(order = 1, override.aes = list(
            colour = "#555555", fill = NA))
    ) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    coord_cartesian(ylim = config$ylim) +
    labs(title = config$labels$title, x = NULL,
        y = if (is.null(config$labels$y)) expression(Log[2] ~ "FC") else config$labels$y) +
    theme_prism(border = TRUE, base_size = 12) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(hjust = 1, face = "plain",
            size = config$comparison_label_size * .pt,
            margin = margin(1, 1, 1.5, 0, unit = "mm")),
        strip.background = element_blank(),
        panel.spacing.y = unit(2, "mm"),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "plain"),
        legend.text = element_text(size = 10, face = "plain"),
        legend.key.spacing.y = unit(1, "mm")
    )

p_strip <- ggplot(blocks) +
    geom_rect(
        aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = order),
        colour = "white", linewidth = 0.3
    ) +
    geom_text(
        aes(x = xmid, y = -0.6, label = order),
        angle = config$order_label_angle, hjust = 1, vjust = 1,
        size = config$order_label_size, colour = "#222222"
    ) +
    scale_fill_manual(values = order_cols, guide = "none") +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0.5, 0, 30, 0, unit = "mm"))

n_panels <- nlevels(d$comparison)
p <- p_main / p_strip +
    plot_layout(heights = c(n_panels, config$strip_height))

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$width, height = config$height)
