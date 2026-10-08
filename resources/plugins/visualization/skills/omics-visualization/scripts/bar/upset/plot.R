#!/usr/bin/env Rscript

# Template-ID: bar-upset
#
# Purpose:
#   Draw an UpSet combination figure: top bars of intersection size
#   (two groups: observed vs expected), a membership matrix with
#   connectors coloured by the number of member sets, and left bars
#   of set size.
#
# Inputs:
#   One row per exclusive combination. Default example:
#     - sets: comma-separated member names
#     - observed: supplied combination count
#     - expected: supplied null count (same units)
#
# Output:
#   A PDF, PNG, or SVG UpSet figure.
#
# Dependencies:
#   ggplot2, readr, ggprism, aplot, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, set_order, labels, order).
#   Edit DATA PREPARATION to change combination order.
#   Edit PLOT only when the three-panel geometry must change.
#   This is one template, not a multi-panel a/b layout.
#   Do not replace it with overlapping Venn circles.
#   UpSetR cannot encode two bar series or n-set matrix colour; keep
#   this ggplot assembly.
#
# Scientific assumptions:
#   Each row is one already-counted combination. This script does not
#   enumerate co-occurrence from raw cells and does not test
#   enrichment. expected is a supplied null, not computed here.
#   Set-size bars sum the supplied observed counts over combinations
#   that contain the set (a display aggregate, not a new assay).
#   Matrix colour is the number of members in the combination.

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
        sets = "sets",
        observed = "observed",
        expected = "expected"
    ),
    set_order = c("LS", "NAc", "BNST", "LH", "mPFC", "BA", "CeA"),
    order_by = "freq",
    decreasing = TRUE,
    labels = list(
        title = "",
        top_y = "Intersection size",
        set_y = "Set size",
        fill = NULL,
        colour = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "aplot", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

sets_col <- config$columns$sets
obs_col <- config$columns$observed
exp_col <- config$columns$expected

df[[obs_col]] <- as.numeric(df[[obs_col]])
df[[exp_col]] <- as.numeric(df[[exp_col]])
ok <- !is.na(df[[sets_col]]) & nzchar(as.character(df[[sets_col]])) &
    is.finite(df[[obs_col]]) & is.finite(df[[exp_col]])
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with invalid sets or counts.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No combinations left to plot.", call. = FALSE)
}

split_sets <- function(x) {
    parts <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
    parts[nzchar(parts)]
}

members <- lapply(df[[sets_col]], split_sets)
n_mem <- vapply(members, length, integer(1))
if (identical(config$order_by, "freq")) {
    o <- order(n_mem, -df[[obs_col]])
    df <- df[o, , drop = FALSE]
    members <- members[o]
    n_mem <- n_mem[o]
}

set_order <- unique(as.character(config$set_order))
seen <- unique(unlist(members, use.names = FALSE))
extra <- setdiff(seen, set_order)
if (length(extra)) {
    set_order <- c(set_order, extra)
}
set_order <- set_order[set_order %in% seen]
n_set <- length(set_order)
n_combo <- nrow(df)
df$combo <- seq_len(n_combo)
df$n_mem <- n_mem

bar_levels <- c("Observed", "Expected")
top <- rbind(
    data.frame(
        x = df$combo - 0.2,
        y = df[[obs_col]],
        group = "Observed",
        stringsAsFactors = FALSE
    ),
    data.frame(
        x = df$combo + 0.2,
        y = df[[exp_col]],
        group = "Expected",
        stringsAsFactors = FALSE
    )
)
top$group <- factor(top$group, levels = bar_levels)

mat <- data.frame(
    X = rep(df$combo, each = n_set),
    Y = rep(seq_along(set_order), times = n_combo),
    present = FALSE,
    n_mem = rep(df$n_mem, each = n_set),
    stringsAsFactors = FALSE
)
for (j in seq_len(n_combo)) {
    hit <- set_order %in% members[[j]]
    idx <- (j - 1L) * n_set + seq_len(n_set)
    mat$present[idx] <- hit
}
n_set_levels <- paste(sort(unique(df$n_mem[df$n_mem > 0L])), "sets")
mat$group <- factor(
    ifelse(mat$present, paste(mat$n_mem, "sets"), "0 sets"),
    levels = c("0 sets", n_set_levels)
)

seg <- data.frame(
    X = df$combo,
    ystart = NA_real_,
    yend = NA_real_,
    group = factor(paste(df$n_mem, "sets"), levels = n_set_levels),
    stringsAsFactors = FALSE
)
for (j in seq_len(n_combo)) {
    idx <- match(members[[j]], set_order)
    idx <- idx[is.finite(idx)]
    if (length(idx)) {
        seg$ystart[[j]] <- min(idx)
        seg$yend[[j]] <- max(idx)
    }
}
seg <- seg[is.finite(seg$ystart), , drop = FALSE]

set_n <- vapply(set_order, function(s) {
    keep <- vapply(members, function(m) s %in% m, logical(1))
    sum(df[[obs_col]][keep])
}, numeric(1))
side <- data.frame(
    Target = factor(set_order, levels = set_order),
    y = seq_along(set_order),
    Cell_num = unname(set_n),
    stringsAsFactors = FALSE
)

safe <- palette_colors("Qualitative.Safe")
bold <- palette_colors("Qualitative.Bold")
algolia <- palette_colors("Brand.Algolia")
fill_top <- c(Observed = algolia[[8]], Expected = safe[[12]])
n_cols <- stats::setNames(bold[seq_along(n_set_levels)], n_set_levels)
col_absent <- safe[[12]]
col_setbar <- algolia[[7]]

legend_compact <- theme(
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(0.32, "cm"),
    legend.spacing.x = unit(4, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt")
)
blank_bg <- theme(
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(2, 4, 2, 2)
) +
    legend_compact

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p_top <- ggplot() +
    geom_col(
        data = top[top$group == "Observed", , drop = FALSE],
        aes(x = x, y = y, fill = group),
        width = 0.5
    ) +
    geom_col(
        data = top[top$group == "Expected", , drop = FALSE],
        aes(x = x, y = y, fill = group),
        width = 0.24
    ) +
    scale_fill_manual(values = fill_top, drop = FALSE, name = config$labels$fill) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_x_continuous(
        limits = c(0.5, n_combo + 0.5),
        expand = c(0, 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = NULL, y = config$labels$top_y) +
    coord_cartesian(clip = "off") +
    theme_prism() +
    blank_bg +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5, 1),
        legend.justification.inside = c(0.5, 1),
        legend.direction = "horizontal",
        plot.margin = margin(2, 4, 2, 2)
    )

p_mat <- ggplot(mat, aes(x = X, y = Y)) +
    geom_segment(
        data = seg,
        aes(x = X, xend = X, y = ystart, yend = yend, colour = group),
        linewidth = 0.85,
        inherit.aes = FALSE
    ) +
    geom_point(
        data = mat[!mat$present, , drop = FALSE],
        colour = col_absent,
        size = 2.7,
        show.legend = FALSE
    ) +
    geom_point(
        data = mat[mat$present, , drop = FALSE],
        aes(colour = group),
        size = 3.9
    ) +
    scale_colour_manual(values = n_cols, drop = FALSE, name = config$labels$colour) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3.2, linewidth = 1))) +
    scale_x_continuous(
        limits = c(0.5, n_combo + 0.5),
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        breaks = seq_along(set_order),
        labels = set_order,
        limits = c(0.5, n_set + 0.5),
        expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_prism() +
    blank_bg +
    theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5, -0.16),
        legend.justification.inside = c(0.5, 1),
        legend.direction = "horizontal",
        plot.margin = margin(2, 4, 36, 2)
    )

p_side <- ggplot(side, aes(x = Cell_num, y = y)) +
    geom_col(fill = col_setbar, width = 0.78, orientation = "y") +
    scale_x_reverse(expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
        breaks = seq_along(set_order),
        labels = set_order,
        limits = c(0.5, n_set + 0.5),
        expand = c(0, 0)
    ) +
    labs(x = config$labels$set_y, y = NULL) +
    theme_prism() +
    blank_bg +
    theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = margin(2, 4, 36, 2)
    )

p_gap <- ggplot() +
    theme_void() +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0)
    )

p <- insert_left(
    insert_top(insert_top(p_mat, p_gap, height = 0.14), p_top, height = 2.05),
    p_side,
    width = 0.22
)
p <- as.patchwork(p) +
    plot_layout(guides = "keep") +
    plot_annotation(
        theme = theme(plot.background = element_rect(fill = NA, colour = NA))
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
suppressWarnings(save_ggplot(p, io$output, width = 10, height = 5.6))
