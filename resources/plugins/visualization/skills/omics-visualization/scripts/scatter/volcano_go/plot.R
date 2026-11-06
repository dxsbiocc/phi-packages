#!/usr/bin/env Rscript

# Template-ID: scatter-volcano-go
#
# Purpose:
#   Draw a volcano and bilateral GO enrichment bars in one ggplot
#   coordinate system. Down terms sit on the left shoulder, up terms
#   on the right. Bar geometry is placed in physical points so the
#   enrichment panels stay aligned with the volcano when the data
#   range changes. This is not a patchwork of two plots.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. deg: gene, log2FC, pvalue
#     2. go:  direction (down|up), term, score
#        term may contain the two-character sequence \n for line breaks.
#        Row order within each direction is the display order (first row
#        is the top bar).
#
# Output:
#   A PDF, PNG, or SVG combined volcano–GO figure.
#
# Dependencies:
#   ggplot2, readr, dplyr, ggrepel
#
# Example:
#   Rscript plot.R deg.tsv go.tsv output.pdf
#
# Agent adaptation:
#   For a new dataset, edit only CONFIG (column names, thresholds, labels).
#   Pass the differential table, then the enrichment table, on the CLI.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   log2FC and p-values are supplied. -log10(p) is a display transform.
#   Point class is a display split at the CONFIG thresholds, not a new test.
#   Enrichment scores are supplied. This script does not run DE or enrichment.
#   Bar length maps each side's score onto that side's |log2FC| > threshold
#   shoulder. Layout follows Nature s41586-026-10459-x Figure 6i;
#   colours come from library palettes (CONFIG palettes / tint).

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

io <- parse_io_args(n_input = 2L, input_names = c("deg", "go"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        gene = "gene",
        log2fc = "log2FC",
        pvalue = "pvalue"
    ),
    go_columns = list(
        direction = "direction",
        term = "term",
        score = "score"
    ),
    thresholds = list(
        log2fc = 1,
        pvalue = 0.05,
        n_labels = 3L
    ),
    labels = list(
        x = "log[2]*\"[FC]\"",
        y = "-log[10]*'[P]'",
        go_y = "GO term description",
        score = "-log[10]*'[P]'"
    ),
    # Library palettes: one hue per direction; zones and bars are
    # white-mixed tints of the same hue.
    palettes = list(
        side = "Diverging.TealRose",
        side_index = c(down = 1L, up = 7L),
        neutral = "Diverging.RdGy",
        ns_index = 7L,
        ns_fc_index = 8L,
        sig = "Diverging.PuOr",
        sig_index = 8L
    ),
    tint = list(
        zone = 0.92,
        bar = 0.55
    )
)

load_packages(c("ggplot2", "readr", "dplyr", "ggrepel"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
deg <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
require_columns(deg, config$columns)
gene_col <- config$columns$gene
fc_col <- config$columns$log2fc
p_col <- config$columns$pvalue

num <- function(v) suppressWarnings(as.numeric(v))
deg[[gene_col]] <- as.character(deg[[gene_col]])
deg[[fc_col]] <- num(deg[[fc_col]])
deg[[p_col]] <- num(deg[[p_col]])
if (all(is.na(deg[[p_col]])) && "padj" %in% names(deg)) {
    deg[[p_col]] <- num(deg[["padj"]])
}
deg <- deg[!is.na(deg[[gene_col]]) & nzchar(deg[[gene_col]]) &
    is.finite(deg[[fc_col]]) & is.finite(deg[[p_col]]), , drop = FALSE]
if (!nrow(deg)) {
    stop("No finite gene / log2FC / pvalue rows.", call. = FALSE)
}
if (any(deg[[p_col]] < 0 | deg[[p_col]] > 1)) {
    stop("pvalue must lie in [0, 1].", call. = FALSE)
}

go <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
require_columns(go, config$go_columns)
dir_col <- config$go_columns$direction
term_col <- config$go_columns$term
score_col <- config$go_columns$score
go[[dir_col]] <- tolower(trimws(as.character(go[[dir_col]])))
go[[term_col]] <- gsub("\\n", "\n", as.character(go[[term_col]]), fixed = TRUE)
go[[score_col]] <- num(go[[score_col]])
if (!all(go[[dir_col]] %in% c("down", "up"))) {
    stop('direction must be "down" or "up".', call. = FALSE)
}
if (any(!is.finite(go[[score_col]]) | go[[score_col]] <= 0)) {
    stop("score must be positive and finite.", call. = FALSE)
}
if (!any(go[[dir_col]] == "down") || !any(go[[dir_col]] == "up")) {
    stop("go table needs at least one down term and one up term.", call. = FALSE)
}

mix_white <- function(col, p) {
    m <- (1 - p) * grDevices::col2rgb(col) / 255 + p
    grDevices::rgb(m[1, ], m[2, ], m[3, ])
}
pc <- config$palettes
side_pal <- palette_colors(pc$side)
neutral_pal <- palette_colors(pc$neutral)
cols <- list(
    deg_down = side_pal[[pc$side_index[["down"]]]],
    deg_up = side_pal[[pc$side_index[["up"]]]],
    ns = neutral_pal[[pc$ns_index]],
    ns_fc = neutral_pal[[pc$ns_fc_index]],
    sig = palette_colors(pc$sig)[[pc$sig_index]]
)
cols$zone_down <- mix_white(cols$deg_down, config$tint$zone)
cols$zone_up <- mix_white(cols$deg_up, config$tint$zone)
cols$bar_down <- mix_white(cols$deg_down, config$tint$bar)
cols$bar_up <- mix_white(cols$deg_up, config$tint$bar)

FC_THR <- config$thresholds$log2fc
P_THR <- config$thresholds$pvalue
Y_HLINE <- -log10(P_THR)
n_lab <- as.integer(config$thresholds$n_labels)

# Physical scale taken from the source figure (pt). Data range may change;
# bar height, pitch, and tick length stay in points.
PT_PER_X <- 32.695
PANEL_W <- 5.802 * PT_PER_X
PANEL_H <- 11.787 * 20.735
PT_TOP <- 4.214 * 20.735

nlp_raw <- -log10(deg[[p_col]])
nlp_raw[deg[[p_col]] == 0] <- NA_real_
if (all(!is.finite(nlp_raw))) {
    stop("pvalue column is all zero.", call. = FALSE)
}
nlp_cap <- max(nlp_raw, na.rm = TRUE) + 0.3
nlp_raw[!is.finite(nlp_raw)] <- nlp_cap

xspan <- max(deg[[fc_col]]) - min(deg[[fc_col]])
if (!is.finite(xspan) || xspan <= 0) {
    xspan <- 1
}
XMIN <- min(deg[[fc_col]]) - 0.04 * xspan
XMAX <- max(deg[[fc_col]]) + 0.04 * xspan
Y_BOT <- min(0, min(nlp_raw)) - 0.30

h <- (PANEL_H - PT_TOP) / (max(nlp_raw) - Y_BOT + 0.25)
w <- PANEL_W / (XMAX - XMIN)
u <- function(pt) pt / h
ux <- function(pt) pt / w

Y_VOL_TOP <- max(nlp_raw) + 0.25
Y_CAP <- Y_VOL_TOP - 0.07
Y_STRIP_TOP <- Y_VOL_TOP + u(15.10)
BAR_H <- u(13.11)
PITCH <- u(14.69)
TICKLEN <- u(2.58)
ASPECT <- h / w

n_rows <- max(sum(go[[dir_col]] == "down"), sum(go[[dir_col]] == "up"))
ROW0 <- Y_STRIP_TOP + u(6.59)
Y_TOP <- ROW0 + (n_rows - 1) * PITCH + BAR_H / 2 + u(6.85)

brks_x <- pretty(c(XMIN, XMAX), n = 6)
brks_x <- brks_x[brks_x >= XMIN & brks_x <= XMAX]
brks_y <- pretty(c(Y_BOT, Y_VOL_TOP), n = 5)
brks_y <- brks_y[brks_y >= Y_BOT & brks_y <= Y_VOL_TOP]

df <- deg %>%
    transmute(
        gene = .data[[gene_col]],
        log2FC = pmax(pmin(.data[[fc_col]], XMAX - 0.02), XMIN + 0.02),
        neg_log10_p = pmin(ifelse(.data[[p_col]] == 0, nlp_cap, -log10(.data[[p_col]])), Y_CAP),
        group = dplyr::case_when(
            neg_log10_p >= Y_HLINE & abs(log2FC) >= FC_THR ~ "deg",
            neg_log10_p >= Y_HLINE ~ "sig",
            abs(log2FC) >= FC_THR ~ "ns_fc",
            TRUE ~ "ns"
        )
    )

lab_deg <- df %>%
    filter(group == "deg") %>%
    mutate(side = log2FC < 0)
lab_L <- lab_deg %>%
    filter(side) %>%
    slice_max(abs(log2FC), n = n_lab, with_ties = FALSE)
lab_R <- lab_deg %>%
    filter(!side) %>%
    slice_max(abs(log2FC), n = n_lab, with_ties = FALSE)

Lmax <- max(go[[score_col]][go[[dir_col]] == "down"])
Rmax <- max(go[[score_col]][go[[dir_col]] == "up"])
if ((-FC_THR) <= XMIN || XMAX <= FC_THR) {
    stop("log2FC range does not extend past both fold-change thresholds.", call. = FALSE)
}
SL <- ((-FC_THR) - XMIN) / Lmax
SR <- (XMAX - FC_THR) / Rmax

go <- go %>%
    group_by(.data[[dir_col]]) %>%
    mutate(
        row = dplyr::row_number(),
        yc = ROW0 + (n_rows - row) * PITCH,
        ymin = yc - BAR_H / 2,
        ymax = yc + BAR_H / 2,
        xmin = if_else(.data[[dir_col]] == "down", XMIN, FC_THR),
        xmax = if_else(
            .data[[dir_col]] == "down",
            XMIN + .data[[score_col]] * SL,
            FC_THR + .data[[score_col]] * SR
        ),
        fill = if_else(
            .data[[dir_col]] == "down",
            cols$bar_down,
            cols$bar_up
        )
    ) %>%
    ungroup()

brks_L <- pretty(c(0, Lmax), n = 4)
brks_L <- brks_L[brks_L >= 0 & brks_L <= Lmax + 1e-8]
brks_R <- pretty(c(0, Rmax), n = 4)
brks_R <- brks_R[brks_R >= 0 & brks_R <= Rmax + 1e-8]
bar_ticks <- data.frame(
    x = c(XMIN + brks_L * SL, FC_THR + brks_R * SR),
    y = Y_STRIP_TOP - u(6.76),
    label = as.character(c(brks_L, brks_R)),
    stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    annotate(
        "rect",
        xmin = c(XMIN, FC_THR), xmax = c(-FC_THR, XMAX),
        ymin = Y_HLINE, ymax = Y_STRIP_TOP,
        fill = c(cols$zone_down, cols$zone_up)
    ) +
    geom_rect(
        data = go,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
        colour = NA
    ) +
    scale_fill_identity() +
    annotate(
        "segment",
        x = c(XMIN, FC_THR), xend = c(XMIN, FC_THR),
        y = Y_STRIP_TOP, yend = Y_TOP, linewidth = 0.088
    ) +
    annotate(
        "segment",
        x = c(XMIN, FC_THR), xend = c(-FC_THR, XMAX),
        y = Y_STRIP_TOP, yend = Y_STRIP_TOP, linewidth = 0.088
    ) +
    annotate(
        "segment",
        x = bar_ticks$x, xend = bar_ticks$x,
        y = Y_STRIP_TOP, yend = Y_STRIP_TOP - TICKLEN, linewidth = 0.088
    ) +
    geom_text(
        data = bar_ticks,
        aes(x = x, y = y, label = label),
        size = 2.21, colour = "black"
    ) +
    annotate(
        "text",
        x = c((XMIN - FC_THR) / 2, (FC_THR + XMAX) / 2),
        y = Y_VOL_TOP + u(0.73),
        label = config$labels$score, parse = TRUE,
        size = 2.64, colour = "black"
    ) +
    annotate(
        "text",
        x = XMIN - ux(12.5), y = (Y_BOT + Y_VOL_TOP) / 2,
        label = config$labels$y, parse = TRUE, angle = 90,
        size = 2.64, colour = "black"
    ) +
    annotate(
        "text",
        x = c(XMIN - ux(4.81), FC_THR - ux(4.90)),
        y = ROW0 + (n_rows - 1) * PITCH / 2,
        label = config$labels$go_y, angle = 90,
        size = 2.21, colour = "black"
    ) +
    geom_text(
        data = go,
        aes(
            x = if_else(.data[[dir_col]] == "down", XMIN + ux(1.55), FC_THR + ux(1.59)),
            y = yc,
            label = .data[[term_col]]
        ),
        hjust = 0, size = 2.21, lineheight = 0.95, colour = "black"
    ) +
    annotate(
        "segment", x = XMIN, xend = XMIN,
        y = Y_BOT, yend = Y_VOL_TOP, linewidth = 0.088
    ) +
    annotate(
        "segment", x = XMIN, xend = XMAX,
        y = Y_BOT, yend = Y_BOT, linewidth = 0.088
    ) +
    geom_point(
        data = filter(df, group == "ns"),
        aes(log2FC, neg_log10_p),
        colour = cols$ns, size = 1.05
    ) +
    geom_point(
        data = filter(df, group == "ns_fc"),
        aes(log2FC, neg_log10_p),
        colour = cols$ns_fc, size = 1.05
    ) +
    geom_point(
        data = filter(df, group == "sig"),
        aes(log2FC, neg_log10_p),
        colour = cols$sig, size = 1.05
    ) +
    geom_point(
        data = filter(df, group == "deg", log2FC < 0),
        aes(log2FC, neg_log10_p),
        colour = cols$deg_down, size = 1.05
    ) +
    geom_point(
        data = filter(df, group == "deg", log2FC > 0),
        aes(log2FC, neg_log10_p),
        colour = cols$deg_up, size = 1.05
    ) +
    annotate(
        "segment",
        x = c(FC_THR, -FC_THR, XMIN),
        xend = c(FC_THR, -FC_THR, XMAX),
        y = c(Y_BOT, Y_BOT, Y_HLINE),
        yend = c(Y_VOL_TOP, Y_VOL_TOP, Y_HLINE),
        linewidth = 0.176, linetype = "22"
    ) +
    geom_text_repel(
        data = lab_L,
        aes(x = log2FC, y = neg_log10_p, label = gene),
        size = 2.21, fontface = "italic",
        max.overlaps = Inf,
        xlim = c(XMIN - ux(4.9), -FC_THR - ux(7.2)),
        ylim = c(Y_HLINE + u(2), Y_VOL_TOP - u(4.7)),
        box.padding = 0.38, point.padding = 0.22,
        segment.colour = "grey35", segment.size = 0.14,
        min.segment.length = 0, seed = 2023
    ) +
    geom_text_repel(
        data = lab_R,
        aes(x = log2FC, y = neg_log10_p, label = gene),
        size = 2.21, fontface = "italic",
        max.overlaps = Inf,
        xlim = c(FC_THR + ux(7.2), XMAX + ux(9.8)),
        ylim = c(Y_HLINE + u(2), Y_VOL_TOP - u(4.7)),
        box.padding = 0.38, point.padding = 0.22,
        segment.colour = "grey35", segment.size = 0.14,
        min.segment.length = 0, seed = 2023
    ) +
    scale_x_continuous(breaks = brks_x, expand = expansion(0)) +
    scale_y_continuous(breaks = brks_y, expand = expansion(0)) +
    coord_fixed(
        ratio = ASPECT, xlim = c(XMIN, XMAX), ylim = c(Y_BOT, Y_TOP),
        expand = FALSE, clip = "off"
    ) +
    labs(x = parse(text = config$labels$x)[[1]], y = NULL) +
    theme_classic() +
    theme(
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.088, colour = "black"),
        axis.ticks.length = unit(2.58, "pt"),
        axis.text = element_text(size = 6.3, colour = "black"),
        axis.title = element_text(size = 7.5, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 2)),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 18, 8, 18, "pt")
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 3.4, height = 4.15)
