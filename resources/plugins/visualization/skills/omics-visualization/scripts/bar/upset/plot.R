#!/usr/bin/env Rscript

# Template-ID: bar-upset
#
# Purpose:
#   Draw an UpSet combination figure. style = "dual" (default): top
#   bars of intersection size (observed vs expected), membership
#   matrix coloured by n-sets, and left set-size bars.
#   style = "grouped": single count bars coloured by a supplied track
#   group (e.g. T1–T7), count/percent labels, group-coloured matrix
#   dots and connectors, alternating row stripes, bottom group labels
#   (preview_grouped.png).
#
# Inputs:
#   dual (example.tsv): sets, observed, expected
#   grouped (example_grouped.tsv): sets, count, percent, group
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
#   Edit CONFIG for columns, set_order, group colours, and style.
#   Do not enumerate combinations from raw membership. Do not add an
#   id for grouped colouring alone. Venn circles stay bar-venn;
#   flower diagrams stay bar-flower.
#
# Scientific assumptions:
#   Each row is one already-counted combination. expected / percent /
#   group are supplied. Set-size bars (dual) aggregate supplied
#   observed counts for display only.

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
    # "dual": observed vs expected + left set sizes (default).
    # "grouped": track-coloured bars / matrix (T1–T7 style).
    style = "dual",
    columns = list(
        sets = "sets",
        observed = "observed",
        expected = "expected",
        count = "count",
        percent = "percent",
        group = "group"
    ),
    set_order = c("PFC", "NAc", "BNST", "LH", "mPFC", "BA", "CeA"),
    group_order = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "Other"),
    order_by = "freq",           # "file" | "freq"
    decreasing = TRUE,
    # Soft bar fills then saturated matrix / label colours (catalog refs).
    palettes = list(
        group_soft = c(
            "Qualitative.Pastel:2",
            "Qualitative.Pastel:5",
            "Qualitative.Pastel:7",
            "Qualitative.Pastel:4",
            "Qualitative.Pastel:1",
            "Qualitative.Pastel:6",
            "Qualitative.Pastel:3",
            "Qualitative.Safe:12"
        ),
        group_hard = c(
            "Qualitative.Bold:4",
            "Qualitative.Bold:6",
            "Qualitative.Bold:5",
            "Qualitative.Bold:1",
            "Qualitative.Bold:8",
            "Qualitative.Safe:5",
            "Qualitative.Bold:7",
            "Qualitative.Safe:12"
        )
    ),
    grouped = list(
        bar_width = 0.72,
        point_size = 3.4,
        absent_colour = "grey82",
        stripe_colour = "grey94",
        label_size = 2.6,
        percent_size = 2.3,
        show_percent = TRUE,
        show_bottom_groups = TRUE,
        size = list(width = 12.5, height = 7.2)
    ),
    labels = list(
        title = "",
        top_y = "Intersection size",
        set_y = "Set size",
        grouped_y = "Number of NDDs",
        fill = NULL,
        colour = NULL
    )
)

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
load_packages(c("ggplot2", "readr", "ggprism", "aplot", "patchwork"))

pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols[[1]]
}

split_sets <- function(x) {
    parts <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
    parts[nzchar(parts)]
}

# -----------------------------------------------------------------------------
# GROUPED STYLE
# -----------------------------------------------------------------------------
if (identical(config$style, "grouped")) {
    df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
    sets_col <- config$columns$sets
    cnt_col <- config$columns$count
    pct_col <- config$columns$percent
    g_col <- config$columns$group
    require_columns(df, list(sets_col, cnt_col, g_col))
    df[[cnt_col]] <- as.numeric(df[[cnt_col]])
    df[[g_col]] <- as.character(df[[g_col]])
    if (pct_col %in% names(df)) {
        df[[pct_col]] <- as.numeric(df[[pct_col]])
    } else {
        tot <- sum(df[[cnt_col]], na.rm = TRUE)
        df[[pct_col]] <- if (tot > 0) round(100 * df[[cnt_col]] / tot) else NA_real_
    }
    ok <- !is.na(df[[sets_col]]) & nzchar(as.character(df[[sets_col]])) &
        is.finite(df[[cnt_col]])
    df <- df[ok, , drop = FALSE]
    if (!nrow(df)) stop("No combinations left to plot.", call. = FALSE)

    members <- lapply(df[[sets_col]], split_sets)
    n_mem <- vapply(members, length, integer(1))
    if (identical(config$order_by, "freq")) {
        o <- order(n_mem, -df[[cnt_col]])
        df <- df[o, , drop = FALSE]
        members <- members[o]
        n_mem <- n_mem[o]
    }

    # Preserve file / group_order block order.
    g_lv <- unique(as.character(df[[g_col]]))
    if (!is.null(config$group_order)) {
        g_lv <- c(
            intersect(as.character(config$group_order), g_lv),
            setdiff(g_lv, as.character(config$group_order))
        )
    }
    df[[g_col]] <- factor(df[[g_col]], levels = g_lv)

    set_order <- unique(as.character(config$set_order))
    seen <- unique(unlist(members, use.names = FALSE))
    set_order <- c(intersect(set_order, seen), setdiff(seen, set_order))
    n_set <- length(set_order)
    n_combo <- nrow(df)
    df$combo <- seq_len(n_combo)
    df$n_mem <- n_mem

    soft <- vapply(seq_along(g_lv), function(i) {
        refs <- config$palettes$group_soft
        pal_ref(refs[[((i - 1L) %% length(refs)) + 1L]])
    }, character(1))
    hard <- vapply(seq_along(g_lv), function(i) {
        refs <- config$palettes$group_hard
        pal_ref(refs[[((i - 1L) %% length(refs)) + 1L]])
    }, character(1))
    names(soft) <- g_lv
    names(hard) <- g_lv
    # Prefer near-black for "Other"
    if ("Other" %in% g_lv) {
        soft[["Other"]] <- "#4d4d4d"
        hard[["Other"]] <- "#222222"
    }

    mat <- data.frame(
        X = rep(df$combo, each = n_set),
        Y = rep(seq_along(set_order), times = n_combo),
        present = FALSE,
        group = rep(as.character(df[[g_col]]), each = n_set),
        stringsAsFactors = FALSE
    )
    for (j in seq_len(n_combo)) {
        hit <- set_order %in% members[[j]]
        idx <- (j - 1L) * n_set + seq_len(n_set)
        mat$present[idx] <- hit
    }
    mat$group <- factor(mat$group, levels = g_lv)

    seg <- data.frame(
        X = df$combo,
        ystart = NA_real_,
        yend = NA_real_,
        group = df[[g_col]],
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

    gc <- config$grouped
    stripe <- data.frame(
        ymin = seq_along(set_order) - 0.5,
        ymax = seq_along(set_order) + 0.5,
        fill = ifelse(seq_along(set_order) %% 2L == 1L, gc$stripe_colour, "white"),
        stringsAsFactors = FALSE
    )

    # Bottom group spans
    grp_lab <- NULL
    if (isTRUE(gc$show_bottom_groups)) {
        spans <- lapply(g_lv, function(g) {
            ix <- which(as.character(df[[g_col]]) == g)
            if (!length(ix)) return(NULL)
            data.frame(
                group = g,
                xmin = min(ix) - 0.45,
                xmax = max(ix) + 0.45,
                x = mean(range(ix)),
                stringsAsFactors = FALSE
            )
        })
        grp_lab <- do.call(rbind, spans)
        if (!is.null(grp_lab)) {
            grp_lab$group <- factor(grp_lab$group, levels = g_lv)
        }
    }

    blank_bg <- theme(
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "none"
    )

    y_top <- max(df[[cnt_col]], na.rm = TRUE)
    p_top <- ggplot(df, aes(x = .data$combo, y = .data[[cnt_col]], fill = .data[[g_col]])) +
        geom_col(width = gc$bar_width, colour = NA) +
        scale_fill_manual(values = soft, guide = "none") +
        scale_x_continuous(limits = c(0.5, n_combo + 0.5), expand = c(0, 0)) +
        scale_y_continuous(
            expand = expansion(mult = c(0, 0.22)),
            limits = c(0, y_top * 1.28)
        ) +
        labs(x = NULL, y = config$labels$grouped_y) +
        theme_classic(base_size = 12) +
        blank_bg +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.text = element_text(colour = "black"),
            plot.margin = margin(6, 8, 0, 4)
        )
    if (isTRUE(gc$show_percent)) {
        # Percent sits just above the bar; count sits above the percent.
        p_top <- p_top +
            geom_text(
                aes(label = paste0(.data[[pct_col]], "%")),
                vjust = -0.25, size = gc$percent_size, colour = "grey35"
            ) +
            geom_text(
                aes(label = .data[[cnt_col]]),
                vjust = -1.45, size = gc$label_size, colour = "grey15"
            )
    } else {
        p_top <- p_top + geom_text(
            aes(label = .data[[cnt_col]]),
            vjust = -0.35, size = gc$label_size, colour = "grey15"
        )
    }

    # Layer order: stripes → absent grey dots → connectors → present dots
    # so combination lines sit above the grey matrix points.
    p_mat <- ggplot() +
        geom_rect(
            data = stripe,
            aes(xmin = 0.5, xmax = n_combo + 0.5, ymin = .data$ymin, ymax = .data$ymax),
            fill = stripe$fill, colour = NA, inherit.aes = FALSE
        ) +
        geom_point(
            data = mat[!mat$present, , drop = FALSE],
            aes(x = .data$X, y = .data$Y),
            colour = gc$absent_colour, size = gc$point_size * 0.85
        ) +
        geom_segment(
            data = seg,
            aes(
                x = .data$X, xend = .data$X,
                y = .data$ystart, yend = .data$yend,
                colour = .data$group
            ),
            linewidth = 0.9, inherit.aes = FALSE
        ) +
        geom_point(
            data = mat[mat$present, , drop = FALSE],
            aes(x = .data$X, y = .data$Y, colour = .data$group),
            size = gc$point_size
        ) +
        scale_colour_manual(values = hard, guide = "none") +
        scale_x_continuous(limits = c(0.5, n_combo + 0.5), expand = c(0, 0)) +
        scale_y_continuous(
            breaks = seq_along(set_order),
            labels = set_order,
            limits = c(0.5, n_set + 0.5),
            expand = c(0, 0),
            trans = "reverse"
        ) +
        labs(x = NULL, y = NULL) +
        theme_classic(base_size = 11) +
        blank_bg +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(colour = "grey70", linewidth = 0.45),
            axis.ticks.y = element_line(colour = "grey70", linewidth = 0.4),
            axis.text.y = element_text(colour = "black", size = 9),
            plot.margin = margin(2, 8, 18, 4)
        )

    if (!is.null(grp_lab) && nrow(grp_lab)) {
        # With y reversed, larger numeric y is toward the visual bottom.
        p_mat <- p_mat +
            geom_label(
                data = grp_lab,
                aes(
                    x = .data$x, y = n_set + 0.85,
                    label = .data$group, fill = .data$group
                ),
                inherit.aes = FALSE,
                colour = "grey15",
                size = 3.0,
                fontface = "bold",
                linewidth = 0,
                label.padding = unit(0.18, "lines")
            ) +
            scale_fill_manual(values = soft, guide = "none") +
            coord_cartesian(ylim = c(n_set + 1.35, 0.4), clip = "off")
    } else {
        p_mat <- p_mat + coord_cartesian(clip = "off")
    }

    p <- insert_top(p_mat, p_top, height = 1.55)
    p <- as.patchwork(p) +
        plot_annotation(
            theme = theme(
                plot.background = element_rect(fill = "white", colour = NA)
            )
        )

    suppressWarnings(save_ggplot(
        p, io$output,
        width = gc$size$width,
        height = gc$size$height
    ))
    quit(save = "no")
}

# -----------------------------------------------------------------------------
# DUAL STYLE (default)
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, list(
    config$columns$sets,
    config$columns$observed,
    config$columns$expected
))

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

members <- lapply(df[[sets_col]], split_sets)
n_mem <- vapply(members, length, integer(1))
if (identical(config$order_by, "freq") || identical(config$order_by, TRUE)) {
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

# Default dual style historically ordered by frequency.
if (!identical(config$order_by, "file")) {
    # already reordered above when freq
}

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
