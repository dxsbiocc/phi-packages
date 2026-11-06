#!/usr/bin/env Rscript

# Template-ID: bar-butterfly-stack
#
# Purpose:
#   Draw a diverging stacked horizontal bar chart: two opposing sides
#   (e.g. TRAEs left, AEs irrespective of causality right) share a
#   category axis; each side stacks supplied severity grades from the
#   centre outward, with end totals and dual legends (preview.png).
#
# Inputs:
#   One row per term-by-side-by-grade. Default example:
#     - term, side, grade, percent
#   percent is non-negative; left/right assignment is by side levels.
#
# Output:
#   A PDF, PNG, or SVG diverging stacked butterfly bar chart.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggnewscale, scales, grid
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, side/grade order, palettes, x limits,
#   labels, and size. Extra terms or grades are extra rows/levels,
#   not a new id. A non-stacked two-group butterfly stays
#   bar-butterfly; tissue communication bars stay bar-communication.
#
# Scientific assumptions:
#   percent, side, grade, and term are supplied display values. This
#   script does not count patients, assign causality, grade events,
#   or select a top-N list. Negating the left side is a display
#   transform only.

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
        term = "term",
        side = "side",
        grade = "grade",
        value = "percent"
    ),
    # First side is drawn to the left of zero.
    side_order = c("TRAE", "AE"),
    side_labels = c(
        TRAE = "TRAEs",
        AE = "AEs irrespective of causality"
    ),
    # Innermost (nearest zero) -> outermost.
    grade_order = c("Grade 1", "Grade 2", "Grade 3", "Grade 4"),
    palettes = list(
        left = "Qualitative.Bold",
        right = "Qualitative.Safe"
    ),
    # 1-based picks into each palette (paper-like contrast).
    palette_index = list(
        left = c(3, 5, 2, 10),
        right = c(6, 8, 1, 9)
    ),
    x_limits = c(-40, 40),
    x_breaks = seq(-40, 40, by = 5),
    bar_height = 0.72,
    label_digits = 1,
    labels = list(
        title = NULL,
        x = "Percentage of patients with TRAEs and AEs irrespective of causality (n = 84)",
        y = NULL,
        caption = paste(
            "Top 10 TRAEs versus AEs irrespective of causality by severity.",
            "a Grade 4 TRAEs include thrombocytopenia (n = 1) and ALT increased (n = 1).",
            "b Includes MedDRA preferred terms of rash and maculopapular rash.",
            "c Includes MedDRA preferred terms of pruritus and generalized pruritus.",
            "At each level of patient summarization, a patient was counted once,",
            "if the patient reported one or more events.",
            "AE, adverse event; ALT, alanine aminotransferase;",
            "AST, aspartate aminotransferase; SAE, serious adverse event.",
            "Data cutoff: 1 December 2021."
        )
    ),
    size = list(width = 13.5, height = 7.2)
)

load_packages(c("ggplot2", "readr", "ggprism", "ggnewscale", "scales", "grid"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
term_col <- config$columns$term
side_col <- config$columns$side
grade_col <- config$columns$grade
val_col <- config$columns$value

df[[val_col]] <- as.numeric(df[[val_col]])
df <- df[is.finite(df[[val_col]]) & df[[val_col]] >= 0, , drop = FALSE]
if (!nrow(df)) stop("No finite non-negative percent rows.", call. = FALSE)

side_lv <- config$side_order
side_lv <- c(
    intersect(side_lv, unique(as.character(df[[side_col]]))),
    setdiff(unique(as.character(df[[side_col]])), side_lv)
)
if (length(side_lv) != 2L) {
    stop("Exactly two side levels are required (left / right).", call. = FALSE)
}
grade_lv <- config$grade_order
grade_lv <- c(
    intersect(grade_lv, unique(as.character(df[[grade_col]]))),
    setdiff(unique(as.character(df[[grade_col]])), grade_lv)
)

df[[side_col]] <- factor(as.character(df[[side_col]]), levels = side_lv)
df[[grade_col]] <- factor(as.character(df[[grade_col]]), levels = grade_lv)

# Term order: descending total on the left side (first side_order entry).
left_side <- side_lv[[1]]
right_side <- side_lv[[2]]
tot_left <- tapply(
    df[[val_col]][df[[side_col]] == left_side],
    as.character(df[[term_col]][df[[side_col]] == left_side]),
    sum
)
term_lv <- names(sort(tot_left, decreasing = TRUE))
extra_terms <- setdiff(unique(as.character(df[[term_col]])), term_lv)
term_lv <- c(term_lv, extra_terms)
df[[term_col]] <- factor(as.character(df[[term_col]]), levels = rev(term_lv))

# Stack from centre outward: left negative, right positive.
stack_side <- function(d, sign) {
    d <- d[order(d[[term_col]], d[[grade_col]]), , drop = FALSE]
    pieces <- lapply(
        split(d, d[[term_col]], drop = TRUE),
        function(g) {
            g <- g[order(g[[grade_col]]), , drop = FALSE]
            w <- g[[val_col]]
            if (sign < 0) {
                xmax <- -c(0, cumsum(w)[-length(w)])
                xmin <- -cumsum(w)
            } else {
                xmin <- c(0, cumsum(w)[-length(w)])
                xmax <- cumsum(w)
            }
            data.frame(
                term = g[[term_col]],
                grade = g[[grade_col]],
                xmin = xmin,
                xmax = xmax,
                ymin = as.numeric(g[[term_col]]) - config$bar_height / 2,
                ymax = as.numeric(g[[term_col]]) + config$bar_height / 2,
                stringsAsFactors = FALSE
            )
        }
    )
    out <- do.call(rbind, pieces)
    rownames(out) <- NULL
    out
}

df_left <- df[df[[side_col]] == left_side, , drop = FALSE]
df_right <- df[df[[side_col]] == right_side, , drop = FALSE]
seg_l <- stack_side(df_left, -1)
seg_r <- stack_side(df_right, 1)
seg_l$grade <- factor(as.character(seg_l$grade), levels = grade_lv)
seg_r$grade <- factor(as.character(seg_r$grade), levels = grade_lv)

tot <- aggregate(
    df[[val_col]],
    by = list(term = df[[term_col]], side = df[[side_col]]),
    FUN = sum
)
names(tot)[3] <- "total"
tot$xend <- ifelse(tot$side == left_side, -tot$total, tot$total)
tot$y <- as.numeric(tot$term)
tot$hjust <- ifelse(tot$side == left_side, 1.15, -0.15)
tot$label <- paste0(
    formatC(tot$total, format = "f", digits = config$label_digits),
    "%"
)

pick_cols <- function(pal_id, idx, n) {
    base <- palette_colors(pal_id)
    idx <- as.integer(idx)
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= length(base)]
    cols <- base[idx]
    if (length(cols) < n) {
        cols <- expand_palette(cols, n)
    } else {
        cols <- cols[seq_len(n)]
    }
    cols
}

n_g <- length(grade_lv)
left_cols <- pick_cols(
    config$palettes$left, config$palette_index$left, n_g
)
right_cols <- pick_cols(
    config$palettes$right, config$palette_index$right, n_g
)
names(left_cols) <- grade_lv
names(right_cols) <- grade_lv

left_title <- if (left_side %in% names(config$side_labels)) {
    unname(config$side_labels[[left_side]])
} else {
    left_side
}
right_title <- if (right_side %in% names(config$side_labels)) {
    unname(config$side_labels[[right_side]])
} else {
    right_side
}

x_breaks <- config$x_breaks
x_labs <- paste0(abs(x_breaks), "%")
y_num <- seq_along(levels(df[[term_col]]))
y_labs <- levels(df[[term_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
extract_legend <- function(pl) {
    g <- ggplotGrob(pl)
    idx <- which(g$layout$name == "guide-box-bottom")
    if (!length(idx)) {
        idx <- which(grepl("^guide-box", g$layout$name))
    }
    for (i in idx) {
        gr <- g$grobs[[i]]
        if (!inherits(gr, "zeroGrob")) {
            return(gr)
        }
    }
    stop("Could not extract legend grob.", call. = FALSE)
}

legend_theme <- theme_prism(base_size = 11) +
    theme(
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = rel(0.85)),
        legend.text = element_text(size = rel(0.8)),
        legend.key.width = unit(0.9, "lines"),
        legend.key.height = unit(0.7, "lines"),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.box.background = element_rect(fill = "white", colour = NA)
    )

leg_data <- data.frame(
    grade = factor(grade_lv, levels = grade_lv),
    x = seq_along(grade_lv),
    y = 1,
    stringsAsFactors = FALSE
)
p_leg_l <- ggplot(leg_data, aes(x = x, y = y, fill = grade)) +
    geom_tile() +
    scale_fill_manual(
        name = left_title,
        values = left_cols,
        breaks = grade_lv,
        drop = FALSE,
        guide = guide_legend(nrow = 2, byrow = TRUE, title.position = "top")
    ) +
    legend_theme
p_leg_r <- ggplot(leg_data, aes(x = x, y = y, fill = grade)) +
    geom_tile() +
    scale_fill_manual(
        name = right_title,
        values = right_cols,
        breaks = grade_lv,
        drop = FALSE,
        guide = guide_legend(nrow = 2, byrow = TRUE, title.position = "top")
    ) +
    legend_theme
leg_l <- extract_legend(p_leg_l)
leg_r <- extract_legend(p_leg_r)

x_lim <- config$x_limits

p <- ggplot() +
    geom_vline(
        xintercept = config$x_breaks,
        colour = "grey88",
        linewidth = 0.35
    ) +
    geom_vline(
        xintercept = 0,
        colour = "grey35",
        linewidth = 0.55
    ) +
    geom_rect(
        data = seg_l,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = grade),
        colour = NA
    ) +
    scale_fill_manual(
        name = left_title,
        values = left_cols,
        breaks = grade_lv,
        drop = FALSE,
        guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    geom_rect(
        data = seg_r,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = grade),
        colour = NA
    ) +
    scale_fill_manual(
        name = right_title,
        values = right_cols,
        breaks = grade_lv,
        drop = FALSE,
        guide = "none"
    ) +
    geom_text(
        data = tot,
        aes(x = xend, y = y, label = label, hjust = hjust),
        size = 2.7,
        colour = "grey20"
    ) +
    scale_x_continuous(
        limits = x_lim,
        breaks = x_breaks,
        labels = x_labs,
        expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
        breaks = y_num,
        labels = y_labs,
        expand = expansion(add = c(0.35, 0.35))
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        caption = config$labels$caption
    ) +
    coord_cartesian(clip = "off") +
    theme_prism(base_size = 11) +
    theme(
        panel.grid = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = rel(0.85), hjust = 1),
        axis.text.x = element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.85), margin = margin(t = 8)),
        legend.position = "none",
        plot.caption = element_text(
            hjust = 0, size = rel(0.62), colour = "grey30",
            margin = margin(t = 8), lineheight = 1.15
        ),
        plot.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.4),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 18, 8, 10)
    )

# Pin legends into the panel npc blanks beside short bottom bars.
gt <- ggplotGrob(p)
panel_i <- which(gt$layout$name == "panel")[[1]]
leg_l_pinned <- grobTree(
    leg_l,
    vp = viewport(
        x = unit(0.20, "npc"), y = unit(0.04, "npc"),
        width = unit(0.30, "npc"), height = unit(0.13, "npc"),
        just = c("center", "bottom")
    )
)
leg_r_pinned <- grobTree(
    leg_r,
    vp = viewport(
        x = unit(0.80, "npc"), y = unit(0.04, "npc"),
        width = unit(0.30, "npc"), height = unit(0.13, "npc"),
        just = c("center", "bottom")
    )
)
gt$grobs[[panel_i]] <- grobTree(
    gt$grobs[[panel_i]],
    leg_l_pinned,
    leg_r_pinned
)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    gt,
    io$output,
    width = config$size$width,
    height = config$size$height
)
