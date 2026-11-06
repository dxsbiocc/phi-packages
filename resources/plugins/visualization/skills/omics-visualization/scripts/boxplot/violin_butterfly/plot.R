#!/usr/bin/env Rscript

# Template-ID: boxplot-violin-butterfly
#
# Purpose:
#   Mirrored horizontal violins: one row per group, the left panel
#   shows one score (reversed axis), the right panel another, with
#   the group names in a shared centre column coloured by a supplied
#   class, nested white boxes, a dashed reference line per side and
#   coloured side titles.
#
# Inputs:
#   example.tsv - one row per observation (cell): group, side, score
#                 (side has exactly the two CONFIG sides)
#   groups.tsv  - optional sidecar beside the input: group, class
#                 (label colour; missing class = "Other")
#
# Output:
#   A PDF, PNG, or SVG butterfly violin plot.
#
# Dependencies:
#   ggplot2, readr, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, which side is left, side titles,
#   reference lines, shared x range, row order, box on/off, notch,
#   class colours, row palette). Extra groups are extra rows. Two
#   scores per group is this id; one score per group is
#   boxplot-violin. Do not compute M1 / M2 or module scores here.
#   Do not add an id for notched nested boxes alone.
#
# Scientific assumptions:
#   Scores and the reference values are supplied. Violins estimate
#   the density of the supplied values per group and side. Row
#   order is display only (file order or a median sort). Notches
#   are a geom_boxplot display of the median CI, not a test.

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
    columns = list(group = "group", side = "side", score = "score"),
    sidecars = list(groups = "groups.tsv"),   # NULL = all labels black
    sides = c(left = "M1", right = "M2"),
    side_titles = c(left = "M1 Feature", right = "M2 Feature"),
    reference = c(left = 0.25, right = 0.25),  # NULL = no dashed line
    x_limits = c(0, 1.5),         # NULL = shared range from both sides
    order = "file",              # "file", "left_median", "right_median"
    show_box = TRUE,
    notch = TRUE,                # nested box median CI notch
    notchwidth = 0.5,
    box_width = 0.18,
    axis_title = "Score",
    # "Family.Name:index" catalog references; unlisted classes are black
    class_colors = c(M1 = "Qualitative.Vivid:1", M2 = "Qualitative.Vivid:8"),
    side_colors = c(left = "Qualitative.Vivid:1", right = "Qualitative.Vivid:8"),
    row_palette = "Qualitative.Pastel",
    reference_color = "Diverging.RdBu:1",
    label_size = 4.2,
    size = list(width = 9, height = 6.8)
)

load_packages(c("ggplot2", "readr", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pref <- function(ref) {
    if (!grepl(".", ref, fixed = TRUE)) return(ref)
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

cc <- config$columns
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, cc)
df <- data.frame(group = as.character(df[[cc$group]]),
    side = as.character(df[[cc$side]]),
    score = as.numeric(df[[cc$score]]), stringsAsFactors = FALSE)
if (any(!is.finite(df$score))) {
    stop("score must be finite.", call. = FALSE)
}
sides <- config$sides
if (!all(df$side %in% sides)) {
    stop(sprintf("side must be one of: %s",
        paste(sides, collapse = ", ")), call. = FALSE)
}

grp <- unique(df$group)
if (config$order %in% c("left_median", "right_median")) {
    s <- if (config$order == "left_median") sides[["left"]] else
        sides[["right"]]
    d <- df[df$side == s, , drop = FALSE]
    md <- tapply(d$score, d$group, stats::median)
    grp <- names(sort(md[grp]))
}
df$group <- factor(df$group, levels = rev(grp))

cls <- stats::setNames(rep("Other", length(grp)), grp)
fn <- config$sidecars$groups
if (!is.null(fn)) {
    path <- file.path(dirname(io$input), fn)
    if (file.exists(path)) {
        gt <- as.data.frame(read_table_auto(path), stringsAsFactors = FALSE)
        require_columns(gt, list(group = "group", class = "class"))
        hit <- as.character(gt$group) %in% grp
        cls[as.character(gt$group)[hit]] <- as.character(gt$class)[hit]
    }
}
lab_col <- vapply(cls, function(k) {
    if (k %in% names(config$class_colors)) pref(config$class_colors[[k]])
    else "black"
}, character(1))

row_fill <- stats::setNames(palette_colors(config$row_palette,
    n = length(grp)), grp)
xl <- if (is.null(config$x_limits)) range(c(0, df$score)) else
    as.numeric(config$x_limits)
brk <- pretty(xl, 3)
brk <- brk[brk >= xl[[1]] & brk <= xl[[2]]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
th <- theme_classic(base_size = 12) +
    theme(
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5,
            size = 12, colour = "black"),
        plot.caption = element_text(hjust = 0.5, size = 15,
            margin = margin(t = 8)),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"
    )

side_panel <- function(key) {
    d <- df[df$side == sides[[key]], , drop = FALSE]
    p <- ggplot(d, aes(x = score, y = group, fill = group)) +
        geom_violin(colour = "grey20", linewidth = 0.35, scale = "width",
            width = 0.85, trim = TRUE)
    if (isTRUE(config$show_box)) {
        p <- p + geom_boxplot(
            width = config$box_width,
            notch = isTRUE(config$notch),
            notchwidth = config$notchwidth,
            fill = "white",
            colour = "grey20",
            linewidth = 0.35,
            outlier.shape = NA
        )
    }
    ref <- config$reference
    if (!is.null(ref) && !is.null(ref[[key]])) {
        p <- p + geom_vline(xintercept = ref[[key]], linetype = "dashed",
            colour = pref(config$reference_color), linewidth = 0.7)
    }
    xs <- if (key == "left") {
        scale_x_continuous(
            transform = "reverse",
            limits = xl,
            breaks = brk,
            position = "top",
            expand = expansion(mult = c(0.02, 0.02))
        )
    } else {
        scale_x_continuous(
            limits = xl,
            breaks = brk,
            position = "top",
            expand = expansion(mult = c(0.02, 0.02))
        )
    }
    p + xs +
        scale_y_discrete(position = if (key == "left") "right" else "left",
            drop = FALSE) +
        scale_fill_manual(values = row_fill) +
        labs(caption = config$side_titles[[key]]) +
        th +
        theme(plot.caption = element_text(colour = pref(
            config$side_colors[[key]])))
}

lab_df <- data.frame(group = factor(grp, levels = rev(grp)), lab = grp,
    col = unname(lab_col[grp]))
mid <- ggplot(lab_df, aes(x = 0, y = group, label = lab, colour = col)) +
    geom_text(size = config$label_size) +
    scale_colour_identity() +
    scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE) +
    labs(title = config$axis_title, caption = " ") +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 15,
            margin = margin(b = 14)),
        plot.caption = element_text(size = 15, margin = margin(t = 8)))

p <- patchwork::wrap_plots(side_panel("left"), mid, side_panel("right"),
    nrow = 1, widths = c(1, 0.62, 1))

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width,
    height = config$size$height)
