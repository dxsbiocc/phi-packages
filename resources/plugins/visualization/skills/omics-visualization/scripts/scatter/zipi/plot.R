#!/usr/bin/env Rscript

# Template-ID: scatter-zipi
#
# Purpose:
#   MetaNet-style Zi–Pi topological-role scatter. Four shaded role
#   regions (Peripheral / Connector / Module hub / Network hub) with
#   threshold lines at Zi and Pi cutoffs. Points are coloured by a
#   supplied group (omics layer, taxon class, …), not by role.
#
# Inputs:
#   A long table (tsv/csv). Default example:
#     - node: node id / label
#     - zi: within-module degree z-score
#     - pi: participation coefficient
#     - group: point colour (e.g. Microbiome / Metabolome / Transcriptome)
#     - role: optional; else derived from thresholds (display only)
#     - show_label: optional 1/0; else label non-Peripheral points
#
# Output:
#   A PDF, PNG, or SVG Zi–Pi role scatter.
#
# Dependencies:
#   ggplot2, readr, ggrepel
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, thresholds, region fills, labels).
#
# Scientific assumptions:
#   zi and pi are supplied. This script does not run module detection,
#   zp_analyse, or MetaNet. Role from thresholds is display-only.

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
        node = "node",
        zi = "zi",
        pi = "pi",
        group = "group",
        role = "role",
        show_label = "show_label"
    ),
    labels = list(
        title = "Topological roles",
        x = "Participation Coefficient",
        y = "Within-module connectivity",
        role = "Topological roles",
        group = NULL
    ),
    thresholds = list(
        zi = 2.5,
        pi = 0.62
    ),
    # MetaNet zp_plot-style region fills (background).
    role_fills = c(
        Peripheral = "#F5F0E6",
        Connector = "#B8D4E8",
        "Module hub" = "#C5DDB0",
        "Network hub" = "#E8A0A0"
    ),
    role_levels = c("Peripheral", "Module hub", "Connector", "Network hub"),
    role_label_text = c(
        Peripheral = "Peripherals",
        Connector = "Connectors",
        "Module hub" = "Module hubs",
        "Network hub" = "Network hubs"
    ),
    palettes = list(
        group = "Qualitative.Safe"
    ),
    group_colours = c(
        Microbiome = "#6BAED6",
        Genes = "#6BAED6",
        Metabolome = "#FD8D3C",
        Transcriptome = "#E7298A",
        Bacteria = "#6BAED6",
        Fungi = "#FD8D3C"
    ),
    point_size = 2.2,
    label_size = 2.5,
    # Expand plot limits slightly past data / cuts.
    pi_max = 1.0,
    zi_pad = 0.35
)

load_packages(c("ggplot2", "readr", "ggrepel"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns[c("node", "zi", "pi")])
node_col <- config$columns$node
zi_col <- config$columns$zi
pi_col <- config$columns$pi
group_col <- config$columns$group
role_col <- config$columns$role
lab_col <- config$columns$show_label

df[[node_col]] <- as.character(df[[node_col]])
df[[zi_col]] <- as.numeric(df[[zi_col]])
df[[pi_col]] <- as.numeric(df[[pi_col]])
ok <- is.finite(df[[zi_col]]) & is.finite(df[[pi_col]]) &
    !is.na(df[[node_col]]) & nzchar(df[[node_col]])
df <- df[ok, , drop = FALSE]
if (!nrow(df)) {
    stop("No finite Zi–Pi rows to plot.", call. = FALSE)
}

zi_cut <- as.numeric(config$thresholds$zi)
pi_cut <- as.numeric(config$thresholds$pi)
derive_role <- function(zi, pi) {
    ifelse(
        zi >= zi_cut & pi >= pi_cut, "Network hub",
        ifelse(
            zi >= zi_cut & pi < pi_cut, "Module hub",
            ifelse(zi < zi_cut & pi >= pi_cut, "Connector", "Peripheral")
        )
    )
}

if (!is.null(role_col) && role_col %in% names(df)) {
    df[[role_col]] <- as.character(df[[role_col]])
    miss <- is.na(df[[role_col]]) | !nzchar(df[[role_col]])
    df[[role_col]][miss] <- derive_role(df[[zi_col]][miss], df[[pi_col]][miss])
} else {
    role_col <- ".role"
    df[[role_col]] <- derive_role(df[[zi_col]], df[[pi_col]])
}
df[[role_col]] <- factor(df[[role_col]], levels = config$role_levels)

has_group <- !is.null(group_col) && group_col %in% names(df)
if (!has_group && "module" %in% names(df)) {
    group_col <- "module"
    has_group <- TRUE
}
if (has_group) {
    df[[group_col]] <- as.character(df[[group_col]])
    df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
} else {
    group_col <- ".group"
    df[[group_col]] <- factor("Node")
}

df$`.label` <- ""
if (!is.null(lab_col) && lab_col %in% names(df)) {
    show <- as.integer(df[[lab_col]])
    show[is.na(show)] <- 0L
    df$`.label`[show > 0L] <- df[[node_col]][show > 0L]
} else {
    df$`.label`[df[[role_col]] != "Peripheral"] <-
        df[[node_col]][df[[role_col]] != "Peripheral"]
}

group_levels <- levels(df[[group_col]])
gcols <- config$group_colours
gcols <- gcols[names(gcols) %in% group_levels]
miss_g <- setdiff(group_levels, names(gcols))
if (length(miss_g)) {
    extra <- palette_colors(config$palettes$group, n = length(miss_g))
    names(extra) <- miss_g
    gcols <- c(gcols, extra)
}
gcols <- gcols[group_levels]

# Axis limits.
x_min <- 0
x_max <- max(as.numeric(config$pi_max), max(df[[pi_col]], na.rm = TRUE) * 1.02)
y_min <- min(df[[zi_col]], na.rm = TRUE) - as.numeric(config$zi_pad)
y_max <- max(df[[zi_col]], zi_cut, na.rm = TRUE) + as.numeric(config$zi_pad)

# Background rectangles (fixed fills — not mapped, so they cannot stain
# the point colour legend keys).
role_fills <- config$role_fills[config$role_levels]
rects <- data.frame(
    xmin = c(x_min, x_min, pi_cut, pi_cut),
    xmax = c(pi_cut, pi_cut, x_max, x_max),
    ymin = c(y_min, zi_cut, y_min, zi_cut),
    ymax = c(zi_cut, y_max, zi_cut, y_max),
    role = factor(
        c("Peripheral", "Module hub", "Connector", "Network hub"),
        levels = config$role_levels
    ),
    fill = unname(role_fills[c("Peripheral", "Module hub", "Connector", "Network hub")]),
    stringsAsFactors = FALSE
)

# Invisible points that only drive the role (square) legend.
role_text <- config$role_label_text
role_leg <- data.frame(
    role = factor(config$role_levels, levels = config$role_levels),
    x = NA_real_,
    y = NA_real_,
    stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
lab_df <- df[nzchar(df$`.label`), , drop = FALSE]
show_group_leg <- isTRUE(nlevels(df[[group_col]]) > 1L)

p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[pi_col]], y = .data[[zi_col]])) +
    ggplot2::geom_rect(
        data = rects,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = rects$fill,
        colour = NA,
        alpha = 0.95,
        show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = zi_cut, colour = "grey25", linewidth = 0.45) +
    ggplot2::geom_vline(xintercept = pi_cut, colour = "grey25", linewidth = 0.45) +
    ggplot2::geom_point(
        data = role_leg,
        ggplot2::aes(x = x, y = y, fill = role),
        shape = 22,
        size = 4.5,
        colour = "grey50",
        stroke = 0.3,
        na.rm = TRUE,
        inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
        ggplot2::aes(colour = .data[[group_col]]),
        size = as.numeric(config$point_size),
        alpha = 0.9,
        show.legend = show_group_leg
    ) +
    ggplot2::scale_fill_manual(
        name = config$labels$role,
        values = role_fills,
        breaks = config$role_levels,
        labels = unname(role_text[config$role_levels]),
        guide = ggplot2::guide_legend(
            order = 1,
            override.aes = list(shape = 22, size = 4.5, colour = "grey60", alpha = 1)
        )
    ) +
    ggplot2::scale_colour_manual(
        name = if (!is.null(config$labels$group)) config$labels$group else waiver(),
        values = gcols,
        guide = if (show_group_leg) {
            ggplot2::guide_legend(
                order = 2,
                override.aes = list(
                    shape = 16, size = 3.2, fill = NA, alpha = 1, stroke = 0
                )
            )
        } else {
            "none"
        }
    ) +
    ggplot2::coord_cartesian(
        xlim = c(x_min, x_max),
        ylim = c(y_min, y_max),
        expand = FALSE,
        clip = "off"
    ) +
    ggplot2::labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(
            size = 14, face = "bold", hjust = 0.5, colour = "#222222"
        ),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "grey40", linewidth = 0.6),
        axis.title = ggplot2::element_text(size = 11, colour = "#333333"),
        axis.text = ggplot2::element_text(size = 9, colour = "#444444"),
        legend.title = ggplot2::element_text(size = 10, face = "bold"),
        legend.text = ggplot2::element_text(size = 9),
        # Opaque white keys — PNG uses transparent page bg; NA keys look black
        # in dark UIs.
        legend.key = ggplot2::element_rect(fill = "white", colour = NA),
        legend.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(12, 14, 10, 10)
    )

if (nrow(lab_df)) {
    p <- p +
        ggrepel::geom_text_repel(
            data = lab_df,
            ggplot2::aes(label = .data$`.label`),
            size = as.numeric(config$label_size),
            colour = "grey15",
            box.padding = 0.28,
            point.padding = 0.18,
            min.segment.length = 0.05,
            segment.colour = "grey55",
            segment.size = 0.25,
            max.overlaps = Inf,
            seed = 1,
            show.legend = FALSE
        )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8.2, height = 6.4)
