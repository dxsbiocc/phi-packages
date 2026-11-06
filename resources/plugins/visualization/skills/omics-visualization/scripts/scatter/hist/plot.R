#!/usr/bin/env Rscript

# Template-ID: scatter-hist
#
# Purpose:
#   Draw a 1D histogram or frequency polygon of a numeric value,
#   optionally coloured / filled by a group (stack or identity),
#   faceted (same column or a separate facet column), horizontal,
#   or as a mirrored pair of series. Facet + identity overlay with
#   inset coloured panel titles is CONFIG (preview_facet_overlay).
#
# Inputs:
#   Default (example.tsv): carat / price / cut / depth from diamonds.
#   Mirror (example_mirror.tsv): var1, var2.
#   Overlay (example_overlay.tsv): type, value.
#   Facet overlay (example_facet_overlay.tsv): photons, series, prepared.
#
# Output:
#   A PDF, PNG, or SVG histogram / freqpoly figure.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for geom, bins, position, orientation, facet, and
#   mirror columns. Separate facet vs fill group, inset titles, and
#   inward ticks are CONFIG. Grouped density histograms with rugs stay
#   scatter-hist-rug. 1D kernel densities stay scatter-density1d.
#   Ridgelines stay line-ridge.
#
# Scientific assumptions:
#   Each row is one observation. Bin counts / densities are display
#   summaries of the supplied values, not a test.

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
    # "single": one series (optional group fill/colour).
    # "mirror": two numeric columns, one up / one down.
    layout = "single",
    geom = "histogram",          # "histogram" | "freqpoly"
    columns = list(
        x = "carat",
        group = NULL,            # fill / colour series
        facet = NULL,            # optional panel column (else facet by group)
        a = "var1",              # mirror only
        b = "var2"
    ),
    labels = list(
        title = NULL,
        x = "carat",
        y = NULL,
        a = "variable1",
        b = "variable2",
        # Named character vector: facet level -> display title (Unicode OK).
        # NULL uses the facet level string as the inset title.
        facet = NULL
    ),
    bins = 30,
    binwidth = NULL,             # overrides bins when set
    position = "stack",          # "stack" | "identity" | "dodge"
    y_stat = "count",            # "count" | "density"
    orientation = "x",           # "x" | "y"
    facet = FALSE,
    facet_scales = "fixed",      # "fixed" | "free_y" | "free"
    facet_nrow = NULL,
    facet_ncol = NULL,
    # "strip": default ggplot strips. "inset": coloured title inside panel.
    facet_title = "strip",
    tick_direction = "out",      # "out" | "in"
    show_legend = NULL,          # NULL = auto (off when faceted)
    alpha = 0.85,
    colour = "white",           # bar border
    linewidth = 0.25,
    group_order = NULL,          # explicit level order for the group
    facet_order = NULL,          # explicit level order for the facet
    palettes = list(
        group = "Qualitative.Bold",
        mirror_a = "Qualitative.Safe:7",
        mirror_b = "Qualitative.Safe:5"
    ),
    base_size = 13,
    size = list(width = 7.5, height = 5)
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

order_levels <- function(vals, preferred) {
    lv <- unique(as.character(vals))
    if (is.null(preferred)) return(lv)
    c(intersect(as.character(preferred), lv), setdiff(lv, as.character(preferred)))
}

df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
if (identical(config$layout, "mirror")) {
    require_columns(df, list(config$columns$a, config$columns$b))
    a <- as.numeric(df[[config$columns$a]])
    b <- as.numeric(df[[config$columns$b]])
    a <- a[is.finite(a)]
    b <- b[is.finite(b)]
    col_a <- pal_ref(config$palettes$mirror_a)
    col_b <- pal_ref(config$palettes$mirror_b)
    hist_args <- list(
        colour = config$colour, alpha = config$alpha,
        linewidth = config$linewidth
    )
    if (is.null(config$binwidth)) {
        hist_args$bins <- config$bins
    } else {
        hist_args$binwidth <- config$binwidth
    }
    p <- ggplot() +
        do.call(geom_histogram, c(list(
            mapping = aes(x = a, y = after_stat(density)),
            fill = col_a
        ), hist_args)) +
        do.call(geom_histogram, c(list(
            mapping = aes(x = b, y = -after_stat(density)),
            fill = col_b
        ), hist_args)) +
        annotate(
            "label", x = Inf, y = Inf, label = config$labels$a,
            hjust = 1.1, vjust = 1.4, colour = col_a, size = 4,
            label.size = 0
        ) +
        annotate(
            "label", x = Inf, y = -Inf, label = config$labels$b,
            hjust = 1.1, vjust = -0.4, colour = col_b, size = 4,
            label.size = 0
        ) +
        geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
        labs(title = config$labels$title, x = config$labels$x, y = config$labels$y)
} else {
    x_col <- config$columns$x
    g_col <- config$columns$group
    f_col <- config$columns$facet
    require_columns(df, list(x_col))
    df[[x_col]] <- as.numeric(df[[x_col]])
    df <- df[is.finite(df[[x_col]]), , drop = FALSE]
    if (!nrow(df)) stop("No finite values.", call. = FALSE)
    use_group <- !is.null(g_col) && g_col %in% names(df)
    use_facet_col <- !is.null(f_col) && f_col %in% names(df)
    do_facet <- isTRUE(config$facet) && (use_facet_col || use_group)
    if (use_group) {
        df[[g_col]] <- factor(df[[g_col]], levels = order_levels(df[[g_col]], config$group_order))
        cols <- expand_palette(
            palette_colors(config$palettes$group), nlevels(df[[g_col]])
        )
        names(cols) <- levels(df[[g_col]])
    }
    if (use_facet_col) {
        df[[f_col]] <- factor(df[[f_col]], levels = order_levels(df[[f_col]], config$facet_order))
    }

    dens <- identical(config$y_stat, "density")
    horiz <- identical(config$orientation, "y")
    common <- list(alpha = config$alpha)
    if (is.null(config$binwidth)) {
        common$bins <- config$bins
    } else {
        common$binwidth <- config$binwidth
    }

    if (horiz) {
        p <- ggplot(df, aes(y = .data[[x_col]]))
        if (identical(config$geom, "freqpoly")) {
            mapping <- if (dens) {
                if (use_group) {
                    aes(x = after_stat(density), colour = .data[[g_col]])
                } else {
                    aes(x = after_stat(density))
                }
            } else if (use_group) {
                aes(x = after_stat(count), colour = .data[[g_col]])
            } else {
                aes(x = after_stat(count))
            }
            p <- p + do.call(geom_freqpoly, c(list(
                mapping = mapping, linewidth = 0.7,
                colour = if (use_group) NULL else pal_ref("Qualitative.Bold:3")
            ), common[!names(common) %in% "alpha"]))
        } else {
            mapping <- if (dens) {
                if (use_group) {
                    aes(x = after_stat(density), fill = .data[[g_col]])
                } else {
                    aes(x = after_stat(density))
                }
            } else if (use_group) {
                aes(x = after_stat(count), fill = .data[[g_col]])
            } else {
                aes(x = after_stat(count))
            }
            args <- c(list(
                mapping = mapping, colour = config$colour,
                linewidth = config$linewidth,
                position = if (use_group) config$position else "stack"
            ), common)
            if (!use_group) args$fill <- pal_ref("Qualitative.Bold:3")
            p <- p + do.call(geom_histogram, args)
        }
    } else {
        p <- ggplot(df, aes(x = .data[[x_col]]))
        if (identical(config$geom, "freqpoly")) {
            mapping <- if (dens) {
                if (use_group) {
                    aes(y = after_stat(density), colour = .data[[g_col]])
                } else {
                    aes(y = after_stat(density))
                }
            } else if (use_group) {
                aes(y = after_stat(count), colour = .data[[g_col]])
            } else {
                aes(y = after_stat(count))
            }
            fp <- list(mapping = mapping, linewidth = 0.7)
            if (!use_group) fp$colour <- pal_ref("Qualitative.Bold:3")
            p <- p + do.call(geom_freqpoly, c(fp, common[!names(common) %in% "alpha"]))
        } else {
            mapping <- if (dens) {
                if (use_group) {
                    aes(y = after_stat(density), fill = .data[[g_col]])
                } else {
                    aes(y = after_stat(density))
                }
            } else if (use_group) {
                aes(y = after_stat(count), fill = .data[[g_col]])
            } else {
                aes(y = after_stat(count))
            }
            args <- c(list(
                mapping = mapping, colour = config$colour,
                linewidth = config$linewidth,
                position = if (use_group) config$position else "stack"
            ), common)
            if (!use_group) args$fill <- pal_ref("Qualitative.Bold:3")
            p <- p + do.call(geom_histogram, args)
        }
    }

    if (use_group) {
        if (identical(config$geom, "freqpoly")) {
            p <- p + scale_colour_manual(values = cols, name = NULL)
        } else {
            p <- p + scale_fill_manual(values = cols, name = NULL)
        }
    }

    facet_var <- if (use_facet_col) f_col else if (use_group) g_col else NULL
    if (do_facet && !is.null(facet_var)) {
        inset <- identical(config$facet_title, "inset")
        p <- p + facet_wrap(
            stats::as.formula(paste("~", facet_var)),
            scales = config$facet_scales,
            nrow = config$facet_nrow,
            ncol = config$facet_ncol
        )
        if (inset) {
            lv <- levels(df[[facet_var]])
            lab_map <- config$labels$facet
            title_txt <- vapply(lv, function(k) {
                if (!is.null(lab_map) && !is.null(lab_map[[k]])) {
                    as.character(lab_map[[k]])
                } else {
                    as.character(k)
                }
            }, character(1))
            # Colour title by matching fill series when levels align; else black.
            title_col <- if (use_group && all(lv %in% names(cols))) {
                unname(cols[lv])
            } else {
                rep("grey15", length(lv))
            }
            lab_df <- data.frame(
                facet = factor(lv, levels = lv),
                label = title_txt,
                colour = title_col,
                stringsAsFactors = FALSE
            )
            names(lab_df)[1] <- facet_var
            lab_df$x <- mean(range(df[[x_col]], na.rm = TRUE))
            lab_df$y <- Inf
            for (i in seq_len(nrow(lab_df))) {
                p <- p + geom_text(
                    data = lab_df[i, , drop = FALSE],
                    aes(x = .data$x, y = .data$y, label = .data$label),
                    colour = lab_df$colour[[i]],
                    inherit.aes = FALSE,
                    vjust = 1.6,
                    size = 3.8,
                    fontface = "bold"
                )
            }
        }
    }
    p <- p + labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    )
}

tick_len <- if (identical(config$tick_direction, "in")) {
    grid::unit(-0.15, "cm")
} else {
    grid::unit(0.15, "cm")
}

p <- p +
    theme_bw(base_size = config$base_size) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
        axis.text = element_text(colour = "black"),
        axis.ticks.length = tick_len,
        axis.ticks = element_line(colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black")
    )

# Inset titles hide strips.
if (!identical(config$layout, "mirror") &&
    isTRUE(config$facet) &&
    identical(config$facet_title, "inset")) {
    p <- p + theme(
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = grid::unit(0.6, "lines")
    )
}

legend_off <- if (!is.null(config$show_legend)) {
    !isTRUE(config$show_legend)
} else {
    isTRUE(config$facet) && !identical(config$layout, "mirror")
}
if (legend_off) {
    p <- p + theme(legend.position = "none")
}
if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
    p <- p + theme(plot.title = element_blank())
}
if (is.null(config$labels$y) || !nzchar(as.character(config$labels$y))) {
    p <- p + theme(axis.title.y = element_blank())
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width, height = config$size$height)
