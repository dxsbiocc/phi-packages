#!/usr/bin/env Rscript

# Template-ID: scatter-ternary
#
# Purpose:
#   Draw a ternary (simplex) scatter: each point is one observation
#   whose position is the relative contribution of three supplied
#   non-negative components.
#
# Inputs:
#   One row per observation. Default example:
#     - EPI, TE, PrE: three lineage (or state / compartment) scores
#     - Group: discrete membership (vertex-only, pairwise, or all three)
#   Optional size column encodes a supplied magnitude (abundance).
#   A numeric colour column (age, score) uses a sequential fill.
#
# Output:
#   A PDF, PNG, or SVG ternary scatter plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggtern, scales
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (component columns, colour,
#   size, vertex labels). Discrete fill is mapped to the three
#   vertices (and their pairwise / interior combinations), not
#   palette order. Continuous colour, and size on or off, is CONFIG,
#   not a new id. Edit DATA PREPARATION to change group order.
#   Edit PLOT only when the glyph must change.
#
# Scientific assumptions:
#   The three components are supplied. This script does not compute
#   lineage scores, NMF, relative abundance, enrichment, or
#   absorption probabilities. Components are closed to a composition
#   (a + b + c = 1) for display. Position is relative contribution,
#   not a Cartesian x-y embedding.

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
        a = "EPI",
        b = "TE",
        c = "PrE",
        colour = "Group",
        size = NULL
    ),
    fill_by_vertex = TRUE,
    labels = list(
        title = "",
        L = "EPI",
        T = "TE",
        R = "PrE",
        colour = "Group",
        size = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggtern", "scales"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
need <- config$columns[c("a", "b", "c")]
require_columns(df, need)

a_col <- config$columns$a
b_col <- config$columns$b
c_col <- config$columns$c
col_col <- config$columns$colour
size_col <- config$columns$size

df[[a_col]] <- as.numeric(df[[a_col]])
df[[b_col]] <- as.numeric(df[[b_col]])
df[[c_col]] <- as.numeric(df[[c_col]])

ok <- is.finite(df[[a_col]]) & is.finite(df[[b_col]]) & is.finite(df[[c_col]]) &
    df[[a_col]] >= 0 & df[[b_col]] >= 0 & df[[c_col]] >= 0
comp_sum <- df[[a_col]] + df[[b_col]] + df[[c_col]]
ok <- ok & is.finite(comp_sum) & comp_sum > 0
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with invalid ternary components.")
    df <- df[ok, , drop = FALSE]
    comp_sum <- comp_sum[ok]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}
df[[a_col]] <- df[[a_col]] / comp_sum
df[[b_col]] <- df[[b_col]] / comp_sum
df[[c_col]] <- df[[c_col]] / comp_sum

colour_mode <- "none"
if (!is.null(col_col) && col_col %in% names(df)) {
    if (is.numeric(df[[col_col]]) && !is.factor(df[[col_col]])) {
        colour_mode <- "continuous"
        df[[col_col]] <- as.numeric(df[[col_col]])
    } else {
        colour_mode <- "discrete"
        df[[col_col]] <- factor(
            df[[col_col]],
            levels = unique(as.character(df[[col_col]]))
        )
    }
} else {
    col_col <- NULL
}

if (!is.null(size_col) && size_col %in% names(df)) {
    df$size_enc <- as.numeric(df[[size_col]])
    size_ok <- is.finite(df$size_enc) & df$size_enc >= 0
    if (any(!size_ok)) {
        message("Dropped ", sum(!size_ok), " row(s) with invalid size.")
        df <- df[size_ok, , drop = FALSE]
    }
} else {
    size_col <- NULL
}

outline <- palette_colors("Brand.Algolia")[[10]]
fill_continuous <- palette_colors("Quantitative.BluGrn")
alg <- palette_colors("Brand.Algolia")
safe <- palette_colors("Qualitative.Safe")
# L / T / R vertices, then LT / LR / TR edges, then interior.
vertex_cols <- c(L = safe[[12]], T = alg[[1]], R = alg[[4]])
pair_cols <- c(LT = alg[[8]], LR = alg[[5]], TR = alg[[2]])
all_col <- alg[[7]]

tokenise <- function(x) {
    tolower(trimws(gsub("[^a-z0-9]+", " ", tolower(x))))
}

has_token <- function(hay, needle) {
    if (!nzchar(needle)) {
        return(FALSE)
    }
    grepl(paste0("(^| )", needle, "( |$)"), hay)
}

group_role <- function(name, lab_L, lab_T, lab_R) {
    key <- tokenise(name)
    hits <- c(
        L = has_token(key, tokenise(lab_L)),
        T = has_token(key, tokenise(lab_T)),
        R = has_token(key, tokenise(lab_R))
    )
    n_hit <- sum(hits)
    if (grepl("\\ball\\b|\\bthree\\b|\\bmixed\\b|\\binterior\\b", key) &&
        n_hit != 1L) {
        return("ALL")
    }
    if (n_hit == 3L) {
        return("ALL")
    }
    if (n_hit == 2L) {
        return(paste0(names(hits)[hits], collapse = ""))
    }
    if (n_hit == 1L) {
        return(names(hits)[hits])
    }
    NA_character_
}

composition_role <- function(a, b, c) {
    m <- c(L = mean(a), T = mean(b), R = mean(c))
    o <- sort(m, decreasing = TRUE)
    if (o[[1]] >= 0.55 && (o[[1]] - o[[2]]) >= 0.12) {
        return(names(o)[[1]])
    }
    if (o[[3]] <= 0.18) {
        keep <- names(o)[1:2]
        keep <- keep[order(match(keep, c("L", "T", "R")))]
        return(paste0(keep, collapse = ""))
    }
    "ALL"
}

role_fill <- function(role) {
    if (identical(role, "ALL") || is.na(role)) {
        return(all_col)
    }
    if (role %in% names(vertex_cols)) {
        return(unname(vertex_cols[[role]]))
    }
    unname(pair_cols[[role]])
}

discrete_fills <- NULL
if (identical(colour_mode, "discrete")) {
    lv <- levels(df[[col_col]])
    discrete_fills <- vapply(lv, function(g) {
        role <- group_role(
            g, config$labels$L, config$labels$T, config$labels$R
        )
        if (isTRUE(config$fill_by_vertex)) {
            if (is.na(role)) {
                sub <- df[df[[col_col]] == g, , drop = FALSE]
                role <- composition_role(
                    sub[[a_col]], sub[[b_col]], sub[[c_col]]
                )
            }
            role_fill(role)
        } else {
            NA_character_
        }
    }, character(1))
    names(discrete_fills) <- lv
    if (!isTRUE(config$fill_by_vertex) || anyNA(discrete_fills)) {
        discrete_fills <- palette_colors("Brand.Algolia", n = length(lv))
        names(discrete_fills) <- lv
    }
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
if (identical(colour_mode, "none") && is.null(size_col)) {
    p <- ggtern(
        df,
        aes(x = .data[[a_col]], y = .data[[b_col]], z = .data[[c_col]])
    ) +
        geom_point(
            shape = 21,
            size = 2.2,
            fill = palette_colors("Brand.Algolia")[[1]],
            colour = outline,
            stroke = 0.25,
            alpha = 0.88
        )
} else if (identical(colour_mode, "none")) {
    p <- ggtern(
        df,
        aes(
            x = .data[[a_col]],
            y = .data[[b_col]],
            z = .data[[c_col]],
            size = size_enc
        )
    ) +
        geom_point(
            shape = 21,
            fill = palette_colors("Brand.Algolia")[[1]],
            colour = outline,
            stroke = 0.25,
            alpha = 0.88
        )
} else if (is.null(size_col)) {
    p <- ggtern(
        df,
        aes(
            x = .data[[a_col]],
            y = .data[[b_col]],
            z = .data[[c_col]],
            fill = .data[[col_col]]
        )
    ) +
        geom_point(
            shape = 21,
            size = 2.2,
            colour = outline,
            stroke = 0.25,
            alpha = 0.88
        )
} else {
    p <- ggtern(
        df,
        aes(
            x = .data[[a_col]],
            y = .data[[b_col]],
            z = .data[[c_col]],
            fill = .data[[col_col]],
            size = size_enc
        )
    ) +
        geom_point(
            shape = 21,
            colour = outline,
            stroke = 0.25,
            alpha = 0.88
        )
}

p <- p + geom_mask()

if (identical(colour_mode, "discrete")) {
    p <- p + scale_fill_manual(
        name = config$labels$colour,
        values = discrete_fills,
        drop = FALSE
    )
} else if (identical(colour_mode, "continuous")) {
    p <- p + scale_fill_gradientn(
        name = config$labels$colour,
        colours = fill_continuous
    )
}

if (!is.null(size_col)) {
    p <- p + scale_size_continuous(
        name = config$labels$size,
        range = c(1.2, 8)
    )
}

p <- p +
    labs(
        title = config$labels$title,
        x = config$labels$L,
        y = config$labels$T,
        z = config$labels$R
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        tern.plot.background = element_blank(),
        tern.panel.background = element_rect(fill = NA, colour = NA),
        tern.panel.grid.major = element_line(
            colour = safe[[12]], linewidth = 0.3, linetype = "22"
        ),
        tern.panel.grid.minor = element_blank(),
        tern.panel.grid.major.T = element_line(
            colour = safe[[12]], linewidth = 0.3, linetype = "22"
        ),
        tern.panel.grid.major.L = element_line(
            colour = safe[[12]], linewidth = 0.3, linetype = "22"
        ),
        tern.panel.grid.major.R = element_line(
            colour = safe[[12]], linewidth = 0.3, linetype = "22"
        ),
        tern.panel.grid.minor.T = element_blank(),
        tern.panel.grid.minor.L = element_blank(),
        tern.panel.grid.minor.R = element_blank(),
        panel.grid = element_blank(),
        tern.axis.arrow.show = FALSE,
        tern.axis.arrow = element_blank(),
        tern.axis.arrow.text = element_blank(),
        tern.axis.line = element_line(colour = outline, linewidth = 0.55),
        tern.axis.line.L = element_line(colour = vertex_cols[["L"]]),
        tern.axis.line.T = element_line(colour = vertex_cols[["T"]]),
        tern.axis.line.R = element_line(colour = vertex_cols[["R"]]),
        tern.axis.title.L = element_text(colour = vertex_cols[["L"]]),
        tern.axis.title.T = element_text(colour = vertex_cols[["T"]]),
        tern.axis.title.R = element_text(colour = vertex_cols[["R"]])
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 6.4, height = 5.8)
