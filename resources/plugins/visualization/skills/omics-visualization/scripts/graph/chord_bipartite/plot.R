#!/usr/bin/env Rscript

# Template-ID: graph-chord-bipartite
#
# Purpose:
#   Draw a bipartite circular chord (phenotype ↔ locus / gene):
#   side-A nodes on one arc with horizontal labels, side-B nodes on the
#   opposite arc with radial labels, coloured ribbons, and degree
#   numbers inside node circles (preview.png). Typical use: coloc
#   phenotype–HF-locus sharing with PP H4 above a supplied cutoff.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. nodes: id, side (A|B or phenotype|locus), optional label, colour
#     2. links: from, to (node ids); optional colour (defaults to from)
#
# Output:
#   A PDF, PNG, or SVG bipartite chord plot.
#
# Dependencies:
#   ggplot2, readr, ggforce, scales, grid
#
# Example:
#   Rscript plot.R nodes.tsv links.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, arc spans, colours, sizes, and labels.
#   Extra phenotypes or loci are extra rows, not a new id. A
#   flow-weighted Circos chord stays graph-chord. Equal-spacing
#   homogeneous circular networks stay graph-circular. Disease-sector
#   bars with gene labels stay bar-circos-links. Genomic rearrangement
#   ribbons stay ideogram-circos. Do not run coloc or GWAS here.
#
# Scientific assumptions:
#   Nodes, sides, and links are supplied. Degree is counted from the
#   link table for display. Ribbon colour defaults to the from-node
#   colour. This script does not compute PP H4 or filter by threshold.

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

io <- parse_io_args(n_input = 2L, input_names = c("nodes", "links"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        id = "id",
        side = "side",
        label = "label",
        colour = "colour",
        from = "from",
        to = "to",
        link_colour = "colour"
    ),
    # Map raw side values → A (left/phenotype) or B (right/locus).
    side_map = c(
        A = "A", phenotype = "A", trait = "A",
        B = "B", locus = "B", gene = "B"
    ),
    # Arc spans in degrees (0 = east, counter-clockwise). Left ≈ phenotypes.
    arc_A = c(95, 265),
    arc_B = c(-75, 75),
    radius = 1,
    node_size = 9.5,
    ribbon_alpha = 0.45,
    ribbon_linewidth = 1.15,
    label_size_A = 3.0,
    label_size_B = 2.7,
    label_colour_A = "grey10",
    label_colour_B = "#C0392B",
    degree_size = 2.6,
    palettes = list(node = "Qualitative.Bold"),
    labels = list(title = NULL),
    size = list(width = 9.5, height = 9.5)
)

load_packages(c("ggplot2", "readr", "ggforce", "scales", "grid"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
nodes <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
links <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)

id_col <- config$columns$id
side_col <- config$columns$side
lab_col <- config$columns$label
ncol_n <- config$columns$colour
from_col <- config$columns$from
to_col <- config$columns$to
lcol <- config$columns$link_colour

require_columns(nodes, list(id = id_col, side = side_col))
require_columns(links, list(from = from_col, to = to_col))

nodes[[id_col]] <- as.character(nodes[[id_col]])
nodes[[side_col]] <- as.character(nodes[[side_col]])
if (is.null(lab_col) || !lab_col %in% names(nodes)) {
    nodes$.label <- nodes[[id_col]]
} else {
    nodes$.label <- as.character(nodes[[lab_col]])
    nodes$.label[is.na(nodes$.label) | !nzchar(nodes$.label)] <- nodes[[id_col]][
        is.na(nodes$.label) | !nzchar(nodes$.label)
    ]
}

sm <- config$side_map
side_raw <- nodes[[side_col]]
side_norm <- unname(sm[side_raw])
side_norm[is.na(side_norm)] <- ifelse(
    tolower(side_raw[is.na(side_norm)]) %in% c("a", "phenotype", "trait"),
    "A", "B"
)
nodes$.side <- side_norm

links[[from_col]] <- as.character(links[[from_col]])
links[[to_col]] <- as.character(links[[to_col]])
ok_l <- links[[from_col]] %in% nodes[[id_col]] &
    links[[to_col]] %in% nodes[[id_col]] &
    links[[from_col]] != links[[to_col]]
links <- links[ok_l, , drop = FALSE]
if (!nrow(links)) stop("No valid links after id matching.", call. = FALSE)

# Degree from undirected endpoint counts.
deg <- table(c(links[[from_col]], links[[to_col]]))
nodes$.degree <- as.integer(deg[nodes[[id_col]]])
nodes$.degree[is.na(nodes$.degree)] <- 0L
nodes <- nodes[nodes$.degree > 0L | nodes[[id_col]] %in% c(links[[from_col]], links[[to_col]]), , drop = FALSE]

# Colours: supplied column, else palette by id within side A (phenotypes).
if (!is.null(ncol_n) && ncol_n %in% names(nodes)) {
    nodes$.colour <- as.character(nodes[[ncol_n]])
} else {
    nodes$.colour <- NA_character_
}
a_ids <- nodes[[id_col]][nodes$.side == "A"]
miss_a <- a_ids[is.na(nodes$.colour[match(a_ids, nodes[[id_col]])]) |
    !nzchar(nodes$.colour[match(a_ids, nodes[[id_col]])])]
if (length(miss_a)) {
    cols_a <- expand_palette(palette_colors(config$palettes$node), length(unique(miss_a)))
    names(cols_a) <- unique(miss_a)
    for (nm in names(cols_a)) {
        nodes$.colour[nodes[[id_col]] == nm] <- cols_a[[nm]]
    }
}
# Side B inherits a muted default if missing; ribbons use from-node colour.
nodes$.colour[nodes$.side == "B" & (is.na(nodes$.colour) | !nzchar(nodes$.colour))] <- "#4A90A4"

# Place nodes on arcs.
place_arc <- function(ids, arc_deg, radius) {
    n <- length(ids)
    if (!n) {
        return(data.frame(id = character(), x = numeric(), y = numeric(),
                          angle = numeric(), stringsAsFactors = FALSE))
    }
    a0 <- arc_deg[[1]] * pi / 180
    a1 <- arc_deg[[2]] * pi / 180
    # Even spacing inclusive of ends when n > 1.
    ang <- if (n == 1L) {
        (a0 + a1) / 2
    } else {
        seq(a0, a1, length.out = n)
    }
    data.frame(
        id = ids,
        angle = ang,
        x = radius * cos(ang),
        y = radius * sin(ang),
        stringsAsFactors = FALSE
    )
}

ord_A <- nodes[nodes$.side == "A", , drop = FALSE]
ord_B <- nodes[nodes$.side == "B", , drop = FALSE]
# Stable order: as in input table.
pos_A <- place_arc(ord_A[[id_col]], config$arc_A, config$radius)
pos_B <- place_arc(ord_B[[id_col]], config$arc_B, config$radius)
pos <- rbind(pos_A, pos_B)
nodes <- merge(nodes, pos, by.x = id_col, by.y = "id", all.x = TRUE, sort = FALSE)

# Ribbon paths: quadratic bezier via mid control near centre.
links$.col <- if (!is.null(lcol) && lcol %in% names(links)) {
    as.character(links[[lcol]])
} else {
    NA_character_
}
from_col_map <- setNames(nodes$.colour, nodes[[id_col]])
links$.col[is.na(links$.col) | !nzchar(links$.col)] <- unname(
    from_col_map[links[[from_col]]]
)
fx <- setNames(nodes$x, nodes[[id_col]])
fy <- setNames(nodes$y, nodes[[id_col]])
links$x <- unname(fx[links[[from_col]]])
links$y <- unname(fy[links[[from_col]]])
links$xend <- unname(fx[links[[to_col]]])
links$yend <- unname(fy[links[[to_col]]])
# Control point: pull toward origin so chords bow inward.
links$cx <- 0.15 * (links$x + links$xend)
links$cy <- 0.15 * (links$y + links$yend)

# Build bezier long data for geom_bezier (ggforce).
bez <- do.call(rbind, lapply(seq_len(nrow(links)), function(i) {
    data.frame(
        x = c(links$x[[i]], links$cx[[i]], links$xend[[i]]),
        y = c(links$y[[i]], links$cy[[i]], links$yend[[i]]),
        group = i,
        colour = links$.col[[i]],
        stringsAsFactors = FALSE
    )
}))

# Label anchors outside the ring.
lab_r <- config$radius * 1.12
nodes$lx <- lab_r * cos(nodes$angle)
nodes$ly <- lab_r * sin(nodes$angle)
# Side A: keep text horizontal; Side B: rotate radially.
nodes$lab_angle <- nodes$angle * 180 / pi
# Flip labels on the left half so they read outward.
flip <- cos(nodes$angle) < 0
nodes$lab_angle[nodes$.side == "B" & flip] <- nodes$lab_angle[nodes$.side == "B" & flip] + 180
nodes$hjust <- 0.5
nodes$hjust[nodes$.side == "A"] <- ifelse(cos(nodes$angle[nodes$.side == "A"]) >= 0, 0, 1)
nodes$hjust[nodes$.side == "B"] <- ifelse(flip[nodes$.side == "B"], 1, 0)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    ggforce::geom_bezier(
        data = bez,
        aes(x = x, y = y, group = group, colour = colour),
        alpha = config$ribbon_alpha,
        linewidth = config$ribbon_linewidth,
        show.legend = FALSE
    ) +
    geom_point(
        data = nodes,
        aes(x = x, y = y, fill = .colour),
        shape = 21,
        size = config$node_size,
        colour = "grey20",
        stroke = 0.35,
        show.legend = FALSE
    ) +
    geom_text(
        data = nodes,
        aes(x = x, y = y, label = .degree),
        size = config$degree_size,
        colour = "white",
        fontface = "bold",
        show.legend = FALSE
    ) +
    # Side A labels (horizontal).
    geom_text(
        data = nodes[nodes$.side == "A", , drop = FALSE],
        aes(x = lx, y = ly, label = .label, hjust = hjust),
        size = config$label_size_A,
        colour = config$label_colour_A,
        show.legend = FALSE
    ) +
    # Side B labels (radial).
    geom_text(
        data = nodes[nodes$.side == "B", , drop = FALSE],
        aes(
            x = lx, y = ly, label = .label,
            angle = lab_angle, hjust = hjust
        ),
        size = config$label_size_B,
        colour = config$label_colour_B,
        show.legend = FALSE
    ) +
    scale_colour_identity() +
    scale_fill_identity() +
    coord_equal(clip = "off") +
    labs(title = config$labels$title) +
    theme_void(base_size = 12) +
    theme(
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(12, 12, 12, 12)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = config$size$width,
    height = config$size$height
)
