# Shared helpers for omics-visualization R plot scripts.
# Plot scripts locate this file by walking up from --file
# until lib/common.R exists.

load_packages <- function(pkgs) {
    if (missing(pkgs) || !length(pkgs)) {
        stop("load_packages() requires a character vector of package names.",
            call. = FALSE
        )
    }
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0) {
        stop(
            paste(
                "Missing R package(s):", paste(missing, collapse = ", "),
                "\nInstall them before running this script."
            ),
            call. = FALSE
        )
    }
    suppressPackageStartupMessages({
        for (pkg in pkgs) {
            library(pkg, character.only = TRUE)
        }
    })
    invisible(pkgs)
}

script_basename <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (!length(file_arg)) {
        return("plot.R")
    }
    basename(sub("^--file=", "", file_arg[[1]]))
}

parse_io_args <- function(n_input = 1L, input_names = NULL) {
    args <- commandArgs(trailingOnly = TRUE)
    script <- script_basename()
    if (is.null(input_names)) {
        placeholders <- rep("input", n_input)
    } else if (length(input_names) != n_input) {
        stop("input_names must have length n_input.", call. = FALSE)
    } else {
        placeholders <- input_names
    }
    input_ph <- paste(sprintf("<%s>", placeholders), collapse = " ")
    usage <- sprintf("Usage: Rscript %s %s <output>", script, input_ph)

    if (length(args) == 0 || args[[1]] %in% c("-h", "--help")) {
        cat(usage, "\n", sep = "")
        quit(save = "no", status = 0)
    }
    if (length(args) != n_input + 1L) {
        stop(usage, call. = FALSE)
    }

    output <- args[[n_input + 1L]]
    # ggplot2/ggprism may open the default device during build; in Rscript
    # that creates a stray Rplots.pdf in the working directory.
    if (!interactive() && identical(grDevices::dev.cur(), 1L)) {
        grDevices::pdf(nullfile())
    }
    if (n_input == 1L) {
        list(input = args[[1]], output = output)
    } else {
        list(input = args[seq_len(n_input)], output = output)
    }
}

read_table_auto <- function(path) {
    if (is.null(path) || !nzchar(path)) {
        stop("Input path is empty.", call. = FALSE)
    }
    if (!file.exists(path)) {
        stop(paste("Input file not found:", path), call. = FALSE)
    }
    ext <- tolower(tools::file_ext(path))
    delim <- if (ext %in% c("tsv", "tab")) {
        "\t"
    } else if (identical(ext, "csv")) {
        ","
    } else {
        stop(
            paste(
                "Unsupported table format:", path,
                "\nUse .tsv, .tab, or .csv."
            ),
            call. = FALSE
        )
    }
    readr::read_delim(
        path,
        delim = delim,
        show_col_types = FALSE,
        progress = FALSE
    )
}

.common_r_path <- local({
    found <- NA_character_
    for (i in rev(seq_len(sys.nframe()))) {
        ofile <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
        if (is.null(ofile) || !nzchar(ofile)) {
            next
        }
        np <- normalizePath(ofile, mustWork = FALSE)
        if (grepl("common\\.R$", np) && file.exists(np)) {
            found <- np
            break
        }
    }
    found
})

.palette_json_path <- local({
    env_root <- Sys.getenv("OMICS_VISUALIZATION_SKILL_ROOT", unset = NA_character_)
    candidates <- character()
    if (!is.na(env_root) && nzchar(env_root)) {
        candidates <- c(candidates, file.path(
            normalizePath(env_root, mustWork = FALSE),
            "references", "palettes", "colors.json"
        ))
    }
    if (!is.na(.common_r_path)) {
        candidates <- c(candidates, file.path(
            dirname(dirname(dirname(.common_r_path))),
            "references", "palettes", "colors.json"
        ))
    }
    for (candidate in unique(candidates)) {
        if (file.exists(candidate)) {
            return(candidate)
        }
    }
    if (length(candidates)) candidates[[1]] else NA_character_
})

.hex_key <- function(x) {
    toupper(sub("^#", "", as.character(x)))
}

.hex_luminance <- function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    as.numeric(0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ])
}

.hex_too_light <- function(hex) {
    .hex_luminance(hex) > 0.90
}

load_palette_library <- function() {
    path <- .palette_json_path
    if (is.na(path) || !file.exists(path)) {
        stop("Cannot find generated references/palettes/colors.json", call. = FALSE)
    }
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Missing R package(s): jsonlite", call. = FALSE)
    }
    jsonlite::fromJSON(path, simplifyVector = TRUE)
}

# Discrete colors from colors.json. Prefer Qualitative, then Brand.
# Artwork and Concept are excluded unless a plot asks for them by id.
palette_pool <- function() {
    lib <- load_palette_library()
    preferred <- list(
        Qualitative = c(
            "Safe", "Bold", "Prism", "Vivid", "Dark2", "Paired", "Antique"
        ),
        Brand = c(
            "Algolia", "Asana", "Aetna", "Codepen", "Trello", "Yo",
            "Zendesk", "Gospel", "Socialbro", "Evaneos"
        )
    )
    out <- character()
    for (fam in names(preferred)) {
        if (is.null(lib[[fam]])) {
            next
        }
        names_in <- names(lib[[fam]])
        order <- unique(c(preferred[[fam]], names_in))
        order <- order[order %in% names_in]
        for (nm in order) {
            cols <- unname(as.character(lib[[fam]][[nm]]))
            cols <- cols[nzchar(cols) & !is.na(cols)]
            cols <- cols[!.hex_too_light(cols)]
            out <- c(out, cols)
        }
    }
    out[!duplicated(.hex_key(out))]
}

# Look up one named palette, e.g. palette_colors("Qualitative.Safe", n = 8).
# If n is longer than that palette, unused hex from Qualitative then Brand
# are appended. Hex is never interpolated.
palette_colors <- function(id, n = NULL) {
    lib <- load_palette_library()
    parts <- strsplit(id, ".", fixed = TRUE)[[1]]
    if (length(parts) < 2L) {
        stop(
            "Palette id must look like Family.Name, e.g. Qualitative.Safe",
            call. = FALSE
        )
    }
    family <- parts[[1]]
    name <- paste(parts[-1], collapse = ".")
    if (is.null(lib[[family]]) || is.null(lib[[family]][[name]])) {
        stop(paste("Unknown palette:", id), call. = FALSE)
    }
    base <- unname(as.character(lib[[family]][[name]]))
    if (is.null(n)) {
        return(base)
    }
    expand_palette(base, n)
}

expand_palette <- function(colors, n) {
    colors <- unname(as.character(colors))
    if (!length(colors)) {
        stop("expand_palette() requires at least one color.", call. = FALSE)
    }
    n <- as.integer(n)
    if (!is.finite(n) || n < 1L) {
        return(colors)
    }
    if (length(colors) >= n) {
        return(colors[seq_len(n)])
    }
    extra <- palette_pool()
    extra <- extra[!.hex_key(extra) %in% .hex_key(colors)]
    combined <- c(colors, extra)
    if (length(combined) >= n) {
        return(combined[seq_len(n)])
    }
    stop(
        paste(
            "Need", n, "colors, but colors.json Qualitative + Brand only",
            "supply", length(combined), "after dropping Artwork, Concept,",
            "and near-white hex."
        ),
        call. = FALSE
    )
}

# ggplot2 4.x dropped ggplot2:::check_linewidth. ggideogram 0.1.0 still
# calls it from GeomIdeogram$draw_panel. Replace that method; a locked
# ggplot2 namespace cannot grow a new binding.
patch_ggideogram_ggplot2 <- function() {
    if (!requireNamespace("ggideogram", quietly = TRUE)) {
        return(invisible(FALSE))
    }
    geom <- ggideogram:::GeomIdeogram
    geom$draw_panel <- function(self, data, panel_params, coord,
                                radius = grid::unit(0.1, "npc"),
                                chrom.col = NULL, chrom.lwd = NULL,
                                chrom.lty = NULL, lineend = "butt",
                                linejoin = "mitre") {
        if (!coord$is_linear()) {
            cli::cli_warning("Ideogram geom expects linear x and y scales")
        } else {
            coords <- coord$transform(data, panel_params)
            ggplot2:::ggname(
                "geom_ideogram",
                ggideogram:::ideogramGrob(
                    coords, radius, chrom.col, chrom.lwd, chrom.lty,
                    linejoin, lineend
                )
            )
        }
    }
    invisible(TRUE)
}

require_columns <- function(df, columns) {
    required <- unique(unlist(columns, use.names = FALSE))
    missing <- setdiff(required, names(df))
    if (length(missing) > 0) {
        stop(
            paste(
                "Missing required column(s):", paste(missing, collapse = ", "),
                "\nAvailable columns:", paste(names(df), collapse = ", ")
            ),
            call. = FALSE
        )
    }
    invisible(df)
}

save_ggplot <- function(plot, path, width = 8, height = 6, dpi = 300) {
    fmt <- tolower(tools::file_ext(path))
    if (!fmt %in% c("pdf", "png", "svg")) {
        stop(sprintf("Unsupported output format: %s", fmt), call. = FALSE)
    }
    out_dir <- dirname(path)
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (inherits(plot, "ggplot")) {
        ggplot2::ggsave(
            filename = path, plot = plot,
            width = width, height = height, units = "in", dpi = dpi,
            bg = if (identical(fmt, "pdf")) "white" else "transparent"
        )
    } else {
        if (identical(fmt, "png")) {
            grDevices::png(
                path, width = width, height = height, units = "in",
                res = dpi, bg = "transparent"
            )
        } else if (identical(fmt, "pdf")) {
            grDevices::pdf(path, width = width, height = height)
        } else {
            grDevices::svg(path, width = width, height = height, bg = "transparent")
        }
        tryCatch(
            grid::grid.draw(plot),
            finally = grDevices::dev.off()
        )
    }
    message("Saved: ", normalizePath(path, mustWork = TRUE))
    invisible(path)
}
