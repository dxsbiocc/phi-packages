#!/usr/bin/env Rscript

# Top-variable-gene selection and row z-scoring adapted from bioSkills's
# differential-expression/de-visualization heatmap.R example and SKILL.md's
# own "row-scaling trap" guidance (z-score per gene via scale='row', not
# raw values). Styling (diverging red/blue palette for the z-score fill)
# adapted from resources/skills/omics-visualization's heatmap/cluster_basic
# template's intent -- implemented with pheatmap, already proven in this
# family's shared container, rather than that template's own
# ComplexHeatmap dependency. See
# resources/palettes/palettes.yaml
# (Diverging.RdBu) for where the palette comes from.

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
parse_args <- function(x) {
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text = x, what = 'character', quiet = TRUE))
    args_vals <- lapply(args_vals, function(z) { length(z) <- 2; z })
    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[!is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE) {
    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))
    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }
    read.delim(file, sep = separator, header = header, row.names = row.names, check.names = check.names)
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', 'de', '$task.ext.prefix'),
    matrix_file    = '$matrix',
    sample_file    = '$samplesheet',
    gene_id_col    = 'gene_id',
    sample_id_col  = 'sample',
    annotate_by    = 'condition',
    top_n          = 50,
    scale_rows     = TRUE
)
opt_types <- lapply(opt, class)

args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)) {
    if (!ao %in% names(opt)) stop(paste("Invalid option:", ao))
    opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
}

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(pheatmap)

################################################
################################################
## READ INPUTS                                ##
################################################
################################################

mat.table              <- read_delim_flexible(opt\$matrix_file, header = TRUE, row.names = NULL, check.names = FALSE)
rownames(mat.table)     <- mat.table[[opt\$gene_id_col]]
mat.table[[opt\$gene_id_col]] <- NULL
mat                     <- as.matrix(mat.table)

sample.sheet <- read_delim_flexible(opt\$sample_file, header = TRUE, check.names = FALSE)
if (!opt\$sample_id_col %in% colnames(sample.sheet)) {
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

missing_samples <- setdiff(colnames(mat), rownames(sample.sheet))
if (length(missing_samples) > 0) {
    stop(paste("Samples in matrix missing from sample sheet:", paste(missing_samples, collapse = ",")))
}
sample.sheet <- sample.sheet[colnames(mat), , drop = FALSE]

if (!opt\$annotate_by %in% colnames(sample.sheet)) {
    stop(paste0("Column '", opt\$annotate_by, "' not found in sample sheet. Available columns: ", paste(colnames(sample.sheet), collapse = ", ")))
}
annotation_col <- sample.sheet[, opt\$annotate_by, drop = FALSE]

################################################
################################################
## SELECT TOP-VARIABLE GENES                  ##
################################################
################################################

rv     <- apply(mat, 1, stats::var) # base R only -- no matrixStats dependency
select <- order(rv, decreasing = TRUE)[seq_len(min(opt\$top_n, length(rv)))]
mat    <- mat[select, , drop = FALSE]

################################################
################################################
## PLOT                                       ##
################################################
################################################

# Diverging red/blue hexes from Phi's Diverging.RdBu
# catalog entry (resources/palettes/palettes.yaml) -- for signed, zero-centered
# values, which row z-scores are.
DIVERGING_RDBU <- rev(c("#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac"))
heatmap_colors <- colorRampPalette(DIVERGING_RDBU)(100)

PlotFile <- paste0(opt\$output_prefix, ".heatmap.pdf")

pheatmap(
    mat,
    scale                    = if (opt\$scale_rows) "row" else "none",
    annotation_col           = annotation_col,
    show_rownames            = nrow(mat) <= 60,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    color                    = heatmap_colors,
    border_color             = NA,
    main                     = paste("Top", opt\$top_n, "variable genes"),
    filename                 = PlotFile,
    width                    = 8,
    height                   = max(6, min(20, nrow(mat) * 0.15))
)

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']], R.version[['minor']], sep = ".")

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-pheatmap:', as.character(packageVersion('pheatmap')))
    ),
    'versions.yml')
