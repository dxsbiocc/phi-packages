#!/usr/bin/env Rscript

# New local script (no nf-core module to vendor). Statistical approach --
# TMM normalization (calcNormFactors), filterByExpr, and the
# quasi-likelihood F-test (glmQLFit/glmQLFTest, the modern edgeR default
# over the legacy exact test) -- follows bioSkills's differential-
# expression/edger-basics skill. Parameter handling (cell-means design,
# formula/comparison override, parse_args/read_delim_flexible helpers)
# mirrors ../../deseq2/templates/deseq2_differential.R for consistency
# within this family, swapping DESeq2's model for edgeR's.

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

#' Turn "null" or empty strings into actual NULL
nullify <- function(x) {
    if (is.character(x) && (tolower(x) == "null" || x == "")) NULL else x
}

#' Check for a non-empty, non-whitespace string
is_valid_string <- function(input) {
    !is.null(input) && nzchar(trimws(input))
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

opt <- list(
    # count_file/sample_file are the only true process inputs (see
    # ../main.nf); everything else here has no channel variable of its own
    # -- it comes solely through the task.ext.args parse below, the same
    # way ../../qc/templates/deseq2_qc.R's options do.
    output_prefix     = ifelse('$task.ext.prefix' == 'null', 'de', '$task.ext.prefix'),
    count_file         = '$counts',
    sample_file        = '$samplesheet',
    contrast_variable  = 'null',
    reference_level    = 'null',
    target_level       = 'null',
    formula            = 'null',
    contrast_string    = 'null',
    sample_id_col      = 'sample',
    gene_id_col        = 'gene_id'
)

keys <- c("formula", "contrast_string", "contrast_variable", "reference_level", "target_level")
opt[keys] <- lapply(opt[keys], nullify)

args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)) {
    if (!ao %in% names(opt)) stop(paste("Invalid option:", ao))
    opt[[ao]] <- args_opt[[ao]]
}
opt[keys] <- lapply(opt[keys], nullify)

if (is_valid_string(opt\$formula)) {
    required_opts <- c("output_prefix", "contrast_string")
} else {
    required_opts <- c("contrast_variable", "reference_level", "target_level", "output_prefix")
}
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string))]
if (length(missing) > 0) stop(paste("Missing required options:", paste(missing, collapse = ", ")))

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(edgeR)
library(limma) # for makeContrasts()

################################################
################################################
## READ IN COUNTS AND SAMPLE METADATA         ##
################################################
################################################

count.table            <- read_delim_flexible(opt\$count_file, header = TRUE, row.names = NULL, check.names = FALSE)
rownames(count.table)   <- count.table[[opt\$gene_id_col]]
count.table[[opt\$gene_id_col]] <- NULL

sample.sheet <- read_delim_flexible(opt\$sample_file, header = TRUE, check.names = FALSE)
if (!opt\$sample_id_col %in% colnames(sample.sheet)) {
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}
sample.sheet <- sample.sheet[!duplicated(sample.sheet[[opt\$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

missing_samples <- setdiff(rownames(sample.sheet), colnames(count.table))
if (length(missing_samples) > 0) {
    stop(paste("Samples in sample sheet missing from count table:", paste(missing_samples, collapse = ",")))
}
count.table <- count.table[, rownames(sample.sheet), drop = FALSE]

################################################
################################################
## BUILD THE DESIGN                           ##
################################################
################################################

if (is_valid_string(opt\$formula)) {
    message("Using user-specified formula: ", opt\$formula)
    model <- opt\$formula
} else {
    contrast_variable <- make.names(opt\$contrast_variable)
    if (!contrast_variable %in% colnames(sample.sheet)) {
        stop(paste0("Chosen contrast variable '", contrast_variable, "' not in sample sheet"))
    } else if (any(!c(opt\$reference_level, opt\$target_level) %in% sample.sheet[[contrast_variable]])) {
        stop(paste("Please choose reference and target levels present in the", contrast_variable, "column"))
    }
    sample.sheet[[contrast_variable]] <- as.factor(sample.sheet[[contrast_variable]])
    # Cell-means (no intercept) design so the contrast below is explicit --
    # same rationale as ../../deseq2/templates/deseq2_differential.R.
    model <- paste0("~ 0 + ", contrast_variable)
}
message("Final design formula: ", model)

design <- model.matrix(as.formula(model), data = sample.sheet)
colnames(design) <- make.names(colnames(design))

################################################
################################################
## RUN EDGER (TMM + quasi-likelihood F-test)  ##
################################################
################################################

dge <- DGEList(counts = round(as.matrix(count.table)))
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

if (is_valid_string(opt\$contrast_string)) {
    if (opt\$contrast_string %in% colnames(design)) {
        qlf <- glmQLFTest(fit, coef = opt\$contrast_string)
    } else {
        contrast_vec <- limma::makeContrasts(contrasts = opt\$contrast_string, levels = colnames(design))
        qlf <- glmQLFTest(fit, contrast = contrast_vec)
    }
} else {
    contrast_expr <- paste0(contrast_variable, opt\$target_level, " - ", contrast_variable, opt\$reference_level)
    contrast_vec  <- limma::makeContrasts(contrasts = contrast_expr, levels = colnames(design))
    qlf <- glmQLFTest(fit, contrast = contrast_vec)
}

################################################
################################################
## OUTPUTS                                    ##
################################################
################################################

if (!is_valid_string(opt\$contrast_string)) {
    opt\$contrast_string <- paste(opt\$target_level, opt\$reference_level, sep = "_vs_")
}
cat("Saving results for ", opt\$contrast_string, " ...\n", sep = "")

res <- topTags(qlf, n = Inf, sort.by = "none")\$table
res <- data.frame(gene_id = rownames(res), res, check.names = FALSE)
write.table(
    res,
    file = paste0(opt\$output_prefix, ".edger.results.tsv"),
    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)

pdf(file = paste0(opt\$output_prefix, ".edger.bcv.pdf"))
plotBCV(dge)
dev.off()

norm_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
norm_counts <- data.frame(gene_id = rownames(norm_counts), norm_counts, check.names = FALSE)
write.table(
    norm_counts,
    file = paste0(opt\$output_prefix, ".normalised_counts.tsv"),
    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)

write(model, file = paste0(opt\$output_prefix, ".edger.model.txt"))

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version      <- paste(R.version[['major']], R.version[['minor']], sep = ".")
edger.version  <- as.character(packageVersion('edgeR'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-edger:', edger.version)
    ),
    'versions.yml')
