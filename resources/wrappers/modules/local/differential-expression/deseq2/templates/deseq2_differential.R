#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string, and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid, non-empty, non-whitespace string; FALSE otherwise.
#' @examples
#' is_valid_string("Hello World") # Returns TRUE
#' is_valid_string("   ")         # Returns FALSE
#' is_valid_string(NULL)          # Returns FALSE

is_valid_string <- function(input) {
  !is.null(input) && nzchar(trimws(input))
}

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
  args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
  args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

  parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE, colClasses = NA){

  ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

  if (ext == "tsv" || ext == "txt") {
    separator <- "\\t"
  } else if (ext == "csv") {
    separator <- ","
  } else {
    stop(paste("Unknown separator for", ext))
  }

  read.delim(
    file,
    sep = separator,
    header = header,
    row.names = row.names,
    check.names = check.names,
    colClasses = colClasses
  )
}

#
#' Turn “null” or empty strings into actual NULL
#'
#' @param x Input option
#'
#' @return NULL or x
#'
nullify <- function(x) {
  if (is.character(x) && (tolower(x) == "null" || x == "")) NULL else x
}

#' Run DESeq2 differential expression results with shared settings
#'
#' A thin wrapper around DESeq2::results() that injects the common parameters
#' (lfcThreshold, altHypothesis, independentFiltering, alpha, pAdjustMethod, minmu)
#' from the global `opt` list and the `dds` object. Use the `...` to supply the
#' branch-specific argument (`name = …` or `contrast = …`).
#'
#' @param ... Additional arguments passed to `results()`, e.g. `name = opt\$contrast_string`
#'   or `contrast = c(variable, target, reference)`.
#' @return A `DESeqResults` object containing the differential expression test results.
#'
run_results <- function(...) {
  results(
    dds,
    lfcThreshold         = opt\$lfc_threshold,
    altHypothesis        = opt\$alt_hypothesis,
    independentFiltering = opt\$independent_filtering,
    alpha                = opt\$alpha,
    pAdjustMethod        = opt\$p_adjust_method,
    minmu                = opt\$minmu,
    ...
  )
}
#' Shrink log2 fold changes via ASHR with shared settings
#'
#' A thin wrapper around DESeq2::lfcShrink() that injects the common parameters
#' (type = "ashr") and the `dds` object. Use the `...` to supply the branch-
#' specific argument (`coef = …` or `contrast = …`).
#'
#' @param ... Additional arguments passed to `lfcShrink()`, e.g. `coef = opt\$contrast_string`
#'   or `contrast = c(variable, target, reference)`.
#' @return A `DESeqResults` object with shrunken log2 fold changes.
run_shrink <- function(...) {
  lfcShrink(
    dds,
    type = 'ashr',
    ...
  )
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
  output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
  count_file = '$counts',
  sample_file = '$samplesheet',
  contrast_variable = '$contrast_variable',
  reference_level = '$reference',
  target_level = '$target',
  contrast_string    = "$comparison",             # Optional full (complex) contrast expression comparison
  formula            = "$formula",              # User-specified formula
  blocking_variables = NULL,
  study_design = 'independent',
  donor_col = NULL,
  control_genes_file = '$control_genes_file',
  transcript_lengths_file = '$transcript_lengths_file',
  sizefactors_from_controls = FALSE,
  gene_id_col = "gene_id",
  sample_id_col = "experiment_accession",
  subset_to_contrast_samples = FALSE,
  exclude_samples_col = NULL,
  exclude_samples_values = NULL,
  test = "Wald",
  fit_type = "parametric",
  sf_type = 'ratio',
  min_replicates_for_replace = 7,
  use_t = FALSE,
  lfc_threshold = 0,
  alt_hypothesis = 'greaterAbs',
  independent_filtering = TRUE,
  p_adjust_method = 'BH',
  alpha = 0.1,
  # Post-hoc |log2FoldChange| cutoff. padj uses alpha (also DESeq2's independent-filtering target).
  log2fc_threshold = 1,
  minmu = 0.5,
  vs_method = 'vst', # 'rlog', 'vst', or 'rlog,vst'
  shrink_lfc = TRUE,
  cores = 1,
  vs_blind = TRUE,
  vst_nsub = 1000,
  round_digits = NULL,
  seed = NULL
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
  if (! ao %in% names(opt)){
    stop(paste("Invalid option:", ao))
  }else{

    # Preserve classes from defaults where possible
    if (! is.null(opt[[ao]])){
      args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
    opt[[ao]] <- args_opt[[ao]]
  }
}


if ( ! is.null(opt\$round_digits)){
  opt\$round_digits <- as.numeric(opt\$round_digits)
}

# If there is no option supplied, convert string "null" to NULL
keys <- c("formula", "contrast_string", "contrast_variable", "reference_level", "target_level", "seed", "blocking_variable", "transcript_lengths_file", "donor_col")
opt[keys] <- lapply(opt[keys], nullify)

if ( ! is.null(opt\$seed)){
  set.seed(opt\$seed)
}

# Validate design mode and threshold ranges (also enforced by the wrapper
# preflight; repeated here so bypassing the wrapper still fails safely).
if (! opt\$study_design %in% c('independent', 'paired')) {
  stop(paste0("study_design must be 'independent' or 'paired' (got '", opt\$study_design, "')"), call. = FALSE)
}
if (! is.finite(opt\$alpha) || opt\$alpha <= 0 || opt\$alpha > 1) {
  stop(paste0("padj_threshold (alpha) must be finite and in (0, 1] (got '", opt\$alpha, "')"), call. = FALSE)
}
if (! is.finite(opt\$log2fc_threshold) || opt\$log2fc_threshold < 0) {
  stop(paste0("log2fc_threshold must be finite and >= 0 (got '", opt\$log2fc_threshold, "')"), call. = FALSE)
}

if (opt\$study_design == 'paired') {
  if (is_valid_string(opt\$formula) || is_valid_string(opt\$contrast_string)) {
    stop("study_design='paired' supports only the simple contrast_variable/reference_level/target_level model; remove formula/comparison or set study_design='independent'", call. = FALSE)
  }
  if (! is_valid_string(opt\$donor_col)) {
    stop("study_design='paired' requires donor_col naming the biological-donor column of the sample sheet; there is no implicit sample-as-donor mapping", call. = FALSE)
  }
  if (! grepl('^[A-Za-z.][A-Za-z0-9._]*\$', opt\$donor_col)) {
    stop(paste0("donor_col must be a plain column identifier (letters, digits, '.', '_'; got '", opt\$donor_col, "')"), call. = FALSE)
  }
}

# Check if required parameters have been provided
if (is_valid_string(opt\$formula)) {
  required_opts <- c('output_prefix', 'contrast_string')
} else {
  required_opts <- c('contrast_variable', 'reference_level', 'target_level', 'output_prefix')
}
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]

if (length(missing) > 0){
  stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('count_file', 'sample_file')){
  if (! is_valid_string(opt[[file_input]])) {
    stop(paste("Please provide", file_input), call. = FALSE)
  }

  if (! file.exists(opt[[file_input]])){
    stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
  }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(DESeq2)
library(BiocParallel)
library(limma)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

count.table <-
  read_delim_flexible(
    file = opt\$count_file,
    header = TRUE,
    row.names = opt\$gene_id_col,
    colClasses = setNames("character", opt\$gene_id_col),
    check.names = FALSE
  )
# Keep identifiers such as 001 distinct from 1, including donor and group labels.
id_columns <- opt\$sample_id_col
if (!is_valid_string(opt\$formula)) id_columns <- c(id_columns, opt\$contrast_variable)
if (opt\$study_design == 'paired') id_columns <- c(id_columns, opt\$donor_col)
sample.sheet <- read_delim_flexible(file = opt\$sample_file,
  colClasses = setNames(rep("character", length(id_columns)), id_columns))

# Deal with spaces that may be in sample column
opt\$sample_id_col <- make.names(opt\$sample_id_col)

if (! opt\$sample_id_col %in% colnames(sample.sheet)){
  stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names. Paired designs instead require unique sample ids.

if (opt\$study_design == 'paired') {
  paired_ids <- sample.sheet[[opt\$sample_id_col]]
  if (any(is.na(paired_ids) | !nzchar(trimws(as.character(paired_ids))))) {
    stop(paste0("Paired design requires a non-empty sample id in every sample sheet row (column '", opt\$sample_id_col, "')"), call. = FALSE)
  }
  paired_dups <- unique(paired_ids[duplicated(paired_ids)])
  if (length(paired_dups) > 0) {
    stop(paste0("Paired design requires unique sample ids; duplicated: ", paste(paired_dups, collapse = ', ')), call. = FALSE)
  }
}

sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt\$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

# Check that all samples specified in the input sheet are present in the counts
# table. Assuming they are, subset and sort the count table to match the sample
# sheet

missing_samples <-
  sample.sheet[!rownames(sample.sheet) %in% colnames(count.table), opt\$sample_id_col]

if (length(missing_samples) > 0) {
  stop(paste(
    length(missing_samples),
    'specified samples missing from count table:',
    paste(missing_samples, collapse = ',')
  ))
} else{
  # Save any non-count data, will gene metadata etc we might need later
  noncount.table <-
    count.table[, !colnames(count.table) %in% rownames(sample.sheet), drop = FALSE]
  count.table <- count.table[, rownames(sample.sheet)]
}

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################
if (! is_valid_string(opt\$formula)) {
  contrast_variable <- make.names(opt\$contrast_variable)
  blocking.vars <- c()

  if (!contrast_variable %in% colnames(sample.sheet)) {
    stop(
      paste0(
        'Chosen contrast variable \"',
        contrast_variable,
        '\" not in sample sheet'
      )
    )
  } else if (any(!c(opt\$reference_level, opt\$target_level) %in% sample.sheet[[contrast_variable]])) {
    stop(
      paste(
        'Please choose reference and treatment levels that are present in the',
        contrast_variable,
        'column of the sample sheet'
      )
    )
  } else if (is_valid_string(opt\$blocking_variables)) {
    blocking.vars = make.names(unlist(strsplit(opt\$blocking_variables, split = ';')))
    if (!all(blocking.vars %in% colnames(sample.sheet))) {
      missing_block <- paste(blocking.vars[! blocking.vars %in% colnames(sample.sheet)], collapse = ',')
      stop(
        paste(
          'Blocking variables', missing_block,
          'do not correspond to sample sheet columns.'
        )
      )
    }
  }

  # Optionally, subset to only the samples involved in the contrast
  if (opt\$subset_to_contrast_samples){
    sample_selector <- sample.sheet[[contrast_variable]] %in% c(opt\$target_level, opt\$reference_level)
    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    count.table <- count.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
  }

}
# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((is_valid_string(opt\$exclude_samples_col)) && (is_valid_string(opt\$exclude_samples_values))){
  exclude_values = unlist(strsplit(opt\$exclude_samples_values, split = ';'))

  if (! opt\$exclude_samples_col %in% colnames(sample.sheet)){
    stop(paste(opt\$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
  }

  print(paste0('Excluding samples with values of ', opt\$exclude_samples_values, ' in ', opt\$exclude_samples_col))
  sample_selector <- ! sample.sheet[[opt\$exclude_samples_col]] %in% exclude_values

  selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
  count.table <- count.table[, selected_samples]
  sample.sheet <- sample.sheet[selected_samples, ]
}

# Paired design: validate the donor mapping before any DESeq2 fitting and block
# on donor. Runs after subsetting/exclusion so the checks see the final samples.
if (opt\$study_design == 'paired') {
  donor_col <- make.names(opt\$donor_col)
  if (! donor_col %in% colnames(sample.sheet)) {
    stop(paste0("donor_col '", opt\$donor_col, "' is not a column of the sample sheet"), call. = FALSE)
  }
  if (donor_col == contrast_variable) {
    stop("donor_col must differ from contrast_variable", call. = FALSE)
  }

  # The paired contrast uses only the two contrast groups.
  paired_selector <- sample.sheet[[contrast_variable]] %in% c(opt\$reference_level, opt\$target_level)
  sample.sheet <- sample.sheet[paired_selector, , drop = FALSE]
  count.table <- count.table[, rownames(sample.sheet), drop = FALSE]

  donor_values <- trimws(as.character(sample.sheet[[donor_col]]))
  if (any(is.na(donor_values) | !nzchar(donor_values) | tolower(donor_values) %in% c('na', 'null'))) {
    stop(paste0("Incomplete donor mapping: every contrast sample needs a donor value in '", opt\$donor_col, "'"), call. = FALSE)
  }
  pair_table <- table(donor_values, as.character(sample.sheet[[contrast_variable]]))
  if (! all(c(opt\$reference_level, opt\$target_level) %in% colnames(pair_table))) {
    stop(paste0("Paired design needs samples in both the '", opt\$reference_level, "' and '", opt\$target_level, "' groups"), call. = FALSE)
  }
  if (any(pair_table != 1)) {
    bad_donors <- rownames(pair_table)[apply(pair_table, 1, function(r) any(r != 1))]
    stop(paste0(
      "Paired design requires exactly one distinct biological sample per donor per group, with every donor in both groups; unbalanced donors: ",
      paste(bad_donors, collapse = ', ')
    ), call. = FALSE)
  }
  if (nrow(pair_table) < 2) {
    stop(paste0("Paired design requires at least two donors (found ", nrow(pair_table), ")"), call. = FALSE)
  }

  sample.sheet[[donor_col]] <- donor_values
  blocking.vars <- unique(c(donor_col, blocking.vars))
}

# Now specify the model. Use cell-means style so we can be explicit with the
# contrasts

if (is_valid_string(opt\$formula)) {
  message("Using user-specified formula: ", opt\$formula)
  model  <- opt\$formula
}  else {
  # Make sure all the appropriate variables are factors

  for (v in c(blocking.vars, contrast_variable)) {
    sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
  }

  if (opt\$study_design == 'paired') {
    # Estimate the between-condition change within donors.
    model <- paste('~', paste(c(blocking.vars, contrast_variable), collapse = ' + '))
  } else {
    model <- '~ 0'

    if (is_valid_string(opt\$blocking_variables)) {
      model <- paste(model, paste(blocking.vars, collapse = ' + '), sep=' + ')
    }

    # Variable of interest goes last, see
    # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs

    model <- paste(model, contrast_variable, sep = ' + ')
  }
}
message("Final design formula: ", model)

################################################
################################################
## Run DESeq2 processes                       ##
################################################
################################################

if (opt\$control_genes_file != ''){
  control_genes <- readLines(opt\$control_genes_file)
  if (! opt\$sizefactors_from_controls){
    count.table <- count.table[setdiff(rownames(count.table), control_genes),]
  }
}

dds <- DESeqDataSetFromMatrix(
  countData = round(count.table),
  colData = sample.sheet,
  design = as.formula(model)
)

# Build in transcript lengths. Copying what tximport does here:
# https://github.com/thelovelab/DESeq2/blob/6947d5bc629015fb8ffb2453a91b71665a164483/R/AllClasses.R#L409

if (!is.null(opt\$transcript_lengths_file)){
  lengths <-
    read_delim_flexible(
      file = opt\$transcript_lengths_file,
      header = TRUE,
      row.names = opt\$gene_id_col,
      check.names = FALSE
    )
  lengths <- lengths[rownames(count.table), colnames(count.table)]
  dimnames(lengths) <- dimnames(dds)
  assays(dds)[["avgTxLength"]] <- lengths
}

if (opt\$control_genes_file != '' && opt\$sizefactors_from_controls){
  print(paste('Estimating size factors using', length(control_genes), 'control genes'))
  dds <- estimateSizeFactors(dds, controlGenes=rownames(count.table) %in% control_genes)
}

dds <- DESeq(
  dds,
  test = opt\$test,
  fitType = opt\$fit_type,
  minReplicatesForReplace = opt\$min_replicates_for_replace,
  useT = opt\$use_t,
  sfType = opt\$sf_type,
  parallel=TRUE, BPPARAM=MulticoreParam(opt\$cores)
)

if (!is.null(opt\$contrast_string)) {
  coef_names <- resultsNames(dds)
  if (opt\$contrast_string %in% coef_names) {
    # Direct coefficient name
    comp.results <- run_results(name = opt\$contrast_string)
    if (opt\$shrink_lfc) {
      comp.results <- run_shrink(res = comp.results, coef = opt\$contrast_string)
    }
  } else {
    # Parse as limma-style contrast expression
    design_mat <- model.matrix(as.formula(model), data = as.data.frame(colData(dds)))
    colnames(design_mat) <- make.names(colnames(design_mat))
    numeric_contrast <- as.numeric(
      limma::makeContrasts(contrasts = opt\$contrast_string, levels = colnames(design_mat))
    )

    # Run DESeq2 results with numeric contrast
    comp.results <- run_results(contrast = numeric_contrast)
    if (opt\$shrink_lfc) {
      comp.results <- run_shrink(res = comp.results, contrast = numeric_contrast)
    }
  }
} else {
  contrast_var_tg_ref <- c(contrast_variable,
                           opt\$target_level,
                           opt\$reference_level)

  comp.results <- run_results(contrast = contrast_var_tg_ref)
  if (opt\$shrink_lfc) {
    comp.results <- run_shrink(res = comp.results, contrast = contrast_var_tg_ref)
  }
}

# See https://support.bioconductor.org/p/97676/

if (!is.null(opt\$transcript_lengths_file)){
  size_factors = estimateSizeFactorsForMatrix(counts(dds) / assays(dds)[["avgTxLength"]])
} else {
  size_factors = sizeFactors(dds)
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

# if no explicit contrast_string was given, synthesize a two‐level default
if (! is_valid_string(opt\$contrast_string)) {
  opt\$contrast_string <- paste(opt\$target_level, opt\$reference_level, sep = "_vs_")
}

cat("Saving results for ", opt\$contrast_string, " ...\n", sep = "")

# Differential expression table- note very limited rounding for consistency of
# results

if (! is.null(opt\$round_digits)){
  comp.results <- apply(data.frame(comp.results), 2, function(x) round(x, opt\$round_digits))
}
comp.results <- `colnames<-`(
  data.frame(
    gene_id = rownames(comp.results),
    comp.results,
    check.names = FALSE
  ),
  c(opt\$gene_id_col, colnames(comp.results))  # Setting all column names
)

# QC contract: reject malformed/NA identifiers, duplicates, and non-NA invalid
# numeric values before anything is published.
qc_gene_ids <- as.character(comp.results[[opt\$gene_id_col]])
if (any(is.na(qc_gene_ids)) || any(!nzchar(trimws(qc_gene_ids))) || any(grepl('[\t\r\n]', qc_gene_ids))) {
  stop("Results contain missing, empty, or multi-line/tab-containing gene identifiers", call. = FALSE)
}
if (anyDuplicated(qc_gene_ids)) {
  stop("Results contain duplicated gene identifiers", call. = FALSE)
}
qc_padj <- comp.results[['padj']]
qc_log2fc <- comp.results[['log2FoldChange']]
if (any(!is.na(qc_padj) & (!is.finite(qc_padj) | qc_padj < 0 | qc_padj > 1))) {
  stop("Results contain non-NA padj values outside [0, 1] or non-finite", call. = FALSE)
}
if (any(!is.na(qc_log2fc) & !is.finite(qc_log2fc))) {
  stop("Results contain non-NA log2FoldChange values that are not finite", call. = FALSE)
}

results_file <- paste(opt\$output_prefix, 'deseq2.results.tsv', sep = '.')
write.table(
  comp.results,
  file = results_file,
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE
)

# Differentially expressed genes: |log2FoldChange| >= log2fc_threshold and padj < alpha.
# Rows with NA fold change or padj are dropped. Up/down split on the sign of
# log2FoldChange (exact zeros stay only in the combined table).
# Compute the screen from exactly the values readers will see in the TSV.
published <- read.delim(results_file, check.names = FALSE, quote = "", colClasses = "character")
log2fc <- as.numeric(published[['log2FoldChange']])
padj <- as.numeric(published[['padj']])
de_mask <- !is.na(log2fc) & !is.na(padj) & abs(log2fc) >= opt\$log2fc_threshold & padj < opt\$alpha
de_genes <- comp.results[de_mask, , drop = FALSE]
filtered_tables <- list(
  list(de_genes, 'deseq2.de_genes_filtered.tsv'),
  list(de_genes[de_genes[['log2FoldChange']] > 0, , drop = FALSE], 'deseq2.de_genes_filtered_up.tsv'),
  list(de_genes[de_genes[['log2FoldChange']] < 0, , drop = FALSE], 'deseq2.de_genes_filtered_down.tsv')
)
for (de_table in filtered_tables) {
  write.table(
    de_table[[1]],
    file = paste(opt\$output_prefix, de_table[[2]], sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
  )
}

# QC sidecar: one header + one record, counts against the published rows with
# the exact thresholds used above. Only rows with finite (non-NA) padj and
# log2FoldChange are evaluable; all-NA results are 'invalid'.
qc_evaluated <- sum(!is.na(padj) & !is.na(log2fc))
qc_hits <- sum(de_mask)
qc_disposition <- if (qc_evaluated == 0) 'invalid' else if (qc_hits > 0) 'eligible' else 'skip'

if (is_valid_string(opt\$formula)) {
  qc_contrast <- paste(opt\$formula, opt\$contrast_string, sep = ' | ')
} else {
  qc_contrast <- paste(opt\$target_level, opt\$reference_level, sep = ' vs ')
}
# Single-line, tab-free, bounded strings for the sidecar.
qc_field <- function(x, max_chars = 200) {
  if (length(x) != 1 || is.na(x) || !nzchar(trimws(x)) || nchar(x) > max_chars || grepl('[[:cntrl:]]', x)) {
    stop("Contrast or gene identifier column cannot be represented by the QC contract.")
  }
  as.character(x)
}

qc_df <- data.frame(
  schemaVersion = 1L,
  disposition = qc_disposition,
  contrast = qc_field(qc_contrast),
  geneIdColumn = qc_field(opt\$gene_id_col),
  padjThreshold = sprintf("%.17g", opt\$alpha),
  log2fcThreshold = sprintf("%.17g", opt\$log2fc_threshold),
  evaluated = qc_evaluated,
  hits = qc_hits,
  resultMd5 = unname(tools::md5sum(results_file)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.table(
  qc_df,
  file = paste(opt\$output_prefix, 'deseq2.qc.tsv', sep = '.'),
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE,
  fileEncoding = 'UTF-8'
)

# Dispersion plot

png(
  file = paste(opt\$output_prefix, 'deseq2.dispersion.png', sep = '.'),
  width = 600,
  height = 600
)
plotDispEsts(dds)
dev.off()

pdf(
  file = paste(opt\$output_prefix, 'deseq2.dispersion.pdf', sep = '.')
)
plotDispEsts(dds)
dev.off()

# R object for other processes to use

saveRDS(dds, file = paste(opt\$output_prefix, 'dds.rld.rds', sep = '.'))

# Size factors

sf_df = data.frame(
  sample = names(size_factors),
  data.frame(size_factors, check.names = FALSE),
  check.names = FALSE
)
colnames(sf_df) <- c('sample', 'sizeFactor')
write.table(
  sf_df,
  file = paste(opt\$output_prefix, 'deseq2.sizefactors.tsv', sep = '.'),
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE
)

# Write specified matrices

normalised_matrix <- counts(dds, normalized = TRUE)
if (! is.null(opt\$round_digits)){
  normalised_matrix <- apply(normalised_matrix, 2, function(x) round(x, opt\$round_digits))
}
normalised_matrix <- `colnames<-`(
  data.frame(
    gene_id = rownames(counts(dds)),  # First column with row names from counts(dds)
    normalised_matrix,                # Other columns
    check.names = FALSE
  ),
  c(opt\$gene_id_col, colnames(normalised_matrix))  # Setting all column names
)

write.table(
  normalised_matrix,
  file = paste(opt\$output_prefix, 'normalised_counts.tsv', sep = '.'),
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE
)

# Note very limited rounding for consistency of results

for (vs_method_name in strsplit(opt\$vs_method, ',')){
  if (vs_method_name == 'vst'){
    vs_mat <- assay(vst(dds, blind = opt\$vs_blind, nsub = opt\$vst_nsub))
  }else if (vs_method_name == 'rlog'){
    vs_mat <- assay(rlog(dds, blind = opt\$vs_blind, fitType = opt\$fit_type))
  }

  if (! is.null(opt\$round_digits)){
    vs_mat <- apply(vs_mat, 2, function(x) round(x, opt\$round_digits))
  }

  vs_mat <- `colnames<-`(
    data.frame(
      gene_id = rownames(counts(dds)),  # First column with row names from counts(dds)
      vs_mat,                           # Other columns from vs_mat
      check.names = FALSE
    ),
    c(opt\$gene_id_col, colnames(vs_mat))  # Setting all column names
  )

  write.table(
    vs_mat,
    file = paste(opt\$output_prefix, vs_method_name, 'tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
  )
}

# Save model to file

write(model, file=paste(opt\$output_prefix, 'deseq2.model.txt', sep = '.'))

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(opt\$output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
deseq2.version <- as.character(packageVersion('DESeq2'))

writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r.version),
    paste('    bioconductor-deseq2:', deseq2.version)
  ),
  'versions.yml')

################################################
################################################
################################################
################################################
