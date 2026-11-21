# A peak x sample count matrix, the alternative entry point when counting already
# happened upstream (featureCounts, or a DiffBind dba.count export). Rows are
# peaks, columns are samples, the first column holds peak ids and the header row
# holds sample names.
path <- if (!is.null(count_path) && nzchar(count_path)) count_path else data_path
if (!file.exists(path)) stop("count table not found: ", path, call. = FALSE)

separator <- if (grepl("[.]csv$", tolower(path))) "," else "\t"
frame <- read.delim(path, sep = separator, row.names = 1, check.names = FALSE)
count_matrix <- as.matrix(frame)
storage.mode(count_matrix) <- "numeric"
if (anyNA(count_matrix)) {
  stop("count table has non-numeric entries; expected integer read counts", call. = FALSE)
}

cat(sprintf("%d peaks x %d samples\n", nrow(count_matrix), ncol(count_matrix)))
cat("samples:", paste(colnames(count_matrix), collapse = ", "), "\n")
cat("library sizes:", paste(format(colSums(count_matrix), big.mark = ","), collapse = ", "), "\n")
analysis_steps <- c(analysis_steps, "load-counts")