# Three metrics decide most QC calls:
#   nFeature_RNA  genes detected per cell  (low = empty/damaged droplet)
#   nCount_RNA    UMIs per cell            (low = same, high = possible doublet)
#   percent.mt    mitochondrial fraction   (high = stressed/dying cell)
# Cytoplasmic RNA leaks out of a torn membrane while mitochondria stay inside,
# so a high mitochondrial fraction marks damage, not biology. Ribosomal and
# hemoglobin fractions are context: they vary by cell type and by contamination.
# mt_pattern comes from the config cell: "^MT-" for human, "^mt-" for mouse.
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
qc_scored <- "percent.mt"

ribo_ok <- tryCatch({
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
  TRUE
}, error = function(e) FALSE)
if (ribo_ok) qc_scored <- c(qc_scored, "percent.ribo")

hb_ok <- tryCatch({
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = hb_pattern)
  TRUE
}, error = function(e) FALSE)
if (hb_ok) qc_scored <- c(qc_scored, "percent.hb")

qc <- obj[[]][, c("nFeature_RNA", "nCount_RNA", intersect(qc_scored, colnames(obj[[]]))), drop = FALSE]
summary(qc)

# The quantiles are what a threshold should be based on, not a round number.
sapply(qc, function(x) round(quantile(x, c(0.01, 0.05, 0.5, 0.95, 0.99), na.rm = TRUE), 3))