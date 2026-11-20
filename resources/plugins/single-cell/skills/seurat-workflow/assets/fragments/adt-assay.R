# A CITE-seq object carries an "ADT" (or "HTO") assay beside RNA. If yours
# arrived as a separate matrix, attach it here.
if ("ADT" %in% names(obj@assays)) {
  DefaultAssay(obj) <- "ADT"
  # CLR centres each cell's protein counts on a log ratio; margin = 2 centres
  # across features, which is what the scNotebooks module uses.
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)

  # Protein and RNA evidence for the same marker should agree.
  candidates <- grep("^adt_", colnames(obj[[]]), value = TRUE)
  if (length(candidates) > 0) {
    show_plot(VlnPlot(obj, features = candidates[1], pt.size = 0), "adt-violin")
  } else {
    message("No adt_* metadata columns; inspect colnames(obj[[]]) for your naming.")
  }
  DefaultAssay(obj) <- "RNA"
} else {
  message("No ADT assay in this object - this notebook's RNA half still applies.")
}