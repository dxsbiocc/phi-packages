# Scaling gives every gene mean 0 and variance 1, so highly expressed genes do
# not dominate PCA purely because their numbers are larger. vars.to.regress
# removes a covariate's linear effect first (use sparingly: it can also remove
# biology that happens to correlate with the covariate).
if (DefaultAssay(obj) == "SCT") {
  message("SCT assay already carries scaled residuals; skipping ScaleData.")
} else {
  obj <- if (length(regress_vars) > 0) {
    ScaleData(obj, vars.to.regress = regress_vars, verbose = FALSE)
  } else {
    ScaleData(obj, verbose = FALSE)
  }
  analysis_steps <- c(analysis_steps, "scale")
}

# Confirm: scaled values centre around zero and range symmetrically.
obj[[DefaultAssay(obj)]]$scale.data[1:5, 1:3]