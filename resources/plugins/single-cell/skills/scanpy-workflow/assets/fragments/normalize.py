# Library size differs by orders of magnitude between droplets, so counts are
# not comparable until each cell sums to the same total. normalize_total rescales
# each cell to scale_factor counts; log1p then stabilizes the variance so a few
# highly expressed genes do not dominate every distance.
sc.pp.normalize_total(adata, target_sum=scale_factor)
sc.pp.log1p(adata)
# Snapshot before scaling: marker tests and panel scoring must read un-scaled
# log expression, and counts keeps the raw UMIs for HVG and pseudobulk.
adata.layers["lognorm"] = adata.X.copy()
analysis_steps.append("normalize")
print(f"normalized to {scale_factor:g} counts per cell and log1p-transformed; "
      "counts keeps raw UMIs, lognorm keeps the pre-scale snapshot")