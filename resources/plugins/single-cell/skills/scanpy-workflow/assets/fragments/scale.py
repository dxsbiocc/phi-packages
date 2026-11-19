# PCA weights every gene equally unless each is scaled to unit variance; without
# scaling, a handful of abundant genes define the axes. Clipping at 10 stops one
# outlier cell from stretching an axis. regress_vars (optional) removes known
# technical covariates first — it is expensive and can remove biology, so the
# default is empty.
if regress_vars:
    missing = [v for v in regress_vars if v not in adata.obs.columns]
    if missing:
        raise ValueError("regress_vars not in obs: " + ", ".join(missing))
    sc.pp.regress_out(adata, regress_vars)
    analysis_steps.append("regress")
sc.pp.scale(adata, max_value=10)
analysis_steps.append("scale")
print("scaled to unit variance, clipped at 10")