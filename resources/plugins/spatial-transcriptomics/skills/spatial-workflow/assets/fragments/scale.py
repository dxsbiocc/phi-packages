# PCA weights every gene equally unless each is scaled to unit variance; without
# scaling, a handful of abundant genes define the axes. Clipping at 10 stops one
# outlier unit from stretching an axis.
sc.pp.scale(adata, max_value=10)
analysis_steps.append("scale")
print("scaled to unit variance, clipped at 10")