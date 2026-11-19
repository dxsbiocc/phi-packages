# Moving objects between the R and Python routes

The plugin ships two private environments: `plugin:seurat` (R 4.5 + Seurat 6)
and `plugin:scanpy` (Python 3.13 + scanpy 1.13). Neither ships the dedicated
converters (`zellkonverter`, `SeuratDisk`), and plugin environments are
read-only — so conversion happens through formats both routes already read,
or through a project environment when full fidelity is required.

## What survives which bridge

| Bridge | Counts | obs metadata | Layers / reductions | Needs |
| --- | --- | --- | --- | --- |
| 10x-style directory | yes | via sidecar CSV | no | nothing extra |
| pseudobulk tables | aggregated | sample sheet | no | nothing extra (DE handoff) |
| `.h5ad` ↔ `.rds` converters | yes | yes | yes | project env with zellkonverter / SeuratDisk |

Prefer the 10x-style bridge for switching routes mid-analysis: both
`sc_inspect` and `sc_py_inspect` read it natively, and the counts layer — the
one thing QC and pseudobulk must not lose — travels intact.

## R → Python (Seurat object to a 10x-style directory)

Run in the Seurat route's environment (Matrix ships with it):

```r
library(Seurat); library(Matrix)
obj <- readRDS("analyzed.rds")
counts <- LayerData(obj, layer = "counts")   # genes x cells
dir.create("bridge", showWarnings = FALSE)
writeMM(counts, "bridge/matrix.mtx")
writeLines(colnames(counts), "bridge/barcodes.tsv")
write.table(
  data.frame(rownames(counts), rownames(counts), "Gene Expression"),
  "bridge/features.tsv", sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = FALSE
)
write.csv(obj[[]], "bridge/obs.csv")         # metadata sidecar
```

Then in the scanpy route:

```python
import scanpy as sc, pandas as pd
adata = sc.read_10x_mtx("bridge", var_names="gene_symbols")
obs = pd.read_csv("bridge/obs.csv", index_col=0)
adata.obs = obs.loc[adata.obs_names]
adata.layers["counts"] = adata.X.copy()
```

## Python → R (AnnData to a 10x-style directory)

```python
from scipy.io import mmwrite
import pandas as pd
counts = adata.layers["counts"]              # cells x genes
mmwrite("bridge/matrix.mtx", counts.T.tocsr())  # 10x layout: genes x cells
pd.Series(adata.var_names).to_csv("bridge/barcodes.tsv", index=False, header=False)
# features.tsv: gene_id, gene_name, type — same three columns as above
pd.DataFrame({
    "id": adata.var_names, "name": adata.var_names,
    "type": "Gene Expression",
}).to_csv("bridge/features.tsv", sep="\t", index=False, header=False)
adata.obs.to_csv("bridge/obs.csv")
```

Then `sc_inspect` reads the directory; attach the sidecar with
`obj <- AddMetaData(obj, read.csv("bridge/obs.csv", row.names = 1))`.

## Full fidelity (when layers and reductions must survive)

Request a project environment (`env_request`) with `zellkonverter` (R side,
reads/writes `.h5ad` via SingleCellExperiment) or `SeuratDisk` (R side,
`SaveH5Seurat` / `Convert`). The plugin environments themselves stay
untouched. This is the only bridge that keeps normalized layers, embeddings,
and graphs; the 10x-style bridge deliberately carries counts + metadata only.

## Shared handoff that needs no conversion

`sc_pseudobulk` and `sc_py_pseudobulk` write the same `gene_id`-first counts
table and `sample` sheet, so a project can run QC/clustering on one route and
the sample-level differential expression on the wrappers without any object
conversion at all.