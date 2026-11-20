# OSTA 章节对照（maintainers）

[OSTA](https://bioconductor.org/books/OSTA/)（Orchestrating Spatial Transcriptomics Analysis
with Bioconductor，Crowell et al., bioRxiv 2025, doi:10.1101/2025.11.20.688607）是 R/Bioconductor
的空间分析书，v1.2.2 共 40 章。本插件是 **Python 单路**（scanpy + squidpy + spatialdata），
把 OSTA 的分析步骤与原理翻译到 Python 生态；R 包名只作概念锚点，不是依赖。许可 CC BY 4.0
（待最终核对 DESCRIPTION）。知识原创蒸馏，不复制上游文本/图/数据。

## 背景（1–8）

| OSTA 章 | R 锚点 | 本插件实现 |
| --- | --- | --- |
| 4 Infrastructure | SpatialExperiment / SingleCellExperiment | AnnData（`obsm["spatial"]` 坐标、`uns["spatial"]` 图像）；多样本/多尺度用 spatialdata `SpatialData` |
| 6 Importing data | TENxVisiumData, SpatialExperiment 读取 | `load` 片段 + `spatialdata-io`（visium/visium_hd）；Stereo-seq GEM/GEF 自研 reader |
| 8 Python interoperability | reticulate 调 Python | 本插件本身即 Python 路，无需互操作层 |

## 测序型平台（9–17）

| OSTA 章 | R 锚点 | 本插件实现 |
| --- | --- | --- |
| 10 Reads to counts | Cell Ranger / spaceranger | wrapper 树（不在 notebook） |
| 11 Quality control | per-spot QC, 组织覆盖 | `qc-metrics`/`qc-plots`/`filter`（空间坐标上的分布与组织图） |
| 12 Intermediate processing | normalize/HVG/PCA/cluster | `normalize`/`hvg`/`scale`/`pca`/`cluster` |
| 13 Deconvolution | RCTD(spacexr), cell2location, SPOTlight, stereoscope, Tangram | `deconvolve`（scipy NNLS，需 scRNA 参考；cell2location/tangram 是项目环境缺口） |
| 14–15 Visium DLPFC / CRC workflow | 完整工作流 | presets（Visium 标准链） |
| 16–17 Visium HD bin / seg | square bin / 细胞分割 | `load` 识别 HD bin 与 seg 输出，分析单元分别为 bin/细胞 |

## 成像型平台（18–26）— v0.1 不覆盖

Xenium/CosMx/MERFISH 不在 v0.1（仅 Visium + Stereo-seq）。但其 **邻域分析（22）** 与
**细胞通讯（23）** 的概念被测序型复用：`nhood`、`communication`。分割（19）对应
Stereo-seq cellbin——分割由 SAW 上游完成，插件只读输出。亚细胞分析（24）、图像深度学习
（33–34）是缺口。

## 平台无关分析（27–34）

| OSTA 章 | R 锚点 | 本插件实现 |
| --- | --- | --- |
| 27 Normalization | logNormCounts, sctransform | `normalize`（library-size + log1p） |
| 28 Dimensionality reduction | PCA, UMAP | `pca`/`cluster`（UMAP） |
| 29 Clustering & annotation | Leiden, marker panel | `cluster`/`markers`/`annotate`/`assign-labels` |
| 30 Feature selection & testing | nnSVG, SPARK, SpatialDE（SVG） | `svg`（squidpy `spatial_autocorr`：Moran's I / Geary's C） |
| 31 Feature sets & signatures | panel scoring, GSVA | `annotate`（panel-mean-minus-rest） |
| 32 Spatial statistics | 邻域富集, 共现, Moran | `nhood`（squidpy `nhood_enrichment`/`co_occurrence`） |
| 33–34 Image analysis / deep learning | 图像分割/深度学习 | 缺口（项目环境） |

## 多样本分析（35–37）

| OSTA 章 | R 锚点 | 本插件实现 |
| --- | --- | --- |
| 35 Differential spatial patterns | 跨条件 SVG 差异 | `spatial-diff`（条件间空间模式差异） |
| 36 Differential colocalization | 跨条件共定位差异 | `spatial-diff`（共定位对比） |
| 37 Structure-based analysis | 组织结构比较 | 部分由 `nhood` 覆盖；其余缺口 |

## 跨平台分析（38–40）— 缺口

Registration（38）、Imputation（39）、Xenium×Visium（40）不在 v0.1，列为后续/项目环境。

## 与 single-cell 插件的分工

OSCA（单细胞）的 QC/normalize/cluster/marker/pseudobulk 语义已在 single-cell 插件实现；
本插件复用其阈值语义与 panel 打分规则，新增的是**空间维度**：坐标、邻接图、SVG、邻域、
deconvolution、通讯、跨样本空间模式。deconvolution 的 scRNA 参考来自 single-cell 插件。