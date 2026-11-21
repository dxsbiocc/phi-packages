# single-cell-atac

单细胞染色质可及性（scATAC）分析插件。**Python 单路**：snapatac2 + scanpy，私有锁定环境
`plugin:scatac`，工具前缀 `scatac_`。主引擎是 **snapatac2**（Rust 后端、file-backed AnnData、
plotly 图），scanpy 只用于 gene activity 矩阵的下游。分析单元是**细胞（barcode）**，不是 spot/bin。

| 输入起点 | 来源 | 说明 |
| --- | --- | --- |
| `fragments.tsv.gz` | 10x Cell Ranger ATAC | 首选，snapatac2 直接导入并算 QC/TSSe |
| BAM | 比对输出 | 先 `make_fragment_file` 转 fragment 再导入 |
| `.h5ad` | 已有 snapatac2/scanpy 对象 | backed 或内存态 |
| peak×cell 矩阵 | csv/tsv/h5/10x | 已带特征，跳过 tile 矩阵 |

**原始 reads 比对 / Cell Ranger ATAC 不在这里做**，那是 nf-core wrapper 的事。

## 设计原则

与 single-cell / spatial 插件同一套契约：

- **插件只做确定性填充/渲染**。意图解析、步骤选择、阈值决策归 agent（规则写在
  [`agents/SingleCellAtac.md`](agents/SingleCellAtac.md) 与技能的 SKILL.md）。
- **片段库是唯一分析源**。`skills/scatac-workflow/assets/fragments/` 的教学代码同时驱动
  脚本工具（薄适配器 + JSON 报告层）与 notebook 渲染器，两条形态共用一份代码。
- **私有锁定环境，不在里面装包**。缺包/缺基因组参考由 `scatac_env_probe` 当成事实报告，
  下一步是项目环境（env_request）或指引下载到项目目录。
- **notebook 是真代码格 + 每步 Markdown 教学**，渲染进用户工作目录，不是拷贝模板目录。

## snapatac2 的两个硬约束

- **backed AnnData**：对象背在 `.h5ad` 文件上自动同步，**进程退出前必须 `close()`**，否则
  HDF5 损坏。脚本链每步落盘并 close；notebook 末尾有 close 格。
- **plotly 图，不是 matplotlib**：`snap.pl.*` 返回 plotly Figure，脚本链存成自包含 `.html`；
  gene activity 下游走 scanpy 时才用 matplotlib。

## 分析链

`scatac_inspect` → `scatac_qc`（片段数/重复率/线粒体/TSSe + 阈值分布）→ `scatac_preprocess`
（tile 矩阵/特征选择/scrublet 双细胞/spectral 降维/Leiden 聚类/UMAP）→ `scatac_peaks`
（per-cluster MACS3 peak + 合并 peak 集 + peak-count 矩阵）→ `scatac_markers`（per-cluster
marker regions，可选 pairwise 差异可及性 `run_diff`）→ `scatac_annotate`（gene activity +
panel/scRNA 参考/显式 mapping 标签）→ `scatac_motif`（marker 区域 motif 富集）→
`scatac_integrate`（多样本批次整合，spectral 嵌入上 harmony/scanorama/mnc）。
`scatac_notebook` 是平行的教学路径；`scatac_env_probe` 在承诺重依赖能力前探测。

## 基因组参考

- **chrom_sizes 自带**：snapatac2 内置 hg38/mm10/hg19/mm39，离线可用。
- **blacklist 内联**：ENCODE hg38/mm10 blacklist 以 base64-gzip 内联在
  `assets/fragments/blacklist.py`，离线解码，不依赖外部文件。
- **基因注释（GFF3 ~50MB）与基因组 FASTA（~845MB）是懒下载**，不随环境分发：由用户提供
  `input-path`，或一次性下载到项目目录；缺失时报缺口，绝不在只读插件环境里下载大参考。
- **motif 库**（JASPAR/cis_bp）是 snapatac2 的下载资源，不随插件分发。

## 知识层

`skills/scatac-workflow/references/` 按功能组织：

- [`analysis-playbook.md`](skills/scatac-workflow/references/analysis-playbook.md) — agent
  的操作手册：每步调用哪个工具、读哪些 JSON 字段与图、怎么定阈值、红线是什么。
- [`red-line-cheatsheet.md`](skills/scatac-workflow/references/red-line-cheatsheet.md) —
  12 条硬红线 × 常见误读速查表，运行时对照。
- [`step-notes/`](skills/scatac-workflow/references/step-notes/) — 每步的目的/原理/解读/陷阱，
  即 notebook 的教学格。
- [`maintainers/implementation-notes.md`](skills/scatac-workflow/references/maintainers/implementation-notes.md)
  — 实现说明与溯源（片段库架构、backed AnnData 约束、snapatac2 API 映射、JSON 契约）。

知识原创蒸馏自 [SnapATAC2](https://github.com/kaizhang/SnapATAC2)（MIT）文档与
[sc-best-practices](https://www.sc-best-practices.org/)（Apache 2.0）；ENCODE blacklist 来源
doi:10.1038/s41598-019-45839-z。不抄上游文本/图/数据。

## 环境

`environments/scatac/environment.yml` 声明依赖（纯 conda，无 pip），`locks/` 是三平台
（darwin-arm64、darwin-x64、linux-x64）的锁定文件。求解走 tuna 镜像，再重写回官方 host。
macs3 作为 snapatac2 依赖随环境分发。

## 缺口（由 `scatac_env_probe` 报事实，走项目环境）

- **chromVAR/ArchR/Signac**：R/Bioconductor 路，未实现。插件的 TF 层是 motif 富集（序列预测），
  **不是 chromVAR 逐细胞 TF 偏差**，也不要把 motif 命中当实测 TF 结合。
- **episcanpy**：虽在 bioconda，但运行时 `from numpy import NaN` 与 numpy 2 不兼容，不随环境分发。
- **scvi-tools / magic-impute**：pip-only，锁定环境不支持 pip。
- **TF-gene GRN**：未实现（snapatac2 的 inject 扩展，不在本插件）。

## 与其他插件的关系

- **single-cell**：label transfer / gene-activity 注释对照的 scRNA 参考来自 SingleCell agent，
  走共享出口（counts+obs 交换 / pseudobulk 表）。gene activity 是可及性**代理**，不是实测 RNA。
- **wrappers**：原始 reads 到 fragment（Cell Ranger ATAC、比对）是 nf-core wrapper 的事，不是 notebook 格。
- **chromatin-profiling**（规划中的姊妹插件）：bulk ChIP/CUT&Tag 是另一插件（R 路、纯下游，
  从已 call 的 peak/counts 起步），与本插件按分析范式分开。