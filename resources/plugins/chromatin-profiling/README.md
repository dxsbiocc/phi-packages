# chromatin-profiling

染色质 profiling（bulk **ChIP-seq** 与 **CUT&Tag/CUT&RUN**）下游分析插件。**R 单路**：
Bioconductor，私有锁定环境 `plugin:chip`，工具前缀 `chip_`。主引擎是 **DiffBind**
（occupancy + affinity 分析）+ **DESeq2/edgeR**（差异结合）+ **csaw**（broad mark 窗口法）+
**ChIPseeker**（peak 注释）。分析单元是 **peak/样本**，不是细胞。**纯下游**：从已 call 的
peak/counts/indexed BAM/bigWig 出发，绝不碰 FASTQ。

| 输入起点 | 来源 | 说明 |
| --- | --- | --- |
| BED/narrowPeak/broadPeak | MACS2（ChIP）/ SEACR（CUT&Tag） | 首选，per-sample peak 集 |
| DiffBind sample sheet CSV | 用户整理 | 驱动 DBA 对象与 contrast（Peaks/bamReads 列） |
| peak×sample count 表 | featureCounts / dba.count 导出 | 跳过 counting，直接 normalize/pca/diff |
| indexed BAM | 比对输出 | FRiP / counting / csaw（需 .bai） |
| bigWig | coverage track | profile / 可视化 |
| 已存 `.rds`（DBA/GRanges） | 本插件 save | 续跑 diff/annotate |

**比对 / peak calling 不在这里做**：ChIP 归 nf-core/chipseq wrapper；CUT&Tag 上游 SEACR
无 wrapper，是缺口，插件直接吃用户 SEACR 产出的 peak。

## 设计原则

与 single-cell / spatial / single-cell-atac 插件同一套契约：

- **插件只做确定性填充/渲染**。意图解析、步骤选择、阈值决策归 agent（规则写在
  [`agents/ChromatinProfiling.md`](agents/ChromatinProfiling.md) 与技能的 SKILL.md）。
- **片段库是唯一分析源**。`skills/chromatin-workflow/assets/fragments/` 的教学代码同时驱动
  脚本工具（薄适配器 + JSON 报告层）与 notebook 渲染器，两条形态共用一份代码。
- **私有锁定环境，不在里面装包**。缺包/缺基因组参考由 `chip_env_probe` 当成事实报告，
  下一步是项目环境（env_request）或指引下载到项目目录。
- **notebook 是真代码格 + 每步 Markdown 教学**，渲染进用户工作目录，不是拷贝模板目录。

## R 路的三个硬约束

- **DiffBind 需 sample sheet + indexed BAM**：`dba(sampleSheet=)` 从 CSV 建 DBA；counting/FRiP
  需 `.bai` 索引（`samtools index` 是上游）。相对路径相对 sheet 目录解析。
- **csaw 在 Apple Silicon 是缺口**：bioconda 无 osx-arm64 构建，csaw 只进 darwin-x64/linux-x64。
  darwin-arm64 上 broad mark 走 DiffBind 替代或换平台/项目环境，由 `chip_env_probe` 如实报告。
- **图是 PNG，不是 plotly**：ggplot2（PCA/MA）、ChIPseeker（特征分布/TSS 距离/profile/heatmap）、
  ComplexHeatmap（UpSet），ragg 渲染。

## 分析链

`chip_inspect` → `chip_qc`（FRiP + 重复性重叠分数 + blacklist 过滤）→ `chip_consensus`
（disjoin/reduce 建 consensus peakset，min_occupancy 决定 region universe）→ `chip_count`
（summarizeOverlaps 计数 + 归一化）→ `chip_diff`（DiffBind DESeq2/edgeR 用于 TF/narrow/CUT&Tag；
csaw 窗口法用于 broad mark）→ `chip_annotate`（ChIPseeker 基因组特征分布/nearest gene/TSS 距离 +
coverage profile/heatmap）→ `chip_visualize`（样本 PCA + peak 重叠 UpSet）。
`chip_notebook` 是平行的教学路径；`chip_env_probe` 在承诺重依赖能力前探测。

## 基因组参考（精简策略）

- **blacklist 内联**：ENCODE hg38（636 区域）/mm10（3435 区域）以 base64-gzip 内联在
  `assets/fragments/blacklist-filter.R`，离线解码；`genome` 选版本，`blacklist` 覆盖其他版本，
  两者都无则跳过并说明。blacklist 基因组版本特异——hg38 区域在 mm10 上是错的。
- **目标基因组 TxDb 来自用户**：ChIPseeker 注释用用户 GTF/GFF3 经 `makeTxDbFromGFF` 建 TxDb，
  或用户 TxDb `.sqlite`（`loadDb`）；**插件不 bundle 目标基因组 TxDb**，缺失报缺口。
  ChIPseeker conda 配方附带的 hg19 TxDb 只是构建期依赖，**绝不作目标基因组的默认注释源**。
- **gene symbol 需用户 org.db**（annoDb）；无则注释仍给 gene id/特征/TSS 距离，报缺口。
- 绝不在只读插件环境里装包或下载大参考。

## 知识层

`skills/chromatin-workflow/references/` 按功能组织：

- [`analysis-playbook.md`](skills/chromatin-workflow/references/analysis-playbook.md) — agent
  的操作手册：每步调用哪个工具、读哪些 JSON 字段与图、怎么定阈值、红线是什么。
- [`red-line-cheatsheet.md`](skills/chromatin-workflow/references/red-line-cheatsheet.md) —
  12 条硬红线 × 常见误读速查表，运行时对照。
- [`step-notes/`](skills/chromatin-workflow/references/step-notes/) — 每步的目的/原理/解读/陷阱，
  即 notebook 的教学格（19 步）。
- [`maintainers/implementation-notes.md`](skills/chromatin-workflow/references/maintainers/implementation-notes.md)
  — 实现说明与溯源（片段库架构、DiffBind/csaw/ChIPseeker API 映射、精简注释与 blacklist 自包含、
  JSON 契约）。

知识原创蒸馏自 [DiffBind](https://bioconductor.org/packages/DiffBind)、
[csaw](https://bioconductor.org/packages/csaw)、[ChIPseeker](https://bioconductor.org/packages/ChIPseeker)
vignette（Artistic-2.0）与 ENCODE ChIP-seq 指南（doi:10.1101/gr.136184.111）；blacklist 来源
doi:10.1038/s41598-019-45839-z。不抄上游文本/图/数据。

## 环境

`environments/chip/environment.yml` 声明依赖（纯 conda + R 源码包，无 pip），`locks/` 是三平台
（darwin-arm64、darwin-x64、linux-x64）的锁定文件。r-base pin 4.5（DiffBind 3.20 / ChIPseeker 1.46 /
csaw 1.44 / DESeq2 1.50 / edgeR 4.10，与 single-cell seurat 路同线）。csaw 是平台条件依赖
（CondaPlatformDependency，排除 darwin-arm64）。GenomeInfoDbData 走 GitHub 镜像 pin（codeload 可达）。
求解走 tuna 镜像，再重写回官方 host。

## 缺口（由 `chip_env_probe` 报事实，走项目环境）

- **csaw（darwin-arm64）**：bioconda 无 osx-arm64 构建；broad mark 在 Apple Silicon 走 DiffBind
  替代或换 linux-x64/darwin-x64/项目环境。
- **IDR**：正式重复性检验需 `idr` pip 包，锁定环境不支持 pip；插件给透明重叠分数代理，不冒充 IDR。
- **profileplyr**：仅 r44，与 r-base 4.5 主线不兼容，不纳入。
- **motif 富集（HOMER/MEME）**：不在 R 下游栈，是缺口；motif 命中是序列预测，不是实测 TF 结合。
- **org.db gene symbol**：需用户 org.db，缺失报 gene id。
- **CUT&Tag 上游 SEACR**：无 wrapper，插件吃用户 SEACR peak；比对/peak calling 是上游缺口。

## 与其他插件的关系

- **single-cell-atac**（姊妹插件）：scATAC 是另一插件（Python/snapatac2，细胞级可及性），
  与本插件按分析范式分开——本插件是 bulk、peak/样本级、R 路。
- **wrappers**：原始 reads 到 BAM/peak（比对、MACS2）是 nf-core/chipseq wrapper 的事，不是 notebook 格。
- **single-cell**：如需 scRNA 对照（peak 注释到表达），参考对象走 SingleCell agent 共享出口。