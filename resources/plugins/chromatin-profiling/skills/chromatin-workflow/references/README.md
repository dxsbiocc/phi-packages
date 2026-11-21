# references — chromatin-workflow

按功能组织，不按来源堆放。两类读者，两个入口。

## Agent — 怎么跑、怎么解读

- [`analysis-playbook.md`](analysis-playbook.md) — 操作手册。每个阶段是一个决策过程：
  **调用**哪个工具、**读**哪些 JSON 字段与图、从读到的数**决定**参数、遵守**红线**。
  阈值起点针对 ChIP/CUT&Tag 数据，不是默认值；先读分布再适配。引擎是 R/Bioconductor：
  DiffBind 3.20（occupancy+affinity）、DESeq2/edgeR（差异）、csaw 1.44（broad mark 窗口法）、
  ChIPseeker 1.46（注释）、ComplexHeatmap（UpSet）、ggplot2（PCA/MA）。**纯下游**：
  从 called peaks/count 表/indexed BAM/bigWig 出发，绝不碰 FASTQ。
- [`red-line-cheatsheet.md`](red-line-cheatsheet.md) — §1–§9 硬红线 × 常见误读速查表，
  运行时对照用；完整决策过程仍以 analysis-playbook.md 为准。
- [`step-notes/`](step-notes/) — 每步一篇短笔记（目的/原理/解读/陷阱），按 step id 命名，
  与 `assets/fragments/<step>.R` 配对。这些就是 notebook 渲染器嵌入的教学 Markdown 格。

## Maintainers — 实现与溯源

- [`maintainers/implementation-notes.md`](maintainers/implementation-notes.md) — 每步计算什么、
  片段库/适配器/报告层怎么分工、DiffBind/csaw/ChIPseeker API 映射、精简注释与 blacklist
  自包含约束、JSON 契约。

## 来源与许可

知识原创蒸馏自下列开放来源，不复制上游文本/图/数据：

| 来源 | URL | 许可 |
| --- | --- | --- |
| DiffBind vignette（Stark & Brown） | https://bioconductor.org/packages/DiffBind | Artistic-2.0 |
| csaw vignette（Lun & Smyth） | https://bioconductor.org/packages/csaw | Artistic-2.0 |
| ChIPseeker vignette（Yu et al.） | https://bioconductor.org/packages/ChIPseeker | Artistic-2.0 |
| ENCODE ChIP-seq 指南（Landt et al.） | doi:10.1101/gr.136184.111 | Genome Res 22:1813 (2012) |
| ENCODE Blacklist（Amemiya, Kundaje, Boyle） | doi:10.1038/s41598-019-45839-z | Sci Rep 9:9354 (2019)；Boyle-Lab/Blacklist ENCODE v2 |
| nf-core/chipseq（上游 wrapper） | https://github.com/nf-core/chipseq | MIT |

ENCODE hg38/mm10 blacklist 以 base64-gzip 内联在 `assets/fragments/blacklist-filter.R`，
脚本工具与渲染出的 notebook 都能离线过滤 artifact 区域。目标基因组注释（GTF/GFF3 或 TxDb）
与 org.db 由用户经 `input-path` 提供或走项目环境，绝不为目标基因组捆绑、绝不装进只读插件环境。
代码为本插件原创，运行在私有 `plugin:chip` 环境。