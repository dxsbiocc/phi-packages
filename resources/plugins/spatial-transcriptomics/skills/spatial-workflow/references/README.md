# references — spatial-workflow

按功能组织，不按来源堆放。两类读者，两个入口。

## Agent — 怎么跑、怎么解读

- [`analysis-playbook.md`](analysis-playbook.md) — 操作手册。每个阶段是一个决策过程：
  **调用**哪个工具、**读**哪些 JSON 字段与图、从读到的数**决定**参数、遵守**红线**。
  阈值起点针对 Visium/Stereo-seq 数据，不是默认值；先读分布再适配。
- [`red-line-cheatsheet.md`](red-line-cheatsheet.md) — §1–§10 硬红线 × 常见误读速查表，
  运行时对照用；完整决策过程仍以 analysis-playbook.md 为准。
- [`step-notes/`](step-notes/) — 每步一篇短笔记（目的/原理/解读/陷阱），按 step id 命名，
  与 `assets/fragments/<step>.py` 配对。这些就是 notebook 渲染器嵌入的教学 Markdown 格。

## Maintainers — 实现与溯源

- [`maintainers/implementation-notes.md`](maintainers/implementation-notes.md) — 每步计算什么、
  片段库/适配器/报告层怎么分工、JSON 契约。
- [`maintainers/osta-map.md`](maintainers/osta-map.md) — OSTA（Bioconductor）章节到本插件
  Python 实现的对照，R 概念如何翻译到 scanpy/squidpy。
- [`maintainers/stereo-seq-methods.md`](maintainers/stereo-seq-methods.md) — Stereo-seq
  数据格式（GEM/GEF、square bin/cellbin）与 7 篇应用论文的方法学要点。

## 来源与许可

知识原创蒸馏自下列开放教材，不复制上游文本/图/数据：

| 来源 | URL | 许可 |
| --- | --- | --- |
| OSTA (Orchestrating Spatial Transcriptomics Analysis with Bioconductor) | https://bioconductor.org/books/OSTA/ | CC BY 4.0（待最终核对） |
| Single-cell best practices（空间章节） | https://www.sc-best-practices.org/ | Apache 2.0；Heumos et al., Nat Rev Genet 2023 |

Stereo-seq 方法学参照用户提供的 7 篇论文（Dev Cell 2022、Zool Res 2024、Cell Discov 2024、
Nat Genet 2024、Cell Rep 2024×2、Cell Metab 2024），见 maintainers/stereo-seq-methods.md。