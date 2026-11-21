# references — scatac-workflow

按功能组织，不按来源堆放。两类读者，两个入口。

## Agent — 怎么跑、怎么解读

- [`analysis-playbook.md`](analysis-playbook.md) — 操作手册。每个阶段是一个决策过程：
  **调用**哪个工具、**读**哪些 JSON 字段与图、从读到的数**决定**参数、遵守**红线**。
  阈值起点针对 scATAC 数据，不是默认值；先读分布再适配。引擎是 snapatac2 2.10
  （file-backed，Rust），scanpy 1.12 只用于 gene-activity 下游与 matplotlib 面板。
- [`red-line-cheatsheet.md`](red-line-cheatsheet.md) — §1–§10 硬红线 × 常见误读速查表，
  运行时对照用；完整决策过程仍以 analysis-playbook.md 为准。
- [`step-notes/`](step-notes/) — 每步一篇短笔记（目的/原理/解读/陷阱），按 step id 命名，
  与 `assets/fragments/<step>.py` 配对。这些就是 notebook 渲染器嵌入的教学 Markdown 格。

## Maintainers — 实现与溯源

- [`maintainers/implementation-notes.md`](maintainers/implementation-notes.md) — 每步计算什么、
  片段库/适配器/报告层怎么分工、snapatac2 API 映射、backed 对象约束、JSON 契约。

## 来源与许可

知识原创蒸馏自下列开放来源，不复制上游文本/图/数据：

| 来源 | URL | 许可 |
| --- | --- | --- |
| SnapATAC2 文档与教程（Zhang et al.） | https://github.com/kaizhang/SnapATAC2 | MIT |
| Single-cell best practices（染色质可及性章节） | https://www.sc-best-practices.org/ | Apache 2.0；Heumos et al., Nat Rev Genet 2023 |
| ENCODE Blacklist（Amemiya, Kundaje, Boyle） | doi:10.1038/s41598-019-45839-z | Sci Rep 9, 9354 (2019)；Boyle-Lab/Blacklist ENCODE v2 |

ENCODE hg38/mm10 blacklist 以 base64-gzip 内联在 `assets/fragments/blacklist.py`，
脚本工具与渲染出的 notebook 都能离线过滤 artifact 区域。snapatac2 离线内置
hg38/mm10/hg19/mm39 的 chrom_sizes；基因注释（GFF3）与基因组 FASTA 是懒下载，
由用户经 `input-path` 提供或一次性下载到项目目录，绝不装进只读插件环境。
代码为本插件原创，运行在私有 `plugin:scatac` 环境。