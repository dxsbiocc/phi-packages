# spatial-transcriptomics

空间转录组分析插件。**Python 单路**：scanpy + squidpy + spatialdata，私有锁定环境
`plugin:spatial`，工具前缀 `st_py_`。覆盖两大测序型平台：

| 平台 | 输入 | 分析单元 |
| --- | --- | --- |
| 10x Visium | Cell Ranger matrix 目录/.h5 + 空间坐标 + 组织图像 | spot（含多个细胞） |
| 10x Visium HD | square bin（8/16µm）或 seg（细胞分割）输出 | bin 或细胞 |
| 华大 Stereo-seq | GEM/GEF（SAW 输出） | square bin（binN）或 cellbin（细胞） |

## 设计原则

与 single-cell 插件同一套契约：

- **插件只做确定性填充/渲染**。意图解析、步骤选择、阈值决策归 agent（规则写在
  [`agents/SpatialTranscriptomics.md`](agents/SpatialTranscriptomics.md) 与技能的 SKILL.md）。
- **片段库是唯一分析源**。`skills/spatial-workflow/assets/fragments/` 的教学代码同时驱动
  脚本工具（薄适配器 + JSON 报告层）与 notebook 渲染器，两条形态共用一份代码。
- **私有锁定环境，不在里面装包**。缺包由 `st_py_env_probe` 当成事实报告，下一步是项目环境。
- **notebook 是真代码格 + 每步 Markdown 教学**，渲染进用户工作目录，不是拷贝模板目录。

## 分析链

`st_py_inspect` → `st_py_qc` → `st_py_preprocess` → 空间图 → `st_py_markers` →
`st_py_annotate` → `st_py_svg`（空间可变基因）→ `st_py_neighborhood`（邻域富集/共现）→
`st_py_deconvolve`（需 scRNA 参考）→ `st_py_communication`（细胞通讯）→
`st_py_spatial_diff`（跨样本差异空间模式）。`st_py_notebook` 是平行的教学路径。

## 知识层

`skills/spatial-workflow/references/` 按功能组织：

- [`analysis-playbook.md`](skills/spatial-workflow/references/analysis-playbook.md) — agent
  的操作手册：每步调用哪个工具、读哪些 JSON 字段与图、怎么定阈值、红线是什么。
- [`step-notes/`](skills/spatial-workflow/references/step-notes/) — 每步的目的/原理/解读/陷阱，
  即 notebook 的教学格。
- [`maintainers/`](skills/spatial-workflow/references/maintainers/) — 实现说明与溯源
  （OSTA 章节对照、Stereo-seq 方法学）。

知识原创蒸馏自 [OSTA](https://bioconductor.org/books/OSTA/)（Bioconductor，CC BY 4.0）与
[sc-best-practices](https://www.sc-best-practices.org/) 空间章节，方法学参照 7 篇 Stereo-seq
高分应用论文；不抄上游文本/图/数据。

## 环境

`environments/spatial/environment.yml` 声明依赖，`locks/` 是三平台（darwin-arm64、
darwin-x64、linux-x64）的锁定文件。求解走 tuna 镜像，再重写回官方 host。

## 与其他插件的关系

- **single-cell**：deconvolution / label transfer 的 scRNA 参考来自 SingleCell agent，
  走共享出口（counts+obs 交换 / pseudobulk 表）。SingleCell 把空间列为缺口，本插件补上。
- **wrappers**：原始 reads 到 counts（Cell Ranger、SAW、比对）是 wrapper 的事，不是 notebook 格。