# Changelog — spatial-transcriptomics

## 0.2.1 (2026-11-24)

- 在报告输出层统一为所有工具产出的图像与结果表写入 `1.0.0` 版 `.phi-artifact.json` 描述旁注并附带 `artifacts` 展示列表（上限 8 项）。

## 0.2.0 (2026-11-23)

纯 conda-forge 零新依赖的空间分析升级（保持 `plugin:spatial` 锁环境不变）：

- **扩展 `st_py_spatial`（`nhood` 片段）**：在原有邻域富集与共现基础上，补齐 `squidpy.gr.interaction_matrix`（归一化接触矩阵）、`squidpy.gr.centrality_scores`（空间图度/聚类/接近中心性）与 `squidpy.gr.ripley`（Ripley's L/F/G 点模式聚集检验）及对应图形输出。
- **新增 `st_py_niche`（`niche` 片段）**：实现 UTAG 风格的空间邻域聚合（支持 `composition` 邻域细胞组成/去卷积模式与 `expression` 空间平滑表达模式）+ Leiden 空间微环境（Niche）划分，支持在提供 `group_column` 时自动计算跨条件/疾病阶段 Niche 比例动态。
- **扩展 `st_py_communication`（`communication` 片段）**：新增可选 `stratify_by` 参数（如 `niche`），支持在全局配受体通讯之外按空间微环境分层执行 `squidpy.gr.ligrec`。
- **新增 `st_py_trajectory`（`trajectory` 片段）**：使用 `scanpy` 内置 PAGA + Diffusion Pseudotime（DPT）将表达流形拟时序映射回组织空间坐标，并支持通过 `gene_sets_path` 在未缩放 `lognorm` 层运行 `sc.tl.score_genes` 空间基因程序打分。
- **Notebook 与教学注释同步**：`assets/steps.json`、Module 4/6 预设、`references/step-notes/` 及端到端运行时测试全部同步更新。

## 0.1.0 (unreleased)

首个版本：Python 单路空间转录组分析插件，镜像 single-cell 的片段库契约。

- 平台覆盖：10x Visium（含 Visium HD 的 square bin 与细胞分割 seg）与华大 Stereo-seq
  （GEM/GEF 读取、square binning 自实现、cellbin 读 SAW 输出）。
- 片段库唯一分析源：`assets/fragments/` 教学代码 + `assets/steps.json` presets/参数默认，
  脚本工具（薄适配器 + JSON 报告层）与 notebook 渲染器共用同一份代码。
- 分析链：inspect/load/bin → QC → normalize/HVG/PCA/cluster/UMAP → 组织坐标空间图 →
  marker/注释 → 空间可变基因（SVG，Moran's I/Geary's C）→ 邻域富集/共现 →
  deconvolution（scRNA 参考）→ 细胞通讯 → 跨样本差异空间模式/共定位。
- 知识层：`references/analysis-playbook.md`（调用→读→决定→红线）与
  `references/step-notes/`（每步目的/原理/解读/陷阱，即 notebook 教学格），
  原创蒸馏自 OSTA（Bioconductor，CC BY 4.0）与 sc-best-practices 空间章节，
  方法学参照 7 篇 Stereo-seq 高分应用论文；不抄上游文本/图/数据。
- 私有锁定环境 `plugin:spatial`（scanpy + squidpy + spatialdata + 全量能力依赖），
  tuna 镜像求解后重写回官方 host；环境内不装包，缺口由 `st_py_env_probe` 报告。
- notebook 渲染器复用 single-cell 的语言无关机器（标题/provenance 改由 steps.json 驱动），
  真代码格 + 每步 Markdown 教学，`self_check` 交付前逐格执行。
- 与 single-cell 插件交接：deconvolution/label transfer 的 scRNA 参考走共享出口
  （counts+obs 交换 / pseudobulk 表）。