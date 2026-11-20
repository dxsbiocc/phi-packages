# Changelog — spatial-transcriptomics

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