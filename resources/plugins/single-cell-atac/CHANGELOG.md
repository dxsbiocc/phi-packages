# Changelog — single-cell-atac

## 0.1.3 (2026-11-25)

- 为 `scatac-workflow` 全部脚本工具与 `scatac_notebook` 补齐 `output` JSON Schema 声明（`step-result.json`、`notebook-result.json`）。

## 0.1.2 (2026-11-25)

- 将 `scripts/_common.py` 中对 `numpy` 的导入改为按需容错加载，使 `scatac_env_probe` 与 `scatac_notebook` 在纯标准库 Python 环境下也可独立运行自检。

## 0.1.1 (2026-11-24)

- 在报告输出层自动为产出的图表（含交互式 Plotly `.html` 报告及 `.png`/`.svg`）与结果表生成 `1.0.0` 版 `.phi-artifact.json` 旁注并返回 `artifacts` 列表（上限 8 项）。

## 0.1.0 (2026-11-20)

首个版本：Python 单路单细胞 ATAC 分析插件，镜像 single-cell / spatial 的片段库契约。

- 输入覆盖：10x Cell Ranger ATAC `fragments.tsv.gz`（首选）、BAM（先转 fragment）、
  已有 `.h5ad`、peak×cell 矩阵（csv/tsv/h5/10x）。分析单元是细胞（barcode），不是 spot/bin。
  原始 reads 比对 / Cell Ranger ATAC 不做（nf-core wrapper）。
- 主引擎 snapatac2 2.10（Rust 后端、file-backed AnnData、plotly 图），scanpy 1.12 只用于
  gene activity 下游。两个硬约束写进 agent 与红线：backed 对象退出前必须 `close()`；
  `snap.pl.*` 是 plotly Figure（存自包含 `.html`），不是 matplotlib Axes。
- 片段库唯一分析源：`assets/fragments/`（16 步教学代码）+ `assets/steps.json`（48 参数、
  16 步依赖、6 个 preset），脚本工具（薄适配器 + JSON 报告层）与 notebook 渲染器共用同一份代码。
- 分析链：inspect → QC（n_fragment/frac_dup/frac_mito/**TSSe** + 阈值分布）→ preprocess
  （tile 矩阵/特征选择带 ENCODE blacklist/scrublet 双细胞/spectral 降维/kNN/Leiden/UMAP）→
  per-cluster MACS3 peak calling + 合并 peak 集 + peak-count 矩阵 → marker_regions
  （cluster-vs-rest，可选 pairwise `diff_test` 差异可及性，`run_diff` 因 O(clusters²) 而 opt-in）→
  gene activity 注释（panel 打分 / scRNA 参考 kNN label transfer / 显式 mapping）→
  motif 富集（TF 调控 motif，序列预测）→ 多样本整合（harmony/scanorama/mnc）。
- 基因组参考：chrom_sizes 由 snapatac2 自带（hg38/mm10/hg19/mm39，离线）；ENCODE blacklist
  （hg38/mm10）base64-gzip 内联在 `assets/fragments/blacklist.py`；基因注释（GFF3 ~50MB）与
  FASTA（~845MB）是懒下载，由用户 input-path 或下载到项目目录，缺失报缺口，不在只读环境里下载。
- 知识层：`references/analysis-playbook.md`（调用→读→决定→红线）、
  `references/red-line-cheatsheet.md`（12 条硬红线速查）、`references/step-notes/`（每步
  目的/原理/解读/陷阱，即 notebook 教学格）、`references/maintainers/implementation-notes.md`
  （架构/backed 约束/API 映射/JSON 契约）。原创蒸馏自 SnapATAC2（MIT）与 sc-best-practices
  （Apache 2.0），blacklist 来源 doi:10.1038/s41598-019-45839-z；不抄上游文本/图/数据。
- 私有锁定环境 `plugin:scatac`（纯 conda，无 pip；macs3 作为 snapatac2 依赖随环境），
  三平台 lock（darwin-arm64/darwin-x64/linux-x64），tuna 镜像求解后重写回官方 host；
  环境内不装包，缺口由 `scatac_env_probe` 报告（唯一声明 output schema 的工具）。
- notebook 渲染器复用 single-cell 的语言无关机器（RENDERER_VERSION=3，标题/provenance 由
  steps.json 驱动），真代码格 + 每步 Markdown 教学，`self_check` 交付前逐格执行。
- 缺口（报事实，走项目环境）：chromVAR/ArchR/Signac（R 路）、episcanpy（numpy 2 不兼容）、
  scvi-tools/magic-impute（pip-only）、TF-gene GRN（未实现）。TF 层是 motif 富集，
  不是 chromVAR 逐细胞偏差；gene activity 是代理，不是实测 RNA。
- 与 single-cell 插件交接：label transfer / gene-activity 注释对照的 scRNA 参考走共享出口
  （counts+obs 交换 / pseudobulk 表）。bulk ChIP/CUT&Tag 归规划中的姊妹插件 chromatin-profiling。