# Changelog

## 0.1.0 (2026-11-19)

首个版本：宏基因组两条路线的下游脚本链 + 知识层，不做 notebook 渲染器。

- 新插件 `metagenomics`，`toolPrefix: mg`，智能体 `Metagenomics`，私有环境
  `metagenomics`（Python 3.13 + pandas/numpy/scipy/matplotlib-base，三平台
  显式 lock：darwin-arm64 83 包 / darwin-x64 83 包 / linux-x64 67 包）。
- 新技能 `profiling-workflow`（7 个工具）：`mg_inspect`（报告类型识别）、
  `mg_merge`（Kraken2/Bracken/MetaPhlAn/MPA → taxa × sample 矩阵，计数与相对
  丰度不混）、`mg_rarefaction`（计数矩阵稀释曲线：序贯超几何抽样、均值±sd、
  未平台期标记）、`mg_diversity`（α 五指数 + Bray-Curtis + PCoA + 种子置换
  PERMANOVA）、`mg_differential`（Mann-Whitney U/Kruskal-Wallis + BH FDR +
  火山图）、`mg_function`（HUMAnN 表归一/regroup/stratified 拆分）、
  `mg_plot`（组成 barplot/热图，热图可选行/列层次聚类排序）。
- 新技能 `mag-workflow`（6 个工具）：`mg_assembly_stats`（QUAST N50/L50）、
  `mg_bin_qc`（CheckM/CheckM2 → MIMAG 分级 + 散点图，caveat 注明严格 HQ 还需
  rRNA/tRNA）、`mg_bin_filter`（阈值过滤复制 + manifest）、`mg_mag_table`
  （质量 + GTDB-Tk + CoverM + dRep 左连接 + 门水平图；CoverM 低
  covered_fraction 掩码默认 0.3 + 可选 genome_size 加权 per-sample 相对丰度，
  支持单/多样本列）、`mg_mag_function`（prokka/prodigal(.gz)/eggNOG 注释汇总：
  基因计数、MIMAG 结构证据 5S/16S/23S+tRNA≥18（预测级）、bin×KO/bin×COG
  矩阵 + COG 堆叠柱图，--bin_map 覆盖 bin 归属）、`mg_mag_tree`（纯 Python
  Newick 解析 GTDB classify.tree：MAG 间 patristic 距离矩阵 + 最近参考表，
  --summary 精确识别用户 tip，无系统发育库依赖）。
- read-based 链增强：`mg_merge --rank T` 保留 MetaPhlAn/MPA SGB（菌株级）行
  并返回 sgb_species 映射；功能差异链显式接通（mg_function 归一/regroup →
  mg_differential values=relative / mg_diversity / mg_plot）；`mg_inspect`
  增识别 prokka .txt/.tsv、prodigal .faa(.gz)、eggNOG .emapper.annotations、
  Newick 树，ready_steps 指向对应 mag 工具。
- 知识层：两份 analysis-playbook（调用 → 读 → 决定 → 红线）、13 份
  step-notes、两份 maintainers/source-map（四来源对照 + 工具引文 + 许可：
  EasyMetagenome GPL-3.0 只借鉴流程；sib-swiss GPL-3.0 只借鉴思路；
  cambiotraining CC BY 4.0 署名；全部原创 prose）。
- 分工：重计算全走 wrapper（模块 wrapper + nf-core/mag 整流程）；dRep/CoverM
  是 mag 5.5.0 之后的独立 wrapper 步骤；数据库下载前确认磁盘与用户同意。
- 产物附 artifact 描述符：9 个产文件工具为每个图/表写 `.phi-artifact.json`
  （artifact contract 1.0.0，含 provenance/尺寸/行列数），引擎以卡片呈现；
  mag 工具保留 `outputs` 语义键供链式调用，并新增 `artifacts` 列表。
- 配套：nf-core/mag starter wrapper 从 5 参数扩充为覆盖组装/分箱/质量/分类
  主要选项（默认值进 params.json，去除 REVIEW 标记）。
- 来源：EasyMetagenome v1.24（iMeta 2025, DOI 10.1002/imt2.70001）、本地
  Snakemake Workflow Meta 流程、sib-swiss 与 cambiotraining 两个宏基因组
  课程、vendored nf-core/mag 5.5.0 schema。