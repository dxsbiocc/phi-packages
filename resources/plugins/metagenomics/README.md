# metagenomics（宏基因组分析）· 0.1.0

两条宏基因组路线的**下游分析链 + 知识层**：read-based 物种/功能 profiling
（Kraken2/Bracken、MetaPhlAn/HUMAnN 报告 → 矩阵 → 稀释曲线/多样性/差异/图表，
含 SGB 菌株级视图与功能差异链）与组装/MAG 路线（QUAST 统计 → CheckM/CheckM2
质量分级 → 过滤 → MAG 汇总表（CoverM 可靠性掩码 + 相对丰度）→ prokka/prodigal/
eggNOG 注释汇总 → GTDB 树距离）。
重计算步骤（质控、分类、组装、分箱、CheckM、GTDB-Tk、dRep、CoverM、注释）
全部走现有 wrapper；插件脚本只做确定性下游统计，跑在一个私有锁定 Python 环境里。

Two metagenomics routes as downstream script chains plus a knowledge layer:
read-based profiling (reports → matrix → rarefaction/diversity/differential/
figures, with SGB strain-level views and functional differential chains) and
assembly/MAGs (QUAST stats → quality tiers → filtering → joined MAG table
with CoverM masking/normalization → annotation summaries → GTDB tree
distances).
All heavy steps run through existing wrappers; plugin scripts do
deterministic downstream statistics in one private locked Python environment.

## 组成

| 组件 | 作用 |
|--|--|
| Agent `Metagenomics` | 路线选择（read-based vs 组装/MAG）、wrapper/插件分工、数据库磁盘与下载确认、阈值决策规则 |
| Skill `profiling-workflow` | read-based 下游：7 个 `mg_*` 工具（inspect/merge/rarefaction/diversity/differential/function/plot） |
| Skill `mag-workflow` | 组装/MAG 下游：6 个 `mg_*` 工具（assembly_stats/bin_qc/bin_filter/mag_table/mag_function/mag_tree） |
| 环境 `metagenomics` | Python 3.13 + pandas/numpy/scipy/matplotlib-base，三平台 lock |

两个技能共用 `plugin:metagenomics`。spec 与三平台 lock 在
`environments/metagenomics/`。v0.1.0 只有脚本链，没有 notebook 渲染器。

## 路线（agent 的职责）

- 已有 Kraken2/Bracken/MetaPhlAn/HUMAnN 报告 → read-based 下游链。
- 已有 contigs/bins/CheckM 表，或目标是基因组解析（MAG）→ 组装/MAG 链。
- 只有原始 reads → 先走 wrapper 上游（fastp → kraken2/bracken 或
  metaphlan/humann3；或 megahit/spades → 分箱 → checkm2），插件从报告/表开始。
- 16S/ITS 扩增子 → ampliseq wrapper 流程，两条链都不接。
- 病毒组、泛基因组、菌株追踪 → 知识层指路（genomad/checkv、bowtie2 等
  wrapper），不在本插件实现。

## 工具

read-based（前缀 `mg`）：

| 工具 | 用途 |
|--|--|
| `mg_inspect` | 识别输入类型（Kraken2/Bracken/MetaPhlAn/MPA/HUMAnN/矩阵/samplesheet，也认 QUAST、bin 质量表、prokka/prodigal/eggNOG 输出与 Newick 树），报告样本名与可做的后续步骤 |
| `mg_merge` | 同类型报告合并成 taxa × sample 矩阵（指定 rank；Kraken2/Bracken 计数、MetaPhlAn 相对丰度；rank T = SGB 菌株级，附 sgb_species 映射），返回 `value_kind` |
| `mg_rarefaction` | 计数矩阵稀释曲线（种子超几何抽样、observed richness 均值±sd、分组着色、未平台期样本标记），测序深度充分性检查 |
| `mg_diversity` | α（observed/Chao1/Shannon/Simpson/InvSimpson）、Bray-Curtis、PCoA、PERMANOVA（种子置换），箱线图与排序图 |
| `mg_differential` | 两组 Mann-Whitney U（log2FC + 效应量 r）/多组 Kruskal-Wallis + BH FDR，火山图；矩阵可以是物种或功能群（mg_function 输出，values=relative） |
| `mg_function` | HUMAnN 风格功能表：列归一（relative/CPM）、按 mapping regroup（UniRef90→KO 等）、stratified/unstratified 拆分；归一/regroup 后接 mg_differential/mg_diversity/mg_plot |
| `mg_plot` | 组成 stacked barplot（top N + Other，可组均值）、taxa 热图（可选 log10 与行/列层次聚类排序） |

组装/MAG（同前缀）：

| 工具 | 用途 |
|--|--|
| `mg_assembly_stats` | QUAST report.tsv → 每个组装的 contigs/总长/N50/L50/GC%（取 >=0 bp 行） |
| `mg_bin_qc` | CheckM2/CheckM 质量表 → MIMAG 分级（HQ ≥90/≤5、MQ ≥50/≤10，阈值可传），分级表 + 完整度-污染散点图 |
| `mg_bin_filter` | 按阈值把通过的 bin FASTA 复制成干净集合 + manifest（默认 MQ 底线 50/10） |
| `mg_mag_table` | 质量表左连接 GTDB-Tk summary（拆 d__…s__ 列）、CoverM 丰度、dRep Cdb 簇 → mag-table.tsv + 门水平柱状图；CoverM 低 covered_fraction 掩码（默认 0.3）+ 可选 per-sample 相对丰度（genome_size 加权） |
| `mg_mag_function` | prokka 目录/prodigal .faa(.gz)/eggNOG .emapper.annotations → 每 bin 基因计数、MIMAG 结构证据（5S/16S/23S + tRNA≥18，预测级）、bin×KO 与 bin×COG 矩阵 + COG 堆叠柱图 |
| `mg_mag_tree` | GTDB classify.tree（纯 Python Newick 解析）→ MAG 间 patristic 距离矩阵 + 每个 MAG 最近 GTDB 参考及距离（--summary 精确识别用户 tip） |

## 知识层

每技能一份 `references/analysis-playbook.md`（agent 决策手册：调用 → 读 JSON →
决定 → 红线），每步一份 `references/step-notes/<step>.md`（目的/原理/结果解读/
陷阱），`references/README.md` 按功能索引，`maintainers/source-map.md` 记录
来源对照、许可与引文。全部原创 prose，不复制上游文本。

## 与 wrapper 的分工

- 上游与重计算：fastp、kneaddata、kraken2、bracken、metaphlan、humann3、
  megahit、spades、metabat2、maxbin2、concoct、dastool、checkm/checkm2、
  gtdbtk、quast、prodigal/prokka/eggnogmapper、abricate/rgi 等模块 wrapper，
  以及整条 **nf-core/mag 5.5.0** 工作流 wrapper（本次一并从 5 参数 starter
  扩充）。mag 自带 prodigal/prokka 基因预测；eggNOG 功能注释是 mag 之后的
  独立 wrapper 步骤。
- **dRep 与 CoverM 不在 mag 5.5.0 里**——作为 mag 之后的独立 wrapper 步骤，
  playbook 写明衔接。
- 数据库体积是硬约束：GTDB ~80+ GB、Kraken2 库数十 GB、MetaPhlAn/HUMAnN
  数到十几 GB；任何触发下载的 wrapper 运行前先确认磁盘与用户同意。

## 来源与许可

- **EasyMetagenome v1.24 / EasyMicrobiome**（Liu 实验室，GPL-3.0；iMeta 2025,
  DOI 10.1002/imt2.70001）：流程分段、QC 参数、dRep 阈值惯例。只借鉴流程知识。
- **Snakemake Workflow 仓库**（本地）：工具配对与输入输出形状核对。未复制代码。
- **sib-swiss/metagenomics-data-analysis-training**（GPL-3.0）：决策循环结构与
  解读注意事项。**只借鉴思路，未复制文本或代码。**
- **cambiotraining/metagenomics**（CC BY 4.0）：组装侧实操细节（CheckM 两步
  用法、bowtie2 去宿主、MIMAG 质量讨论）。署名见两个 source-map.md。
- **nf-core/mag 5.5.0**（vendored wrapper）：参数名与默认值对照其
  nextflow_schema.json。

统计实现为原创：Shannon/Simpson/Chao1、Bray-Curtis、PCoA（Gower 中心化）、
PERMANOVA（Anderson 2001）、Mann-Whitney/Kruskal-Wallis（scipy）、BH FDR、
稀释曲线（序贯超几何抽样）、热图聚类（average-linkage，scipy）。
完整引文清单在
[source-map.md](skills/profiling-workflow/references/maintainers/source-map.md)。

## 边界

- 不在受管理环境里装包；环境只有 pandas/numpy/scipy/matplotlib。
- 不做 LEfSe LDA/cladogram（lefse + graphlan wrapper）、Krona（krona
  wrapper）、ANCOM-BC、网络图——缺口如实报告，不假装能跑。
- MIMAG HQ 分级只基于完整度/污染；严格 HQ 还需 rRNA/tRNA 证据——
  `mg_mag_function` 从 prokka 注释给出预测级结构判据（mimag_structural_pass），
  非 barrnap/tRNAscan-SE 确证；工具 caveat 与报告都要写明级别。
- 不做树渲染（环境无系统发育库）；mg_mag_tree 只给距离矩阵与最近参考。
- 不做病毒组/泛基因组分析；16S 归 ampliseq。