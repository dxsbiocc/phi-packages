# ChIP/CUT&Tag 硬红线 × 常见误读速查表

运行时对照用。完整决策过程（调用 → 读 → 决定 → 红线）以 [`analysis-playbook.md`](analysis-playbook.md) 为准；
本表把每条红线压成"红线 / 为什么 / 常见误读 / 正确做法"四列，便于在每一步快速自检。

| # | 硬红线 | 为什么 | 常见误读 | 正确做法 |
| --- | --- | --- | --- | --- |
| 1 | **纯下游**：从 called peaks / counts / BAM / bigWig 出发，不做比对与 peak calling | FASTQ→BAM→peak 是上游；ChIP 归 nf-core/chipseq wrapper，CUT&Tag 归用户 SEACR run | 在插件里跑 bowtie2/MACS2/SEACR，或期望输入 FASTQ | 输入是 BED/narrowPeak/broadPeak、sample sheet、count 表、indexed BAM；上游缺失指引 wrapper 或用户 SEACR（CUT&Tag 无 wrapper 是缺口） |
| 2 | **只读环境**：缺包/缺参考报事实 + 项目环境，不装包不下载大参考 | 插件环境预建锁定，装包破坏可复现性与签名 | 在插件里 BiocManager::install / pip 装 idr / 下载全基因组 GTF 进插件环境 | `chip_env_probe` 报缺口 → env_request 项目环境，或用户经 input-path 提供 GTF/TxDb/org.db/blacklist |
| 3 | **精简注释**：目标基因组 TxDb 来自用户；附带 hg19 TxDb 绝不默认注释源 | ChIPseeker conda 配方硬依赖 hg19 TxDb（构建期），它不是目标基因组 | 拿附带的 hg19 TxDb 注释 hg38/mm10 peak；无 org.db 却声称给了基因符号 | 用户 `gtf`（makeTxDbFromGFF）或 `txdb`（loadDb）；无则报缺口；无 org.db 报 gene id 不报符号 |
| 4 | **blacklist 基因组版本特异**：hg38 区域在 mm10 上是错的 | 传错版本静默删错 locus，污染 consensus 与差异分析 | 不传 genome 就默认 hg38 套到小鼠；或干脆不过 blacklist | 传 `genome`（hg38/mm10）用内联 ENCODE v2，或 `blacklist` 用其他版本；两者都无则跳过并说明 |
| 5 | **DiffBind vs csaw 是生物学决策** | TF/narrow mark 信号尖峰适合 peak-based；broad mark 弥散 over kilobases，peak-based 会漏 | 用 DiffBind 一套到底跑 H3K27me3；或在 Apple Silicon 硬等 csaw | TF/narrow（H3K4me3/H3K27ac）与 CUT&Tag→DiffBind（DESeq2/edgeR）；broad（H3K27me3/H3K36me3/H3K9me3）→csaw；csaw 在 darwin-arm64 是缺口，DiffBind 替代或换 linux-x64/darwin-x64 |
| 6 | **FRiP 是信噪比头条，阈值来自分布；IDR 是缺口** | 低 FRiP = 抗体没富集，多数 reads 是背景；IDR 需 idr pip 包，环境未带 | FRiP 阈值凭默认；把 peak 重叠分数当 IDR 正式结论 | 读 `chip_qc` 的 frip 分布再判断（TF ChIP 理想 >1%，broad 更低，CUT&Tag 往往很高）；重叠分数是透明代理，说清不是 IDR |
| 7 | **重复是实验单元，不是 peak** | 条件级声称的重复是生物学样本；peak 是观测维度不是重复 | 把成千上万 peak 当重复做条件差异；每组 1 个样本就声称条件差异 | 每组 ≥ `min_members`（默认 2）生物学重复才做条件级声称；否则只描述样本内模式 |
| 8 | **归一化/探索 count 不进差异检验** | DESeq2/edgeR/csaw 各自内部归一化（median-of-ratios/TMM），外部 CPM 破坏均值-方差假设 | 拿 `chip_count` 的 normalized_counts 喂差异检验 | normalized_counts 只用于 PCA/可视化；差异回到原始 count，引擎内部归一化 |
| 9 | **consensus / min_occupancy 是 region universe 决策** | 一个位点要在几个样本出现才进 consensus，决定检验的区间全集与多重检验负担 | min_occupancy 凭默认一套到底；不过 blacklist 就建 consensus | 先 blacklist 过滤再 consensus；按重复数与生物学问题定 min_occupancy；说清选择与理由 |
| 10 | **注释/motif 是推断/预测，不是实测调控** | ChIPseeker 按位置注释到最近基因；motif 富集是序列层面过代表检验 | 把 peak 注释当实测靶基因调控；把 motif 命中当实测 TF 结合 | 注释是位置推断需结合表达/功能验证；motif 富集（HOMER/MEME）是缺口、不在 R 下游栈，要它报缺口走项目环境 |
| 11 | **DiffBind 需 sample sheet + indexed BAM** | dba(sampleSheet=) 从 sheet 建 DBA；counting/FRiP 需 .bai 索引 | 无 sheet 硬跑 diff；BAM 无 index 就 counting | diff/count/qc 从 DiffBind sample sheet 出发（Peaks/bamReads 列驱动加载，相对路径相对 sheet 目录）；BAM 需 .bai（samtools index 是上游） |
| 12 | **JSON 是审计记录；"没跑"不说成"没问题"** | 每步返回前后计数与缺口是交付凭证 | 自己重算工具已给的数；把 csaw 缺口/注释缺口说成已完成 | 引用工具 JSON 的 peak 数/FRiP/occupancy/差异数/特征分布；明确哪步因缺包或缺参考没执行 |

## 使用方式

- 每跑一步前，扫一遍与本步相关的行（QC→1/6/11；consensus→4/9；count→8/11；diff→5/7/8/11；annotate→3/10；任何重依赖→2/5）。
- 交付前对照第 2 行（缺包/缺参考有没有被装进插件或说成没问题）与第 12 行（"没跑"有没有被说成"没问题"）。
- 阈值与策略永远"先读分布/报告，再决定"，不是先验默认。引用工具 JSON 里的数字，不要自己重算。