# Changelog — chromatin-profiling

## 0.1.1 (2026-11-24)

- `chip_diff` 在产出差异结合结果表时返回 `scientificContext`（`valueSemantics: "de_stats"`），供下游绘图与语义校验读取。
- 报告层统一为各步骤产出的图表与结果表自动写入 `1.0.0` 版 `.phi-artifact.json` 旁注并返回 `artifacts` 展示列表（上限 8 项）。

## 0.1.0 (2026-11-20)

首个版本：R 单路 bulk ChIP-seq / CUT&Tag / CUT&RUN 下游分析插件，镜像 single-cell /
spatial / single-cell-atac 的片段库契约。

- 输入覆盖：已 call 的 BED/narrowPeak/broadPeak（MACS2/SEACR）、DiffBind sample sheet CSV、
  peak×sample count 表、indexed BAM、bigWig、已存 `.rds`（DBA/GRanges）。分析单元是 peak/样本，
  不是细胞。**纯下游**：比对/peak calling 不做（ChIP 归 nf-core/chipseq wrapper；CUT&Tag 上游
  SEACR 无 wrapper 是缺口，插件吃用户 SEACR peak）。
- 主引擎 DiffBind 3.20（occupancy+affinity，DBA 对象/sample sheet/自定义 contrast，内部调
  DESeq2/edgeR）+ csaw 1.44（broad mark 窗口法）+ ChIPseeker 1.46（注释）。三个硬约束写进 agent
  与红线：DiffBind 需 sample sheet + indexed BAM；csaw 在 darwin-arm64 是平台缺口；图是 PNG
  （ggplot2/ChIPseeker/ComplexHeatmap，ragg），不是 plotly。
- 片段库唯一分析源：`assets/fragments/`（19 步教学代码）+ `assets/steps.json`（26 参数、
  19 步依赖、6 个 preset），脚本工具（薄适配器 + JSON 报告层）与 notebook 渲染器共用同一份代码。
- 分析链：inspect → QC（**FRiP** 信噪比头条 + 重复性重叠分数 + ENCODE blacklist 过滤）→
  consensus peakset（disjoin/reduce，min_occupancy 定 region universe）→ read counting
  （summarizeOverlaps Union + 归一化）→ 差异结合（DiffBind DESeq2/edgeR 用于 TF/narrow/CUT&Tag；
  csaw 窗口法用于 broad mark H3K27me3/H3K36me3/H3K9me3）→ ChIPseeker 注释（基因组特征分布/
  nearest gene/TSS 距离 + coverage profile/heatmap）→ 样本 PCA + peak 重叠 UpSet。
- 基因组参考（精简策略）：ENCODE blacklist（hg38 636/mm10 3435 区域）base64-gzip 内联在
  `assets/fragments/blacklist-filter.R`，genome 选版本/blacklist 覆盖，版本特异；目标基因组 TxDb
  来自用户 GTF/GFF3（makeTxDbFromGFF）或 TxDb .sqlite（loadDb），不 bundle，缺失报缺口；
  ChIPseeker 附带的 hg19 TxDb 只是构建期依赖，绝不作目标基因组默认注释源；gene symbol 需用户
  org.db，无则报 gene id。绝不在只读环境里装包或下载大参考。
- 知识层：`references/analysis-playbook.md`（调用→读→决定→红线）、
  `references/red-line-cheatsheet.md`（12 条硬红线速查）、`references/step-notes/`（19 步
  目的/原理/解读/陷阱，即 notebook 教学格）、`references/maintainers/implementation-notes.md`
  （架构/DiffBind/csaw/ChIPseeker API 映射/精简注释与 blacklist 自包含/JSON 契约）。原创蒸馏自
  DiffBind/csaw/ChIPseeker vignette（Artistic-2.0）与 ENCODE ChIP-seq 指南
  （doi:10.1101/gr.136184.111），blacklist 来源 doi:10.1038/s41598-019-45839-z；不抄上游文本/图/数据。
- 私有锁定环境 `plugin:chip`（纯 conda + R 源码包，无 pip；r-base pin 4.5，与 seurat 路同线），
  三平台 lock（darwin-arm64/darwin-x64/linux-x64），csaw 是平台条件依赖（排除 darwin-arm64），
  GenomeInfoDbData 走 GitHub 镜像 pin；tuna 镜像求解后重写回官方 host；环境内不装包，缺口由
  `chip_env_probe` 报告（唯一声明 output schema 的工具）。
- notebook 渲染器复用 single-cell 的语言无关机器（RENDERER_VERSION=3，标题/provenance 由
  steps.json 驱动），真代码格 + 每步 Markdown 教学，`self_check` 交付前逐格执行。
- 缺口（报事实，走项目环境）：csaw（darwin-arm64）、IDR（pip，重叠分数是代理不冒充 IDR）、
  profileplyr（仅 r44）、motif 富集 HOMER/MEME（不在 R 下游栈，序列预测非实测结合）、
  org.db gene symbol、CUT&Tag 上游 SEACR（无 wrapper）。
- 运行验证（2026-11-20，darwin-arm64，从 lock 全新重建 scratch env）：presets 1–6 self-check
  全绿。1（DiffBind+DESeq2）与 3（DiffBind+edgeR）在合成 BAM fixture 上给出差异结合峰（方向
  正确、FDR 极低）；2（csaw）在 Apple Silicon 优雅降级后仍跑通 coverage-profile + annotate；FRiP、
  重复性重叠、共识峰、计数矩阵、ChIPseeker 注释均正确。
- 修复：`normalize` 的 `dba.normalize` 归一化方式用 `normalize=DBA_NORM_LIB`（`method=` 是差异引擎
  参数，误填会 "Invalid method"）；`diff-deseq2`/`diff-edger` 的 `dba.report` 显式 `method=` 与
  `dba.analyze` 一致（否则 edgeR 报 "DESeq2 analysis has not been run"）；`setup` 在非 IRkernel
  下打开空设备，避免 `dba.plotMA` 等副作用在工作目录生成 Rplots.pdf。
- 环境锁排序：三平台 lock 重排为 conda-forge → bioconda 常规 → bioconda 数据包（go.db、
  txdb-hg19 等排最后）。否则 post-link 的 `R CMD INSTALL` 依次缺 yq（conda-forge 未先链接）和
  Biobase/S4Vectors（平台子目录包未先链接），sqlite 装不上。此排序在每次 re-solve 后需重做。

- 与 single-cell-atac 按分析范式分开：本插件是 bulk、peak/样本级、R 路；scATAC 是细胞级、
  Python/snapatac2 路。如需 scRNA 对照走 SingleCell agent 共享出口。