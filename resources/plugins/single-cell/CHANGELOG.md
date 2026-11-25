# Changelog

## 0.7.1 (2026-11-25)

- 为 `seurat-workflow` 与 `scanpy-workflow` 全部脚本工具与 Notebook 渲染工具补齐 `output` JSON Schema 声明（`step-result.json`、`notebook-result.json`）。

## 0.7.0 (2026-11-25)

- Python 路（`scanpy-workflow`）新增配受体细胞通讯推断能力：新增片段 `assets/fragments/communication.py`、教学笔记 `references/step-notes/communication.md`、工具 `sc_py_communication` 及 Notebook Preset 8（基于未缩放 `lognorm` 层实现 CellPhoneDB 多亚基最小表达法则 + 细胞标签置换检验 + BH FDR，内置精选人类/小鼠跨物种配受体复合物起始库并支持用户传入自定义 `interactions_path`，自动产出互作计数热图、气泡图及结果表）。

## 0.6.1 (2026-11-24)

- R 路 `sc_pseudobulk` 与 Python 路 `sc_py_pseudobulk` 在导出样本级聚合计数矩阵时同步返回 `scientificContext`（`valueSemantics: "raw_counts"`、`sampleIds` 与 `donorBySample` 映射），与下游差异表达 Wrapper 的矩阵语义契约对齐。
- R 路与 Python 路报告层统一为所有产出的图表与表格自动生成 `1.0.0` 版 `.phi-artifact.json` 旁注并返回 `artifacts` 展示列表（上限 8 项）。

## 0.6.0 (2026-11-19)

每路只有**一份分析代码**，知识按功能归位：片段库既是 notebook 的细胞，也是
脚本工具的执行体；目的/原理/结果解读/陷阱成为技能内容，agent 据此选阈值、读结果。

- 两路统一为「片段库 = 唯一分析源 + 薄适配器 + 报告层」：`assets/fragments/`
  是教学代码本体；`scripts/_fragments.{R,py}` 在同一全局环境执行同一批片段；
  工具脚本只把 CLI 旗标映射成与 config cell 同名的配置全局，`scripts/reports/*`
  只做 JSON 序列化。片段的教学输出走 stderr（R `sink(stderr())`、Python
  `redirect_stdout`），stdout 只有一份 JSON。重构前后逐字段核对：R 7/7、
  Python 10/10 份报告与旧契约一致（仅输出路径不同）。
- 知识层：30 份 `references/step-notes/<step>.md`（每步一份，目的/原理/结果解读/陷阱，
  即 notebook 讲解格，与片段按步骤 id 配对）；两路各一份 `references/analysis-playbook.md`——agent 决策手册，按阶段
  「调用 → 读（确切 JSON 字段与图）→ 决定（参数起点+调整规则）→ 红线」，原理压缩为
  一行，来源与许可在文末（OSCA CC BY 4.0；sc-best-practices Apache 2.0，Heumos/Schaar/
  Lance et al., Nat Rev Genet 2023, DOI 10.1038/s41576-023-00586-w；singlecellcourse.org
  是 Cambridge/Hemberg 课程 GPL-3，不是 OSCA）。`references/README.md` 按功能索引
  （agent 决策 / 每步解读 / 跨路交接 vs `maintainers/` 溯源）。全部原创 prose，不复制上游文本；两个 SKILL.md 增加指引。
- 目录按功能重排：assets 只放渲染源材料（`fragments/` 代码去掉冗余语言层平铺、
  `steps.json` 上移到 assets 顶层）；Markdown 全部归 references 分组（`step-notes/`、
  `maintainers/`、playbook、跨路交接）。命名统一 kebab-case 且揭示功能：
  `r_interop.md`→`cross-route-handoff.md`、`pipeline-notes.md`→
  `maintainers/implementation-notes.md`。渲染器与加载器改读新路径，两副本仍字节相同；
  新增布局测试（steps.json 每步必须有代码片段与 step-note）。
- `scnotebooks-curriculum` 独立技能折进 `seurat-workflow/references/maintainers/`
  （curriculum + module-map 两份对照文档），插件内技能 3→2、仓库插件内技能 7→6、
  技能总数 23→22；agent、README、SOURCE.json 与测试同步。
- 许可更正：NBIS workshop-scRNAseq 的许可以仓库 LICENSE 为准 **CC-BY-4.0**（站点
  页脚标 GPL-3，与仓库不符）。steps.json attribution 与 7 个预设 intro、
  nbis-porting.md、SOURCE.json、设计文档、README 统一更正。
- 渲染器补 bool/list 参数类型（`flag_bool`/`flag_list`，`resolution` 默认 `[0.8]`），
  两技能渲染器仍字节相同（仓库测试断言）。
- 修 `show_plot` 在 Rscript 下留 `Rplots.pdf`：内联 `print(plot)` 只在 IRkernel
  命名空间已加载时执行，脚本/self-check 不再在用户工作目录落默认 PDF 设备产物；
  `adt-assay`/`clonotypes` 的裸图表达式改走 `show_plot`（此前 ADT 小提琴图在
  notebook 里根本不显示）。
- 验证：两路 JSON parity 全一致；R 预设 m4/m5/m8/m9 与 Python 预设 p1–p7
  self-check 全过；仓库 35 项测试通过；Phi 校验器 `check:plugins` 3 插件 0 错、
  `check:skills` 22 技能 0 错。

## 0.5.0 (2026-11-19)

插件内维护**两套平行工具链**：R/Seurat 6 与 Python/scanpy 1.13，各带一个私有锁定
环境；R 还是 Python 是**用户的路**，由 agent 按倾向选，渲染器与脚本只做确定性执行。

- 新增 `environments/scanpy/`：`environment.yml`（Python 3.13 + scanpy 1.13 +
  ipykernel）与三平台显式 lock。scanpy 1.13 要求 python >=3.13；harmonypy、scrublet
  走 bioconda；leidenalg + python-igraph（Leiden 用 `flavor="igraph"`）；scikit-misc
  供 seurat_v3 高变基因。bbknn 在 bioconda 无 py313 build，作为缺口由
  `sc_py_env_probe` 报告，不在环境里装包。
- 新增技能 `scanpy-workflow`（10 个 `sc_py_*` 工具）：inspect / qc（Scrublet）/
  preprocess / integrate（Harmony、ComBat）/ markers / annotate / pseudobulk /
  trajectory（PAGA + 扩散伪时间）/ env_probe / notebook。脚本链与 R 路同构：
  每步读一个 `.h5ad` 写一个 `.h5ad`，参数记在返回 JSON。
- 片段库 `assets/fragments/python/` + `steps.json`：预设 1–7 对应 NBIS
  workshop-scRNAseq scanpy labs（QC、降维、整合、聚类+marker、样本级差异、注释、
  轨迹）。上游 GPL-3.0，只作参考：片段代码原创，intro 与 SOURCE.json 注明出处。
- 渲染器两技能**字节相同**（仓库测试断言）：`--language python` 走同一套
  steps/config/inject/merge/self-check 逻辑；kernel 名 `plugin-single-cell-scanpy`。
- scanpy 1.13 / anndata 0.13 API 迁移：accessor 字符串（`obs.<col>`、`obsm.X_pca`、
  `layers.counts`、`use=`/`out=`/`rng=`）；`.raw` 已弃用，改为归一化后立刻留
  `layers["lognorm"]` 快照（scale 前、全基因），markers/annotate/ComBat 读它；
  `sc.io.*` 读入、显式 `diffmap`/`dendrogram`、`importlib.metadata` 取版本。
  迁移前后结果核对一致（注释打分、聚类成员；igraph flavor 仅交换簇编号）。
- agent `SingleCell` 重写：**不绑定环境**（技能自带环境优先于 agent 绑定，见契约），
  只构建所选路的环境；选路规则（显式倾向优先、`.rds`→R、`.h5ad`→Python、仅 10x
  先问）、不混路、跨路只走 pseudobulk 或显式互转；两路能力缺口表入 agent 文档。
  两个 SKILL.md 的 Rules 首条镜像同一规则。
- 吸收独立技能 `resources/skills/scanpy`（K-Dense BSD-3）进插件后删除：同一能力
  只有一个家。独立技能 17→16、插件内技能 6→7，README 与 SOURCE.json 同步。
- 修正 self-check 弹原生图窗：子进程跑在 Jupyter 之外，macOS 默认 macosx backend
  加 scanpy `autoshow=True` 会每张图开一个阻塞窗口。现在 self-check 用
  `MPLBACKEND=Agg` + `autoshow=False` 前缀；notebook 内联展示不受影响。
- 修正 `batch-compare` 片段：整合后 UMAP 此前未计算就画 after 图（KeyError）；
  现在片段内先算整合后 neighbors+UMAP 再对比。
- 验证：Python 预设 1–7 self-check 全过（8/13/15/16/17/18/19 格）；R 路回归
  m4/m5/m8/m9 全过（17/18/18/18 格）；Python 脚本链 inspect→qc→preprocess→
  integrate（Harmony 与 ComBat）→markers→annotate→pseudobulk→trajectory 全过；
  仓库 35 项测试通过；Phi 校验器 `check:plugins` 3 插件 0 错、`check:skills`
  23 技能 0 错。

## 0.4.0 (2026-11-19)

插件改为维护**私有环境**，不再借用共享的 `phi:r@1`；Seurat 升到最新的 **6.0.0**。
脚本链与 notebook 共用同一套锁定版本。

- 新增 `environments/seurat/`：`environment.yml`（R 4.5 + Seurat 6.0.0 + IRkernel）
  与三平台显式 lock（darwin-arm64 / darwin-x64 / linux-x64）。Seurat 6.0.0、
  RcppThread 2.4.0、GenomeInfoDbData 无可用 conda 构建，按 `sourcePackages` 源码
  pin（CRAN 精确版本 / GitHub 40 位 sha + sha256）；其余依赖走 conda，原生编译用
  `compilers`。r-base 锁 4.5（SeuratObject 5.4.0 只有 r44/r45 构建）。
- `phi-package.yaml` 声明 `environments: { seurat: ... }`；技能与 agent 的环境引用
  由 `phi:r@1` 改为 `plugin:seurat`。
- notebook kernel 由 `phi-r` 改为 `plugin-single-cell-seurat`（渲染器默认与
  `steps.json` 同步）。Phi 主程序扩展 managed-kernel 同步：已启用插件中自带
  `r-irkernel`/`ipykernel` 的私有环境会同步出 `plugin-<pluginId>-<envName>` 内核，
  环境就绪后可选，禁用/卸载时清理（Phi `managed-kernels.ts`，16 项测试通过）。
- 文案与对照表（SKILL.md、README、module-map、pipeline-notes、env_probe 等）由
  phi-r 改述为插件环境；修正 harmony 实为缺包（旧文档误称经 Seurat imports 提供）。
- 修正残留的旧工具名 `seurat_*` → `sc_*`（toolPrefix 为 `sc`；`seurat_clusters`
  等 Seurat 列名不动）。
- 验证：本地按 lock 构建 darwin-arm64 环境并源码安装三个包，`packageVersion("Seurat")`
  == 6.0.0；模块 4/5/8/9 与 steps+inject+merge+config 组合渲染 self-check 全通过，
  Python 渲染仍用 `phi-python`；`env_probe` 报告 Seurat 6.0.0、core/blocked 模块与
  文档一致；仓库 34 项测试通过。

## 0.3.0 (2026-11-19)

`sc_notebook` 从固定模块生成器重构为**片段库 + 渲染器**：插件提供模板代码，按
agent 传入的步骤与参数渲染成一份填好占位的分析 notebook，落在用户工作目录。
意图解析、步骤选择、阈值判断是 agent 的职责，规则写在 SKILL.md；渲染器只做
确定性的填充与组装，不从数据推断任何决策。

- 19 段分析代码从 `build_notebook.py` 的内嵌字符串抽成 `assets/fragments/r/`
  片段库（代码 + 讲解各一份），`steps.json` 声明步骤依赖与参数默认值。
- 按传入步骤组装：`steps` 数组列出步骤，声明的依赖自动补全（如 `cluster,markers`
  自动带上质控到 PCA 的前置链）；模块 4/5/8/9 保留为预设。
- 单一参数源：所有阈值、路径、列名渲染进首格 config cell，分析格引用变量；
  每个值标注来源（default / config / flag）。支持 `config` JSON 覆盖与 `mt_pattern`
  （鼠源 `^mt-`）、`seed` 参数。
- `inject` 把外部 R/Python 脚本原样注入为代码格；`merge` 合并其他 notebook 的
  全部或部分 cell（`path.ipynb[0:6]`）。一个 notebook 一个 kernel，语言不一致
  明确报错。`language=python` 支持纯注入/合并渲染。
- 可移植性：数据与输出路径改为相对 notebook 目录；`self_check` 在输出目录真跑
  每格（原来在临时目录、绝对路径）。
- 渲染溯源写进 notebook metadata：步骤、参数及来源、片段 sha256、注入/合并来源、
  数据指纹、渲染时间。
- 验证：合成 400 细胞 × 300 基因 10x 数据，模块 4/5/8/9 与 steps+inject+merge+config
  组合渲染全部 self-check 通过；语言冲突、未知步骤、未知 config 键均明确报错。
- 分工规则入 SKILL.md：「渲染器填充，agent 决策」「同一次质控读数传给两条路」，
  与脚本链「阈值是参数不是猜测」的既有规则一致。
- 插件 id 由 `single-cell-training` 改为 `single-cell`：定位早已不止教学移植，
  与标题「单细胞分析」、智能体 `SingleCell`、前缀 `sc` 及兄弟插件的名词风格一致。
  课程对照（`scnotebooks-curriculum`）仍是组件之一。

## 0.2.0 (2026-11-16)

补上 notebook 路径：从「只有脚本」变成「notebook 探索 + 脚本复现」两条路。

- 新工具 `sc_notebook`：在用户项目里生成可运行的 R notebook（模块 4/5/8/9），
  每格是真实 Seurat 代码并配讲解，用 Phi 的 notebook 工具逐格执行、按需修改。
  模块 4 覆盖质控到 pseudobulk 全流程；5 增加批次对比；8 增加蛋白 assay；9 增加
  克隆型连接。
- `self_check` 参数会把每格真的跑一遍，任一格失败就报错——保证交付给学员的
  notebook 是能跑的，而不是看起来能跑。
- 生成时修正三个实测踩到的坑：`Read10X` 需要 `.gz`；`min.features` 会静默删光
  细胞后抛 `LogMap` 报错（现在先查分布再明确拒绝）；面板基因全缺时 `DotPlot`
  崩溃（现在只画存在的基因，并在一个可用基因都没有时停下来说明标识符类型不对）。
- 原始测序数据处理不做：比对和 Cell Ranger 属于 wrapper 树，notebook 从计数矩阵
  或 `.rds` 开始。

## 0.1.0 (2026-11-16)

首个版本，把 scNotebooks 单细胞课程里能执行的分析步骤移植成 Phi 插件。

- 新插件 `single-cell-training`，`toolPrefix: sc`，智能体 `SingleCell`。
- 新技能 `seurat-workflow`（7 个工具）：`sc_inspect`、`sc_qc`、`sc_preprocess`、
  `sc_markers`、`sc_annotate`、`sc_pseudobulk`、`sc_env_probe`。全部跑在 `phi:r@1`。
- 新技能 `scnotebooks-curriculum`：13 个模块的对照表，只有说明文档。
- `sc_env_probe` 逐个模块报告缺哪个包、哪些部分跑不了：`phi-r` 有 Seurat、
  SingleCellExperiment、ComplexHeatmap，没有 monocle3、liana、semla、ArchR、
  SoupX、scDblFinder 等。跑不了的模块（06、07、10、11）在文档和工具返回里都标出来，
  不在环境里装包。
- 不复制 notebook 原文、图片或示例数据；上游是 CC BY 4.0
  （Rojas-Hidalgo et al., *Nat Genet*, 2026），在 README 与 SOURCE.json 中注明。
