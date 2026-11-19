# single-cell（单细胞分析）· 0.6.0

两套平行的单细胞 RNA-seq 工具链：**R/Seurat 6** 与 **Python/scanpy 1.13**，各带一个
私有锁定环境，按用户倾向选路。附带 [scNotebooks](https://github.com/integrativebioinformatics/scNotebooks)
13 个教学模块的对照说明；Python 路的 notebook 预设改编自
[NBIS workshop-scRNAseq](https://nbisweden.github.io/workshop-scRNAseq/home_contents.html)
的 scanpy labs。

Two parallel single-cell toolkits (R/Seurat 6 and Python/scanpy 1.13), each with its
own locked private environment, routed by user preference. Includes a map of the 13
scNotebooks teaching modules; Python notebook presets adapt the NBIS scanpy labs.

## 组成

| 组件 | 作用 |
|--|--|
| Agent `SingleCell` | 角色、边界与 **R/Python 选路规则**；不绑定环境，只构建所选路的环境 |
| Skill `seurat-workflow` | R 路可执行流程：8 个 `sc_*` 工具（含 notebook 渲染） |
| Skill `scanpy-workflow` | Python 路可执行流程：10 个 `sc_py_*` 工具（含 notebook 渲染） |
| 参考 `scnotebooks-curriculum` | 13 个模块的对照表，折在 `seurat-workflow/references/maintainers/`，只有说明文档 |
| 环境 `seurat` | R 4.5 + Seurat 6.0.0 + IRkernel，三平台 lock |
| 环境 `scanpy` | Python 3.13 + scanpy 1.13 + ipykernel，三平台 lock |

两个技能分别跑在 `plugin:seurat` 与 `plugin:scanpy`，不依赖共享 base 环境。spec 与
三平台 lock 在 `environments/<name>/`。notebook kernel 名为
`plugin-single-cell-seurat` / `plugin-single-cell-scanpy`，由 Phi 主程序在环境就绪后同步。

## 选路（agent 的职责，不是渲染器的）

- 用户明说 R/Seurat 或 Python/scanpy → 走那一路，全程不混。
- 只给文件：`.rds` → R 路；`.h5ad` → Python 路；只有 10x 矩阵目录 → 问一句倾向。
- 一路做完再换路属于新任务；跨路交接只走两个接口：`pseudobulk` 计数表（交给
  `deseq2/limma/edger` wrapper）或显式的对象互转请求。
- 规则全文在 `agents/SingleCell.md` 与两个 SKILL.md 的 Rules 首条。

## 工具

R 路（前缀 `sc`，环境 `plugin:seurat`）：

| 工具 | 用途 |
|--|--|
| `sc_inspect` | 读 `.rds`、10x 矩阵目录或计数表，报告 shape、assay/layer、降维与 metadata 列，列出已能做的步骤 |
| `sc_qc` | 线粒体/核糖体/血红蛋白比例、双细胞检测、按阈值过滤，报告每道过滤删了多少细胞，输出小提琴图 |
| `sc_preprocess` | 归一化（LogNormalize 或 SCTransform）、高变基因、scale、PCA、近邻图、多分辨率聚类、UMAP/t-SNE |
| `sc_markers` | `FindAllMarkers`，每簇一个 CSV + marker 热图 |
| `sc_annotate` | marker 面板打分标注或显式 mapping；输出 dot plot |
| `sc_pseudobulk` | 样本（或样本×分组）聚合计数表与样本表，列名按 DE wrapper 约定，返回 `wrapper_hint` |
| `sc_env_probe` | 逐模块报告缺包与跑不了的部分 |
| `sc_notebook` | 渲染 R notebook：预设模块 4/5/8/9 或 steps 数组，依赖自动补全，config cell 集中参数，可注入/合并 |

Python 路（前缀 `sc_py`，环境 `plugin:scanpy`）：

| 工具 | 用途 |
|--|--|
| `sc_py_inspect` | 读 `.h5ad`、10x `.h5`/矩阵目录或计数表，报告 shape、layers（含 lognorm 快照）、obsm 与 obs 列 |
| `sc_py_qc` | QC 指标、Scrublet 双细胞、按阈值过滤，报告每道过滤删了多少细胞，输出分布图 |
| `sc_py_preprocess` | normalize+log1p（留 lognorm 快照）、seurat_v3 高变基因、scale、PCA、Leiden、UMAP/t-SNE |
| `sc_py_integrate` | Harmony 或 ComBat 批次整合，整合前后 UMAP 对比，整合后重聚类 |
| `sc_py_markers` | `rank_genes_groups`，每簇一个 CSV + marker 热图 |
| `sc_py_annotate` | marker 面板打分标注或显式 mapping；输出 dot plot |
| `sc_py_pseudobulk` | 与 R 路同构的样本级计数表（表达过滤按单元计数，语义一致），返回 `wrapper_hint` |
| `sc_py_trajectory` | PAGA 抽象图 + 扩散伪时间（DPT），指定根细胞 |
| `sc_py_env_probe` | 报告 scanpy 生态缺包（bbknn、scVI、CellTypist 等）与跑不了的部分 |
| `sc_py_notebook` | 渲染 Python notebook：预设 NBIS labs 1–7 或 steps 数组，其余能力与 R 路渲染器相同 |

两路的 `build_notebook.py` 是**同一份字节相同的渲染器**（仓库测试断言），语言由
`--language` / `language` 参数决定；一个 notebook 一个 kernel，语言不一致明确报错。

## 两条路径（每路都有）

**边做边看（学习、探索）**：渲染 notebook，再用 Phi 的 `notebook_run_cell` 逐格跑。

```
sc_notebook(module=4, data="tenx_dir/", self_check=true)          # R
sc_py_notebook(module=3, data="batched.h5ad", self_check=true)    # Python
notebook_read → notebook_run_cell → 用户改 config cell → 只重跑受影响的格
```

模板是渲染的源材料，不是拷贝的对象：片段库（`assets/fragments/`）按 agent 传入的
步骤与参数组装成填好占位的 notebook，落在用户工作目录。选哪些步骤、阈值取多少，
是 agent 依据用户问题和 QC 分布做的判断（规则在 SKILL.md）；渲染器只做确定性填充。
所有阈值、路径、列名集中在首格 config cell；渲染溯源（步骤、参数来源、片段哈希、
数据指纹）记在 notebook metadata。`self_check` 在输出目录真跑每格，headless
（Agg + autoshow=False），不会弹原生图窗。

**一次跑完（可复现）**：

```
inspect → qc → preprocess → [integrate] → markers → annotate → pseudobulk
```

每一步读一个对象写一个对象（`.rds` / `.h5ad`），可单独重跑，参数记录在返回 JSON 里。
Python 路在 markers 之后可接 `trajectory`。两条路、两种载体的结论应当一致；
不一致就是其中一条有 bug。

**原始数据处理都不做**：比对、Cell Ranger 属于 wrapper 树。notebook 与脚本链都从
计数矩阵或已有对象开始。

## 单一分析源与知识层（0.6.0）

每路只有**一份分析代码**：`assets/fragments/`（单语言平铺；预设与参数默认值在
`assets/steps.json`）。notebook 的细胞就是片段
原文；脚本工具是薄适配器——`scripts/_fragments.{R,py}` 在同一全局环境里跑同一批
片段，适配器只把 CLI 旗标映射成与 config cell 同名的配置全局，`scripts/reports/*`
只做 JSON 序列化。片段的教学输出走 stderr，stdout 只有一份 JSON（R 用
`sink(stderr())`，Python 用 `redirect_stdout`）。重构前后逐字段核对：R 7/7、
Python 10/10 份报告与旧契约一致（仅输出路径不同）。

知识按功能归位：每步一份 `references/step-notes/<step>.md`（目的/原理/结果解读/陷阱，
即 notebook 的讲解格，与片段按步骤 id 配对）；两路各有一份
[analysis-playbook.md](skills/seurat-workflow/references/analysis-playbook.md)——agent 的
决策手册，按阶段「调用 → 读（确切 JSON 字段与图）→ 决定（参数起点+调整规则）→ 红线」，
原理压缩成一行，来源与许可在文末（OSCA CC BY 4.0、sc-best-practices Apache 2.0、
Cambridge 课程 GPL-3、NBIS CC-BY-4.0）。`references/README.md` 按功能索引：agent 决策
（playbook、step-notes、跨路交接）与维护者溯源/对照（`maintainers/`）分开。全部为原创 prose，不复制上游文本。

## 与上游课程的关系

- **scNotebooks**（CC BY 4.0，Rojas-Hidalgo et al., *Nat Genet*, 2026）：课程不是流水线
  （Colab 装包、`/content/` 写死路径、状态留在 kernel）。插件把可执行部分改写成脚本与
  片段，跑不了的部分说清楚。能跑：模块 04 主线、05 预处理半、08/09 表达量部分；
  跑不了：06（monocle3）、07（liana）、10（semla）、11（ArchR）——R 环境没有这些包；
  03（Cell Ranger）归 wrapper 树、12（SCAPE-APA）归独立工具。详见
  [模块对照表](skills/seurat-workflow/references/maintainers/scnotebooks-module-map.md)。
- **NBIS workshop-scRNAseq**（CC-BY-4.0，以仓库 LICENSE 为准；站点页脚标 GPL-3 与仓库不符）：scanpy labs 只作**参考**，不分发原文。
  片段代码为本插件原创，语义对齐 labs 1–7（QC、降维、整合、聚类+marker、样本级差异、
  注释、轨迹）；每个预设的 intro 与 SOURCE.json 注明出处与许可。移植决策见
  [NBIS 移植说明](skills/scanpy-workflow/references/maintainers/nbis-porting.md)。

## 边界

- 不改上游仓库内容，不分发它们的 notebook、图片或示例数据。
- 不在受管理环境里装包。缺包就报缺包：R 路缺 monocle3/liana/semla/ArchR/scDblFinder；
  Python 路缺 bbknn（bioconda 无 py313 build）、scVI/CellTypist/Scanorama/gseapy/空间组。
  `sc_env_probe` / `sc_py_env_probe` 逐模块报出来。
- Python 路能做轨迹（PAGA+DPT）与双细胞（Scrublet）；R 路对应模块仍标缺包，不假装能跑。
- 不做空间转录组与 scATAC-seq。

## 为什么不用 Nextflow wrapper

单细胞流程做成 wrapper 不是不行，但不合适，原因是执行模式不同：

- **wrapper 走 `-profile conda`，运行时解依赖**；**插件环境是锁 + 只读**（客户端从不解
  依赖）。Seurat/scanpy 大版本之间对象结构不兼容（本插件锁 Seurat 6.0.0 / scanpy 1.13），
  交给求解器换天换机器可能解出不同版本。
- **粒度不匹配**：单细胞分析是人在回路（看完分布定阈值、看完 marker 再注释），每步之间
  要读结果改参数；脚本工具一次调用做完一步并返回结构化结果。
- **该用 wrapper 的一步确实用了**：两路 `pseudobulk` 都把对象变成 `deseq2/limma/edger`
  适配器能吃的输入并返回 `wrapper_hint`；R 路接口已用真实 DESeq2 wrapper 实测跑通。

详细论证见 [结合设计](../../../docs/design/scnotebooks-integration.zh-CN.md) 第 8、10 节。