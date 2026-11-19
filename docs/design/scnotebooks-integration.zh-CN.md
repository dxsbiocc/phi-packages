# scNotebooks 与 Phi 的结合

日期：2026-11-16（2026-11-19 更新：插件 id 由 `single-cell-training` 改为 `single-cell`；`sc_notebook` 重构为片段库 + 渲染器，见插件 CHANGELOG 0.3.0。本文第 9 节记录的是 0.2.0 的生成器方案。2026-11-19 再更新：插件改为私有环境 `plugin:seurat`（R 4.5 + Seurat 6.0.0 + IRkernel），不再借用 `phi:r@1`，notebook kernel 为 `plugin-single-cell-seurat`，见 CHANGELOG 0.4.0；第 3 节的环境核查是 0.1.0 时针对 phi-r 的记录，缺包结论对私有环境同样成立。）
上游：[integrativebioinformatics/scNotebooks](https://github.com/integrativebioinformatics/scNotebooks)，commit `55229bb496f5b6b96eae266ac867e5068d600464`
落地：[`resources/plugins/single-cell/`](../../resources/plugins/single-cell/)

## 1. 上游是什么

13 个单细胞与空间组学教学 notebook，每种语言一套（English / Español / Português / 中文），共 52 个文件。

| 事实 | 数值 |
| --- | --- |
| notebook 数 | 13 × 4 语言 |
| 英文代码单元 | 510 |
| 语言之间代码是否相同 | 13/13 完全相同，翻译只改文字 |
| 未压缩 notebook 体积 | ~37 MB（Module 01 一个就 14 MB，几乎全是 base64 图片） |
| 去掉 base64 与输出后 | 1.12 MB |
| 运行期装包调用 | 46 处 `install.packages`/`BiocManager::install`/`remotes::install_github` |
| Colab 专用产物 | 45 处 apt/bspm 初始化、7 处 `/content/`、38 处 `shell_call`、53 处 `!` 壳命令 |
| 许可证 | CC BY 4.0（Rojas-Hidalgo et al., *Nat Genet*, 2026） |
| 附带的资源目录 | `scNotebooks-Resources/`，Git LFS 指针 + 一个 147 MB 的 tarball |

**关键结论：这是课程，不是流水线。** 代码单元是给人一课一课读的，旁边就是讲解。每个 notebook 自己在 Colab 里装包，路径写死 `/content/`，状态留在 kernel 变量里（`sampleA`、`sc_merged`、`SpatialData`）。

## 2. 为什么不能直接转换

Phi 的三种内容类型都要求「装进去就能跑」，而上面三件事恰好都违反这一点：

| 上游特征 | 与 Phi 的冲突 |
| --- | --- |
| 运行期装包 | 受管理环境创建后**只读**；在环境里 `install.packages` 是破坏可复现性，不是功能 |
| `/content/` 路径 | 引擎按项目目录解析路径 |
| kernel 内状态 | 工具调用没有 kernel |

所以 13 个 notebook 不能机械转成 13 个 skill，也不能把 `.ipynb` 塞进包里——那只是把不能跑的东西换个地方放。

## 3. 环境的真实情况（这一步决定了方案）

查 `Phi/resources/runtime/environments/phi-r/environment.yml`：

| `phi-r` 直接提供 | `phi-r` 没有 |
| --- | --- |
| Seurat、SeuratObject、SingleCellExperiment、ComplexHeatmap、jsonlite、ggplot2、patchwork | monocle3、liana、SoupX、DropletUtils、semla、ArchR、scDblFinder、pheatmap、clusterProfiler、org.Hs.eg.db、SingleR、celldex、miloR、clustree、presto、hdf5r |

`Matrix`、`uwot` 作为 Seurat 的依赖间接存在，但不作为独立声明（`harmony` 不是 Seurat 依赖，实际缺包）。

对缺的包查了 conda 是否有货，因为「加进 phi-r」只有在包存在时才是一个选项：

- 有：`r-monocle3`（bioconda 0.2.0）、`r-archr`（bioconda 1.0.3，仅 osx-arm64）、`r-soupx`、`r-clustree`、`r-harmony`、`bioconductor-milor`、`r-presto`、`r-hdf5r`。
- 没有：`r-liana`、`r-semla`、`r-singlet`、`r-seuratwrappers`、`bioconductor-dropletutils`、`r-cellchat`——这些要走 `sourcePackages` 钉住 commit。

这个事实直接否掉了「把 13 个模块都做成能跑」的想法，也让方案从「移植课程」变成「移植能跑的部分 + 把跑不了的部分变成工具能回答的问题」。

## 4. 落地：一个插件，多个技能

```
resources/plugins/single-cell/              toolPrefix: sc
  phi-package.yaml                          environments: seurat / scanpy（两个私有锁定环境；技能用 plugin:<name>）
  agents/SingleCell.md                      角色、边界与 R/Python 选路（不绑定环境，见 §10）
  skills/seurat-workflow/                   8 个脚本工具（含 notebook 渲染），可执行
    scripts/_common.R                       共用参数/IO/JSON 输出
    scripts/project_inspect.R               sc_inspect
    scripts/qc_filter.R                     sc_qc
    scripts/preprocess_cluster.R            sc_preprocess
    scripts/find_markers.R                  sc_markers
    scripts/annotate_clusters.R             sc_annotate
    scripts/pseudobulk.R                    sc_pseudobulk
    scripts/env_probe.R                     sc_env_probe
    schemas/env-probe.json                  输出校验
  skills/scanpy-workflow/                   10 个脚本工具（含 notebook 渲染），0.5.0 新增，见 §10
  skills/scnotebooks-curriculum/            只有说明，无脚本
    references/module-map.md                13 个模块逐条对照
```

为什么是插件而不是 13 个独立 skill：同一套环境、同一个角色、同一条链。拆成 13 个会把读对象和解析参数的代码复制 13 份，目录里也会多出 13 个近乎相同的条目。为什么是两个技能而不是一个：课程对照只有文档、没有脚本，单独放才能让人只问「模块 6 讲什么」时不必拉进一整套 R 流程。0.5.0 起执行技能再按工具链分成 `seurat-workflow` 与 `scanpy-workflow`：环境不同、对象格式不同，但角色与出口相同（§10）。

### 工具链

| 工具 | 替代上游哪一段 | 关键差异 |
| --- | --- | --- |
| `sc_inspect` | Module 04「Data input」、Module 09 读矩阵 | 一次报告 shape/assay/layer/metadata，并列出这个对象已能做哪几步；`.h5` 和 `.h5ad` 明确拒绝并给替代路径 |
| `sc_qc` | Module 04 质控、双细胞、SoupX 段 | 阈值变参数，报告每道过滤删了多少细胞和指标分布，存活率低于一半给警告 |
| `sc_preprocess` | Module 04 归一化→聚类全链 | 一次确定性 pass，记录方差解释、簇大小、参数；重复 `--resolution` 得到多列 |
| `sc_markers` | Module 04 差异表达 | 每簇一个 CSV，整表挂到对象上供下一步复用 |
| `sc_annotate` | Module 04 注释（SingleR 与手工面板两条路） | 面板是 JSON，打分是「面板均值减其余均值」，面板大小可比；mapping 必须覆盖所有簇 |
| `sc_pseudobulk` | 上游**没有**这一步 | 按样本×分组聚合，产出差异表达 wrapper 直接能吃的计数表 + 样本表 |
| `sc_env_probe` | 上游的安装单元 | 逐模块算出缺哪个包、核心能否跑；把「装什么」变成可查询的事实 |
| `sc_notebook` | 模块 04/05/08/09 的可执行部分 | 渲染器：片段库 + steps/config 组装成 notebook，self-check 真跑每格（§9）；0.5.0 起与 Python 路共用同一份渲染器（§10） |

`sc_pseudobulk` 是唯一一个上游没有、但补上后才完整的一步：notebook 全程单样本描述细胞类型，从不做样本间比较，而这正是差异表达和审稿人要看的东西。

### 比的是行为，不只是代码

| 上游 | 这里 | 原因 |
| --- | --- | --- |
| 看图定阈值 | 参数 + 分布报告 | 流水线要能在没人看图时复现 |
| 总是跑 scDblFinder | `--detect-doublets`，默认关，尽力而为 | 包不在 phi-r；静默跳过比明确说没跑更糟 |
| 面板写在代码单元 | JSON 面板文件 | 可评审、面板大小可比 |
| 手动 subset 删细胞 | 每道过滤的删除计数 | 可审计，阈值过严会警告 |
| 单样本分析 | 样本×分组 pseudobulk + design 报告 | 接上差异表达 |

## 5. 模块对照结论

| 结论 | 模块 |
| --- | --- |
| 核心能跑 | 04 主线、05 的预处理一半、08、09 的表达量部分 |
| 只能读 | 01（数据库部分由已有 connector 覆盖）、02（绘图归 Visualization 插件）、13 |
| 不属于 R 流程 | 03（Cell Ranger，归 wrapper 树）、12（SCAPE-APA，Python/BAM/GFF） |
| 缺包跑不了 | 04 的 Milo/SoupX/GO/scDblFinder 段、06（monocle3）、07（liana）、10（semla）、11（ArchR） |

这份划分不是写在文档里的断言，而是 `sc_env_probe` 每次调用现算出来的。

## 6. 验证

- `node --test tests/*.test.mjs`：34 项全过（含 README 计数与内容边界）。
- Phi 校验器对 `phi-packages` 源：`check-plugins` 3 个插件 0 错 0 警告；`check-skills` 23 个技能 0 错 0 警告。
- 端到端（本机 R 4.4.2 + Seurat 5.1，合成数据 400 细胞 / 300 基因）：inspect → qc（400→396）→ preprocess（2 簇）→ markers（41 个，2 簇）→ annotate（2 个标签）→ pseudobulk（4 个样本，300 基因），写出 6 张图、2 张 marker 表、计数表与样本表。
- 参数按引擎实际拼法验证：下划线名、布尔不带值、数组重复传、绝对路径。
- 失败路径验证：缺列、mapping 不全、缺包，都返回 `{"error": ...}` 且给出可执行的下一步。

## 7. 已知边界

- 上游 notebook、图片、示例数据一概不复制；课程正文仍在原仓库（CC BY 4.0，保留署名）。
- `phi-r` 里没有的包不装。想跑 06/07/10/11 需要项目环境（`env_request`）加包，`r-liana`、`r-semla`、`r-singlet` 还得走 `sourcePackages`。
- 本机 Seurat 5.1 / ggplot2 4.0 下 `VlnPlot(ncol=)` 会因 patchwork 的 `&` 算子报错，脚本改成逐个子图 `wrap_plots` 规避；`phi-r` 也是 ggplot2 4.x，所以这不是本机特例。
- 双细胞检测、SingleR、`.h5` 读取都依赖缺包，实际是「可发现但当前不可用」，返回里写明。

## 8. 为什么是脚本工具而不是 Nextflow wrapper

这是设计评审时被问到的问题，答案不是偏好，是两种执行模式的区别。

**wrapper 用的是 `-profile conda`，Phi 的环境是「锁 + 只读」。** 两者由 Phi 的两段不同代码负责：

| | wrapper（`conda-profile.ts`） | skill / 脚本工具（`phi-r`） |
| --- | --- | --- |
| 依赖从哪来 | 每个 Nextflow 进程的 `conda "${moduleDir}/environment.yml"` | 预先按 explicit lock 建好的 prefix |
| 谁解依赖 | Nextflow 运行时解，用自带 micromamba（`conda.useMicromamba = true`） | 客户端**从不解**，CI 出锁 |
| 可复现性 | 每次运行依赖求解结果 | 锁定了精确 URL + md5 |
| 读取方式 | 只读 prefix，`R_LIBS_USER` 指向 prefix | 同左 |

`phi-r` 的契约写着「clients never solve」，而 wrapper 的 conda profile 是**刻意保留的例外**（源码注释原话：*"That is the one intended exception"*）。所以把 Seurat 流程做成 wrapper，等于把一个需要固定锁的分析链交给运行时求解器——换平台、换一天，conda 可能解出不同版本的 Seurat，而 Seurat 5 与 4 的对象结构不兼容。

**wrapper 树里没有 Seurat。** `resources/wrappers/` 有 148 个 nf-core 模块，单细胞侧只有 `scanpy`（filter/pca/scrublet/hashsolo）、`anndata`、`cellbender`——全是 Python。`grep -ic seurat` = 0。也就是说 scNotebooks 的 R 内容在 wrapper 树里没有现成的落脚点，要做也是从零写 13 个 `main.nf`。

**粒度不匹配。** wrapper 的调用方式是「发现 → inspect → 跑一次 → 拿输出目录」。单细胞分析是**人在回路**的：看完 QC 分布再定阈值，看完 marker 再标细胞类型，看完高变基因再定 PC 数。每步之间要读 JSON、改参数、重跑。为每一步付一次 Nextflow 启动 + 依赖求解（实测本机 Nextflow 冷启 1 秒，但首次 conda 求解要几分钟），而且中间状态要落盘成 `.rds` 才能传递。脚本工具一次调用做完一步、直接返回结构化结果，正好对上。

**真正该用 wrapper 的地方，确实用了。** 单细胞做完要答「两个条件有没有差异」，那一步是典型批处理、无人在回路、输入是两张表——正是 wrapper 的形状。仓库里已经有 `deseq2 / limma / edger / timeseries` 四个适配器。所以 `sc_pseudobulk` 的职责不是重新实现统计，而是**把 Seurat 对象变成 wrapper 能吃的输入**。

### 这个接口最初是错的，已经修好

第一版 `sc_pseudobulk` 写出 `gene` + `S1_control` 列、样本表用 `column/sample/group`。对照 wrapper 的 fixture 才发现在两处对不上：

| | 第一版（错） | wrapper 实际读 | 现在 |
| --- | --- | --- | --- |
| 计数表首列名 | `gene` | `gene_id`（`--gene_id_col` 默认值） | `gene_id` |
| 样本表连接列 | `column` | `sample`（`--sample_id_col` 默认值） | `sample`，另留 `source_sample` |

修完后**实测跑通**：用 micromamba conda profile 跑 `deseq2/main.nf`，读我产出的两张表，DESeq2 正常算出 size factor（1.008 / 1.017 / 1.013 / 0.995）、写出 `deseq2.results.tsv`（121 行）。最后一个 VST 作图步骤因为合成数据只有 120 个基因而报 `less than 'nsub' rows`——这是测试数据太小的预期结果，不是接口问题。

`sc_pseudobulk` 现在返回 `wrapper_hint`，把这些参数值直接交给 wrapper：`counts`、`samplesheet`、`sample_id_col: sample`、`gene_id_col: gene_id`、`contrast_variable`。没有分组列时该字段为 `null`——因为没有可比的条件，不能给一个会误导的参数。

### 如果确实想走 wrapper

那么正确的形状不是「把 Seurat 流程改成 wrapper」，而是：给 wrapper 树补一个 `modules/local/seurat/<step>/`，每个模块自己的 `environment.yml` 把 Seurat 和需要的包钉到具体版本，配 `wrapper/` 适配器。代价与收益：
- 收益：HPC/远程执行、Nextflow 的缓存与重试、运行记录。
- 代价：依赖变成运行时求解（失去锁的确定性）；每步要落盘；人在回路的迭代变慢；13 个 `main.nf` + 适配器 + 冒烟测试的维护量。

对**教学与探索**（scNotebooks 的本职）不合适；对**已定型的批量重跑**（比如每周对 200 个样本跑固定流程）合适。这两件事可以并存：脚本工具确定参数，定型后再固化成 wrapper。

## 9. Jupyter/notebook 模式：为什么之前没做，现在补上了

这一节记录一次设计失误和它的修正。

### 我没有用 notebook 模式的真实原因

不是不需要，而是我一开始把问题搞错了对象：我在**核实能分发什么**，而用户问的是**能用什么**。

查过的事实：
- `PackageType = 'skill' | 'plugin' | 'wrapper' | 'mcp'`——**没有 `notebook` 类型**。
- skill 契约第 23 行：*"Nothing is executed from `references/` or `assets/`"*。所以把 `.ipynb`
  放进包里只是占地方，不会被执行。
- notebook 工具（`notebook_list/read/insert_cell/update_cell/delete_cell/run_cell/save`）
  是引擎内置 core tools，固定扫描**用户项目目录**，不接受包提供 notebook。

由此我得出「不能做成 notebook 模式」——**这一步推理是错的**。正确的结论是：**不能"分发"
notebook，但完全可以在用户项目里"生成"notebook**。而 `phi-r` 已经装了 `r-irkernel`，
`phi-jupyter` 提供 server，`notebook_run_cell` 可用。能力一直在那里。

这个失误的性质值得记下来：我把"内容分发单位的限制"误当成"产品能力的限制"，于是在一个
本质是 notebook 的课程上，交付了纯脚本，丢掉了教学载体。

### 现在的方案：生成，不分发

`sc_notebook` 工具（`scripts/build_notebook.py`）在用户项目里写出一个 `.ipynb`：

- **每格是真实的 Seurat 代码**，不是 shell 调脚本。学员看到 `FindVariableFeatures(...)`
  本身，可以改参数。
- **每格上面一段 markdown 讲为什么**，按你的要求：讲解 + 可运行 cell。
- **按需逐格运行**，用 Phi 现成的 notebook 工具。改完阈值只重跑那一格。
- **`self_check`** 生成时把每格真的执行一遍，任一格失败就整体失败。交付给学员的
  notebook 必须是跑得通的，而不是看起来能跑。

模块 4/5/8/9，都通过 self-check（16–17 个代码格，实测全过）。

### 两条路的分工

| 场景 | 用什么 | 特征 |
| --- | --- | --- |
| 边做边看、改参数、教学 | `sc_notebook` + `notebook_run_cell` | 真代码在内联格子里，逐格可控 |
| 一次跑完、可复现、能审计 | `sc_*` 脚本工具 | 每步 `.rds`，参数记在 JSON |
| 原始测序数据 → 计数矩阵 | wrapper 树（Cell Ranger 等） | 批处理，无人在回路 |
| 样本间差异表达 | `sc_pseudobulk` → `deseq2` wrapper | 已实测跑通 |

两条路应当给出一致结论。不一致就是其中一条有 bug——这本身是个有用的交叉验证。

### 生成过程中实测踩到的三个坑（已写进生成代码）

1. `Read10X` 只接受 `.gz`（`matrix.mtx.gz` / `features.tsv.gz` / `barcodes.tsv.gz`）。
2. `CreateSeuratObject(min.features = N)` 若把所有细胞都删掉，Seurat 不报"细胞为空"，
   而是在后面抛 `invalid class "LogMap" object: Rownames must be supplied`。现在先算
   每细胞检测基因数并打印分布，阈值会删光时明确拒绝并说明可能原因（子集矩阵、
   feature ID 类型不对、文库过浅）。
3. 面板基因在对象里一个都没有时，`DotPlot` 直接报错；更糟的是打分全为 `NA` 却仍然
   输出一个"最高分"标签——等于凭空编细胞类型。现在只画存在的基因，一个可用基因都没有
   就停下来，提示检查标识符类型（人 `MT-` / 鼠 `mt-` / Ensembl ID）。

第 3 条是这次最值得记的：一个"不报错但输出错误答案"的分支，比崩溃危险得多。

## 10. 双工具链：R/Seurat 与 Python/scanpy（0.5.0）

用户决定：插件内维护两套工具链、两个独立私有环境，按用户倾向 R 还是 Python 选路；
NBIS workshop-scRNAseq 的 scanpy labs 作为 Python 路 notebook 的模板参考。

### 为什么是插件内两套，而不是两个插件

同一个角色、同一条边界、同一个 pseudobulk 出口。拆成两个插件会把 agent 与课程对照
复制两份，且把「选 R 还是 Python」变成**安装时**的决定；而它应该是**每次任务**的
决定——同一个人这周做 Seurat 项目、下周接 scanpy 的 `.h5ad`。环境则必须分开：
scanpy 1.13 要求 python >=3.13，与 R 环境无任何共享依赖，锁也各自三平台。
bbknn 在 bioconda 没有 py313 build：报缺口，不降 python、不在受管理环境装包
（§7 的边界不变）。

### 选路归 agent，环境归技能

契约顺序是「技能自带环境 → agent 会话环境 → 报错」。据此 `SingleCell` **不声明
environment**：选哪路只构建哪路的环境，agent 的 bash 落在 host（只做文件操作）。
选路规则写在 agent 与两个 SKILL.md 的 Rules 首条（"R or Python is the user's route,
not yours."）：显式倾向优先；只给文件时 `.rds`→R、`.h5ad`→Python；只有 10x 矩阵先问；
全程不混路；跨路交接只走 pseudobulk 计数表或显式互转请求。渲染器与脚本工具不参与
选路——它们只执行传进来的语言与步骤。

### 共享渲染器，字节相同

两个技能的 `scripts/build_notebook.py` 是同一份文件（仓库测试断言字节相等），
语言由 `steps.json` 的 `language`/`kernels` 与 `--language` 决定；一个 notebook 一个
kernel，注入/合并语言不一致明确报错。预设：R 路模块 4/5/8/9，Python 路 NBIS labs 1–7
（QC、降维、整合、聚类+marker、样本级差异、注释、轨迹）。片段库各语言独立
（当时为 `fragments/r/`、`fragments/python/`，0.6.0 去掉冗余语言层平铺，见 §11），语义对齐：同一次质控读数传给两条路、
阈值是参数不是猜测。

self-check 的一次实测教训：验证子进程跑在 Jupyter **之外**，macOS 默认 macosx
backend 加 scanpy `autoshow=True` 会每张图开一个阻塞的原生窗口（用户真的看到了弹窗）。
修正：self-check 子进程 `MPLBACKEND=Agg` + `autoshow=False` 前缀；notebook 内联展示
不受影响。值得记的一般结论：**验证代码的执行环境不等于交付环境**，验证必须自己
声明 backend，不能继承宿主的交互式默认。

### scanpy 1.13 / anndata 0.13 的语义迁移

- accessor 字符串取代裸列名与旧关键字：`harmony_integrate(..., "obs.<batch>")`、
  `hvg(batch="obs.<batch>")`、`neighbors(use="obsm.X_pca", out=..., rng=...)`、
  `leiden(out="obs.<col>", flavor="igraph")`、`pca(mask="var.highly_variable")`。
  ComBat 仍收裸列名（实测），不统一改成 accessor。
- `.raw` 在 anndata 0.13 弃用：改为归一化+log1p 后立刻留 `layers["lognorm"]`
  快照（scale 前、全基因），markers/annotate/ComBat 读它；inspect 报 `snapshot`。
  语义与 Seurat 路「RNA assay 保留 log 表达、scale 只进 PCA」对齐。
- 显式 `diffmap` 再 `dpt`、热图/dotplot 前显式 `dendrogram`、`sc.io.*` 读入、
  `importlib.metadata` 取版本：把弃用警告变成显式步骤，而不是留着默认行为漂变。
- 核对而非假设：迁移前后注释打分与聚类成员一致（igraph flavor 只交换簇编号）；
  pseudobulk 表达过滤两路都按**单元**（样本×分组的列）计数而非细胞数，与
  `pseudobulk.R` 同语义——合成混淆数据上曾因此「没有基因过过滤」，是参数问题不是 bug。

### NBIS 作为模板源与许可

NBIS labs 许可以仓库 LICENSE 为准：CC-BY-4.0（站点页脚标 GPL-3，与仓库不符，
0.6.0 已统一更正）。只作参考，不分发原文。片段代码为本插件原创，语义对齐 labs；
每个预设 intro 与 SOURCE.json 注明出处与许可。与 scNotebooks（CC BY 4.0）同一原则：
课程不是流水线（§2），插件改写可执行部分、说清楚跑不了的部分。Python 路补上后，
轨迹（PAGA+DPT）与双细胞（Scrublet）在 Python 路能跑；R 路对应模块仍标缺包。

### 吸收独立 scanpy 技能

原 `resources/skills/scanpy`（K-Dense BSD-3）与插件 Python 路能力重叠：吸收进
`scanpy-workflow` 后删除独立技能。同一能力只有一个家；独立技能 17→16、插件内技能
6→7，README、SOURCE.json 与计数测试同步。

### 验证

Python 预设 1–7 self-check 全过（8/13/15/16/17/18/19 格）；R 路回归 m4/m5/m8/m9
全过（17/18/18/18 格）；Python 脚本链 inspect→qc→preprocess→integrate（Harmony 与
ComBat）→markers→annotate→pseudobulk→trajectory 全过；仓库 35 项测试通过；
Phi 校验器 `check:plugins --packages` 3 插件 0 错、`check:skills --packages` 23 技能 0 错。

## 11. 单一分析源与知识层（0.6.0）

0.5.0 之后每路仍有两份分析代码：片段库（notebook 用）与脚本工具（JSON 链用）。
语义对齐靠人工核对维持，改一处忘另一处只是时间问题。0.6.0 把两路都统一成
**一份分析源**，并把「为什么这样做、结果怎么读」写进技能，让 agent 对插件的使用、
分析、解读有依据，而不是只看见旗标。

### 架构：片段库 + 薄适配器 + 报告层

```
assets/steps.json              预设与参数默认值（渲染配置）
assets/fragments/              唯一分析代码（教学代码 + 记账全局；单语言技能平铺，无冗余语言层）
references/step-notes/         每步讲解 md（notebook 说明格同源，与片段按步骤 id 配对）
scripts/_fragments.{R,py}      在同一全局环境执行同一批片段
scripts/<tool>.{R,py}          薄适配器：CLI 旗标 → 与 config cell 同名的配置全局
scripts/reports/*              只做 JSON 序列化（脚本专属）
```

- notebook 的细胞就是片段原文；脚本工具不再有自己的分析实现。适配器只做三件事：
  映射旗标、在片段之间放必要的守卫（如整合必须有 batch 列）、决定跑哪些片段。
- 记账全局（figure_paths、qc_filters、markers、pb_info……）在 setup 片段里初始化，
  报告层只读它们。R 用 `sys.source` 进 globalenv，Python 用 `exec(compile(...))`
  进共享 namespace，两路同构。
- stdout 纪律：片段的教学输出（cat/print）走 stderr——R 适配器 `sink(stderr())`，
  Python `run_fragments` 用 `contextlib.redirect_stdout(sys.stderr)`；stdout 永远
  只有一份 JSON。这是「脚本返回结构化结果」契约的一部分。
- notebook 与脚本的默认值差异保留在适配器覆盖里（如 detect_doublets notebook
  默认开、脚本默认关；min_genes notebook 200、脚本不设），不改片段语义。
- 验证方式是**逐字段 parity**：重构前先跑一遍旧脚本存 JSON 基线，重构后同输入
  再跑，归一化输出路径后深比较。R 7/7、Python 10/10 一致。契约没有变，变的只是
  代码的家。

### 知识层：step-notes + analysis-playbook.md

- 每步一份 `references/step-notes/<step>.md`（两路共 30 个）：目的 / 原理 / 结果解读 / 陷阱，
  短 prose，渲染 notebook 时作为该步的说明细胞，agent 读片段时按步骤 id 读它。
  Markdown 全部归 references 按功能分组（agent 决策 / 每步解读 / 跨路交接 /
  `maintainers/` 溯源），assets 只放代码与渲染配置；目录按功能重排，命名统一
  kebab-case 且揭示功能（`cross-route-handoff.md`、`maintainers/implementation-notes.md`）。
- 两路各一份 `references/analysis-playbook.md`（初名 analysis-principles.md，
  原理essay；按「教 agent 决策」的标准重写并更名）：按阶段「调用 → 读（确切 JSON
  字段与图）→ 决定（参数起点+调整规则）→ 红线」。操作细节以代码契约为准，例如
  `sc_qc` 首读要放宽 `mt_max`（适配器默认 10，会先 gate 掉你想读的分布）、阈值起点
  通过 distribution/MAD 交叉核对、R 路脚本链没有整合工具（batch-compare 只诊断）。
  原理压缩成一行「why」，来源与许可移到文末；`references/README.md` 按受众索引
  （agent 决策 vs 维护者溯源/对照）。
- 来源核实过许可：OSCA 在 bioconductor.org/books/release（CC BY 4.0）；
  sc-best-practices.org 是 Apache 2.0（引 Heumos, Schaar, Lance et al., Nat Rev
  Genet 2023, DOI 10.1038/s41576-023-00586-w）；singlecellcourse.org 是
  Cambridge/Hemberg 课程（GPL-3），**不是 OSCA**，引用时不能混。
- 全部原创 prose：蒸馏观点，不复制上游文本/图/数据——与 §2 的原则一致。

### 课程折叠与许可更正

- `scnotebooks-curriculum` 独立技能折进 `seurat-workflow/references/maintainers/`（curriculum
  与 module-map 两份文档）。课程对照是 R 路教学的一部分，不值得一个独立技能；
  插件内技能 3→2、仓库插件内技能 7→6、技能总数 23→22。§4 的目录树是 0.5.0 布局。
- NBIS 许可以仓库 LICENSE 为准：**CC-BY-4.0**。站点页脚标 GPL-3 与仓库不符，
  0.5.0 沿用了页脚的错误；0.6.0 在 steps.json attribution、7 个预设 intro、
  nbis-porting.md、SOURCE.json、本文 §10 与 README 统一更正。

### 顺手修掉的两个显示 bug

- `show_plot` 的内联 `print(plot)` 在 Rscript 下会打开默认 PDF 设备，在用户工作
  目录留 `Rplots.pdf`。现在只在 IRkernel 命名空间已加载时内联打印（核实过
  IRkernel 并不把 `jupyter.in_kernel` 置 TRUE，不能用它判断）；PNG 落盘不受影响。
- `adt-assay` 的 VlnPlot 在块中间，notebook 里从不显示；`clonotypes` 的 DimPlot
  靠细胞末表达式自动打印。两者都改走 `show_plot`，显示与落盘一致。

### 验证

两路 JSON parity 全一致（仅输出路径不同）；R 预设 m4/m5/m8/m9 与 Python 预设
p1–p7 self-check 全过；仓库 35 项测试通过；Phi 校验器 `check:plugins --packages`
3 插件 0 错、`check:skills --packages` 22 技能 0 错。
