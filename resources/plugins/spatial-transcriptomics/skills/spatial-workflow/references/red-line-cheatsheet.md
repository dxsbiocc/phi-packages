# 红线速查表（agent 运行时对照）

把 [`analysis-playbook.md`](analysis-playbook.md) §1–§10 的**硬红线**与**常见误读**汇成一表，
供运行时快速对照。本表是摘要，完整的 Call/Read/Decide 决策过程与 *why* 仍以
analysis-playbook.md 为准；每步的教学短笔记在 [`step-notes/`](step-notes/)。

## 贯穿全局（任何阶段都适用）

| 硬红线 | 常见误读 |
| --- | --- |
| 表达 of record = `layers["lognorm"]`，marker/panel/SVG/通讯绝不用 scaled 值 | 在 z-score 上跑检验/打分 |
| 空间统计需 `obsm["spatial"]` 坐标 + 邻接图；无坐标只能做非空间分析 | 在表达 PCA/UMAP 上建空间图 |
| `{"skipped": reason}` 的缺席 ≠ 干净阴性结果 | 把 skipped 当"无显著" |
| 模型法缺口（cell2location/tangram/stereopy/liana/scvi）→ 项目环境，绝不只读环境内装 | 在插件环境 pip install |
| 每个阈值是 Visium/Stereo-seq 起点而非默认；先读分布再定，并说明选了何值为何 | 默认值原样到输出 = 没做决策 |

## 分阶段

| 阶段（工具） | 硬红线 | 常见误读 |
| --- | --- | --- |
| **§1 inspect**（`st_py_inspect`） | 平台/分析单元先命名；Stereo-seq bin 主动选；无 counts 层不做 QC/seurat_v3/去卷积；`sample/group_column` 选供体/条件列 | 接受 `bin_auto`；load 时静默过滤；把 cluster 列当 group 列 |
| **§2 QC**（`st_py_qc`） | `kept_fraction`<0.5 或零单元 → 停下重读；`pct_counts_mt` 缺席/中位0 → regex/标识符错先修；QC 只在 raw counts | 不看分布接受默认阈值；never-fire 的 mt gate 当通过；pooled violin 掩盖边缘低质量 |
| **§3 preprocess+图**（`st_py_preprocess`） | 空间图建在坐标非表达 PCA；绝不在 2D embedding 聚类；一 cluster 一 sample = 批次伪影；`neighbor_method` 必须明示 | counts flavor HVG 用在 log 值；dims 只取 elbow 漏稀有域；`grid` 用于不规则 bin；只跑单一 resolution |
| **§4 markers**（`st_py_markers`） | 校正 p 描述性（double dipping）非推断；检验不在 scaled；`avg_log2FC` 非可移植效应量；marker = 高检出 AND 其余低检出 | 把 marker p 当显著性；只看大 FC 不看特异性；标识符错致 marker 全消失 |
| **§5 annotate**（`st_py_annotate`） | 零可用 panel 基因 = 标识符不匹配；报告每标签命中数 + unassigned；空间标签 = 区域/优势类型非纯细胞；粗标签不传播成亚型 | 窄 margin 险胜当鉴定；强制贴标签隐藏 unassigned；通用 panel 当组织特异 |
| **§6 SVG/邻域**（`st_py_spatial`） | 自相关/富集完全依赖图+标签；SVG p 解析（正态假设）是筛查非机制；深度/批次基因也高分；**富集是邻接非相互作用** | SVG 显著当机制；正 z-score 当信号；不看空间图信排名；小样本弱 null 仍信 z |
| **§7 去卷积**（`st_py_deconvolve`） | 比例是 NNLS **估计非细胞计数**；需 scRNA 参考；无共享基因 = 标识符不匹配先修；参考粒度/组织匹配界定答案 | 比例当细胞数；以为缺失类型可恢复；小 bin/cellbin 强行去卷积 |
| **§8 通讯**（`st_py_communication`） | **筛查：共表达非信号证明**；ligrec 读 lognorm（`use_raw=False`）；继承图+标签；符号大小写不敏感但需验证物种 pairs；starter 不穷尽 | hit 当信号证明；不确认配体/受体域相邻；人符号资源当已在小鼠验证 |
| **§9 跨样本差异**（`st_py_spatial`+`group_column`） | **SAMPLE 是实验单元非 spot**；<2 重复/组无条件级声明；描述性筛查非差异显著性检验；差异可由 composition/depth 驱动；<2 条件 skipped | spot 当重复（伪重复）；空间差异直接归因调控；不看 composition |
| **§10 notebook**（`st_py_notebook`） | 渲染器填充、你决策；每个 config 值带 source tag，`[default]` = 没决策；交付前 `self_check=true`；script 与 notebook 同阈值 | 把默认当已决策；两路阈值不一致 |

## 一句话记忆

**counts/lognorm 是地基，坐标+图是空间前提；p 值多为描述/筛查非推断；标签是区域假说、
比例是估计、通讯是邻接共表达；实验单元是 sample；skipped ≠ 阴性，缺口走项目环境。**