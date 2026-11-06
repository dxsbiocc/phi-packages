[English](README.md) · **中文**

# 组学可视化

这是一个通用的 Agent Skill：把**已经算好的组学结果**做成可编辑的科研图。它是一套规范模板目录，不是统计引擎，也不是凭记忆现编图的地方。

技能的目标，是让 Agent **像严谨的分析者那样选图**：先说明为什么这张图能帮助展示数据，再检查数据分布，说清主张，匹配数据形态；只有预览会影响选择时才展示真实预览；最后改一份已知脚本，而不是从零即兴写 ggplot2。默认润色遵循 Nature-style 哲学：主张优先、克制、紧凑，并且在最终尺寸下可读。

给 Agent 执行的规则在 [SKILL.md](SKILL.md)。本页面向使用者和维护者。

## 能做什么

- 用两级目录（家族 → 模板）推荐图形，依据 `use_when` / `avoid_when`，而不是猜文件名。选模板是 Agent 的事；工具只负责查找、计数、复制和运行。
- 用户只想看示例时，用 `viz_catalog` 按别名或家族查找已安装的 `preview.png`，不模拟数据，也不重新绘制。
- 修改已有图时，编辑项目里的脚本并用原输入重绘。用户给参考图时，先看图再选模板。
- 选图前用 `viz_profile` 计数：列类型、缺失、范围、标识符重叠。
- 用 `viz_catalog` 查一个家族或精确名称。它不排名，也不猜测。
- 多个候选会改变科学读法时展示匹配的 `preview.png`。候选以编号列表呈现，便于对比。
- 用 `viz_prepare` 把选中的 `plot.R` **复制到你的项目**，改列名、改样式，再用 `viz_render` 输出 PDF / PNG / SVG。
- 使用 OmicsAgent ColorPicker 里的具名配色，不编造 hex。
- 多张图要支撑同一主张时，按多面板规则拼图，并做版式核查。
- 当目录和数据形态明确指向一个模板时直接推进；只有候选会改变科学读法时才展示备选。

它**不会**做差异表达、富集检验或因果推断。这些结果还不存在时，技能应明确指出前提，而不是编造数字。少数模板会从**已有矩阵**做展示用聚类（`tree-dendrogram`、聚类热图），那就是这张图本身，不能当成用户从未做过的上游分析。

## 设计哲学

### 图是论证，不是装饰

标题、视觉编码、面板顺序和不确定性都属于证据。第一个问题是：_读者看完这张图，应该能守住什么主张？_ 图种服从主张。「组学数据」并不默认画火山图。

### 图要服务数据

可视化开始前先声明：这张图要让数据中的什么更容易被看见、比较、验证或质疑。Agent 应先查看行数、取值范围、类别数量、稀疏性、缺失、离群值，以及矩阵/网络/层级/时间/基因组区间结构，再推荐模板。若用户点名的图形会隐藏关键分布，技能应推荐更合适的目录模板，并说明取舍。

### Nature-style 克制

默认风格是紧凑、证据驱动的：视觉词汇尽量少，次要标记保持安静，单位和变换写清楚，不用装饰背景，不用彩虹色板，也不依赖只有放大才能读清的标签。具体投稿规则仍以目标期刊为准；本技能提供的是设计姿态。

### Agent 选图，工具查找

Agent 不得凭训练数据推荐模板，脚本也不得按列名或关键词给模板打分。Agent 读目的、`viz_profile` 给出的事实，以及 `viz_catalog` 里的一个家族，再比较 `use_when`、`avoid_when` 和 `input_shape`。`viz_catalog` 按家族、精确 id 或精确别名（含中文）查找，不做评分或快筛。

[references/template_contracts.json](references/template_contracts.json) 列出每个模板需要的输入角色。`viz_prepare` 用这些角色检查 Agent 提供的 `role=column` 映射，不负责选模板。目录维护用 [scripts/validate_template_contracts.py](scripts/validate_template_contracts.py) 校验模板 id、source 路径、preview 路径和角色名。

用户点名某种图时，用 `alias` 查找。新图按 [references/plots.yaml](references/plots.yaml) → 对应家族 catalog 选择。别名用来理解需求；多个模板会改变读法时，不足以单独定稿。

### 预览用于消歧

使用随模板提供的预览，避免凭想象选图。若一个模板已经明确匹配主张和数据形态，记录 `id` 并直接推进。若多个可行模板会导向不同科学读法，最多列出 4 个，每项含 `id`、标题、一句 `use_when` 和 PNG。

预览 PNG 是**透明底**。深色编辑器会把空 alpha 合成成黑色。那是界面，不是图的背景，不要去「修」。

### 先复制，再改

捆绑脚本是规范库。日常使用应只把选中的 `plot.R` 拷进项目并修改这个副本，`scripts/lib/common.R` 等只读 helper 优先用安装目录的绝对路径 `source()`。只有确实需要一次性修改 helper 时才复制 `common.R`。不要为了渲染把 `references/` 或 palette catalog 拷进项目；配色从安装目录里的 references 选择 palette id。不要改技能安装目录。没有隐藏的绘图配置语言：列名和标签改 `CONFIG`，变形和排序改 `DATA PREPARATION`，只有几何必须变时才改 `PLOT`。

### 一种规范实现

选中的模板决定语言。当前模板是 R（`ggplot2` + `ggprism`）。不要仅为风格统一把已验证的脚本改写成另一种语言，也不要给同一个 `id` 再挂第二套官方实现。

### 诚实优先于好看

不得编造样本量、检验、p 值、校正 p 值、效应量或生物学含义。不得悄悄过滤、填补、截断或抽样。每次实质性变换都要报告变换前后的行数。视觉或统计上的分开，不等于生物学重要。默认用更利于色觉的 `Qualitative.Safe`；除非你明确要求装饰风，否则不用 Artwork / Concept 色板。

### 多面板是主张，不是仪表盘

带编号的拼图从**一个主主张**和必要性检验开始：去掉某面板并不少一个独特推断步骤，就合并、外移或删除。版式来自 [references/layouts.yaml](references/layouts.yaml)。几何要做核查；核查失败时，不得声称对齐已经通过。圆形分类树嵌在外圈富集轨里是**一张**模板（`tree-enrichment-ring`），不是多面板。

## 选图流程

```text
无数据、只要示例
    → viz_catalog（别名或家族）
    → 已安装的 preview.png
    → 用户选择

数据 + 主张
    → 目的
    → viz_profile             （给表计数）
    → plots.yaml              （选定一个家族）
    → viz_catalog family      （读 use_when / avoid_when / input_shape）
    → preview.png 列表        （仅对实质候选展示）
    → palettes.yaml           （只用具名 hex）
    → viz_prepare             （复制 plot.R；检查 role=column）
    → source 安装目录 helper
    → nature-style pass       （需要时做最终尺寸润色）
    → viz_render，检查成品
```

若是多面板需求，先写整图主张和每块图的证据角色，再把每一块必要的图走同一套目录。

## 模板库

**211** 个单图模板，分 **12** 个家族，另有少量拼图版式。部分模板还带 `preview_<name>.png` 的 CONFIG 变体：

| 家族       | 典型用途                                                   |
| ---------- | ---------------------------------------------------------- |
| `scatter`  | 坐标、降维、UMAP 图谱 / 环图 / 标记基因点图、泰森多边形、火山图 / 多组 / 条形码 / 双侧 GO / MA、曼哈顿图 / QQ、一维直方图 / 密度、二维密度、残差、双边误差、放大插图、哑铃图、无箱线蜂群、SVG 点、棒棒糖、三元相图、相关散点矩阵、分组拟合、通路分面、Zi–Pi、潜变量分面 |
| `heatmap`  | 矩阵、聚类、cutree 分块得分、相关、关联点阵 / 气泡 / 对角线分割 / 观察点阵、Mantel 类展示、DE 热图对齐通路富集、通路×分组富集、Circos 注释 / 队列 / 分割环、癌种变异频率、饱和突变能量、oncoprint |
| `bar`      | 大小比较、堆叠、堆叠加甜甜圈、排行富集柱 / 聚类 / 棒棒糖、细胞通讯进出柱、径向甜甜圈 / 径向折线、多轴折柱、Circos 连线、分类 SVG 标记、UpSet、韦恩图、花瓣集合图、泳道图 |
| `boxplot`  | 分布、带箱线蜂群 / 缺口 / 队列箱线、小提琴、雨云图、letter-value 箱线、配对和蝴蝶对比 |
| `line`     | 时间、排序、堆积面积、生存曲线、脊线图、GSEA、差值色带、ROC、螺旋图 |
| `pie`      | 构成比、环形、南丁格尔 / 花瓣饼图、外环柱的环状比例图       |
| `radar`    | 少量命名数值轴，含显著性标记                               |
| `graph`    | 点线网络（给定坐标，或 `graph-stress` 算法布局）、TF–靶点调控网、多组学符号网络、模块共现联图、流量 / 二分弦图 |
| `sankey`   | 守恒流量、平行集合图、等宽蛋白–性状冲积图、富集流向加排行柱 |
| `sunburst` | 环形层级，**圆心角**编码大小，含分类花瓣旭日               |
| `tree`     | 点线树、hclust 树、treemap / 圆堆积 / 冰柱、通路富集环图、环形细胞注释树、daylight / 系统发育环 |
| `ideogram` | 基因结构、染色体 G 带、位点 / 多性状位点、密度填充、基因组窗口覆盖度、蛋白突变棒棒糖、蒲公英图、Circos / 共线性 / 线性 riparian / 亚基因组 Circos、嵌套缩放环图 |

`circular` 以及 `outer = "bar"` / `"point"` 只改 CONFIG，不另开模板。`scatter-svg` / `bar-svg-icon` 的 SVG 文件也是 CONFIG。`graph-force` 使用已有 x/y（`layout = "manual"`）。由边表计算布局用 `graph-stress`。

每个模板带 `plot.R`、示例数据和 `preview.png`。会改变图形的 CONFIG 选项另存为 `preview_<name>.png`，并写在 catalog 的 `variants` 里。`viz_catalog` 按 `ids` 查询时会把这些图作为额外的 `preview_markdown` 返回。用户只想看示例、尚无数据表时，用 `viz_catalog` 按别名或家族查找内置预览，不模拟数据，也不重新绘制。catalog 记录 `id`、别名、语言、状态，以及脚本期望的输入形态。

修改已生成的图时，读取项目中的原绘图脚本和输入数据并直接修订；用户给参考图要求新绘制时，先查看参考图，再以用户数据选择合适模板。

配色见生成的[推荐配色表](references/palettes.yaml)，完整色表在[颜色目录](references/palettes/colors.json)。这些文件从 Phi 共享的 `resources/palettes/` 逐字节同步，让安装后的插件保持自包含；应修改共享源，不要直接修改此副本。默认：分组用 `Qualitative.Safe`，热图用 `Quantitative.BluGrn`，有正负的值用 `Diverging.RdBu`，对齐 OmicsAgent 产品时用 `Brand.Algolia`。

## 目录结构

```text
SKILL.md                         给 Agent 的契约
README.md / README.zh-CN.md      本说明
references/
  template_contracts.json        每个模板需要的输入角色
  plots.yaml                     家族索引
  palettes.yaml                  生成的推荐配色
  palettes/colors.json           生成的完整 hex
  catalog/<family>.yaml          模板条目
  layouts.yaml                   多面板版式
  multipanel-composition.md
  nature-figure-principles.md
  rendered-layout-qa.md
scripts/
  viz.py                         catalog / profile / prepare / render
  viz_catalog.py                 精确目录查找
  viz_profile.py                 表计数和列检查
  viz_templates.py               复制并解析 plot.R
  qa_single_plot.py              单图成品轻量 QA
  validate_template_contracts.py 目录与角色一致性检查
  lib/common.R                   共用读写与保存
  <family>/<template>/plot.R
  layouts/<layout>/compose.R
```

## 怎么用

在兼容 Agent Skill 的环境中，当你已有结果表或矩阵、需要科研图时，调用 **omics-visualization**。带上数据（或清楚的列约定）。Agent 应先声明展示目的，用 `viz_profile` 给表计数，用 `viz_catalog` 查一个家族；当一个模板明显匹配时直接推进；只有可行候选会讲出实质不同的科学故事时，才出示预览让你选。

自己跑某个模板：

```bash
skill_run({ skill: "omics-visualization", script: "bar/basic/plot.R", args: ["scripts/bar/basic/example.tsv", "output.pdf"] })
```

依赖写在各脚本头部。缺包应报告，不要悄悄安装。

## 这不是什么

- 不能代替 DESeq2、limma、clusterProfiler 或任何上游检验。
- 不是通用商业看板，也不是 AI 插画工具。
- 不是批量编写或维护模板库本身的工作流。那是同插件的 `create-template` 技能（`skill://create-template`）。

[切换到 English](README.md)
