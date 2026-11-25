# Changelog

## 0.1.1 (2026-11-24)

- `mb_preprocess` 与 `mb_univariate` 新增 `scientificContext` 语义输出（按处理步骤自动标记 `log_transformed` / `normalized` / `raw_intensity` + `sampleIds`，以及单变量检验表的 `de_stats`），同步更新输出 JSON Schema。

## 0.1.0 (2026-11-20)

首个版本：非靶向 LC-MS 特征表下游脚本链 + 知识层，不做 notebook 渲染器，
不做上游峰提取（仓库暂无 msconvert/xcms/MZmine/MS-DIAL wrapper，报缺口）。

- 新插件 `metabolomics`，`toolPrefix: mb`，智能体 `Metabolomics`，私有环境
  `metabolomics`（Python 3.13 + pandas/numpy/scipy/matplotlib-base，与
  metagenomics 同栈；三平台显式 lock，micromamba 2.1.1 dry-run 求解）。
- 新技能 `lcms-workflow`（9 个工具）：
  - `mb_inspect`：厂商导出识别（xcms/MZmine/MS-DIAL/MetaboAnalyst/通用）、
    metadata/MSP/MGF/映射/质量库/列表分类、RT 单位推断、QC/blank 识别。
  - `mb_qc`：QC RSD（默认 30%）+ 检出率（默认 0.8）过滤，无 QC 样本即失败；
    产过滤矩阵 + 逐特征报告 + RSD 分布图。
  - `mb_preprocess`：固定顺序 impute→normalize→transform；NaN 与 0 均视为
    缺失；归一 tic/median/pqn/is；log 拒绝零/缺失。
  - `mb_multivariate`：PCA（SVD）+ PLS-DA（原创 numpy NIPALS PLS2）、
    uv/pareto 缩放、VIP、R2Y/Q2 cum（LOO，n>50 改 7-fold，折内重缩放）、
    种子标签置换 Q2/R2Y 经验 p；OPLS-DA 报缺口。
  - `mb_univariate`：Welch/Wilcoxon/Kruskal-Wallis + BH FDR + FC（logged
    标志）+ 火山图。
  - `mb_annotate`：13C 同位素、pos/neg 加合物对、H2O/CO2 中性丢失，
    ppm+RT 容差（秒/分自动换算），union-find 伪谱分组（加合物对推导中性
    质量共识）；MSI level 3 caveat。
  - `mb_ms2_match`：MSP（NIST/MoNA/RIKEN 方言）/MGF 解析、前体 ppm 门、
    sqrt 强度 + 贪心碎片配对余弦；MSI level 2 caveat。
  - `mb_pathway`：不内置数据库；ids 模式单侧超几何 + BH；masses 模式
    mummichog 式加合物假设匹配 + 背景 m/z 种子置换经验 p + BH。
  - `mb_plot`：方差 top N 热图，log10 + 平均连接层次聚类（行/列/both）。
- 知识层：analysis-playbook（调用→读→决定→红线）、9 篇 step-notes、
  source-map（book CC-BY-4.0 署名；两篇知乎 403 不可访问已注明；经典方法
  引文均为思路级，全部 prose 与代码原创）。
- 全部产物带 artifact 描述符（contract 1.0.0）。
- 仓库级测试 `tests/test_metabolomics_plugin.py`；登记 SOURCE.json 与
  根 README 计数（插件 5→6，插件内技能 9→10）。