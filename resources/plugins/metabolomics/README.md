# metabolomics（代谢组学分析）· 0.1.0

非靶向 LC-MS 特征表的**下游分析链 + 知识层**：厂商导出识别（xcms/MZmine/
MS-DIAL/MetaboAnalyst）→ QC 样本 RSD/检出率过滤 → 填补/归一化（TIC/median/
PQN/内标）/log 变换 → PCA 与 PLS-DA（VIP + 种子置换 Q2）→ 单变量检验
（Welch/Wilcoxon/Kruskal-Wallis + BH FDR + 火山图）→ 确定性 m/z 注释
（同位素/加合物/中性丢失 → 伪谱）→ MS2 谱库匹配（MSP/MGF 余弦）→ 通路富集
（ID 超几何 / mummichog 式 m/z 置换）→ 聚类热图。
上游（msconvert 转换、xcms/MZmine/MS-DIAL 峰提取与对齐）**不在本插件**，
仓库暂无对应 wrapper——知识层明确报缺口并向用户要导出表。
插件脚本只做确定性下游统计，跑在一个私有锁定 Python 环境里。

Untargeted LC-MS feature-table downstream chain plus a knowledge layer:
vendor export inspection → QC RSD/detection filtering → imputation/
normalization (TIC/median/PQN/IS)/log transform → PCA & PLS-DA (VIP +
seeded permutation Q2) → univariate tests (Welch/Wilcoxon/KW + BH FDR,
volcano) → deterministic m/z annotation (isotopes/adducts/neutral losses →
pseudospectra) → MS2 library matching (MSP/MGF cosine) → pathway enrichment
(ID hypergeometric / mummichog-style m/z permutation) → clustered heatmaps.
Upstream peak picking (msconvert, xcms, MZmine, MS-DIAL) is out of scope and
has no wrapper here; the knowledge layer reports the gap honestly.
All scripts are deterministic downstream statistics in one private locked
Python environment.

## 组成

| 组件 | 作用 |
|--|--|
| Agent `Metabolomics` | 路线判断（特征表下游 vs 上游缺口）、QC/阈值决策、MSI 鉴定分级表述、数据库许可确认 |
| Skill `lcms-workflow` | 9 个 `mb_*` 工具（inspect/qc/preprocess/multivariate/univariate/annotate/ms2_match/pathway/plot） |
| 环境 `metabolomics` | Python 3.13 + pandas/numpy/scipy/matplotlib-base，三平台 lock（与 metagenomics 同栈） |

spec 与三平台 lock 在 `environments/metabolomics/`。v0.1.0 只有脚本链，
没有 notebook 渲染器。

## 边界（agent 的职责）

- 已有特征表（+metadata/QC 标注）→ 本插件下游链。
- 只有原始 `.raw`/`.mzML` → 上游缺口：向用户要 xcms/MZmine/MS-DIAL 导出表，
  不假装能峰提取。
- GC-MS 特征表机械兼容，但保留指数（RI）不支持；NMR、靶向 MRM、脂质种属
  命名不在范围。
- OPLS-DA、QC 漂移校正（crmn/MetNormalizer）、GSimp 填补、CAMERA 完整伪谱、
  分子式解析、SIRIUS/MetFrag/GNPS → 知识层报缺口，不硬造。
- 通路/化合物/谱库**不内置**：KEGG/HMDB 映射、质量数据库、MSP/MGF 库由用户
  提供，先确认许可与版本再使用。

## 知识层

`skills/lcms-workflow/references/`：analysis-playbook（调用→读 JSON→决定→
红线）、9 篇 step-notes、maintainers/source-map。主要来源为
RforMassSpectrometry《Metabolomics Book》（CC-BY-4.0，思路与结构借鉴，
全部文字原创、无代码复制）；用户提供的两篇知乎教程 403 不可访问，已在
source-map 注明。

## 验证

仓库级 `tests/test_metabolomics_plugin.py`（结构 + 功能 + 反例 + 确定性，
合成夹具种子固定）；`node --test tests/*.test.mjs` 覆盖插件内技能契约。