# Protein Design Plugin · 0.2.0

专注从头蛋白设计、蛋白-蛋白结合体（PPI Binder）与蛋白-小分子结合体设计的 Phi 专属领域插件。

## 分工边界

- **重型计算（GPU 推理、批量结构预测与流水线调度）**：统一由 `resources/wrappers/` 中的 Nextflow Wrapper 执行，通过 `wrapper_search` / `wrapper_inspect` / `wrapper_plan_run` 调度：
  - 模块：`local/modules/rfd3`、`local/modules/ligandmpnn`、`local/modules/boltz2`、`local/modules/alphafold3`、`local/modules/protenix`、`local/modules/chai1`、`local/modules/esmfold2`、`local/modules/bindcraft`、`local/modules/foldseek-easysearch`、`local/modules/foldseek-easycluster`、`local/modules/protein-metrics`、`local/modules/design-filter`
  - 子工作流：`local/subworkflows/seq-fold-eval`（支持 `predictor: boltz2 | alphafold3 | protenix | chai1`）、`local/subworkflows/design-postprocess`（支持可选 Foldseek 结构聚类去冗余）
  - 端到端工作流：`local/workflows/protein-binder-design`（支持四选一结构预测器与可选 `FOLDSEEK_EASYCLUSTER` 结构去冗余）
- **领域插件（本目录）**：提供专属专家智能体 `ProteinDesign`、三大使用教程与决策手册 Skill，以及无第三方依赖的轻量确定性 Python 校验/诊断工具（`pd_*` 前缀，复用共享 `phi:python@1` 环境）。

## 内置组件

| 组件 | 路径 | 用途 |
| --- | --- | --- |
| 专家智能体 | `agents/ProteinDesign.md` | 编排蛋白设计 Wrapper、执行烟雾测试 → Pilot → Campaign 推进与迭代优化 |
| 端到端教程 Skill | `skills/protein-design-workflow/` | 工具选型、参数决策手册、分片扩展与多轮迭代策略 |
| 输入准备 Skill | `skills/structure-input-prep/` | Contig/配体语法教程 + `pd_inspect_structure`、`pd_validate_rfd3_spec` |
| 结果评估 Skill | `skills/design-evaluation/` | 置信度/RMSD/界面/亲和力阈值解读 + `pd_funnel_report`、`pd_developability_audit` |

## 模型权重与许可说明

本插件与配套 Wrapper 仅包含原创脚本与配置（MIT），**不分发任何模型权重文件**。运行时需自行提供对应目录：
- `rfd3_weights_dir`：RFdiffusion3 / Foundry 检查点（BSD-3-Clause）
- `mpnn_weights_dir`：LigandMPNN / ProteinMPNN 检查点（MIT）
- `boltz_cache_dir`：Boltz-2 检查点与 CCD 缓存（MIT）
- `af3_model_dir`：AlphaFold 3 模型参数（受 Google DeepMind 非商业条款约束，须由用户自行申请与存放）
