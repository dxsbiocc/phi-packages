# 运行时与内容分发实施计划

日期：2026-09-29（按 [运行时基础设计](../design/phi-runtime-foundation.zh-CN.md) §8 重写；同日补充：R 源码包、删除 docx 和 scvi-tools skill、`phi-nextflow` 纳入 nf-core 和 nf-test）
更新：2026-10-07 — 可视化依赖并入 `phi-r`。
依据：[运行时基础设计](../design/phi-runtime-foundation.zh-CN.md)（已确认）、[内容分发设计](../design/phi-content-distribution-design.zh-CN.md)、[决策记录](../decisions/content-distribution.zh-CN.md)。
英文版：[content-distribution-implementation.md](content-distribution-implementation.md)。

## 1. 原则

1. **由内而外，由小到大。** 步骤顺序由运行时基础设计的分层决定：运行时 → 执行原语 → 使用方 → 绑定 → 插件 → 分发 → 远程。每一步只依赖它内侧已经完成的层。
2. **一步验收通过，才进入下一步。** 每一步都有明确的完成标准；未通过就不开始外层工作。
3. **契约随所属层一起冻结。** 某一层完成时，它对外暴露的契约冻结为 v1，之后只按内容分发设计 §4.4 做增量修改。
4. **主线之外的工作单独列出。** 不属于这条分层主线、也不修改内核的工作（数据访问改造、多智能体编排）列为支线，不打乱主线顺序。
5. 每个子任务 0.5–3 人日，能单独合并、单独回退，合并后 app 照常可用。

## 2. 总览

| 步骤 | 层 | 目标 | 冻结的契约 | 完成标准 | 估时 | beta |
|---|---|---|---|---|---|---|
| 0 | 准备 | 清理和基线，不涉及新架构 | — | 安装包干净；有依赖清单和上下文基线 | 3 天 | 是 |
| 1 | L0 + L1 | 内置 micromamba、`~/.phi/runtime`、按锁文件建环境 | 环境 | 干净账户上按锁文件建出环境；可重复、只读、可回收 | 2 周 | 是 |
| 2 | L2 | 唯一的执行原语与隔离 | 执行（`runInEnvironment`） | 三项隔离测试在本地和 CI 通过 | 1 周 | 是 |
| 3 | L3 | `skill_run` 和脚本工具，skill 脚本在受管理环境中运行 | Skill、`skill_run`、脚本工具 | 没有本机 Python 科学栈的机器上，官方 skill 能运行 | 2.5 周 | 是 |
| 4 | L4 | agent 会话绑定环境；可视化改写为命令行程序加脚本工具，移出引擎 | Agent 定义、产物 | 没有本机 R 的机器上，可视化能出图；引擎中没有可视化代码 | 3.5 周 | 是 |
| 5 | L3 | 其余使用方：Nextflow、Jupyter、MCP stdio；显式选用本机版本；环境面板 | — | 没有本机 nextflow、conda、jupyter 时 wrapper 和 notebook 可用 | 2 周 | 是 |
| 6 | L5 | 插件结构和加载器；可视化以插件形式安装 | 插件 | 插件的安装、运行、升级、卸载闭环 | 1.5 周 | 否 |
| 7 | L6 | 内容分发：包、安装器、目录、按需启用、远程 registry | 包、启用状态、Wrapper、连接器 | 见内容分发设计 | 分批 | 部分 |
| 8 | 远程 | 远程 / HPC 上的运行时与环境 | — | SSH 项目和离线集群上跑通第 3 步的 skill | 2 周 | 否 |

步骤 1–6 完成后，内核固定（命名规范见[运行时基础设计 §10](../design/phi-runtime-foundation.zh-CN.md)）：skill、agent、插件、wrapper、notebook 使用环境的方式都不再变化；步骤 7、8 只是在内核之上分发内容、扩展到远程主机。

估时为单人粗估。建议 beta 覆盖步骤 0–5。

## 3. 各步骤明细

格式：内容 / 主要涉及 / 验收 / 估时。

### 步骤 0 准备

| ID | 内容 | 主要涉及 | 验收 | 估时 |
|---|---|---|---|---|
| 0.1 | 清理 `resources/` 中的运行残留；新增 `scripts/check-resources.mjs`（发现未被 git 跟踪的文件就报错）并接入 `npm run lint` | `resources/`、`scripts/` | `git status --ignored resources` 干净 | 0.5 天 |
| 0.2 | `electron-builder.yml` 的 `files` 排除运行残留 | `electron-builder.yml` | `build:unpack` 产物里没有 `.nextflow` | 0.5 天 |
| 0.3 | 删除 `resources/skills/create-database-connector` 及其测试引用 | `resources/skills/`、`tests/resources.test.ts`、`tests/phi-agents.test.ts` | `npm test` 通过 | 0.5 天 |
| 0.4 | `scripts/eval/` 入库：README、`bun run eval:fetch` | `scripts/eval/` | 一条命令复现 URL 读取评测 | 0.5 天 |
| 0.5 | **依赖清单**：逐个 skill 和可视化脚本，列出 Python 包、R 包、外部命令（已知：`soffice`、pandoc、poppler、tesseract、Node 的 `docx-js`、R） | `docs/runtime/dependency-inventory.md` | 每个带脚本的内容都有清单，作为步骤 3、4 的环境规格输入 | 1 天 |

### 步骤 1 运行时与环境模型（L0 + L1）

新模块放在 `src/main/agent/envs/`。现有的 `src/main/agent/environment/` 是本机工具探测，到步骤 5 再收缩。

| ID | 内容 | 主要涉及 | 验收 | 估时 |
|---|---|---|---|---|
| 1.1 | **冻结环境契约 v1**：`environment.yml` 支持的字段、各平台 explicit 锁文件格式、`host:` 段、`sourcePackages:` 段（R 源码包，见基础设计 §3.5）、`env.json` 元数据 schema、`envId` 计算规则、状态机 | `docs/contracts/environment.schema.json`、`docs/contracts/env-metadata.schema.json` | schema 与样例通过校验 | 1 天 |
| 1.2 | micromamba 打包：`resources/runtime/manifest.json`（版本、各平台 URL 和 sha256）；`scripts/runtime/fetch-micromamba.mjs` 在构建和开发时下载并校验；`extraResources`；`getMicromambaPath()` | `resources/runtime/`、`scripts/runtime/`、`electron-builder.yml`、`envs/paths.ts` | 开发模式和打包产物中都能执行 `micromamba --version` | 2 天 |
| 1.3 | 运行时目录与调用封装：初始化 `~/.phi/runtime/{envs,pkgs,logs,state}`；生成 `mambarc`（conda-forge、bioconda、strict、镜像、代理取自设置）；`micromamba(args)` 封装固定 `MAMBA_ROOT_PREFIX`、`--rc-file`，并清除 `CONDA_*`、`MAMBA_*` 等变量 | `envs/runtime.ts` | 在装有用户 conda 的机器上，封装调用不读取 `~/.condarc`（测试断言） | 1.5 天 |
| 1.4 | 锁文件工具：`scripts/runtime/lock-env.mjs` 用 `micromamba --platform <p> --dry-run` 为三个平台生成 explicit 锁 | `scripts/runtime/` | 一个只含 python 的样例规格生成三份锁 | 1.5 天 |
| 1.5 | `ensureEnvironment(spec, lock)`：计算 `envId` → 已就绪直接返回 → 否则加文件锁防止并发重复构建 → `micromamba create -p … -f lock` → 安装源码包（1.7）→ 捕获激活快照 → 前缀设为只读 → 写 `env.json` → 更新 `state/environments.json`；构建进度事件和日志 | `envs/ensure.ts`、`envs/activation.ts`、`envs/index-store.ts` | 干净账户上建出样例环境；重复调用不重建；向前缀写文件失败；任一步失败时前缀被删除、状态为 `failed` | 3 天 |
| 1.6 | 引用计数、GC、doctor：记录引用方；删除无引用的环境；对比 `micromamba list --json` 与锁文件，发现漂移时重建 | `envs/gc.ts`、`envs/doctor.ts` | 测试覆盖回收、漂移检测和重建 | 1.5 天 |
| 1.7 | 源码包安装（R）：按 `sourcePackages` 下载 CRAN 精确版本或 GitHub 完整提交 sha 的归档 → 校验 sha256 → 缓存到 `~/.phi/runtime/sources/<sha256>` → 按声明顺序 `R CMD INSTALL`（不自动拉依赖）→ 记录到 `env.json`；只使用环境内的编译器 | `envs/source-packages.ts` | 用一个纯 R 的小包做真实安装；覆盖 sha256 不符、缺少依赖、缓存命中不联网 | 2 天 |

### 步骤 2 执行原语（L2）

| ID | 内容 | 主要涉及 | 验收 | 估时 |
|---|---|---|---|---|
| 2.1 | **冻结执行契约 v1**：净化规则（保留 / 清除 / 设置的变量清单、`PATH` 构成）、`runInEnvironment` 的输入输出（`envId`、解释器路径、退出码、输出截断） | `docs/contracts/execution.md` | 规范评审通过 | 0.5 天 |
| 2.2 | 净化规则模块与单元测试 | `envs/sanitize.ts` | 每一条规则都有测试 | 1 天 |
| 2.3 | `environmentVariables(envRef)`、`runInEnvironment(envRef, argv, opts)`：超时、取消、流式输出；宿主依赖所在目录加入 `PATH` | `envs/run.ts` | 单元测试 | 1.5 天 |
| 2.4 | 金丝雀环境和隔离测试：(a) `sys.executable`、`sys.prefix` 位于环境前缀下；(b) 本机有、环境里没有的包导入失败；(c) 用户 site 被禁用；R 环境可用时再检查 `.libPaths()` | `tests/runtime/`、`npm run test:runtime` | 本地全部通过 | 1.5 天 |
| 2.5 | CI：在 macOS 和 Ubuntu 的干净 runner 上跑 `test:runtime` | CI 配置 | CI 通过；这是进入步骤 3 的门禁 | 1 天 |

### 步骤 3 skill 在环境中运行（L3）

| ID | 内容 | 主要涉及 | 验收 | 估时 |
|---|---|---|---|---|
| 3.1 | **冻结 Skill 契约 v1、`skill_run` 契约 v1、脚本工具声明 v1**：frontmatter 的 `environment`（`phi:python@1` / `./environment.yml` / `plugin:<名>` / `project:default`）、`attachTo`、`scripts`（`name`、`run`、`description`、`args`、`approval`、`output`，路径 `format` 为 `input-path` / `project-path`）；带脚本的 skill 必须声明环境；`skill_run` 的参数、返回（含 `envId`、解释器路径）、错误、审批语义；`validateSkill()` | `docs/contracts/skill.schema.json`、`src/main/agent/content/validate-skill.ts` | 样例通过；`npm run lint` 校验所有 skill | 2.5 天 |
| 3.2 | `phi-python` v1：按 0.5 的清单写规格（含 pandoc、poppler、tesseract、nodejs、ipykernel），生成锁文件，放在 `resources/runtime/environments/phi-python/` | `resources/runtime/environments/` | 三个平台的锁文件都能建出环境 | 2 天 |
| 3.3 | `skill_run` 核心工具：解析顺序（skill 声明 → 会话绑定 → `phi-python` 并提示）；环境未就绪时提示用户构建并显示大小和进度，不回落到本机；纳入 ask 模式审批；注册给主 agent 和专家 | `src/main/agent/content/skill-run.ts` | 单元测试加一次真实运行 | 3 天 |
| 3.4 | 脚本工具注册：按 `scripts` 声明生成 `<toolPrefix>_<name>` 工具；参数校验；`input-path` / `project-path` 路径约束；审批级别；经 `runInEnvironment` 执行；输出 JSON 校验；读取 `<文件>.phi-artifact.json` | `src/main/agent/content/script-tools.ts` | 用一个样例 skill 覆盖参数、路径越界、审批、输出校验、错误退出 | 3 天 |
| 3.5 | 迁移 scanpy：声明环境，SKILL.md 中的脚本调用改为 `skill_run` | `resources/skills/scanpy/` | 在没有本机 Python 科学栈的机器上跑通示例 | 1 天 |
| 3.6 | 迁移其余 Python 类 skill（pptx、xlsx、pdf、markitdown、matplotlib、scikit-learn、scvelo、rdkit）；LibreOffice 声明为宿主依赖 | `resources/skills/` | 每个 skill 通过 `validateSkill()`，并经 `skill_run` 跑通一个示例 | 3 天 |

### 步骤 4 agent 绑定环境；可视化移出引擎（L4）

| ID | 内容 | 主要涉及 | 验收 | 估时 |
|---|---|---|---|---|
| 4.1 | **调研**：在 omp 的 `task` / 注册表上给 Wrapper、Visualization 做原型，决定专家委派是否迁移到 omp、`agents/registry.ts` 中重复的部分是否退役；结论写进决策记录 | `src/main/agent/agents/` | 决策记录更新 | 1 周 |
| 4.2 | **冻结 Agent 定义契约 v1**：沿用 omp 字段（`name`、`description`、`tools`、`spawns`、`model`、`thinkingLevel`），加 Phi 的 `environment`、`visibility`、`skills`、`delegationMode`、`delegation`、`fallback`（旧的 `delegation_mode` 作为别名读取）；结构化结果用 `outputSchema` | `docs/contracts/agent.schema.json`、`agents/definition.ts` | 现有 3 个 agent 通过校验 | 2 天 |
| 4.3 | bash 注入扩展：绑定了环境的会话里，在 `tool_call` 事件中把 `environmentVariables` 并入 bash 调用的 `env` 参数 | `src/main/agent/agents/`、`omp/omp-sdk-worker.ts` | 绑定会话里 `which python` 指向环境；主 agent 的 bash 不受影响 | 2 天 |
| 4.4 | 建立插件形态的目录 `resources/plugins/visualization/`：把 `resources/agents/Visualization.md`、`resources/skills/omics-visualization` 移入；把完整的可视化依赖并入官方 `phi-r` 规格和锁文件（159 个 R 脚本，及为纯标准库脚本提供的 Python 3.12）。没有 conda 构建的 7 个直接 R 包——gground、ggideogram、ggcor、linkET、ggsankey、ggsvg、ggmagnify——加上 ggmagnify 的 gridGeometry 依赖，作为 8 个有顺序、锁定提交的 `sourcePackages`；保留 `phi-r` 已有的 GenomeInfoDbData 源码包，并在此解决 ggideogram 与 ggplot2 4.x 的兼容问题。插件不声明私有环境 | `resources/plugins/visualization/`、`resources/runtime/environments/phi-r/` | 三个平台都能建出 `phi-r`；全部模板冒烟渲染通过；IRkernel、Seurat 和 SingleCellExperiment 仍可加载 | 4 天 |
| 4.5 | **冻结产物契约 v1**：`<文件>.phi-artifact.json` 的字段（类型 `figure` / `table` / `structure` / `molecule` / `network` / `report`、标题、来源）；引擎按产物展示 | `docs/contracts/artifact.schema.json`、展示层 | 样例通过；通用产物查看器可以读取产物 | 2 天 |
| 4.6 | 可视化改写为命令行程序 `scripts/viz.py`（子命令 `examples`、`route`、`prepare`、`render`）：`route` 调用现有的 `route_template.py`；`prepare`、`examples` 从 TS 改写为 Python；`render` 在同一环境中调用 `Rscript` 并运行 QA，输出 `figure` 产物 | `resources/plugins/visualization/skills/omics-visualization/scripts/` | 四个子命令各有测试；输出通过 JSON Schema 校验 | 4 天 |
| 4.7 | 在 SKILL.md 中声明 4 个脚本工具（`toolPrefix: viz`，工具名保持 `viz_examples`、`viz_route`、`viz_prepare`、`viz_render`）；Visualization agent 和 `omics-visualization` skill 都声明 `environment: phi:r@1`，原先按工作流过滤工具的逻辑移入 agent 指令；创建会话时确保环境就绪，否则提示构建 | SKILL.md、`Visualization.md`、会话创建 | 没有本机 R 的机器上出图；可视化评测不退化 | 2 天 |
| 4.8 | 删除引擎中的 `src/main/agent/visualization/` 以及 `visualizationToolNamesForWorkflow` 等专用逻辑 | `src/main/agent/` | 引擎中没有可视化代码；`npm test` 通过 | 1 天 |
| 4.9 | `env_request`：agent 请求额外包 → 用户确认 → 以"原规格 + 额外包"求解新的项目环境并生成锁 → 该项目改为绑定新环境 | `src/main/agent/envs/`、项目状态 | 测试覆盖请求、确认、绑定切换 | 2.5 天 |

步骤 6 之前，`resources/plugins/visualization/` 由一段临时的内置加载逻辑挂载（直接注册其中的 agent 和 skill），步骤 6 换成正式的插件加载器。

### 步骤 5 其余使用方（L3）

| ID | 内容 | 主要涉及 | 验收 | 估时 |
|---|---|---|---|---|
| 5.1 | `phi-nextflow` 环境（nextflow + openjdk + nf-core + nf-test；nf-core 和 nf-test 用于编写和测试 wrapper module，约 380 MB），wrapper 执行器默认用它启动 nextflow；保留显式选用本机 nextflow：沿用 `environment.json` 的 `customPaths`，检查版本是否满足 wrapper 最低要求，满足才可用并标注"本机（不受管理）" | `wrappers/composition/executor.ts`、`environment/store.ts` | 本机没有 nextflow 时 wrapper 可运行；版本不足的本机 nextflow 被拒绝 | 2.5 天 |
| 5.2 | `-profile conda` 配置 `conda.useMicromamba = true`、micromamba 路径、`conda.cacheDir` 放在 `~/.phi/runtime` 下 | wrapper 的 `nextflow.config` 生成逻辑 | 没有本机 conda 时 `npm run smoke:wrappers` 通过 | 2 天 |
| 5.3 | `phi-jupyter` 环境运行 Jupyter 服务端；默认内核来自 `phi-python`（ipykernel）和 `phi-r`（irkernel）；同时列出本机已有的内核，标注"本机（不受管理）"，只能显式选择 | `notebook/analysis-kernels.ts` | 本机没有 jupyter 时 notebook 可用；本机内核可选但不作为默认 | 3 天 |
| 5.4 | MCP stdio 服务的启动走 `environmentVariables`（机制就绪；用户自己在 `mcp.json` 里配置的服务保持原样） | MCP 相关代码 | 单元测试 | 1 天 |
| 5.5 | 环境面板改版：受管理环境（状态、大小、引用方、重建、清理）；宿主依赖（Docker、Singularity、LibreOffice）；可显式选用的本机工具（nextflow、Jupyter 内核） | `features/environment/`、`environment/detect.ts` | UI 可用；旧的本机探测字段迁移 | 2.5 天 |

远程项目中的 wrapper 仍然使用远程主机上的 nextflow，放到步骤 8 处理。

### 步骤 6 插件（L5）

| ID | 内容 | 主要涉及 | 验收 | 估时 |
|---|---|---|---|---|
| 6.1 | **冻结插件契约 v1**：目录结构、`phi-package.yaml`（`type: plugin`、`toolPrefix`）、`environments`、`components`、agent 的 `visibility`、预留 `orchestrator` | `docs/contracts/plugin.schema.json` | 样例通过 | 2 天 |
| 6.2 | 插件加载器：从本地目录安装到 `~/.phi/packages/plugin/<id>/<version>/`；确保环境；注册 agent、skill 和脚本工具；升级时先建新环境再切换；卸载时移除引用并回收环境 | `src/main/agent/plugins/`（新建） | 测试覆盖安装、升级、卸载 | 3 天 |
| 6.3 | 给 `resources/plugins/visualization/` 补上不含 `environments` 的 `phi-package.yaml`，改由插件加载器安装；移除步骤 4 的临时加载逻辑。提升内置插件版本，使升级时移除旧 `plugin:viz` 引用；现有 GC 在该环境无引用后回收它 | `resources/`、`src/main/agent/` | 通过插件加载器安装或升级后工作正常；不再有可视化私有环境被引用 | 1 天 |
| 6.4 | 插件页（列表、从本地安装、卸载）；pi 插件页改名为"开发者扩展"，移到高级设置 | `features/plugin/` | UI 可用 | 2 天 |

### 步骤 7 内容分发（L6）

在内核之上搬运内容，细节见内容分发设计。按以下批次推进：

1. **包契约与本地 registry**：冻结包契约 v1；按 git 跟踪文件白名单打包；安装器（校验、staging、原子安装、`.source.json`、升级、卸载、GC），安装时调用步骤 1 的 `ensureEnvironment`。
2. **skill 目录与按需启用**：冻结启用状态契约 v1；内置 skill 按启用状态过滤；Skills 页"从目录添加"和开关，不再改写 skill 文件；从会话记录迁移。
3. **wrapper 与连接器分包**（beta 之后）：冻结 Wrapper、连接器契约 v1；按工具族打包 wrapper 并拼装目录树；MCP 目录由 registry 提供。
4. **远程 registry**（beta 之后）：签名 index、静态存储加 CDN、更新提示、镜像、离线导入、安装包瘦身。

### 步骤 8 远程 / HPC

- 按需下载 linux 版 micromamba 并上传到远程主机；远端运行时目录与本地布局一致。
- 按锁文件在远端建环境；离线集群使用预先下载的包缓存或 conda-pack；锁文件兼容 glibc 2.17。
- 远程 wrapper 使用远端的 `phi-nextflow` 环境。

## 4. 支线（不改内核，不打乱主线）

| 支线 | 前提 | 内容 |
|---|---|---|
| 数据访问改造 | 独立支线 | 本次发布已移除旧数据库工具链；公共 URL 读取方式与下一版本按来源族拆分的 MCP connector 计划见 [data-access-implementation.md](data-access-implementation.md)。 |
| 多智能体编排 | 步骤 6 | 冻结编排契约；在 omp 上补硬性预算、存储工具、人工检查点和运行视图；Co-Scientist 类参考插件 |

## 5. 契约冻结时间表

| 契约 | 冻结于 | 所属层 |
|---|---|---|
| 环境 | 1.1 | L1 |
| 执行 | 2.1 | L2 |
| Skill / `skill_run` / 脚本工具 | 3.1 | L3 |
| Agent 定义 | 4.2 | L4 |
| 产物 | 4.5 | L4 |
| 插件 | 6.1 | L5 |
| 包、启用状态 | 步骤 7 第 1、2 批 | L6 |
| Wrapper、连接器 | 步骤 7 第 3 批 | L6 |
| 核心获取工具 | 数据访问支线 | 支线 |
| 编排 | 编排支线 | 支线 |

## 6. 第一个迭代

步骤 0 全部，然后 1.1 → 1.2 → 1.3 → 1.4 → 1.5 → 1.6 → 1.7（约 2.5 周）。

完成后：安装包干净，有依赖清单；Phi 拥有自己的运行时，能在干净账户上按锁文件建出只读、可复现、可回收的环境，并且完全不读取用户的 conda 配置。
