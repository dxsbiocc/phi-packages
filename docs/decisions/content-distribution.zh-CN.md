# Phi 内容分发决策

日期：2026-09-29
状态：方向已确定。待定问题见设计文档 §16。
设计文档：[phi-content-distribution-design.zh-CN.md](../design/phi-content-distribution-design.zh-CN.md)（英文版：[phi-content-distribution-design.md](../design/phi-content-distribution-design.md)）。
本记录的英文版：[content-distribution.md](content-distribution.md)，两者如有出入以英文版为准。

## 背景

Phi 把所有领域能力都打进安装包（`resources/skills`、`resources/wrappers`、`resources/db-connectors`、`resources/agents`），并对每个用户全量加载。结果是：不发布新版 app 就无法更新内容；每个 skill 的描述都在主 agent 的提示词里；脚本依赖取决于用户 `PATH` 上碰巧有什么。与此同时，领域逻辑也渐渐长进了引擎（`src/main/agent/db/`、`src/main/agent/visualization/`）。2026-09-29 的评测显示，`db_*` 工具链的准确率并不比直接读取 URL 高，耗时却约为两倍（设计文档附录 A）。

## 决策

1. **引擎与内容分离。** 引擎（内核）包含 agent 运行时、通用核心工具、安装器、环境管理、审批、远程执行，以及 `create-wrapper` 这类内置 skill，不包含任何领域逻辑。所有领域相关的东西都是内容。
2. **先定内核契约，再冻结。** 包、Skill、Wrapper、连接器、Agent 定义、插件、环境、核心服务、产物，以及（预留的）编排，都写成成文契约，配 JSON schema、统一的校验器（`phi validate`）和一致性测试。每个契约先由参考实现验证，并在它的第一个使用者动工之前冻结为 v1。此后次版本只允许增量修改；破坏性修改需要 ADR，并留出至少两个 app 版本的弃用窗口。
3. **单元与插件。** 单一工具是独立单元：skill、wrapper 工具族或 MCP 连接器。只有复合能力才做成插件，即包含多种互相依赖的组件类型、需要专属环境或需要专属专家 agent。单元和插件共用一种包格式、一个安装器和一个 registry。
4. **先本地 registry，后远程。** 在安装器、加载器和目录界面稳定之前，registry 由 `resources/` 生成；带签名、托管在静态存储上的远程 registry 放到 beta 之后。
5. **白名单打包。** 包只能由 git 跟踪的文件或 manifest 构建，不能整目录复制；CI 设体积预算。
6. **按类型路由。** 主 agent 只看到专家入口和少量知识型 skill；wrapper 归 `wrapper_agent`，API skill 归对应的专家。只有已安装且已启用的内容才进入运行时。
7. **受管理的运行环境。** Phi 内置 micromamba。官方的脚本型 skill 共用一个由开发者维护、锁定、经 CI 检查的 `phi-python`；用户 skill 通过 dry-run 选择已有环境，或新建按内容寻址、不可变的环境。脚本通过 `skill_run` 运行。
8. **数据访问改用 API skill。** 用核心获取工具（重试、限速、白名单与审计、分页和批量下载、结果查看器）加约 6 个领域 API skill 取代 `db_*` 工具链。MCP 用于需要认证或私有的数据源，以及已有的官方服务。只有当 `api-skill` 组在扩展评测中全面不差于旧工具链时，旧工具链才退役。`create-database-connector` skill 立即删除。
9. **引擎中不再有领域工具。** 领域逻辑写成插件或 skill 里的命令行程序，通过在 SKILL.md 中声明脚本工具获得带类型参数的工具；引擎只负责注册、校验、在对应环境中运行，以及展示类型化产物。可视化第一个迁移，改写为 `scripts/viz.py`，工具名 `viz_examples`、`viz_route`、`viz_prepare`、`viz_render` 保持不变。
10. **多智能体编排基于 omp。** 编排在 v1 中预留扩展点（插件的 `orchestrator` 槽位、agent 的 `visibility`、`orchestration.*` 命名空间）。它通过一层很薄的 Phi 适配器使用 omp 的 `task`、`spawns`、agent 注册表和 `hub`；Phi 补上硬性预算、类型化的运行存储、审批、人工检查点、按包隔离和运行视图。先提供 agent 驱动的编排，声明式和程序式层级以后再加。
11. **pi 插件属于开发者扩展。** 现有的 pi 插件页改名并移到高级设置；界面上的"插件"只指 Phi 插件。
12. **运行时基础先行**（[运行时基础设计](../design/phi-runtime-foundation.zh-CN.md)，2026-09-29 确认）：只有主 agent 的 `bash` 使用本机环境，其余执行都在受管理环境中运行。官方环境按 explicit 锁文件安装，客户端不求解。环境前缀只读，额外依赖生成新的项目环境。LibreOffice、Docker、Singularity 是宿主依赖，只检查不安装。Nextflow 和 Jupyter 默认使用受管理环境；notebook 内核和 Nextflow 允许显式选用本机版本（Nextflow 需通过版本检查），标注为不受管理；内容执行始终使用受管理环境。只打包当前平台的 micromamba，远程主机按需获取。建设顺序由内而外：运行时 → 执行 → 使用方 → 绑定 → 插件 → 分发。
13. **统一命名规范**，定义在运行时基础设计 §10：核心工具为 `<领域>_<动词>`，并有保留前缀；脚本工具为 `<toolPrefix>_<名>`；环境名为 `phi-<用途>`，引用写作 `phi:<名>@<主版本>`、`plugin:<名>`、`project:default`；所有类型的包统一使用 `phi-package.yaml`；字段用 camelCase，与 omp 一致；所有包安装在 `~/.phi/packages/<type>/<id>/<version>/`。

14. **专家委派继续使用 Phi 自己的运行时；omp `task` 用于多智能体编排**（步骤 4.1 调研，2026-09-30，omp 18.1.10）。调研从一个按 Phi 方式构建的父会话发起真实的 omp `task` 子 agent，结论如下：
    - 父会话在进程内传入的 `customTools`（Phi 的 `viz_*`、wrapper、数据库和脚本工具都是这样构建的，并各自持有 `requestHost` 闭包）子 agent **看不到**（"tool not found"）。
    - 通过**扩展**注册的工具（`pi.registerTool`）子 agent **可以**使用，并且在父进程中执行，闭包照常工作；父会话扩展的 `tool_call` 钩子也能看到子 agent 的每次调用。前提是父会话没有设置 `restrictToolNames`，因为只有这时 omp 才转交已加载的扩展。
    - 子 agent 以 `tools.approvalMode: "yolo"` 运行。Phi 的询问模式审批只能借助被转交的扩展钩子继续生效，而且钩子拿不到 Phi 的运行上下文（`agentRunId`、对话里的卡片）。
    - agent 只从 `.omp/agents` 目录和 omp 插件中发现，不读取 Phi 扫描到的定义。
    - omp 没有按 agent 绑定环境的机制。

    因此，把 Phi 的三个专家迁到 `task` 上，意味着所有 Phi 工具都要改为通过扩展注册并按 agent 限定范围，还要重做 `agents/registry.ts` 和界面已经提供的能力：每个对话的并发与次数上限、带对话卡片和后台任务面板的后台运行、中途引导（steer）、汇报协议、回退策略、用量日志、远程项目防护，以及绑定到父会话的审批。这相当于重写，内测风险大，用户几乎感受不到收益。

    决定：内测期间保留 Phi 的委派运行时，`agents/registry.ts` 暂不退役任何部分。Agent 定义契约（4.2）采用 omp 的字段名（`name`、`description`、`tools`、`spawns`、`model`、`thinkingLevel`），保证定义仍能被 omp 加载。多智能体编排（决定 10）使用 omp `task` / `hub`，Phi 工具通过注册工具的扩展提供，编排会话不限制工具名；编排适配层必须为审批补回运行上下文。

## 影响

- 内部 beta roadmap 随之调整：skills、wrappers、连接器改为可从目录安装、启用、停用，取代原来"Skills 和 MCP 页面只读"的规定；DB 连接器原型冻结，等待退役。
- 工作按分层计划推进（[实施计划](../roadmap/content-distribution-implementation.zh-CN.md)），顺序由运行时基础设计推出：运行时 → 执行 → 使用方 → 绑定 → 插件 → 分发 → 远程。步骤 0–5 建议纳入 beta；插件、内容分发、远程 / HPC（步骤 6–8）随后。数据访问和编排是支线。每个契约在所属层完成时冻结。
- omp 调研（步骤 4.1）已完成：专家委派继续使用 Phi 自己的运行时（决定 14），保留 `agents/registry.ts`。
- 内容作者只需面对一个校验器和一种包格式；引擎的修改由一致性测试把关。

## 考虑过的其他方案

- **全部继续打包在 app 里**：最简单，但内容和发版绑死，每加一样东西路由负载就增加一分。
- **一切都做成插件**：单一工具和用户内容也被迫打包；否决，改为"单元加复合插件"。
- **把 db 连接器改写成 MCP 服务**：MCP 只是传输协议；评测表明开销来自查询抽象本身，所以 MCP 只保留给需要认证或协议复杂的数据源。
- **自研编排运行时**：与 omp 已有的能力重复；否决，改为在 omp 之上做适配。
- **所有内容共用一个环境**：无法满足互相冲突的用户 skill；否决，改为 `phi-python` 加按内容寻址的环境。
