# Phi 运行时基础设计：环境、执行、绑定与插件结构

日期：2026-09-29
更新：2026-10-07 — 可视化改用共享的 `phi-r` 环境。
英文版：[phi-runtime-foundation.md](phi-runtime-foundation.md)。
状态：**已确认（2026-09-29，§9 的 6 条决策全部确认）。** 本文先于 [内容分发设计](phi-content-distribution-design.zh-CN.md) 和 [实施计划](../roadmap/content-distribution-implementation.zh-CN.md)：那两份文档里凡是涉及环境、执行、插件结构的内容，都以本文为准。

## 0. 为什么先写这份文档

内容分发（skill、wrapper、连接器、插件的安装）是外层；它依赖几个还没有定下来的底层问题：

1. 内置的 micromamba 放在哪里、怎么调用、和用户自己的 conda 怎么隔离？
2. 环境是什么对象、怎么标识、怎么创建和回收？
3. 一个进程"在环境中运行"具体指什么，怎么保证不串到本机环境？
4. skill、专家 agent、核心工具、notebook、wrapper、MCP 服务，各自怎么使用环境？
5. 谁声明环境、谁绑定环境、冲突时听谁的？
6. 插件的目录结构是什么，插件里的 agent 怎么管理自己的环境？

本文按**由内而外**的顺序回答：L0 运行时布局 → L1 环境模型 → L2 执行原语 → L3 使用方 → L4 绑定规则 → L5 插件结构 → L6 分发。每一层只依赖它内侧的层；**建设也按这个顺序，由小到大**（§8）。

## 1. 现状

| 使用方                  | 现在怎么拿到运行环境                                                                           | 问题                           |
| ----------------------- | ---------------------------------------------------------------------------------------------- | ------------------------------ |
| skill 脚本              | agent 的 `bash` 调本机 `python` / `Rscript`                                                    | 依赖用户本机装了什么；无法复现 |
| 可视化工具 `viz_render` | 在 `PATH` 上找 `Rscript`、`python3`                                                            | 同上；找不到时只能报"缺依赖"   |
| notebook 内核           | 读本机 `jupyter kernelspec`                                                                    | 同上                           |
| wrapper（Nextflow）     | 本机的 `nextflow`；`-profile conda` 用本机 conda                                               | 用户得自己装 nextflow 和 conda |
| 环境面板                | `environment/detect.ts` 探测本机的 micromamba、nextflow、docker、singularity、jupyter、Rscript | 只能"探测"，不能"提供"         |

结论：Phi 目前没有自己的运行时，所有执行都依赖本机。这是所有上层问题的根源。

## 2. L0 运行时布局

### 2.1 micromamba 二进制

- **随 app 打包**，位置：`resources/runtime/micromamba/<platform>-<arch>/micromamba`（`darwin-arm64`、`darwin-x64`、`linux-x64`）。通过 `extraResources` 放在 asar 之外，打包后路径为 `process.resourcesPath/runtime/micromamba/...`；开发模式下使用仓库内同一路径。
- 版本和 sha256 写在 `resources/runtime/manifest.json`。构建脚本从 micromamba 官方发布地址下载并校验，不把二进制提交进 git。
- 只打包当前构建平台对应的那一个。远程 linux 主机需要的 linux 版本，在 P8 阶段按需下载并上传（§7）。
- 引擎用一个函数 `getMicromambaPath()` 定位它；任何代码都不从 `PATH` 上找 micromamba。

### 2.2 与用户 conda 完全隔离

每次调用 micromamba 都显式带上：

- `MAMBA_ROOT_PREFIX=~/.phi/runtime`；
- `--rc-file ~/.phi/runtime/mambarc`，并清除 `CONDARC`、`MAMBARC`、`CONDA_*` 等环境变量，不读取用户的 `~/.condarc`；
- 不执行 `shell init`，不修改用户的 shell 配置，不把任何东西放进用户的 `PATH`。

用户自己的 conda / mamba 不受影响，Phi 也不受用户 conda 配置的影响。

### 2.3 目录布局

```text
~/.phi/runtime/
  mambarc                 # channels（conda-forge、bioconda）、channel_priority: strict、镜像、代理
  pkgs/                   # 包缓存；所有环境从这里硬链接安装
  envs/
    <envId>/              # 环境前缀（prefix），创建后设为只读
      .phi/env.json       # 环境元数据（§3.3）
  state/environments.json # 环境索引：envId → 状态、引用方
  logs/                   # 构建日志
```

`envId = <name>-<hash12>`，`hash` 由"平台 + 锁文件内容"计算（§3.2）。

## 3. L1 环境模型

### 3.1 环境的种类

| 种类     | 例子                                                 | 谁定义                              | 可变性                       |
| -------- | ---------------------------------------------------- | ----------------------------------- | ---------------------------- |
| 基础环境 | `phi-python`、`phi-r`、`phi-nextflow`、`phi-jupyter` | Phi 开发者，随 app 或 registry 发布 | 不可变；新版本是新环境       |
| 包环境   | 插件私有环境；某个 skill 自带的环境                  | 插件或 skill 作者                   | 不可变                       |
| 项目环境 | 某项目额外需要的包                                   | 用户（经 agent 请求、用户确认）     | 不可变；追加依赖会生成新环境 |

`phi-r` 是 R notebook、scanpy R 互操作和内置可视化插件共用的基础
环境。可视化插件不再声明私有包环境。

### 3.2 规格与锁文件

- **规格**：`environment.yml`（name、channels、dependencies；可选 `pip:`）。
- **锁文件**：每个平台一份 **explicit 锁**（`@EXPLICIT` 格式，逐行列出包的 URL 和 md5），由 CI 生成。客户端用 `micromamba create -p <prefix> -f <lock>` 安装，**不在客户端求解**：结果确定，而且只要包缓存里有，就可以离线安装。
- 官方的基础环境和包环境必须带锁文件。用户项目环境没有预先生成的锁，就在本地求解一次，然后把求解结果导出为 explicit 锁保存下来，之后同样按锁安装。
- **宿主依赖**（conda 提供不了的，如 macOS 上的 LibreOffice、Docker、Singularity）在规格的 `host:` 段声明，只检查是否存在、记录路径，不安装。
- **源码包**：conda 上没有构建的 R 包，在规格的 `sourcePackages:` 段声明，建环境时从源码安装，装好后与 conda 部分一起冻结（§3.5）。

### 3.3 环境元数据 `.phi/env.json`

```json
{
  "envId": "phi-r-3f9a1c2b7d10",
  "name": "phi-r",
  "kind": "base",
  "platform": "darwin-arm64",
  "lockSha256": "…",
  "createdAt": "…",
  "micromambaVersion": "2.x",
  "activation": { "set": { "CONDA_PREFIX": "…", "JAVA_HOME": "…" }, "pathPrepend": ["…/bin"] },
  "host": {},
  "sourcePackages": [
    { "language": "r", "name": "ggsankey", "source": "github", "ref": "…", "sha256": "…" }
  ],
  "status": "ready"
}
```

- `activation` 是**激活快照**：建好环境后，用 `micromamba run -p <prefix> env` 捕获一次激活后的变量，与净化后的基础变量做差，得到这个环境需要设置的变量（包括 `activate.d` 脚本设置的，如 `JAVA_HOME`、`GDAL_DATA`、`PROJ_LIB`）。之后每次执行都直接套用快照，不再启动 shell 激活，速度快，而且可以注入到任何进程。
- 状态：`absent` → `building` → `ready`；另有 `failed`、`drifted`（doctor 发现与锁文件不一致）。

### 3.4 生命周期

- **确保存在（ensure）**：给定规格和锁 → 计算 envId → 已就绪就直接返回；否则构建（带进度和日志）：按锁文件安装 conda 部分 → 安装源码包（§3.5）→ 捕获激活快照 → 把前缀设为只读 → 标记 `ready`。任何一步失败，整个前缀删除，状态记为 `failed`，不留下半成品环境。
- **只读**：环境建好后把前缀设为只读。agent 在环境里执行 `pip install` 会直接失败，从机制上保证环境不可变，不依赖提示词约束。
- **引用计数与回收**：`environments.json` 记录每个环境被哪些包、插件、项目引用；没有引用的环境由 GC 删除。包缓存单独按容量清理。
- **修复**：doctor 对比环境和锁文件，发现漂移就重建。

### 3.5 源码包

用于 conda 上没有构建的包，目前只支持 R（来自 CRAN 或 GitHub）。

```yaml
# environment.yml 中 Phi 的扩展段（生成锁文件时会被剥离，不传给 micromamba）
sourcePackages:
  - language: r
    name: ggsankey
    source: github # cran | github
    repo: davidsjoberg/ggsankey # 仅 github
    ref: 5a3b1c… # github 为完整的提交 sha；cran 为精确版本号，如 1.2.3
    sha256: 9f2e… # 源码归档的 sha256
```

规则：

- **完全锁定**：GitHub 只接受完整的提交 sha（不接受分支或标签），CRAN 只接受精确版本（从 CRAN 归档地址下载）；每个归档都必须带 sha256。`sourcePackages` 整段计入 envId 的 hash，所以源码包的任何变化都会生成新环境。
- **不自动拉依赖**：安装时使用 `R CMD INSTALL`（等价于 `install.packages(..., repos = NULL, dependencies = FALSE)`），源码包的依赖必须已经由 conda 部分提供，或者在同一段中排在它前面。缺少依赖会直接报错，不会从网络上补装，保证环境内容完全由规格决定。
- **编译**：需要编译本地代码的包，要求 conda 部分包含相应的编译工具（如 `r-base` 配套的 `compilers`）；构建时只使用环境内的编译器，不使用本机的编译器。
- **缓存与离线**：下载的归档缓存在 `~/.phi/runtime/sources/<sha256>`，校验通过后才会使用；缓存里有就不联网。以后远程 registry 上线后，这些归档可以和包一起提供。
- **冻结**：源码包安装完成后才捕获激活快照、把前缀设为只读，所以之后不能再修改。
- **记录**：已安装的源码包（名称、来源、ref、sha256）写入 `env.json`，doctor 核对时一并检查。

## 4. L2 执行原语

引擎只提供**一个**"在环境中运行"的实现，所有使用方都调用它：

```ts
environmentVariables(envRef, options): Record<string, string>   // 供需要自己启动进程的使用方
runInEnvironment(envRef, argv, { cwd, stdin, timeout, signal }) // 直接运行
```

### 4.1 净化后的基础变量

子进程的环境变量**不是继承 Phi 进程的全部变量**，而是：

- **保留**：`HOME`、`USER`、`LOGNAME`、`TMPDIR`、`TERM`、`LANG`、`LC_*`（缺省时补 UTF-8，复用 `render.ts` 里的 `utf8LocaleEnv`）、`HTTP(S)_PROXY`、`NO_PROXY`、`SSH_AUTH_SOCK`；
- **清除**：`PYTHONPATH`、`PYTHONHOME`、`VIRTUAL_ENV`、`CONDA_*`、`MAMBA_*`、`R_LIBS`、`R_LIBS_USER`、`R_LIBS_SITE`、`R_HOME`、`JAVA_HOME`、`PERL5LIB`、`LD_LIBRARY_PATH`、`DYLD_*`、`NODE_PATH`；
- **设置**：`PYTHONNOUSERSITE=1`、`R_LIBS_USER=`（空）、`MPLBACKEND=Agg`；
- **PATH**：`<环境>/bin` + 激活快照里的其他路径 + 已声明宿主依赖所在的目录 + 系统最小路径（`/usr/bin:/bin:/usr/sbin:/sbin`）。**不包含**用户的 conda、Homebrew、pyenv、`~/.local/bin`。

### 4.2 可验证

- 每次执行在结果里返回 `envId` 和实际解释器路径，并写入运行记录（可复现）。
- 一致性测试（§8 第 2 步）证明：解释器位于环境前缀下；本机装了、但环境里没有的包导入失败；在没有任何系统科学计算栈的干净机器上也能运行。

## 5. L3 使用方

所有使用方都通过 L2 执行，只有两类例外：

- **主 agent 的 bash**：它代表用户在自己的机器上工作，保持使用本机环境，不做绑定。
- **面向用户交互、且用户显式选择本机版本的工具**：notebook 内核和 Nextflow 默认使用受管理环境，但用户可以显式选用本机的内核或 nextflow（§5.2）。

**内容执行**（skill 脚本、脚本工具、插件和专家 agent）**始终**使用受管理环境，不复用本机，即使本机已经装有同名工具。

| 使用方                      | 现在                                           | 以后                                                                                                                                                                                                                                            |
| --------------------------- | ---------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| skill 脚本                  | agent 用 `bash` 调本机解释器                   | 声明为脚本工具（§5.1）或经 `skill_run(skill, script, args)` 调用 → `runInEnvironment(skill 的环境)`                                                                                                                                             |
| 专家 / 插件 agent 的 `bash` | 本机 shell                                     | 会话绑定环境：Phi 扩展在 `tool_call` 事件中改写 bash 调用，把 `environmentVariables(agent 的环境)` 并入 bash 的 `env` 参数（omp 支持返回替换后的输入）。agent 写 `python x.py` 就会用环境里的 python                                            |
| 领域工具（如 `viz_render`） | 写在引擎里的 TS 工具，在 `PATH` 上找 `Rscript` | 移出引擎，改为插件里的命令行程序，以脚本工具（§5.1）的形式注册；工具名不变                                                                                                                                                                      |
| MCP stdio 服务              | —                                              | 用 `environmentVariables` 启动                                                                                                                                                                                                                  |
| notebook 内核               | 本机 kernelspec                                | Jupyter 服务端在 `phi-jupyter` 中运行；默认内核来自受管理的分析环境（`phi-python` 里的 `ipykernel`、`phi-r` 里的 `irkernel`）；本机已有的内核同时列出，供用户显式选择（§5.2）                                                                   |
| wrapper（Nextflow）         | 本机 `nextflow`，本机 conda                    | 默认用 `phi-nextflow` 环境中的 `nextflow` 和 Java；`-profile conda` 设置 `conda.useMicromamba = true`，指向内置 micromamba，`conda.cacheDir` 放在 `~/.phi/runtime` 下；用户可以显式指定本机 nextflow（§5.2）；Docker / Singularity 仍是宿主依赖 |
| 环境面板                    | 探测本机工具                                   | 展示受管理的环境（状态、大小、引用方、宿主依赖）；本机探测只保留宿主依赖（Docker、Singularity、LibreOffice）和用户可显式选用的工具（nextflow、Jupyter 内核）                                                                                    |

### 5.1 脚本工具

引擎里不再有领域工具。skill 和插件通过在 SKILL.md 中**声明脚本工具**，得到带类型参数的工具；引擎按声明自动注册，不需要为任何领域写代码。

```yaml
# SKILL.md frontmatter（片段）。Phi 的扩展字段都放在 phi 块里，见 docs/contracts/skill.md
name: omics-visualization
description: …
phi:
  environment: phi:r@1
  scripts:
    - name: route
      run: [python, ./scripts/viz.py, route]
      description: 分析结果表，列出适合的图表模板候选
      args: # JSON Schema
        type: object
        required: [data, purpose]
        additionalProperties: false
        properties:
          data: { type: string, format: input-path } # 只读，必须位于项目内
          purpose: { type: string }
      approval: read # read | write
      output: ./schemas/route-result.json # JSON Schema
    - name: render
      run: [python, ./scripts/viz.py, render]
      args:
        type: object
        required: [script, output]
        additionalProperties: false
        properties:
          script: { type: string, format: input-path }
          output: { type: string, format: project-path } # 可写，必须位于项目内
      approval: write
```

引擎的职责（全部通用）：

- **注册**：工具名为 `<toolPrefix>_<name>`（§10），注册给声明了该 skill 的 agent；skill 通过 `attachTo` 指向主 agent 时，注册给主 agent。
- **参数**：按 `args` 的 JSON Schema 生成工具参数并校验。
- **路径**：`input-path` 只读、`project-path` 可写，都必须解析到项目目录内，否则拒绝执行。
- **审批**：按 `approval` 纳入现有审批机制。
- **执行**：经 `runInEnvironment` 在 skill 声明的环境中运行 `run` 指定的命令，参数以 `--<名> <值>` 传入。
- **输出**：程序向标准输出写一个 JSON 对象，引擎按 `output` 校验；失败时退出码非零，并输出 `{"error": "..."}`。
- **产物**：程序生成 `<文件名>.phi-artifact.json` 描述产物（类型、标题、来源），引擎据此在界面中展示。

agent 仍然可以在绑定了环境的会话里用 `bash` 直接运行这些程序；脚本工具只是更可靠的调用方式：参数有 schema，路径有约束，审批有级别。

### 5.2 显式选用本机版本

- **notebook 内核**：内核列表同时列出受管理的内核和本机已有的内核（读取本机 kernelspec）。本机内核标注为"本机（不受管理）"，需要用户主动选择，不作为默认。
- **Nextflow**：沿用现有 `~/.phi/environment.json` 中的 `customPaths` 机制，用户在环境面板里指定本机 nextflow 路径。指定时检查版本是否满足 wrapper 的最低要求，不满足则拒绝；满足则标注为"本机（不受管理）"。这主要服务于 HPC 上使用站点自带 nextflow 的情况。
- 本机版本不提供可复现保证；运行记录中写明使用的是本机版本及其路径。

## 6. L4 绑定与解析规则

### 6.1 谁可以声明环境

| 声明位置                      | 写法                                                                                            | 含义                                               |
| ----------------------------- | ----------------------------------------------------------------------------------------------- | -------------------------------------------------- |
| skill frontmatter 的 `phi` 块 | `environment: phi:python@1` / `environment: ./environment.yml` / `environment: plugin:statistics` | 这个 skill 的脚本和脚本工具在哪个环境里跑          |
| agent frontmatter             | `environment: phi:r@1` / `environment: plugin:statistics`                                      | 这个 agent 会话的 bash 和核心工具用哪个环境        |
| 插件 manifest                 | `environments: { statistics: {...} }`                                                          | 插件自带的具名环境，供插件内的 skill 和 agent 引用 |
| 项目                          | `.phi/environment.yml`（可选），引用写作 `project:default`                                      | 项目需要的额外依赖                                 |

`plugin:<name>` 只能引用同一插件内的环境；跨插件引用不允许。

内置可视化插件在 skill 和 agent frontmatter 中都使用官方环境
`phi:r@1`，因此它的 manifest 不写可选的 `environments` 字段。

### 6.2 解析顺序

- **`skill_run(skill, …)`**：skill 自己声明的环境 → 所在 agent 会话绑定的环境 → `phi-python`（同时提示"未声明依赖"）。
- **agent 会话**：agent 声明的环境 → 不绑定（与主 agent 一样用本机）。
- 一个会话里同时存在两个环境时不混用：agent 的 bash 用 agent 的环境；它调用的 skill 若声明了别的环境，经 `skill_run` 在那个环境里单独运行。

### 6.3 需要额外的包怎么办

环境只读，agent 不能往里装包。需要额外依赖时：

1. agent 调用 `env_request(packages, reason)`（预留的核心工具）；
2. 用户在界面上确认；
3. 引擎以"原环境规格 + 额外包"求解出一个新的**项目环境**，导出锁文件，记录到项目状态；
4. 该项目中这个 agent / skill 此后绑定新环境。原环境和其他项目不受影响。

## 7. L5 插件结构

### 7.1 目录

```text
<plugin-id>/
  phi-package.yaml
  agents/
    Visualization.md            # frontmatter: visibility, environment, skills, tools, spawns
  skills/
    omics-visualization/
      SKILL.md                  # frontmatter: environment: phi:r@1, scripts: [...]
      scripts/  references/  assets/
  mcp/                          # 可选：stdio 服务定义
  wrappers/                     # 可选
  orchestrator/                 # 预留
```

### 7.2 可视化插件 manifest

```yaml
id: visualization
type: plugin
version: 1.0.2
toolPrefix: viz # 脚本工具名前缀，见 §10
components:
  agents: [agents/Visualization.md]
  skills: [skills/omics-visualization]
```

`可视化` manifest 不含 `environments` 字段。它的 agent 和 skill
frontmatter 都绑定 `phi:r@1`。

### 7.3 可视化插件如何使用共享环境

1. **安装插件时**：引擎注册 agent、skill 和脚本工具；没有插件私有环境需要构建。
2. **创建 agent 会话时**：读取 `environment: phi:r@1` → ensure 官方环境就绪（没就绪就提示用户构建，而不是静默回落到本机）→ 挂上 bash 注入扩展 → 核心工具通过运行上下文拿到同一个环境 → 插件内的 skill 经 `skill_run` 解析到同一个环境。
3. **运行中**：环境只读；需要额外包走 §6.3。
4. **移除旧私有环境的升级时**：活动版本原子切换后，旧版本的
   `plugin:viz` 引用被移除；现有 GC 会在它无其他引用后回收该环境。
5. **卸载插件时**：删除插件组件；插件不拥有共享的官方 `phi-r` 环境。

agent 本身不"管理"环境，它只**声明**要用哪个环境；创建、激活、隔离、回收都由引擎负责。这样插件作者不写任何环境管理代码，引擎也不需要为某个插件做特殊处理。

## 8. 建设顺序（由内而外，由小到大）

每一步都以上一步完成并通过测试为前提。

| 步骤 | 层      | 内容                                                                                                                                                                           | 完成标准                                                                 |
| ---- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------ |
| 1    | L0 + L1 | 打包 micromamba；`~/.phi/runtime` 布局与 `mambarc`；按 explicit 锁创建环境、激活快照、只读、索引、GC                                                                           | 在干净账户上，用一个只含 python 的锁文件建出环境；重复 ensure 不重建     |
| 2    | L2      | `environmentVariables` / `runInEnvironment`；净化规则；金丝雀测试（解释器位置、反向导入、干净机器）                                                                            | 三项隔离测试全部通过                                                     |
| 3    | L3      | `skill_run`；`phi-python` 首版；迁移一个 skill（如 scanpy）走 `skill_run`                                                                                                      | 该 skill 在没有本机 Python 科学栈的机器上可用                            |
| 4    | L4      | agent frontmatter 的 `environment`；bash 注入扩展；可视化改写为命令行程序 `scripts/viz.py` 并声明为脚本工具，在 `phi-r` 环境中运行；删除引擎中的 `src/main/agent/visualization/` | 可视化在没有本机 R 的机器上出图；引擎中没有可视化代码                    |
| 5    | L3      | 其余使用方：notebook 内核（`phi-jupyter`）、Nextflow（`phi-nextflow` + `conda.useMicromamba`）、MCP stdio；显式选用本机版本（§5.2）                                            | 在没有本机 jupyter / nextflow / conda 的机器上，notebook 和 wrapper 可用 |
| 6    | L5      | 插件目录结构和 manifest；把已经改写好的可视化在本地打包成插件，按 §7.3 安装、运行、卸载                                                                                        | 插件的完整生命周期跑通                                                   |
| 7    | L6      | 内容分发：单元、目录、安装器、registry（见内容分发设计）                                                                                                                       | 按内容分发设计的标准                                                     |
| 8    | 远程    | 远程 / HPC：上传 linux micromamba，按锁文件在远端建环境；离线集群用包缓存或 conda-pack                                                                                         | 在 SSH 项目和离线集群上跑通步骤 3 的 skill                               |

步骤 1–6 就是"内核"：做完之后，skill、agent、插件、wrapper、notebook 使用环境的方式全部固定下来，外层的分发只是在这之上搬运文件。

## 9. 已确认的基础决策（2026-09-29）

1. **只有主 agent 的 bash 使用本机环境**，其余执行一律走受管理环境。
2. **官方环境用 explicit 锁文件，客户端不求解**；只有用户项目环境会在本地求解一次。
3. **环境建好后前缀只读**，额外依赖通过 `env_request` 生成新的项目环境。
4. **LibreOffice、Docker、Singularity 作为宿主依赖**，只检查、不安装。
5. **Nextflow 和 Jupyter 默认使用受管理环境**；notebook 内核和 Nextflow 允许用户显式选用本机版本（Nextflow 需通过版本检查），标注为不受管理。内容执行（skill、脚本工具、插件和专家 agent）始终使用受管理环境。（2026-09-29 修订）
6. micromamba 随 app 打包，**只打包当前平台**，远程所需的版本按需下载。
7. **引擎中不再有领域工具**：领域逻辑写成插件或 skill 里的命令行程序，通过在 SKILL.md 中声明脚本工具获得带类型参数的工具（§5.1）。`viz_*` 是第一个迁移对象，工具名保持不变。
8. **统一命名规范**见 §10。

## 10. 命名规范

| 对象                         | 规则                                                                                                                                                                                                                                   | 例子                                                                                         |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| 核心工具（引擎）             | snake_case，`<领域>_<动词>`                                                                                                                                                                                                            | `skill_run`、`env_request`、`http_fetch`、`wrapper_search`、`agent_status`                   |
| 引擎保留的工具前缀           | 内容不得使用                                                                                                                                                                                                                           | `skill_`、`env_`、`http_`、`wrapper_`、`agent_`、`db_`（已退役但仍保留）、`mcp__`（omp 的 MCP 工具） |
| 脚本工具                     | `<toolPrefix>_<脚本名>`；`toolPrefix` 在包 manifest 中声明，2–12 位小写字母或数字，registry 内唯一，不得与保留前缀冲突                                                                                                                 | `viz_route`、`viz_render`                                                                    |
| 包 id、skill 名、插件 id     | kebab-case                                                                                                                                                                                                                             | `omics-visualization`、`visualization`、`protein-apis`                                       |
| agent 名                     | PascalCase                                                                                                                                                                                                                             | `Visualization`、`Database`、`Wrapper`                                                       |
| 环境名                       | kebab-case；Phi 维护的环境加 `phi-` 前缀，按用途命名；插件内环境不加前缀                                                                                                                                                               | `phi-python`、`phi-r`、`phi-nextflow`、`phi-jupyter`；`statistics`                           |
| 环境引用                     | `<作用域>:<名>[@<主版本>]`，或相对路径                                                                                                                                                                                                 | `phi:python@1`、`phi:r@1`、`plugin:statistics`、`project:default`、`./environment.yml`       |
| envId                        | `<作用域>-<所有者>-<名>-<hash12>`；`phi` 作用域省略所有者                                                                                                                                                                              | `phi-python-3f9a1c2b7d10`、`plugin-reports-statistics-…`、`project-<项目短 id>-default-…` |
| 包清单                       | 所有类型统一为 `phi-package.yaml`，用 `type` 区分                                                                                                                                                                                      | `type: skill` / `wrapper` / `mcp` / `plugin`                                                 |
| frontmatter 和 manifest 字段 | camelCase，与 omp 一致；例外：SKILL.md 的标准字段沿用 Agent Skills 规范（如 `allowed-tools`、`disable-model-invocation`），Phi 的字段一律放在 `phi` 块里；旧的 snake_case 字段（如 `delegation_mode`）作为兼容别名读取，校验器提示迁移 | `thinkingLevel`、`outputSchema`、`attachTo`、`toolPrefix`、`delegationMode`                  |
| 脚本入口                     | `scripts/<toolPrefix>.py`，子命令与脚本工具名一致                                                                                                                                                                                      | `scripts/viz.py route`                                                                       |
| 产物描述文件                 | `<文件名>.phi-artifact.json`                                                                                                                                                                                                           | `fig.png.phi-artifact.json`                                                                  |
| 契约文件                     | `docs/contracts/<名>.schema.json`（kebab-case），文字规范为 `<名>.md`                                                                                                                                                                  | `environment.schema.json`、`env-metadata.schema.json`、`execution.md`                        |
| 用户目录                     | `~/.phi/runtime/`（micromamba 根、环境、包缓存）、`~/.phi/packages/<type>/<id>/<version>/`（所有包，含插件）、`~/.phi/state/`（启用状态等）                                                                                            |                                                                                              |
| 仓库目录                     | `resources/runtime/environments/<环境名>/`、`resources/plugins/<插件 id>/`                                                                                                                                                             | `resources/runtime/environments/phi-python/`                                                 |

实施计划按 §8 编排：[内容分发实施计划](../roadmap/content-distribution-implementation.zh-CN.md)。内容分发设计中环境（§9）和插件结构（§11）两章以本文为准。
