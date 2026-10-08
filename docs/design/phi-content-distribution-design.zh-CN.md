# Phi 内容分发：可安装单元、插件与运行环境

日期：2026-09-29
状态：**设计草案，实施前需负责人确认。** 阶段 0–5 建议纳入内部 beta，阶段 6–9 放到 beta 之后（见 §15）。

英文版：[phi-content-distribution-design.md](phi-content-distribution-design.md)。两份内容保持一致，如有出入以英文版为准。
决策记录：[content-distribution.zh-CN.md](../decisions/content-distribution.zh-CN.md)。

## 1. 概述

Phi 的领域能力（skills、wrappers、数据库连接器、agent 定义、可视化模板）目前都放在安装包的 `resources/` 下，并且对每个用户全量加载。本设计把它们变成**按需分发、按需加载的内容**：

- app 只保留一个小的**引擎**（内核）：agent 运行时、通用核心工具、安装器、环境管理、审批、远程执行，以及少量内置 skill。引擎里不包含任何领域逻辑。
- 引擎为插件、skill、wrapper、连接器、agent 定义、运行环境和核心服务提供一组固定的、**带版本号的契约**（§4.4）。契约先用参考实现验证，然后冻结；此后新能力以内容的形式加入，而不是修改引擎。
- 所有领域相关的内容都变成可安装的**单元**（skill、wrapper、MCP 连接器）；只有复合能力才做成**插件**。
- 单元和插件共用一种包格式、一个安装器、一个 registry。registry 先在本地（由 `resources/` 生成），之后换成带签名的远程 registry，这样内容更新不需要重新发布 app。
- 带脚本的内容在**受管理、有锁文件的 micromamba 环境**中运行，环境由 Phi 自动选择。
- `db_*` 连接器工具链被 **API skill 加一个核心获取工具**取代，依据是附录 A 的评测。

目标按优先级排列：(1) 内容不锁死在安装包里；(2) 主 agent 的路由负载随用户安装的内容增长，而不是随 Phi 知道的全部内容增长；(3) 脚本可复现地运行，用户不需要自己修环境。

## 2. 不在范围内

- 带账号、评分、付费或第三方发布的公开市场。第一版远程 registry 只发布第一方内容。
- 从内容包加载任意可执行的 TypeScript。内容只能是声明式文件，加上可以沙箱化的脚本和 MCP 服务（§11.2）。
- Windows 上的受管理环境（bioconda 在 Windows 上基本不可用）。
- 改写 Nextflow wrapper 的执行方式。Nextflow 继续按 process 管理环境（`conda` / `docker` / `singularity` profile）。

## 3. 现状

| 内容 | 位置 | 如何进入 agent | 问题 |
|---|---|---|---|
| 21 个 skill | `resources/skills/` | 整个目录被追加到运行时的 skill 路径（`appendExistingBundledSkillPaths`） | 每个 skill 的描述都在主提示词里；安装包因此多了 65M（其中 57M 是 `omics-visualization`） |
| 约 150 个 module / subworkflow wrapper | `resources/wrappers/` | `wrapper_agent` + `wrapper_search`；`~/.phi/wrappers/packs/<version>/` 下的 overlay pack 整树替换 | 路由负载低，但只有一个整体 pack；安装时整目录复制，会带上开发残留 |
| 32 个数据库连接器 | `resources/db-connectors/` | Database agent 使用 7 个 `db_*` 工具，背后是 `src/main/agent/db/` 下约 2 万行代码 | 比直接 fetch 更慢、更不可靠（附录 A） |
| 3 个 agent 定义 | `resources/agents/` | 创建会话时扫描 | Visualization 与内置工具和 skill 强耦合 |
| pi 插件 | 通过运行时从 npm / git 安装 | extension、prompt、theme | 属于开发者扩展，不是科研能力 |
| 远程 MCP 连接器 | 通过目录对话框写入 `~/.phi/mcp.json` | MCP 工具（可延迟加载） | 已经是本设计想推广的"需要什么加什么"模式 |

脚本依赖目前没有管理：`viz_render` 直接调用 `PATH` 上的 `Rscript` / `python3`，skill 脚本则用 agent 的 `bash` 碰巧找到的 Python。

## 4. 概念

### 4.1 引擎、单元、插件

- **引擎**：随 app 发布、跟 app 一起升级版本的代码。包括通用核心工具（`wrapper_*`、核心获取工具、`skill_run`、产物展示、文件和 shell 工具）、安装器、环境管理，以及扩展 Phi 自身的**内置 skill**（`create-wrapper`）。内置 skill 不是包，也不出现在任何目录里。
- **单元**：一种类型的一份内容，即一个 skill、一个 wrapper 工具族或一个 MCP 连接器。agent 运行时实际使用的是单元。
- **插件**：把几种不同类型的单元打包在一起的包，这些单元必须一起安装、升级和卸载；插件可以带自己的专家 agent 和专属环境。

插件回答"能力从哪里来"，单元回答"agent 能用什么"。两层并存：官方和共享内容通过分发获得，用户和项目内容直接以单元形式创作。

### 4.2 做成单元还是插件？

满足下列任一条，才做成插件：

1. 包含两种及以上互相依赖的组件类型（agent + skill + 模板，或 skill + wrapper）；
2. 需要与 `phi-python` 冲突的专属环境（§9）；
3. 需要专属的专家 agent（有独立的流程和指令）。

其余情况都是单元。单一工具永远不做成插件。仅仅带脚本，也不足以让一个 skill 变成插件（§9.4）。按现有内容，只有 `omics-visualization`（连同 Visualization agent 和它的模板）符合条件。

### 4.3 运行时内容的归属

内容按类型路由，不论它来自插件还是单独添加：

| 内容 | 归属 |
|---|---|
| 专家 agent | 主 agent，作为一个委派入口 |
| wrapper | `wrapper_agent` 的检索池 |
| API / 知识型 skill | 它声明的专家（例如 Database），没有声明时给主 agent |
| 需要专属环境的脚本型 skill | 所属插件的 agent（环境绑定在该会话上） |
| MCP 连接器 | 所属插件的 agent，没有时给主 agent（工具延迟加载） |

由此带来的结果：主 agent 只看到专家入口和少量知识型 skill，提示词随**已安装的专家数**增长，而不是随工具数增长。在编排运行之外，专家之间不互相调用，只由主 agent 编排；专家之间通过项目文件加简短说明交接。在插件的编排运行内部（§11.4），由插件自己的编排器在引擎强制的限额内协调插件自带的 agent。

名称按包加命名空间：`visualization/omics-visualization`、`gatk4/haplotypecaller`。

### 4.4 内核契约与稳定性

引擎先固定。内容能依赖的一切都是成文的契约，配有 JSON schema、校验器和一致性测试；契约之外的东西，内容一律不能依赖。

| 契约 | 定义的内容 | 起点 |
|---|---|---|
| 包 | `phi-package.yaml`、文件布局、`files.json`、版本、`dependsOn`、`requires`（§5） | 新建 |
| Skill | `SKILL.md` frontmatter（名称、描述、`attachTo` 所属专家、环境、声明的脚本），`references/` / `scripts/` / `assets/` 布局，`skill_run` 语义 | 运行时的 skill 格式加 Phi 扩展字段 |
| Wrapper | `wrapper.yaml`、`params.json`、`main.nf` 入口、输出、profile、按工具族打包和 `dependsOn` | 现有的 wrapper 技术设计 |
| 连接器 | MCP（标准协议）加 Phi 元数据：传输方式、认证类型、分类、stdio 服务的环境 | 现有的 MCP 目录 |
| Agent 定义 | frontmatter（`name`、`description`、`tools`、`spawns`、`model`、`thinkingLevel` 沿用 omp；另加 Phi 的 `visibility`、`skills`、`delegationMode`、`delegation`、`fallback`、环境绑定），结构化结果用 `outputSchema` | 现有的 `resources/agents` 格式，与 omp 的 agent 定义对齐 |
| 插件 | 组件、环境绑定、入口 agent、`requires.coreTools`（§11） | 新建 |
| 环境 | `environment.yml` 加各平台锁文件、运行时注入的变量（§9.4） | 新建 |
| 核心服务 | 内容可以调用的引擎工具：获取、`skill_run`、`wrapper_*`、产物展示、文件 / shell / 审批行为 | 部分已有 |
| 编排 | 构建在 omp 的 `task` / 注册表 / `hub` 之上的多 agent 运行模型：内部 agent、由 manifest 生成的 `spawns`、运行级存储、预算、检查点、可选的编排器组件（§11.4） | v1 中预留扩展点；之后以次版本、增量的方式加入契约，并随参考实现一起冻结（阶段 9） |
| 产物 | 内容交给界面的类型化输出：`table`、`figure`、`structure`、`molecule`、`network`、`report`，每个文件附一个描述它的 JSON 旁注文件 | 由现有的 db 结果查看器推广而来 |

规则：

1. **引擎里不放领域逻辑。** 领域行为都放在内容里。现有的 `viz_*` 工具和 `db_*` 工具链是目前仅有的两处违例，都要移出引擎（§10、§11.3）。界面负责渲染类型化产物，永远不需要插件提供代码。
2. **先验证，再冻结。** 每个契约在定为 v1 之前，都要有一个参考实现：一个 API skill（`protein-apis`）、一个 wrapper 工具族（`samtools`）、一个 MCP 连接器（一个远程目录条目加一个 stdio 服务），以及可视化插件。这一步发现的问题只修改草案，不修改 v1。v1 同时预留编排所需的扩展点（插件的 `orchestrator` 槽位、agent 的 `visibility` 字段、`orchestration.*` 服务命名空间），这样以后加入多 agent 编排只是增量修改。
3. **契约版本。** 每个契约都带 `contractVersion`（semver），包要声明自己面向的契约版本。次版本只允许增量修改（新增可选字段、新的产物类型、新的核心服务）；未知的可选字段一律忽略。
4. **破坏性修改是例外。** 升级主版本需要一份 ADR，并设置至少两个 app 版本的弃用窗口（期间新旧两个版本都被接受），同时提供迁移说明或自动兼容层。
5. **所有人用同一个校验器。** 同一个校验器（`phi validate`）在包构建器、registry 的 CI 和用户创作的内容上运行，作者看到的错误和 Phi 实际执行的检查完全一致。
6. **一致性测试套件。** 每个契约都有一组必须能加载、路由和运行的样例包；这套测试属于 `npm test` 的一部分，是修改引擎的门禁。

## 5. 包模型

### 5.1 Manifest

每个包的根目录都有一个 `phi-package.yaml`：

```yaml
schemaVersion: 1
id: protein-apis                 # 在 registry 内唯一，[a-z0-9-]
type: skill                      # skill | wrapper | mcp | plugin
version: 1.2.0                   # semver
title: 蛋白质数据库 API
summary: UniProt, PDB, AlphaFold and InterPro lookups
minAppVersion: 0.9.0
requires:
  coreTools: [fetch, skill_run]  # 本包依赖的引擎能力
dependsOn:                       # 其他包，安装时解析
  - id: samtools
    version: ">=1.0.0"
environment:                     # 可选，见 §9
  spec: environment/environment.yml
  lock: environment/conda-lock.yml
files: files.json                # 生成的白名单，每个文件带 sha256
```

插件的 manifest 还要列出组件：

```yaml
type: plugin
components:
  agents: [agents/Visualization.md]   # visibility: entry | internal（在 frontmatter 中声明）
  skills: [skills/omics-visualization]
  wrappers: []
  mcp: []
  orchestrator: null                  # 预留：workflow.yaml 或脚本（§11.4）
```

可视化插件刻意不声明私有环境：它的 agent 与 skill frontmatter 都绑定共享的
官方环境 `phi:r@1`。

### 5.2 白名单打包

包只能从白名单构建（git 跟踪的文件，或者 wrapper 现有的 `index.json` digest 列表），不能整目录复制。构建器拒绝 `.nextflow/`、`.nextflow.log*`、`work/`、`results/`、`.DS_Store`、`__pycache__` 以及所有未跟踪的文件。包的体积超过所属类型的预算时，CI 直接失败，除非 manifest 里写明了理由。这一机制取代 `copyWrapperSourceTree` 的整目录复制。

### 5.3 版本与兼容性

- 每个包独立使用 semver。安装前检查 `minAppVersion` 和 `requires.coreTools`；不兼容的包会显示出来，但不可安装。
- 核心工具带 API 版本号。删除工具或改动工具签名属于引擎的破坏性变更，必须升级这个版本号。

## 6. Registry 与安装器

### 6.1 Registry

registry 是一个 `index.json`，每个包列出：`id`、`type`、`version`、`summary`、`sha256`、`size`、`url`、`dependsOn`、`minAppVersion`，以及目录展示用的元数据（分类、预览图）。有两种来源，实现同一个接口：

- **本地 registry**：构建时由 `resources/` 生成。阶段 2–6 使用，所有安装和升级路径都能在不联网的情况下测试。
- **远程 registry**：阶段 7。静态文件放在对象存储加 CDN 上；index 用保存在 CI 中的 ed25519 私钥签名，app 内置公钥，拒绝未签名或签名不符的 index。镜像只要提供同一把密钥签名的内容，就可以使用。

### 6.2 安装流程

1. 解析依赖，展示安装计划（包、体积、环境开销）。
2. 下载（远程）或读取（本地）每个归档，校验 `sha256`。
3. 解压到 `~/.phi/.staging/` 下的临时目录，拒绝绝对路径、`..` 片段，以及指向根目录之外的符号链接。
4. 按 `files.json` 逐个校验文件。
5. 原子地 rename 到目标位置，写入 `.source.json`（`{ registry, id, version, sha256, installedAt, installedBy: user|dependency }`）。
6. 如果声明了环境，就构建或复用环境（§9）。也可以推迟到首次使用时构建，届时要明确显示进度。

只有当没有其他已安装的包依赖它时，才能卸载一个包；仅作为依赖安装的包，在最后一个依赖方卸载后由 GC 回收。升级时先把新版本并排装好，再原子切换。

### 6.3 磁盘布局

```text
~/.phi/
  packages/<type>/<id>/<version>/  # 已安装的包；type 为 skill、wrapper、mcp、plugin
  wrappers/tree/                   # 拼装好的 wrapper 目录树（§8.2）
  skills/<name>/                   # 用户创作的 skill（不是包）
  wrappers/custom/<id>/            # 用户创作的 wrapper
  runtime/                         # micromamba 根目录：envs/、pkgs/、mambarc（运行时基础设计 §2）
  state/enabled.json               # 启用状态，全局加项目级覆盖
  .staging/
<project>/.phi/skills/...          # 项目级单元
```

## 7. 加载、启用与目录界面

- 一个**内容加载器**合并四种来源：内置、已安装的包（单元和插件组件）、用户创作、项目。每一项都记录来源。
- 安装是全局的；启用状态是全局的，支持项目级覆盖（例如单细胞项目只启用相关能力）。
- 只有已启用的内容才进入运行时。主提示词里可以放一行简短说明，列出"可安装但未安装"的能力，这样 agent 会建议"安装 X"，而不是以为这个能力不存在。
- 界面：
  - **Skills / Wrappers / 连接器各页面**列出所有单元并标注来源（`插件：可视化`、`我创建的`、`项目`、`导入`），支持启用和停用，并提供"从目录添加"对话框，交互仿照 `McpConnectorCatalogDialog`。属于插件的单元可以停用，但不能单独卸载。
  - **插件页**只列复合能力，卡片上展示示例产出和环境开销。
  - 现有的 pi 插件页改名为**开发者扩展**，移到高级设置里。界面上的"插件"一词只指 Phi 插件。
- 用户之后可以把自己的一组单元"打包成插件"分享出去。正是两层模型让这条路径成为可能。

## 8. 单元类型

### 8.1 Skills

- 由 `SKILL.md` 加可选的 `references/`、`scripts/`、`assets/` 组成。
- 带脚本的 skill 要声明环境（§9），并在文档里用 `skill_run` 的方式调用脚本。第三方 skill 如果仍然写的是 `python x.py`，通过 `bash` 依然能运行，只是没有环境保障。
- 内置（不打包）的：`create-wrapper`。

### 8.2 Wrappers

- 包的粒度是 `modules/<provider>/<tool>/` 下现有的**工具族**（例如 `samtools`、`gatk4`、`bcftools`）。每个 subworkflow 单独成包，用 `dependsOn` 声明它引用的 module。场景组合（例如"变异检测"）做成只含依赖声明的元包。
- subworkflow 通过相对路径引用 module（`../../../../modules/nf-core/...`）。因此安装后的 wrapper 包会被拼装到 `~/.phi/wrappers/tree/` 下同一棵目录树里，与源码布局保持一致，而不是去改写 include 路径。
- 这修订了 `composition/packs.ts` 的前提：目录树按包合并、按包校验（各自的 `files.json`），不再整体替换、整体校验。按包安装上线后，overlay pack 机制随之废弃。
- `wrapper_search` 只返回已安装且已启用的包里的 wrapper。

### 8.3 MCP 连接器

- 远程 HTTP 连接器沿用现有的目录流程；目录改为由 registry 提供数据，这样新增条目不需要发布 app。
- 本地 stdio 连接器是带环境（§9）的包，启动命令由环境管理器解析。
- MCP 用于三类数据源：需要身份或凭据的源（实验室 LIMS、机构或商业数据库）、直接调用确实很困难的协议，以及已有的官方服务（目录里已有 PubMed、bioRxiv、ClinicalTrials）。公开的 REST 数据库做成 API skill（§10），不做成 MCP。

## 9. 运行环境

> **细节以 [运行时基础设计](phi-runtime-foundation.zh-CN.md) 为准**（运行时布局、explicit 锁文件、只读前缀、激活快照、`runInEnvironment`、绑定规则）。本章与之不一致的地方，以基础设计为准。

### 9.1 环境管理器

- Phi 内置 **micromamba**（单个静态二进制，BSD 许可证），不需要、也不使用用户自己的 conda。
- 环境放在 `~/.phi/runtime/envs/<envId>/`，`envId` 由平台和锁文件内容计算。依赖相同则共享同一个环境。micromamba 从包缓存用硬链接安装，所以每多一个环境，额外占用通常只有几十到几百 MB，而不是完整复制一份。
- 环境**不可变**。依赖一变化就生成新环境，旧环境在没有包引用后被回收。
- 可以配置 channel 镜像（tuna / ustc）和代理。doctor 会比对环境与锁文件，发现漂移时提供重建。

### 9.2 共享基础环境

- `phi-python`：官方 skill 使用的、锁定的科学计算 Python 栈（numpy、pandas、matplotlib、scikit-learn、scanpy、scvelo、rdkit、openpyxl、defusedxml、pillow、markitdown……），包列表来自现有 skill 实际的 import。由 Phi 开发者维护，带版本号，按平台锁定（osx-arm64、osx-64、linux-64）。
- `phi-r`：可选，按需构建，由 R notebook、scanpy R 互操作和可视化插件共用。

### 9.3 为包解析环境

- **官方包**必须能被当前的 `phi-python` / `phi-r` 满足，否则要自带锁文件（即插件）。CI 运行 `micromamba install -n phi-python --file environment.yml --dry-run --json`，结果不是"无变更"就判定失败。客户端永远不会为官方 skill 求解环境。
- **用户添加的 skill** 在添加时解析：
  1. 没有环境文件 → 使用 `phi-python`，并标记"未声明依赖"；
  2. 对 `phi-python` 做 dry-run，结果无变更 → 直接使用；
  3. 对每个已有的受管理环境做 dry-run → 复用第一个无变更的；
  4. 把 `phi-python` 的规格和 skill 的规格合并求解 → 新环境；
  5. 如果合并时冲突，只用 skill 自己的规格求解 → 新环境。
  第 4、5 步要展示体积和耗时估计，并由用户确认。
- 环境不能叠加：Python 的 `site-packages` 和 R 的库路径无法跨 conda 环境叠加，所以"基础环境加额外包"永远是一个重新完整求解的新环境。
- 用户 skill 可以有 `pip:` 段，但它不在锁文件的保证范围内；官方包只用 conda-forge / bioconda。如果以后依赖 PyPI 的包变多，再评估 pixi。

### 9.4 使用环境

| 执行路径 | 机制 |
|---|---|
| skill 脚本 | `skill_run(skill, script, args)` 解析出 skill 的环境，按绝对路径执行脚本 |
| 插件 agent 的 `bash` | 会话创建时把插件环境放到 `PATH` 最前，设置 `CONDA_PREFIX`、`PYTHONNOUSERSITE=1`，清空 `R_LIBS_USER` |
| 执行内容的核心服务（`skill_run`、`wrapper_*`） | 接收一个运行上下文，其中带有环境里解释器的绝对路径，不再去 `PATH` 上查找 |
| stdio MCP 服务 | 用环境里的解释器启动 |
| wrapper | 不变，由 Nextflow 按 process 管理环境 |

包的环境是只读的。项目进行中需要的额外包，装进项目环境（它本身也是一个完整求解的环境，§9.3）。运行记录（`wrappers/reproducibility.ts` 及同类记录）要写入环境 hash。

### 9.5 远程与离线主机

阶段 8。在远程主机上用 micromamba 按同一份锁文件重建环境；对离线集群，在本地或登录节点用 `conda-pack` 打包，再通过现有的远程上传链路传上去。面向这类集群时，linux-64 的锁文件必须兼容 glibc 2.17。

## 10. 数据访问重构

### 10.1 依据

附录 A：在 14 个结构化查询、每题跑 2 遍的评测中，直接读 URL 的准确率不低于 `db_*` 工具链（27/28 和 26/28，对比 25/28），耗时约为一半，工具调用约为三分之一；`db_query` 有 40% 的调用失败。工具链的失败来自连接器覆盖不全，fetch 的失败来自上游错误且没有重试。

### 10.2 目标形态

- **核心获取工具（引擎）**：支持 HTTP GET/POST，带重试和退避、按主机限速、域名白名单和出口审计（复用 `policy-*.ts`）、分页和批量下载（来自 `db_download`）、可选地从凭据存储注入凭据，并把结果登记给 UI 结果查看器（`DbQueryResultPreview` 及各查看器）。
- **API skill（单元）**：约 6 个领域 skill，取代 32 个连接器：`protein-apis`、`genomics-apis`、`chemistry-apis`、`pathway-network-apis`、`clinical-cancer-apis`、`ontology-apis`。`SKILL.md` 保持简短，各数据库的细节放在 `references/<db>.md`（端点、查询语法、字段、示例、常见坑）。现有的连接器 manifest 就是编写这些参考文档的素材。
- **辅助脚本**（只用标准库，通过 `skill_run` 运行）：处理 GraphQL（gnomAD）、SPARQL（UniProt、WikiPathways），以及 Entrez 的批量和限速访问。
- **Database agent** 仍然是主 agent 的唯一入口，工具换成核心获取工具、`skill_run` 和这些 API skill。
- **MCP** 用于需要认证或私有的数据源（§8.3）。

### 10.3 退役门槛

只有当 `api-skill` 组在扩展评测（§14）中全面不差于 `db` 组时，才退役 `db_*` 工具和连接器 manifest 格式（`create-database-connector` skill 已提前删除）。扩展评测要包括批量、分页、GraphQL、SPARQL 和限速类任务。在此之前，两条路径并存。

## 11. 插件

> 插件的目录结构，以及插件 agent 如何使用环境，见 [运行时基础设计](phi-runtime-foundation.zh-CN.md) §6–§7。

### 11.1 形态

插件是 `type: plugin` 的包，组件安装在 `~/.phi/packages/plugin/<id>/<version>/` 下，并按类型路由（§4.3）。一个插件最多有一个**入口** agent，也就是主 agent 能看到的那个；插件还可以带任意数量的**内部** agent，它们只对插件自己的编排器可见（§11.4）。私有环境在 manifest 中声明，组件绑定在 agent 或 skill frontmatter 中声明；可视化改用 `phi:r@1`。

### 11.2 代码边界

插件不携带在 Phi 内部运行的 TypeScript。新的程序化能力只能来自三处：(a) 在 `requires.coreTools` 中引用的核心工具；(b) 在组件解析到的受管理环境中运行的脚本；(c) MCP 服务。这样审核、安全和兼容性才可控。

### 11.3 第一个插件：可视化

`agents/Visualization.md` + `skills/omics-visualization`（脚本、模板、预览图），两者都绑定共享的 `phi:r@1` 环境。现在的 `viz_*` 工具（`src/main/agent/visualization/`，约 1100 行）是放在引擎里的领域逻辑；它们会改写为命令行程序 `scripts/viz.py`，以脚本工具的形式声明（运行时基础设计 §5.1），工具名保持不变，图表以 `figure` 产物的形式返回。这一步在实施计划的步骤 4 完成，早于插件机制；步骤 6 只负责打包。可视化插件同时是插件、官方环境引用、脚本工具、产物四个契约的参考实现。

之后可能的插件：单细胞分析（agent + scanpy/scvi skill + starsolo wrapper + torch 环境）、bulk RNA-seq（agent + skill + wrapper 元包）。

### 11.4 多智能体编排（预留）

有些能力不是一个专家就能完成的，而是一个团队多轮协作的结果。例如类似 AI Co-Scientist 的系统：由一个主管驱动生成、反思、排名（锦标赛式的两两辩论）、演化、相近性分析和元评审等 agent，持续运行数小时，不断积累假设和评分。Phi 应当能运行这类系统，用户也应当能定义自己的系统，而且**不需要为每个系统修改引擎**。

#### 11.4.1 omp 运行时已经提供的能力

Phi 运行在 `@oh-my-pi/pi-coding-agent`（18.1.10）之上，它本身已经带有一套多智能体运行时。Phi 应当在它之上构建，而不是重新实现：

| omp 能力 | 位置 | 对 Phi 的价值 |
|---|---|---|
| `task` 工具 | `src/task/` | 拉起具名 agent；批量模式下多个任务共享一份 `context` 并行执行；每次拉起可指定 `outputSchema`，按 `strict` / `permissive` 校验结构化结果；可选的 worktree 隔离 |
| Agent 定义 | `src/task/agents.ts` 及 discovery | Markdown 加 frontmatter（`name`、`description`、`tools`、`spawns`、`model`、`thinkingLevel`、`blocking`）；来源分为 `bundled` / `user` / `project` |
| 拉起策略 | `task/spawn-policy.ts` | 用 frontmatter 里的 `spawns` 限制一个 agent 能拉起哪些 agent；`task.maxRecursionDepth`（默认 2）限制嵌套深度 |
| Agent 注册表与生命周期 | `src/registry/` | 进程级的 agent 注册表，id 稳定；状态分为 `running` / `idle` / `parked` / `aborted`；parked 的 agent 保留会话记录，**需要时可以恢复** |
| `hub` 工具加 IRC 总线 | `tools/hub/`、`irc/bus.ts` | agent 之间互发消息（`send`、`wait`、`inbox`、`list`），可以唤醒空闲的 agent 或恢复 parked 的 agent；后台作业（`start`、`ps`、`logs`、`stop`） |
| 限额 | `task/executor.ts`、`task/provider-concurrency.ts` | 每个 agent 的软性请求预算，超出时引导收尾；按模型提供商限制并发 |
| 长时间循环 | `src/goals/`、`src/autoresearch/` | 带持久化状态和 token 统计的目标模式；带仪表盘、可续跑的 autoresearch 循环 |

Phi 目前基本绕开了这些能力：专家会话用 `restrictToolNames` 创建，而 Phi 自己的 `agents/registry.ts` 又重新实现了一遍并行、后台、steer 和停止。阶段 1 需要决定是否把 Phi 的专家迁移到 omp 的 `task` 和注册表上（建议迁移，见 §11.4.4）。

#### 11.4.2 契约归 Phi，实现用 omp

内核稳定的原则（§4.4）仍然成立：插件只依赖 Phi 的契约，绝不直接依赖 omp 的内部实现，因为 omp 迭代很快。编排契约和 Agent 定义契约被设计成**一个由 Phi 维护、能很薄地映射到 omp 的子集**，由引擎里的一层适配器完成映射：

- Phi 的 agent frontmatter 在 omp 已有的地方沿用 omp 的字段名（`name`、`description`、`tools`、`spawns`、`model`、`thinkingLevel`），再加上 Phi 自己的字段（`visibility`、`skills`、`environment`、`delegation*`）。
- 结构化结果尽量用 `outputSchema`（JSON Schema），而不是 Phi 专用的解析协议。
- 适配器由一致性测试覆盖；如果某次 omp 升级破坏了它，会在发布前的 CI 里失败。omp 版本保持锁定。

omp 缺少、由 Phi 以引擎服务形式补上的部分：

- **硬性预算**：按每次编排运行限制 token、费用、墙钟时间和累计 agent 运行次数，叠加在 omp 的软性请求预算之上。
- **运行级的类型化存储**：集合及其 JSON schema 在 manifest 中声明（假设、评审、比赛记录），历史只追加，可在界面中浏览。agent 通过存储工具访问，不随意写文件。
- **审批、项目边界、远程项目防护**：对每一个被拉起的 agent 生效（沿用 Phi 现有的扩展）。
- **人工检查点**：在界面上呈现（`ask` / `checkpoint`）。
- **运行视图**：展示 agent 树、消息、存储和预算，支持暂停 / steer / 停止，完成后唤醒主 agent（沿用现有的 `run-host.ts`）。
- **按包隔离**：插件的 agent 只能拉起同一插件内的 agent，以及 `uses` 中列出的公共专家。`spawns` 由 manifest 生成，而不是信任 agent 文件里写的内容。

#### 11.4.3 三种编写层级，同一个运行时

1. **agent 驱动（默认，最先提供）**：入口 agent 担任主管，它的 `spawns` 列出插件的内部 agent。主管使用 `task`（批量、`outputSchema`）、`hub` 消息和存储工具。除了上面那些 Phi 服务，不需要任何新的运行时概念。
2. **声明式** `workflow.yaml`（以后）：顺序、并行扇出 / 汇聚、对存储中的记录逐项处理、循环直到满足条件或预算耗尽、人工检查点。面向想要可重复流程、又不想写代码的用户。
3. **程序式**编排器（以后）：插件环境中的一个脚本，通过本地 JSON-RPC 通道驱动同一套原语，用于确定性的流程（Elo 锦标赛、演化选择）。

三种层级都通过同一个适配器拉起 agent，所以限额、审批、存储和运行视图的行为完全一致。

#### 11.4.4 阶段 1 要做的决定

在冻结 v1 之前，先通过适配器把 Phi 现有的专家（Database、Wrapper、Visualization）在 omp 的 `task` 和注册表上做一个原型。如果原型能保持现有行为（审批路由、远程防护、运行卡片、steer / 停止、唤醒），就退役 `agents/registry.ts` 中重复的部分，专家委派和编排都基于 omp；如果做不到，把差距写进 ADR，v1 继续使用 Phi 自己的注册表。

#### 11.4.5 类似 Co-Scientist 的插件示意

```yaml
id: co-scientist
type: plugin
components:
  agents:
    - agents/CoScientist.md        # visibility: entry；主管
    - agents/Generation.md         # visibility: internal
    - agents/Reflection.md
    - agents/Ranking.md
    - agents/Evolution.md
    - agents/MetaReview.md
  skills: [skills/elo-tournament]  # 只用标准库的脚本，通过 skill_run 运行
  orchestrator: null               # 层级 1：由主管 agent 负责编排
uses: [Database]
store:
  hypotheses: schemas/hypothesis.json
  reviews: schemas/review.json
  matches: schemas/match.json
limits:
  concurrentAgents: 4
  totalAgentRuns: 200
  wallTime: 6h
  tokens: 5000000
```

主管批量拉起生成任务（`outputSchema` 为假设）→ 反思 agent 逐条评审假设 → 排名 agent 两两辩论并记录比赛结果，Elo 更新由 skill 脚本对存储执行 → 演化 agent 改进排名前 k 的假设 → 元评审 agent 写出反馈，作为下一轮的输入 → 预算耗尽或结果收敛时停止 → 输出一份 `report` 产物和排好序的假设列表。以后可以把这个循环改成层级 3 的程序式编排器，agent 本身不需要改动。

## 12. 安全与信任

- 信任层级与 wrapper 的 `trustTier` 保持一致：`builtin`、`official`（签名 registry）、`user`（本地创作）、`imported`（第三方文件）。
- `skill_run` 纳入审批体系：在 `ask` 模式下，`imported` skill 的脚本首次运行需要确认；`official` 内容视为可信。
- 核心获取工具保留现有的 URL 策略和出口审计。
- 包的解压限制在目标路径内（§6.2）；远程 registry 必须验证签名。

## 13. 迁移

- 改动上线后首次启动时，用户会话历史中出现过的 skill 会从本地 registry 自动安装，其余内容收进目录。如果更简单，也可以直接重置 beta 数据（beta roadmap 本来就允许干净重置）。
- `~/.phi/wrappers/installed/` 下的内置条目改为按包安装；`custom` wrapper 移到 `~/.phi/wrappers/custom/`。
- `~/.phi/wrappers/packs/` 下的 overlay pack 失效，由 GC 删除。
- `resources/` 变成本地 registry 的源码树。远程 registry 上线后，安装包只保留引擎内容（内置 skill、引擎所需的 agent 定义、micromamba）。

## 14. 以评测作为门禁

`scripts/eval/` 变成常规回归套件：

- **数据访问**：扩展 `db-vs-fetch.ts`，加入 `api-skill` 组，以及批量、分页、GraphQL、SPARQL、限速类任务；
- **路由**：在给定的已安装内容下，主 agent 是否选对了专家或 skill，缺少能力时是否会建议安装；
- **上下文开销**：用 `npm run usage:report` 统计首轮提示词大小和每次运行的 token。

每个阶段的出口标准（§15）都引用这些数据。凭据允许的话，至少用两个模型跑。

## 15. 阶段

> **实施顺序以 [内容分发实施计划](../roadmap/content-distribution-implementation.zh-CN.md) 为准**（P0–P9，每个契约在第一个使用者动工前冻结）。下表只作为工作内容的逻辑分组保留。

| 阶段 | 范围 | 出口标准 |
|---|---|---|
| 0 | 本文档、`docs/decisions/` 中的 ADR、roadmap 更新 | 负责人确认 |
| 1 | **内核契约 v1**（§4.4）：规范、JSON schema、`phi validate`、一致性样例、每个契约的参考实现；在 omp 的 `task` / 注册表上做专家原型（§11.4.4）；然后冻结 | 所有参考实现都通过校验并能运行；关于 omp 的决定写入 ADR；契约标记为 v1 |
| 2 | 内容卫生（白名单构建器、体积预算）、包格式、本地 registry、安装器、带来源记录的内容加载器、评测套件 | 一个 skill 的安装、升级、卸载端到端跑通；评测套件一条命令跑完 |
| 3 | skills、wrappers、连接器的目录界面；wrapper 按工具族拆包和目录树拼装；启用管理；迁移 | 路由评测不退化；主 agent 首轮提示词有可测量的下降 |
| 4 | micromamba、`phi-python`、环境解析、`skill_run`、环境管理页 | 带脚本的 skill 在干净机器上无需手动配置即可运行 |
| 5 | 核心获取工具、API skill、辅助脚本、扩展评测、退役 db 工具链 | 扩展评测中 `api-skill` ≥ `db`；`src/main/agent/db/` 大幅缩减 |
| 6 | 插件加载器、可视化插件（可视化逻辑移出引擎）、`phi-r`、"开发者扩展"改名 | 可视化插件能安装、在共享 `phi-r` 中运行、干净卸载；引擎里不再有 `viz_*` 代码 |
| 7 | 带签名的远程 registry、更新、镜像、离线导入、安装包瘦身 | 全新安装能从服务器获取内容；离线也能首次启动 |
| 8 | 远程 / HPC 上 skill 和环境的对齐 | skill 脚本能在 SSH 项目和离线集群上运行 |
| 9 | 基于 omp 的编排（§11.4）：由适配器生成的 `spawns`、硬性预算、存储工具、人工检查点、运行视图；层级 1 的类似 Co-Scientist 参考插件；冻结编排契约。层级 2、3 按需跟进 | 参考插件在预算内完成多轮运行，能在 app 重启后恢复（恢复 parked 的 agent），并且可以 steer 和停止 |

依赖关系：0 → 1 → 2 → 3；2 → 4 和 2 → 5（4 与 5 可以并行）；4 → 6；3、5、6 → 7 → 8；6 → 9。阶段 2–8 都是在实现已冻结的契约；如果其中某个阶段需要修改契约，必须按 §4.4 的规则走，而不是直接改引擎打补丁。建议：阶段 0–5 纳入内部 beta，6–9 放到 beta 之后。

## 16. 待定问题（附建议的默认答案）

1. beta 范围：阶段 0–5 是否纳入 beta？*建议：是。*
2. 阶段 4 和 5 是否并行？*建议：是。*
3. 阶段 5 期间，是否保留内置的 db 连接器作为离线兜底？*建议：保留，直到通过退役门槛。*
4. 更新策略：自动更新还是提示更新？*建议：提示更新；对已记录版本的项目不做静默升级。*
5. 是否为协作者提供项目级的包和环境锁？*建议：阶段 2–5 先写进运行记录，项目级锁以后再做。*
6. registry 托管：静态对象存储加 CDN，不要后端？*建议：是；国内用户走国内 CDN。*
7. 签名私钥的保管：CI secret 还是离线签名？*建议：第一方 registry 用 CI secret。*
8. 界面命名：Phi 插件叫"插件"，命名空间分隔符用 `/`。*建议：是。*
9. 是否把 `viz_*` 移出引擎、放进可视化插件？*已确认：是，在实施计划步骤 4 改写为命令行程序加脚本工具声明；引擎只保留产物展示。*
10. 契约的形式：JSON Schema 文件放在 `docs/contracts/` 下，引擎内置校验器？*建议：是；schema 是唯一的事实来源，文字规范链接到 schema。*
11. 编排的编写方式：先提供基于 omp 的 agent 驱动编排（层级 1），声明式和程序式以后再加？*建议：是。程序式编排器加入后，在沙箱方案评审通过之前，只对 `official` 和 `user` 两个信任层级开放。*
12. 是否把 Phi 的专家委派迁移到 omp 的 `task` / 注册表上，并退役 `agents/registry.ts` 中重复的部分？*建议：阶段 1 根据原型结果决定（§11.4.4）。*

## 附录 A：db 工具链与 fetch 的对比评测（2026-09-29）

评测脚本：`scripts/eval/db-vs-fetch.ts`（题目在 `db-vs-fetch-tasks.ts`，评分在 `db-vs-fetch-score.py`）。每一组都作为独立的内存 agent 会话运行，只挂自己那组工具。模型：`cursor/cursor-grok-4.6-fast`（经 `scripts/eval/cursor-bridge.mjs` 调用）；经 Cursor 调用 Claude 被地区限制拦下。14 道题，每题 2 遍；标准答案已对照线上 API 核实。

| 组别 | 正确 | 工具调用/题 | 失败调用/题 | 平均耗时 s | p90 s | 输出 token/题 |
|---|---|---|---|---|---|---|
| `db`（Database 提示词 + 7 个 `db_*` 工具） | 25/28 | 8.0 | 1.2 | 45.7 | 89.3 | 6891 |
| `fetch`（只有 URL 读取） | 27/28 | 3.0 | 0.2 | 23.9 | 44.4 | 2941 |
| `fetch-hints`（URL 读取 + 约 15 行 API 地址） | 26/28 | 2.5 | 0.1 | 20.9 | 37.8 | 2338 |

- `db_query` 79 次调用中失败 32 次（过滤条件或字段格式错误）。发现类工具（`db_search`、`db_routes`、`db_domain`、`db_resolve`）调用 124 次，而取数据的调用只有 79 次。
- `db` 组的失败：PDBe `entry_summary` 缺少分辨率（2 次）；AlphaFold 域缺少全局 pLDDT（1 次）。
- `fetch` 组的失败：AlphaFold 返回 403（2 次）、Ensembl 返回 500（1 次），且没有重试。
- 未覆盖：批量下载、分页、需要凭据的 API、冷门数据库、限速、UI 结果查看器。Cursor 不上报输入 token，所以只比较了输出 token。只用了一个 fast 档模型。
