# Phi Packages

**中文** | [English](README.en.md)

Phi 的独立公开内容仓库，保存可安装领域能力的源码。

## 当前内容

| 内容 | 数量 |
| --- | ---: |
| 独立 skills | 14 |
| MCP connector 定义 | 19 |
| Phi 复合插件 | 1 |
| 插件内 skills | 1 |
| Wrapper 适配器 | 652 |

```text
resources/
  skills/                独立技能、参考资料和脚本
  connectors/            MCP 连接器声明
  plugins/               Visualization 插件及其组件
  wrappers/              Nextflow module、subworkflow、workflow 和适配器
docs/
  contracts/             Phi 内容契约参考快照
  design/                内容分发、运行环境和 wrapper 设计
  decisions/             已记录的分发决策
  roadmap/               内容分发实施计划参考
SOURCE.json              导入来源、版本、范围和排除项
```

每个可安装包独立管理版本；仓库的提交或发布版本不要求所有内容包同步升级。

## 与 Phi 主程序的关系

主程序仓库：[dxsbiocc/phi](https://github.com/dxsbiocc/phi)。

这次导入是当前受 Git 跟踪文件的工作区快照，来源提交及纳入的未提交内容见
[SOURCE.json](SOURCE.json)。Phi 原有资源仍保留，应用仍使用现有加载路径。

本仓库目前提供内容源码，尚未接通 Phi 的在线下载，也尚未生成已签名远程软件源。
后续需要由 Phi 的统一校验器和包构建器生成独立 `.tar.gz` 包、目录索引及签名，
再配置应用端远程软件源。核心运行时代码、安装器、校验器及应用界面留在 Phi 主仓库。

Phi 官方共享环境和应用级色板由主程序提供，本仓库不复制其核心定义。
插件自身运行所需的资源仍随插件保存；`docs/` 中的契约和设计为参考快照，
仍以 Phi 主仓库中的定义为准。
运行时二进制、环境安装目录、账户凭据及运行产物未导入。

## 分发边界

本仓库只保存 skill、connector、插件、wrapper 及它们自身必需的资源。
应用级 `resources/palettes/` 和核心 `resources/runtime/` 不属于内容分发范围。
差异表达 wrapper 的旧本地镜像目录已移除，其依赖迁移说明见
[差异表达依赖](resources/wrappers/modules/local/differential-expression/README.md)。

边界检查：`node --test tests/content-boundaries.test.mjs`。

## 导入边界与第三方内容

- `resources/plugins/office/` 暂未导入：现有 Office skills 许可证明确禁止复制和第三方分发；
  为保持组件完整性，本次排除整个 Office 插件。需要获得相应授权或替换受限组件后再纳入。
- `resources/skills/create-wrapper/` 留在 Phi：它属于引擎内置创作技能。
- 保留已导入文件原有的许可证、引用、作者和上游来源记录。本仓库包含 nf-core 等上游内容；
  导入不改变其原有许可条件。
- 本仓库目前未声明统一的开源许可证。

完整边界见 [SOURCE.json](SOURCE.json)，初始快照的检查结果见 [VALIDATION.md](VALIDATION.md)。

## 可安装 connector

[BioMCP](resources/connectors/biomcp/README.zh-CN.md) 提供独立的
[本地软件源分发包](https://github.com/dxsbiocc/phi-packages/releases/tag/biomcp-v1.0.0)，
可通过 Phi 现有目录安装功能添加，无需改动项目源码。

## 内容约定

- 尽量沿用现有目录布局及内容格式；包清单遵循 [Package contract](docs/contracts/package.md)。
- Wrapper 保留其 include 目标、辅助文件和测试数据，避免只复制适配器造成依赖丢失。
- 修改插件时同时核对组件引用；本地 MCP 服务仍在用户机器的受管理环境中运行。
- 不提交密钥、账户配置、环境安装目录或 Nextflow 运行结果。
