# Phi 的 BioMCP connector

[English](README.md) | **中文**

此 `phi-packages` connector 在独立的 Phi 托管环境中运行官方
[BioMCP](https://github.com/genomoncology/biomcp) 原生程序。

## 版本与平台

- Connector 包版本：**1.1.0**，需要 Phi **1.0.1** 或更新版本。
- 上游包：**biomcp-cli 0.9.1**，使用 [PyPI](https://pypi.org/project/biomcp-cli/0.9.1/) 官方原生 wheel。
- 支持 macOS Apple Silicon、macOS Intel 和 Linux x86-64；Linux 需要 glibc 2.28 或更新版本。
- [官网图标](../../../docs/content-icons.md#biomcp-website-icon)继续由内容包提供。

## 安装与启动

在 Phi 的**连接器目录**中安装或更新 BioMCP。Phi 校验包后先准备其声明的原生
托管环境，完成后才启用 connector。安装阶段下载所选平台约 15–17 MB 的 wheel，
校验固定大小和 SHA-256，然后只提取声明的原生程序。进度和错误统一归属托管安装流程；
准备失败时 connector 保持不可用，可以重新安装。

`environment.yml` 是唯一的程序下载固定信息来源，声明 native backend、
`biomcp` 可执行程序，以及各平台官方 wheel 的地址、大小、校验值和准确 ZIP 成员路径。
`locks/<平台>.txt` 仅包含 `@EXPLICIT`，没有 Conda 依赖。程序安装到 Phi 的托管
运行时前缀中；connector 直接启动 `biomcp serve`。

此包不再包含 Python 启动器，也不依赖 Python 运行时、`pip install`、Conda 包安装、
Git clone 或首次启动时的软件下载安装。准备完成后启动直接复用已安装程序。
在线生物医学查询仍需访问相应上游服务。旧 1.0.x 启动器的缓存保留原样，新包不再使用。

1.1.0 源码需按[签名目录发布流程](../../../docs/publishing.md)重新构建发布，应用才能安装。
之前的独立 [1.0.0 发行包](https://github.com/dxsbiocc/phi-packages/releases/tag/biomcp-v1.0.0)
包含旧启动器，不提供此次托管安装流程。

## 凭据与来源

提供 BioMCP 只读生物医学工具，多数查询无需密钥；需要凭据的上游服务仍需各自的 API key。
参见[上游密钥说明](https://biomcp.org/getting-started/api-keys/)。

wheel 版本、地址、SHA-256 与大小沿用原 1.0.x 固定信息。三个已校验 wheel 均包含
`biomcp_cli-0.9.1.data/scripts/biomcp`，只安装该成员。`SOURCE.json` 记录官方仓库、
PyPI 分发来源、版本与图标来源；上游 MIT 版权及许可保留在
[LICENSE.upstream](LICENSE.upstream)。数据源授权另行适用，参见
[上游数据授权说明](https://biomcp.org/reference/source-licensing/)。

## 验证

将 Phi 与 `phi-packages` 并排放置，使用 Phi 已有的解析器和构建器：

```sh
cd ../Phi
node --import ./scripts/test-loader.mjs --test ../phi-packages/tests/biomcp-native.test.mjs
node --import ./scripts/test-loader.mjs --test tests/envs-application-native.test.ts tests/envs-application-artifacts.test.ts tests/envs-ensure.test.ts
```

内容测试验证原生安装声明、原固定信息、各平台 lock、环境所有权和身份、许可及归档允许列表。
Phi 的共享测试覆盖下载校验、缓存复用、并发、受限 ZIP 成员提取、取消和运行时就绪状态。
原启动器测试的职责转移到这些托管运行时测试，避免保留第二套安装实现。
将 `PHI_BIOMCP_WHEEL_DIR` 指向含三个已校验官方 wheel 的临时目录，可额外离线验证
Phi 原生安装器，无需测试期间联网。
