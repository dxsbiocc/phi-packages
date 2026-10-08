# Phi 的 BioMCP connector

[English](README.md) | **中文**

这是 `phi-packages` 维护的 [BioMCP](https://github.com/genomoncology/biomcp) 适配包，
使用 Phi 现有的包安装和 stdio MCP 机制，不需要修改 Phi 项目源码。

## 版本与支持平台

- Connector 包版本：**1.0.0**。
- 上游 BioMCP：**0.9.1**，下载地址与校验值固定在 `upstream.json`。
- 支持 macOS Apple Silicon、macOS Intel、Linux x86-64。Linux 二进制需要
  glibc 2.28 或更新版本，与上游 wheel 的平台声明一致。
- 使用 Phi 已有的 `phi:python@1` 托管环境。核心环境定义和锁文件留在 Phi，
  本包只用 Python 标准库，没有新增 Python 库依赖。

## 不改源码的安装方式

从[发行页面](https://github.com/dxsbiocc/phi-packages/releases/tag/biomcp-v1.0.0)下载
[biomcp-registry-1.0.0.zip](https://github.com/dxsbiocc/phi-packages/releases/download/biomcp-v1.0.0/biomcp-registry-1.0.0.zip)，
解压后，先到 Phi 的 **Skills** 页打开技能目录，点击带文件夹图标的 **添加** 按钮，使用目录选择器
选择解压目录。这一现有入口会登记所有包类型的软件源。此源只有 connector，
因此技能列表为空是正常现象；关闭对话框，再打开 **连接器** 目录并安装 BioMCP。
解压目录包含 `index.json` 和 `mcp-biomcp-1.0.0.tar.gz`。
当前连接器对话框本身尚无本地目录选择按钮。
如果 Phi 提示 Python 环境尚未就绪，使用现有环境管理功能构建所引用的环境。

这个本地软件源没有签名，会显示为导入来源。本次使用现有的本地目录安装路径，
尚未接入 Phi 的自动远程软件源获取。

## 启动与缓存

启动器根据平台下载官方 PyPI 原生 wheel，检查固定的大小和 SHA-256，
只提取 `biomcp` 可执行程序。它不会执行 `pip install` 或修改托管环境。

缓存位于包目录之外：

- macOS：`~/Library/Caches/Phi/connectors/biomcp/<版本>/<平台>/`。
- Linux：`~/.cache/Phi/connectors/biomcp/<版本>/<平台>/`。

Linux 默认缓存不随 Phi 托管环境的 XDG 缓存目录变化，确保预热结果可以复用。
文件锁避免并发重复下载。每次启动都会校验缓存；程序损坏时从已验证 wheel 修复。
随后直接执行 `biomcp serve`，保留 MCP 标准输入输出流，启动器消息只写入 stderr。

首次启动会下载约 15–17 MB。慢网或离线使用前，可由运行 Phi 的同一用户预热缓存：

```bash
python /path/to/biomcp/server.py --prepare
```

通过 Phi 安装后，脚本位于 `~/.phi/packages/mcp/biomcp/1.0.0/server.py`。
测试或离线准备时可用 `--cache-dir /path/to/biomcp-cache` 指定专用缓存目录。
首次下载需要网络；在线生物医学查询仍需要访问相应上游数据源。

## 功能与凭据

提供 BioMCP 标准只读 MCP 工具，涵盖文献、基因、变异、药物和临床试验检索。
多数查询无需密钥。需要凭据的上游数据源仍需各自的 API key；本包不会扩展
Phi 当前 stdio 包契约的凭据字段，也不会把密钥写进源码。
参见[上游密钥说明](https://biomcp.org/getting-started/api-keys/)。

## 来源与验证

`upstream.json` 记录官方 `biomcp-cli` wheel 的版本、地址、大小及 SHA-256。
上游 BioMCP 的 MIT 版权与许可保留在 [LICENSE.upstream](LICENSE.upstream)，
这不构成本适配器或整个仓库的统一开源许可。
数据源的许可另行适用，参见[上游数据授权说明](https://biomcp.org/reference/source-licensing/)。

离线测试：

```bash
python3 -B -m unittest discover -s tests -p test_biomcp_connector.py
```

另外使用 macOS ARM64 官方程序和适配启动器验证 MCP 初始化及工具发现；
该验证不执行真实生物医学查询，也不需要调用有凭据的数据源。
