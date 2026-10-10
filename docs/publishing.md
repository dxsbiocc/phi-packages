# 发布官方软件源 / Publishing the official catalog

[中文 README](../README.md) · [English README](../README.en.md)

Phi 默认从本仓库的 `catalog-v1` GitHub Release 读取已签名的在线目录。
启动或打开目录只获取 `index.json`、`index.sig.json`、图标和 connector manifest；
点击安装时才下载所选软件包与必要依赖的归档。
断网时只使用通过内置公钥验证的现有缓存，不回退到应用的领域内容源码。

Phi uses this repository's `catalog-v1` GitHub Release as its default signed catalog.
Catalog browsing fetches the index, signature, icons, and small connector manifest sidecars.
Installation downloads only the selected package and its dependency closure.
Offline browsing uses an existing verified cache.

## 准备发布 / Prepare a release

两个仓库并排放置；构建、校验、签名实现使用 Phi 当前版本，不在内容仓库复制安装器。
需要 Phi 的开发依赖和 Node.js 24。私钥保存在仓库外，绝不提交或上传。

Keep `Phi` and `phi-packages` side by side. The content repository calls Phi's current
shared builder and validators. Install Phi's development dependencies and use Node.js 24.
Store the signing private key outside both repositories; never commit or upload it.

```sh
node scripts/prepare-catalog.mjs \
  --phi ../Phi \
  --out /tmp/phi-catalog-20261008 \
  --key ~/.phi/publishing/phi-packages-ed25519.pem
```

输出的 `assets/` 包含扁平文件名的 `.tar.gz`、图标、connector YAML、索引与签名，
可直接上传为 GitHub Release assets。GitHub Release 不支持目录路径，因此先扁平化再签名。
输出目录必须是新的路径，防止混入旧版本文件。

The resulting `assets/` directory contains flat archive, icon, connector YAML, index,
and signature filenames, ready for GitHub Release assets. Flattening happens before signing.
Use a fresh output path to keep earlier release files out of the new catalog.

上传所有 assets 到本仓库 `catalog-v1` release，索引与签名最后上传，避免客户端读取半成品。
客户端验证签名以及每个 sidecar/归档的 SHA-256 和大小，并通过完整缓存代际原子切换。
保留同一包版本的不可变归档；更改内容应提升包版本。

Upload all assets to the repository's `catalog-v1` release, uploading the index and signature
last. Phi validates the signature and every asset's SHA-256 and size before atomic cache activation.
Keep an existing package version's archive immutable; bump the package version when its contents change.

wrapper 默认版本为 `1.0.0`，必须与 Phi 内置 wrapper 树的版本一致，否则在线条目会被当作降级且依赖范围无法满足。
后续可显式使用 `--wrapper-version` 提升该组生成软件包版本；其他软件包使用各自 manifest 的版本。

Wrappers default to `1.0.0`, which must match the wrapper tree bundled with Phi; a lower version looks like a downgrade and breaks dependency ranges.
For later changes, pass an explicit `--wrapper-version`; other packages retain their manifest versions.

## 安装路径 / Installation paths

| 内容 / Content | 安装位置 / Installed location |
| --- | --- |
| Plugin | `~/.phi/packages/plugin/<id>/<version>/` |
| Skill | `~/.phi/packages/skill/<id>/<version>/` |
| Connector | `~/.phi/packages/mcp/<id>/<version>/` |
| Wrapper | `~/.phi/wrappers/tree/` 中保留依赖使用的目录结构 / preserved component tree |
| 官方源缓存 / Official cache | `~/.phi/cache/registries/phi-packages/generations/<index-sha256>/` |

图标始终跟随内容安装到其入口旁，公共目录中的 sidecar 只供安装前预览。
Connector credentials、运行环境和用户配置保留在用户数据目录，不进入 release。

Icons install beside their content entry; catalog sidecars provide previews before installation.
Connector credentials, installed runtimes, and user configuration never enter a release.
