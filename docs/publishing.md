# 发布官方软件源 / Publishing the official registry

[中文 README](../README.md) · [English README](../README.en.md)

Phi 官方软件源由 `phi-site` 官网提供的 HTTP Registry API（`https://phi.omicsagent.cn/api/registry/`）直接分发标准无签名包目录，不再使用 Ed25519 签名或 `index.sig.json`。

客户端浏览目录时请求 `index.json`、图标与 manifest；点击安装时按 `index.json` 中的 `archive`、`sha256` 与 `size` 下载所需 `.tar.gz` 归档及依赖并校验完整性。

The official Phi package registry is served directly by `phi-site` via its HTTP Registry API (`https://phi.omicsagent.cn/api/registry/`) without Ed25519 signing or `index.sig.json`.

Clients fetch `index.json`, icons, and manifests while browsing, then download and SHA-256-verify only the selected `.tar.gz` archives and dependencies upon installation.

## API 端点 / Endpoints (`phi-site`)

| 路径 / Path | 内容 / Content | 缓存策略 / Cache-Control |
| --- | --- | --- |
| `/api/registry/index.json` | 注册表索引（含包版本、`size` 与 `sha256`） | `no-cache` |
| `/api/registry/<name>.tar.gz` | 不可变安装包归档 | `public, max-age=31536000, immutable` |
| `/api/registry/icons/<name>` | 包图标 | `public, max-age=31536000, immutable` |
| `/api/registry/manifests/<name>` | 包 manifest（`text/yaml`） | `public, max-age=31536000, immutable` |

## 从 `phi-site` 构建与部署 / Build and deploy from `phi-site`

生产服务器通过 `~/bin/sync-phi-packages.sh` 自动拉取 `phi-packages` `main` 分支、重新生成官网目录与注册表并重建站点；也可在并排放置的 `Phi`、`phi-packages` 与 `phi-site` 目录下手动执行：

```sh
PHI_DIR=../Phi \
PHI_PACKAGES_DIR=. \
PHI_REGISTRY_DIR=../phi-registry \
../phi-site/scripts/build-registry.sh
```

`phi-site/scripts/build-registry.sh` 会在临时目录中调用构建器生成新的无签名注册表，校验 `index.json` 后通过 `carry-registry-assets.mjs` 为仍在上一代索引中的旧版本保留不可变归档，再原子替换服务目录；若构建失败则保持当前在线服务目录不变。

On the production server, `~/bin/sync-phi-packages.sh` pulls `phi-packages`, rebuilds the catalog and registry via `phi-site/scripts/build-registry.sh`, carries forward previous-generation immutable assets if needed, and atomically swaps the served directory.

## 本地构建与校验 / Local unsigned registry build & verification

本地离线构建标准无签名注册表（需要并排放置的 `Phi` 仓库依赖与 Node.js 24，不需要任何私钥或签名配置）：

```sh
node scripts/prepare-catalog.mjs \
  --phi ../Phi \
  --out /tmp/phi-registry
```

也可在 GitHub Actions 中手动触发 **Verify registry build**（`.github/workflows/publish-catalog.yml`）验证无签名目录的完整构建与 `SHA-256` / 大小一致性（若需要覆盖已有不可变版本，可可选传入 `--wrapper-version`）。

Run `node scripts/prepare-catalog.mjs --phi ../Phi --out /tmp/phi-registry` locally to build the unsigned registry directory without any signing key. No Ed25519 private key or `index.sig.json` file is required.

## 安装路径 / Installation paths

| 内容 / Content | 安装位置 / Installed location |
| --- | --- |
| Plugin | `~/.phi/packages/plugin/<id>/<version>/` |
| Skill | `~/.phi/packages/skill/<id>/<version>/` |
| Connector | `~/.phi/packages/mcp/<id>/<version>/` |
| Wrapper | `~/.phi/wrappers/tree/` 中保留依赖使用的目录结构 / preserved component tree |
| 注册表缓存 / Registry cache | `~/.phi/cache/registries/` |

图标始终跟随内容安装到其入口旁，目录中的 manifest 供安装前预览。Connector credentials、运行环境和用户配置保留在用户数据目录，不进入发布目录。

Icons install beside their content entry; manifest sidecars provide previews before installation. Connector credentials, installed runtimes, and user configuration never enter the registry.
