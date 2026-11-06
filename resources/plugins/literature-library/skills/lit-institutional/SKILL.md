---
name: lit-institutional
description: >-
  Optional campus / SSO authorized full-text pathway. Own account only
  (VPN/proxy, session reuse, SAML host-concern, or user-provided authorized
  file). Explicit --institutional_access required. Credentials stay in the
  user vault — never stored in SQLite. No paywall bypass, no borrowed logins.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: institutional_access
      description: >-
        Mount an already-authorized local PDF, or fetch via LIT_INSTITUTIONAL_PROXY
        / ambient VPN (mockable). Records access_origin + access_method on pub_files.
      run: [python, ./scripts/institutional_access.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          institutional_access:
            type: boolean
            description: Required gate (must be true)
          method:
            type: string
            description: vpn_proxy|session_reuse|saml_browser|user_provided_file
          file_path:
            type: string
            format: input-path
          url:
            type: string
          dest_dir:
            type: string
            format: project-path
          doi:
            type: string
          pmid:
            type: string
          title:
            type: string
          mock_fixture:
            type: string
            format: input-path
          proxy:
            type: string
          no_ocr:
            type: boolean
          ocr_cmd:
            type: string
        required:
          - project_root
          - collection_id
          - institutional_access
          - method
      approval: write
      timeoutSeconds: 300
---

# lit-institutional

校内 / SSO **本人账号**授权全文入口（脚手架）。

## 硬规则

- 必须传 `--institutional_access`
- 只用用户**自己的**校园网 / VPN / 代理 / 已登录会话；**禁止**学长账号或共享凭据
- 凭据只在用户 vault；插件 / SQLite **不存**账号密码
- 全文只收：合法 OA，或用户已获授权的本地文件 / 经本人代理拉下的文件
- `saml_browser`：Playwright/SAML 由宿主 agent 完成；本脚本只记录 method，并提示改用 `user_provided_file`

## 方法

| method | 行为 |
|--|--|
| `user_provided_file` | 挂载已授权本地 PDF（推荐） |
| `vpn_proxy` | 经 `LIT_INSTITUTIONAL_PROXY` 拉 URL → 落盘 → 挂载 |
| `session_reuse` | 假设环境已有 VPN/会话；可直连或同代理 |
| `saml_browser` | 不启动浏览器；返回宿主侧指引 |

## 示例

```bash
# 已下载的授权 PDF
python skills/lit-institutional/scripts/institutional_access.py \
  --project_root P --collection_id C --institutional_access \
  --method user_provided_file --file_path ./papers/foo.pdf --doi 10.1234/x

# 离线 mock fetch
python skills/lit-institutional/scripts/institutional_access.py \
  --project_root P --collection_id C --institutional_access \
  --method vpn_proxy --url https://publisher.example/pdf \
  --mock_fixture ./tests/fixtures/tiny.pdf
```

元数据写入 `pub_files.access_origin` / `access_method`（schema v6）。
