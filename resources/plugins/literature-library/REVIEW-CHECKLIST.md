# literature-library 进仓检查清单

首次打包落地：`resources/plugins/literature-library/`，版本 **0.6.0**（已去掉 `-draft`）。Live `demo-*` 目录不在本包内。

## 评审门禁

- [x] 用户已授权进仓（「解决这些问题」覆盖 deferred 进仓）
- [x] 本树是 phi-packages 中的首次落地，不是本地草案硬阻断
- [x] 未把 hypothesis / scientific-gate / omics 塞进本插件
- [x] 未合并 research-evidence（桥接测试在兄弟插件缺席时跳过）
- [x] **无 embedding** 依赖或代码路径
- [x] TextIn **未**集成为 KB 骨架
- [x] 无商业文献网关、CNKI、Google Scholar 爬取、付费墙绕过
- [x] 插件内不打包 API key，也不内嵌 MCP server（仅可选 MCP client）

## 包元数据

- [x] `id: literature-library`，`toolPrefix: lit`，`type: plugin`
- [x] 无私有 `environments:`；技能用 `phi:python@1`
- [x] `components.skills` 列全 skill 目录（含 `lit-institutional` 与全部 wiki-*）
- [x] semver：`0.6.0`
- [x] `toolPrefix: lit` 与本仓库另一插件 `viz` 及独立 skill 前缀不冲突

## Agent / Skills

- [x] 脚本均在 `skills/<name>/scripts/`；stdout 单 JSON
- [x] live：collection / run-config / ingest / enrich / query / export / db
- [x] search：in-plugin HTTP live（pubmed/…/medrxiv + pmc/doaj/datacite/zenodo；S2 soft-fail）；不内嵌 MCP server
- [x] `lit-expand` live HTTP；默认不写库；`--confirm_write` 入库（0.4.4）
- [x] Wiki：init / ingest / lint / query / suggest-tree / apply-tree（`--confirm`）
- [x] collection create 默认 wiki-init（`--no_wiki` 可关）
- [x] 可选 MCP *client*（`providers_mcp` / `--backend mcp`）；HTTP 仍默认；宿主需提供 `LIT_MCP_CMD`
- [x] research-evidence `literature-search` 可 `--delegate_lit` / `--hits_jsonl`（插件不合并）
- [x] `lit-ingest` 含 `watch_rebind.py` 与 `pdf_extract.py`（OCR soft-fail）
- [x] `lit-institutional`：本人账号脚手架

## 测试

- [x] `tests/test_dedup.py`（无网络）
- [x] `tests/run_smoke.sh`（mock search + bib ingest）
- [x] `tests/run_smoke_v2.sh`（wiki vault + suggest/apply）
- [x] 本仓库没有宿主 `validatePlugin`；按 `docs/contracts/plugin.md` 做了清单校验（summary ≤ 300、skills 列全、脚本声明）
- [ ] 联网冒烟：真实 pubmed HTTP（可选，不在离线套件内）

- [x] lit-enrich：PMID/DOI HTTP + mock；ids ingest 默认 enrich（0.4.2）
- [x] Wiki agent-compile（无 plugin API keys；0.4.0+）

## 仍开口（进仓后可迭代）

- [x] 可选 MCP client（0.6.0）；插件仍不提供 MCP server；宿主 `LIT_MCP_CMD` 待生产接线
- [x] 种子扩展确认写入（`lit-expand` + `--confirm_write`，0.4.4）
- [x] PDF `pdftotext` + 可选 OCR soft-fail（0.6.0）；TextIn 非默认 / 非 KB
- [x] path-watch rebind（0.6.0）
- [x] institutional_access 脚手架（本人账号；0.6.0）；真实 SSO 浸泡仍开放
- 跨项目全局库挂载
- [x] 跨文献 concept 深合并（contrib + agent/mock apply，0.4.5）

- [x] providers_http：openalex/crossref/arxiv/biorxiv/medrxiv live；S2 soft-fail（0.4.3）；+ pmc/doaj/datacite/zenodo（0.4.9）
- [x] `tests/test_providers_http.py` 离线（`demo-search-v043/` 不随包分发）
- [x] lit-expand HTTP + `tests/test_expand.py`（`demo-expand-v044/` 不随包分发）
- [x] concept deep merge + `tests/test_concept_merge.py`（`demo-concept-merge-v045/` 不随包分发）
- [x] human_queue + coverage + version chain；`tests/test_queue_coverage_versions.py`（0.4.7）；**不接**三方商业源
- [x] HTTP version detect (OpenAlex/Crossref, `--mock`) + wiki source repoint；`tests/test_versions_http_wiki.py`（0.4.8）
- [x] protocol freeze + lit-run-close gate + dig queue + preview_fetch；`tests/test_protocol_ops_v050.py`（0.5.0）
- [x] 文件入库打磨 + OA 软验证 + pub_files；`tests/test_file_ingest_v051.py`（0.5.1）
- [x] gold-set recall + OOW + uncontacted + expand saturation；`tests/test_recall_saturation_v052.py`（0.5.2）
- [x] UX polish: related-link resolve, --query, rebind_files, coverage hints, wiki_agent_status；`tests/test_ux_v053.py`（0.5.3）
- [x] deferred 五项：path-watch / PDF extract+OCR / institutional / MCP client / packaging；`tests/test_deferred_v060.py`（0.6.0）

| 项 | 值 |
| --- | --- |
| 插件 id | literature-library |
| toolPrefix | lit |
| 版本 | 0.6.0 |
| 知识层 | SQLite + LLM Wiki（无 embedding） |
