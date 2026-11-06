# literature-library DESIGN（v0.6.0）

首次进入 `phi-packages` 的打包版本（`resources/plugins/literature-library/`）。`-draft` 后缀已去掉。Live `demo-*` 目录不随包分发。biology / scientific claims：N/A。

## 已定决策（Q1–Q12 + V2）

| Q | 决策 |
|--|--|
| Q1 | 独立插件 `literature-library`；生产可走 `phi-literature` MCP client（非 Claude）；**草案 0.4.5=插件内 HTTP CLI（search+expand）+ concept 深合并，不内嵌 MCP server** |
| Q2 | 默认**项目级** SQLite（`library.sqlite`） |
| Q3 | v1 源 = MCP 已支持集（pubmed…biorxiv/medrxiv）；不做 Google Scholar |
| Q4′ | 独立插件含文件入库 + 去重 + 知识库；hypothesis/gate **不**并入 |
| Q5′ | v1 知识库 = 可查询 SQLite（题录+摘要+标签+run）；**不做 embedding** |
| Q6′ | 文件源：PDF 文本 sidecar / .bib / RIS / DOI-PMID 清单 |
| Q7 | **collection（课题）** 为一级容器；run/文件挂在下面 |
| Q8 | 种子扩展可选有限额 → **live HTTP**（`lit-expand`：similar/citations/references；默认预览，`--confirm_write` 才入库）；默认仍人发起 search/文件入库 |
| Q9 | 检索 = SQL + FTS（title/abstract/tags）+ wiki markdown 关键词；**无 embedding** |
| Q10 | V2：**主题树（LLM Wiki）为主**；embedding **砍掉不做** |
| Q11 | 建议枝需确认后写入（`lit-wiki-suggest-tree` → `lit-wiki-apply-tree --confirm`） |
| V2 | 知识结构 = **LLM Wiki 三层**；TextIn 文档树 **不是** KB 骨架 |

另：search **自动** ingest；软冲突 → `needs_review`；**真**强 ID 冲突 → `needs_review`；互补强 ID（空壳 doi + 空壳 pmid + 同篇联合命中 / `#conflict-` 簇且标题相容）→ 自动合并到 DOI 优先规范主键（`merge_ids_cli` 可 dry-run / `--apply` 修复旧库）。

## Why heuristic (v0.2) then LLM (v0.3)

Heuristics enabled offline smoke without API keys. Intended product path is LLM Wiki synthesis; v0.3 flips the default.

## V2：LLM Wiki 层（相对 SQLite）

SQLite 仍管 **文献身份 / 去重 / runs**。Wiki 是并列的**编译知识层**（Karpathy/llmwikis 风格）：

| 层 | 路径 | 可变性 |
|--|--|--|
| raw | `collections/<id>/wiki-vault/raw/` | **不可改**（只追加） |
| wiki | `…/wiki-vault/wiki/` | 可编译 markdown（sources / concepts / entities / index / log） |
| schema | `…/wiki-vault/AGENTS.md` | 约定文件 |

- **Collection ≈ 一个 wiki vault 根**
- **主题树** = `wiki/index.md` 标题目录 + 可选 `theme_nodes` / `theme_pubs`（文献可多挂）
- **无 embedding / 向量库**
- **TextIn 文档树**：单篇 PDF 内部标题解析用（可选、外发、以后再接）；**不**当课题知识库骨架

### 对比 TextIn

| | LLM Wiki | TextIn 文档树 |
|--|--|--|
| 粒度 | 跨文献的概念 / 主题 / 实体页 | 单 PDF 内部标题层级 |
| 对本插件 | 课题知识库主路径 | 最多当 ingest 可选解析器（未集成） |

### Skills（Wiki）

0.6.1：下表是 `lit-wiki` 技能里的工具行为，不再是各自一个 skill 目录。共享实现在 `skills/lit-wiki/scripts/wiki_lib.py` 与 `wiki_llm.py`。检索 / 入库 / 去重 / coverage / 版本链在技能 `literature-library`。

| 工具（原 skill 名） | 行为 |
|--|--|
| `lit_wiki_init` | 建 vault 骨架 |
| `lit_wiki_ingest` | raw 快照 + agent/mock/heuristic 编译 sources（Summary/claims/methods）+ concepts/entities + 空 tags 自动回写 SQLite + index/log |
| `lit_wiki_lint` | 孤儿页 / 缺 index / 断链 / needs_review 横幅 |
| `lit_wiki_query` | wiki 关键词 + SQLite FTS（无向量） |
| `lit_wiki_suggest_tree` | 提案 JSON only |
| `lit_wiki_apply_tree` | 仅 `--confirm` 时写入 |

`lit_collection` create **默认**调用 wiki-init（`no_wiki` 可关）。

## 边界

| 组件 | 职责 |
|--|--|
| **literature-library（本插件）** | 找得到、存得住、不重复、Wiki 编译层、可导出 |
| **phi-literature MCP** | 联网多源检索 / related / fulltext |
| **research-evidence** | 证据表、citation-verify、假设卡、scientific-gate；消费本库 export |
| **paper-repro** | 沉积复现；可读文献库但不负责建库 |

## Schema 概览（SQLite，schema_version=6）

- `collections` — 课题
- `publications` — 规范主键 `doi:…` / `pmid:…` / …；v5 起含 `oa_status` / `oa_json`（候选 OA 链接软验证结果）
- `pub_files` — 授权本地文件路径 + `content_hash` + `ingest_origin` + `access_origin` / `access_method`（v6；不存 PDF 字节）
- `aliases` — 所有规范化 ID → publication
- `collection_pubs` — 多对多 + `source_run_id`
- `runs` / `run_hits` — 检索可追溯（含 `dedupe_action`；v4：`status` / `coverage_path` / `gate_json`）
- `theme_nodes` / `theme_pubs` — 树节点；V2 apply-tree 可写入
- `pubs_fts` — FTS5（失败则 LIKE）
- `publication_versions` — preprint↔journal（`preprint_of` / `published_as` / `same_work`；`needs_confirm`）
- `publications.review_reasons_json` — 可选；人审队列可挖掘的原因码

## Publication primary key（去重 / 规范主键）

`publications.id` 是**自然语言主键**（不是 UUID），按强 ID 优先级选择/升格：

| 优先级 | 形式 | 说明 |
|--|--|--|
| 1 | `doi:10.…` | 最高；bioRxiv/medRxiv 有 DOI（`10.1101/…`）时走 DOI，无独立命名空间 |
| 2 | `pmid:…` | |
| 3 | `pmcid:PMC…` | |
| 4 | `arxiv:…` | |
| 5 | `openalex:W…` | |
| 6 | `soft:<16hex>` | title+first author+year 指纹；**永不自动合并** → `needs_review` |
| 7 | `tmp:<12hex>` | 无任何身份时的占位 |

**身份唯一性**

- `aliases.alias` PRIMARY KEY（所有规范化 ID → 一篇 publication）
- schema v2：`publications.doi/pmid/pmcid/arxiv/openalex` 上 **partial UNIQUE INDEX**（排除 `id` 含 `#conflict-` 的审阅行；NULL 可多条）
- Semantic Scholar / S2 **不进** `STRONG_ID_KEYS`

**合并策略（不变 + 升格补丁）**

- 软指纹冲突 → `needs_review`，不自动合并
- 真强 ID 冲突（同一 draft 指向两篇已有文献）→ `needs_review`，必要时 `#conflict-` 后缀
- 互补多命中簇（doi-only + pmid-only + 联合命中）→ 自动合并到 DOI 优先规范主键
- **Canonical upgrade（0.4.6）**：单篇强 ID 命中后若合并出更高优先级 ID（例：先 `pmid:28622514`，后到 DOI+PMID），`_reassign_publication_id` 把 `publications.id` 升格为 `doi:…`，并重挂 aliases / collection_pubs / run_hits / theme_pubs / FTS；返回 `action=merged` + `previous_id`

## 仍未接线 / 开口

- Phi host 生产 `LIT_MCP_CMD` 桥（插件内已有可选 MCP *client*；HTTP 仍是默认；**不**内嵌 MCP server；不打包 API key）
- TextIn PDF API（明确不做为本版依赖，也不是 KB 骨架）
- Concept merge 质量依赖宿主 agent 填 JSON（脚本不发明合成；`--mock` 仅离线）
- **不接**三方商业文献源（SmartLib / global-biblio-base 等）；不接 CNKI / Google Scholar 爬取；不做 embedding
- 校内通路只有本人账号脚手架（`lit-institutional`）；真实 SSO/代理浸泡不在本包

## Frozen search protocol（0.5.0）

每次 `lit-run-config` 写入 `runs/<run_id>/run.json`（+ `run.yaml`）并冻结：

| 字段 | 含义 |
|--|--|
| `protocol_version` | 当前 `0.5.2` |
| `query` | 检索式 |
| `sources` | 开放源列表（≤8；无商业源） |
| `date_from` / `date_to` | 时间窗 |
| `limit` / `limit_per_source` | 每源上限 |
| `provider_notes` | 各源适配器说明 |
| `frozen_at` / `created_at` | 冻结时间戳 |
| `gold_set` | 已核验强 ID 列表（优先 DOI；pmid/arxiv 亦可）；建议 5–10，&lt;5 → coverage `evidence_weak`（不硬失败） |
| `out_of_window_controls` | 1–2 篇窗外真论文 ID，用于证明日期过滤 |
| `no_commercial_sources` | 恒为 true |
| 金标准核验 | **调用方**应先验 DOI（doi.org）；插件离线仅做格式软检查 |

`lit-search` 执行后追加 `executed_at`、`providers`、`sources_executed`，并写耐久 `hits.jsonl` + `hits.meta.json`。

`lit-run-close`（或 `lit-search --close_run`）写 `coverage.json`/`coverage.md`，按门槛标 `passed` / `needs_attention`：

- `--min_unique_pubs`
- `--max_needs_review` / `--max_version_pending` / `--max_unenriched_shells`
- `--fail_under_coverage` → CLI exit 2

## 环境

仅 `phi:python@1`。无私有 conda。

## 0.3.0 变更

- **Wiki compile (superseded by 0.4 agent default):** introduced structured compile schema + mock/file providers; HTTP path removed in 0.4.
- `lit-db`: complementary strong-ID auto-merge + `merge_ids_cli.py` repair
- `lit-ingest`: post-ingest high-confidence merge pass
- `lit-wiki-lib`: entity fragment suppression (no `l1`/`ifn` from PD-L1/IFN-γ); noisy_entity lint; layilin→LAYN

## 0.4.0 变更

- **Wiki compile default = `agent`**: invoking agent's model fills structured JSON; scripts = prompt+schema+validate+write.
- **No plugin API keys** for normal use. Removed OpenAI-compatible HTTP client / `LIT_WIKI_LLM_API_KEY` requirement.
- Two-step UX: `wiki_compile_request` → agent synthesizes → `wiki_compile_apply` / `--llm-json`.
- Keep `mock` (smoke only) and `heuristic` (explicit emergency fallback).
- Preserve ID-merge + entity fragment guards from 0.2.2 / 0.3.

## 0.4.1 变更

- **Ids-file co-identity:** `doi:` + `pmid:` on consecutive lines, same line, or `DOI\tPMID` parse as **one** publication (DOI-preferred canonical) even as metadata-empty shells.
- Blank-line groups and alternating dual-ID sequences supported; block layout (all DOIs then all PMIDs) is **not** auto-paired.
- `--request-out` alias for wiki compile request path.

## 0.4.2 变更

- **lit-enrich:** ID-only shells → real title/abstract/authors/year/venue via PubMed efetch (PMID) and PubMed/Europe PMC (DOI); richest-merge; `--mock` fixture (Zheng).
- **Ids ingest default:** auto-enrich shells unless `--no_enrich` (optional `--enrich_mock`).
- **wiki compile_request:** warn + `prefer_enrich_first` when abstract/title still shell-like.
- providers_http: `parse_pubmed_xml`, `fetch_by_pmid`, `fetch_by_doi`.

## 0.4.3 变更

- **检索源 HTTP 扩展：** `providers_http.py` 接通 openalex / crossref / arxiv / biorxiv / medrxiv；`semantic_scholar`/`S2` 无 key 时 soft-fail。
- **架构：** 对外 lit-* CLI；插件内 HTTP 适配；**不**在插件内嵌第二套 MCP server；生产可选再接 MCP client。
- `search.py --sources` 可覆盖 run.json；`--mock` 保留。
- bioRxiv/medRxiv：无关键词 API → DOI detail / 日期浏览 / Europe PMC JOURNAL 回退（meta.note 说明）。
- 离线 `tests/test_providers_http.py`；联网 `demo-search-v043/`。

## 0.4.4 变更

- **lit-expand live：** PubMed elink（similar+PMID）、OpenAlex related/referenced/cites、OpenCitations DOI 回退；stdlib HTTP only。
- 默认预览；`--confirm_write` → `collections/<id>/expands/<exp-id>/hits.jsonl` + ingest（合成 runs 行）；`--mock` / `--dry-run` / `--skip_ingest`。
- providers_http 共享 helpers：`expand_related` / `pubmed_neighbors` / `hydrate_openalex_ids` / …
- schema `expand-result.json` 去掉 `stub: const true`；Q8 不再 stub。
- 离线 `tests/test_expand.py`；联网 `demo-expand-v044/`。

## 0.4.5 变更

- **跨文献 concept 深合并：** ingest 写 `wiki/concepts/_contrib/<slug>/<safe_id>.json`；≥2 源时 `concept_merge_pending` + soft stub。
- CLI：`concept_merge_request.py` / `concept_merge_apply.py`（`--mock` 或 `--llm-json`）；重写 Definition + ## Across papers，保留 ## Sources。
- `needs_review` 贡献强制 banner；可选 `canonical_slug`/`merge_into` 重定向旧页。
- compile schema：concepts 可带 `evidence`/`note`（字符串仍兼容）。
- 无 MCP、无新依赖、无插件 API key。离线 `tests/test_concept_merge.py`；`demo-concept-merge-v045/`。

## 0.4.6 变更

- **Publication PK 合同锁定：** `publications.id` = 自然键（`doi:` > `pmid:` > `pmcid:` > `arxiv:` > `openalex:` > `soft:` > `tmp:`）；schema_version → 2；强 ID 列 UNIQUE INDEX。
- **Canonical upgrade 补丁：** 先 PMID（或 PMCID/arXiv）入库、后到更高优先级 ID 时升格主键并重挂依赖；软指纹仍不自动合并。
- 文档：DESIGN「Publication primary key」；README dedup 条；离线 `tests/test_dedup.py` 升格用例。

## 0.4.7 变更

- **lit-review-queue：** 身份人审队列（`soft_fingerprint_conflict` / `conflicting_strong_ids` / `empty_unenriched_shell` / `pending_version_link`）；JSON + 可选 markdown；可按 `collection_id` 过滤；空队列友好。
- **lit-coverage：** 按 `run_id` / `collection_id` 报告 hits vs unique pubs、`dedupe_action` 计数、needs_review、缺强 ID、空壳无摘要。
- **publication_versions + lit-versions：** preprint↔journal 链接；detect / apply / propose / confirm；高置信可合并到期刊 DOI 规范主键（预印本 DOI/arXiv 保留为 aliases）。
- schema_version → **3**；`review_reasons_json`；不接三方商业源。

## 0.4.8 变更

- **HTTP version detect：** `lit-versions --http` 可选查询 OpenAlex（locations 多 DOI）与 Crossref（`relation`：`is-preprint-of` / `has-preprint` 等）；`--mock` + fixtures 离线；网络错误 soft-fail。
- 高置信且库内双侧 DOI 已存在 → 与离线启发式一并 propose / `--auto_merge_high_confidence`（期刊 DOI 主键；预印本 ID 为 aliases）。
- **wiki_repoint：** 主键升格 / 版本合并 / `merge_ids --apply` 后，将 `wiki/sources/<old_safe>.md` 与 `raw/publications/<old_safe>/` 迁到新 slug，更新 index/log 链接，旧路径留 stub；共享 helper 在 `lit-wiki-lib`。
- schema 仍为 **v3**；不接商业源；无 MCP；无新依赖。

## 0.4.9 变更

- **开放检索源扩展（无商业网关）：** `providers_http.py` 新增 **pmc**（NCBI E-utilities `db=pmc` esearch+esummary）、**doaj**（DOAJ `/api/search/articles`）、**datacite**（DataCite REST `/dois`）、**zenodo**（Zenodo/Invenio `/api/records`）。
- 统一 hit 形状；错误 soft-fail（空 hits + `meta.error`）；礼貌 User-Agent / mailto（复用 `OPENALEX_MAILTO` / `CONTACT_EMAIL`）。
- `SOURCE_ALIASES` / `lit-run-config` ALLOWED_SOURCES 同步；单 run 源上限 8（草案，非 MCP）。
- 离线 fixtures + `tests/test_providers_http.py`；不接 CORE/Scopus/WoS/SmartLib/CNKI；S2 仍可选 soft-fail；无 MCP；无新 pip 依赖。

## 0.5.0 变更

- **Protocol freeze + coverage gate：** `lit-run-config` 冻结协议字段；`lit-search` 写 `hits.meta.json`；新技能 **`lit-run-close`**（亦可 `search --close_run`）写 coverage 并标记 `runs.status`。
- **Human dig：** `lit-review-queue --dig` 空 backlog 静默 exit 0；非空带 `suggested_next_action`；`--write_report` 仅非空写入。
- **Version preview_fetch：** 对侧 DOI 不在库时，`--preview_fetch` 仅拉元数据建 shell + `pending_version_link`（不下载 PDF）；`--action confirm` 合并到期刊 DOI。
- schema_version → **4**（`runs.status` / `coverage_path` / `gate_json`）；不接商业源；无 MCP；无新 pip 依赖。


## 0.5.1 变更

- **文件入库打磨：** sidecar / 文件名 DOI·PMID·PMCID·arXiv 提取；`.pdf` 仅走 sidecar 或文件名壳（不解析二进制）。
- **`pub_files`：** 记录授权本地 `file_path` + `content_hash` + `ingest_origin=local_file`；合并/升格时随主键重挂。
- **自动 enrich：** 本地文件空壳（缺 title/abstract）与 ids 一样默认调用 `lit-enrich`（可 `--no_enrich` / `--enrich_mock`）。
- **OA 软验证：** `oa_verify.py` + ingest/enrich `--verify_oa`；状态 `verified_oa|unknown|paywall_or_blocked|not_pdf` 写入 `oa_status`/`oa_json`；**从不**把未验证/付费墙正文当库内全文。
- schema_version → **5**；不接商业源；无 MCP；无新 pip 依赖；无 CNKI/Scholar 爬取。


## Recall evidence + saturation honesty（0.5.2）

Inspired conceptually by open research-skills lit-search practice (gold-set coverage, out-of-window controls, declared-but-uncontacted sources, snowball saturation) — **not** vendored code; kept lightweight (no PRISMA machinery).

| 机制 | 行为 |
|--|--|
| **gold_set recall** | coverage 报告 found/missed + `evidence_weak`（空或 &lt;5） |
| **out_of_window_controls** | 窗外对照若出现在窗内 hit corpus → leakage |
| **declared-but-uncontacted** | 配置源但未接触/失败/跳过必须列出原因；空白节≠无缺口 |
| **snowball saturation** | `lit-expand --max_rounds`；连续 2 轮 new-included &lt; 5% seeds **或**无新种子 → `saturated`；触顶未饱和 → `stop_reason=max_rounds_not_saturated`，**禁止**假称闭包 |

## 0.5.2 变更

- **Run protocol：** `gold_set` / `out_of_window_controls` 冻结进 run.json；`protocol_version` → `0.5.2`。
- **lit-coverage / lit-run-close：** recall_evidence（gold recall、OOW leakage、uncontacted sources）；`evidence_weak`。
- **lit-expand：** `--max_rounds` / `--low_yield_ratio` / `--seeds_per_round`；输出 `saturated` / `stop_reason` / `rounds`。
- schema 仍为 **v5**；不接商业源；无 MCP；无新 pip 依赖；离线 `tests/test_recall_saturation_v052.py`。

## UX polish（0.5.3）

Real acceptance friction → fixes (stdlib only; no commercial sources):

| Pain | Fix |
|--|--|
| Merge `related_concepts` wrote broken links (`nsclc.md`) | Resolve to existing `wiki/concepts/*.md` via exact / slugify / aliases / light fuzzy; drop unresolved → `related_unresolved` |
| CLI `--query` missing | `lit-query` / `lit-wiki-query` accept `--query` as alias of `--q` |
| Local PDF path drift | `rebind_files.py` matches `content_hash`, dry-run default, `--confirm` writes |
| Hard to see next wiki step | `wiki_agent_status.py` lists missing sources + pending merges; optional `--emit_requests` |
| Low gold recall silent | coverage `hints[]` (increase limit / seed ids / expand / check query) when relevant |
| Speculative entities | Soft-drop entity rows whose name/note contains "implied" / "not claimed" |
| Wiki query semantics unclear | Help text documents AND for multi-term queries |

**Default wiki path = agent compile** (request → invoking model JSON → apply). `--compile mock` is offline CI only — do not use mock for live acceptance.


## Landed in 0.6.0

| Item | Design |
|--|--|
| Path watch | `watch_rebind.py` polls `--scan_dir` (stdlib); reuses `rebind_pub_files`; dry-run default |
| PDF extract | `pdftotext` digital text; optional OCR soft-fail (`tesseract` / `LIT_PDF_OCR_CMD`); TextIn not KB / not required |
| Institutional | `lit-institutional` + `--institutional_access`; own VPN/proxy/session/file only; vault credentials never in DB |
| MCP client | optional `--backend mcp` → `providers_mcp` → phi-literature tools; HTTP default; no MCP server in plugin |
| Packaging | version `0.6.0` in this repository (first package landing; `-draft` dropped). `demo-*` trees are not included |

schema_version → **6**. Stdlib-first; no new pip deps.
