# literature-library（文献库）· 0.7.0

Project-level literature library for Phi: search open sources, import your own files, deduplicate,
export citations, and compile saved papers into a wiki. No API keys, no embeddings.

项目级文献库：检索开放来源、导入自己的文件、自动去重、导出参考文献，并把已入库的论文整理成 wiki。

## 组成

| 组件 | 作用 |
|--|--|
| Agent `LiteratureLibrary` | 角色与边界，约 40 行；工作流都在技能里 |
| Skill `manage-library` | 建库、检索、导入、扩展、查询、导出、维护（9 个工具） |
| Skill `compile-wiki` | 把论文写成 wiki：摘要、要点、概念页、主题目录（用 `manage-library` 里的 3 个 wiki 工具） |

`manage-library` 自带全部脚本，`compile-wiki` 只有说明文档，所以两个技能之间没有代码依赖。

## 工具（前缀 `lit_`）

| 工具 | 用途 |
|--|--|
| `lit_status` | 库的概况和下一步建议 |
| `lit_collection` | 新建或列出课题（collection） |
| `lit_search` | 检索**外部**来源并去重入库 |
| `lit_import` | 导入 `.bib`、`.ris`、ID 清单、PDF 或整个文件夹 |
| `lit_expand` | 从一篇种子论文找相似、引用它、被它引用的论文（默认只预览） |
| `lit_query` | 查**库内**的论文和 wiki 页面 |
| `lit_get` | 读一篇论文的完整记录、文件、版本、wiki 页，可带 PDF 文本 |
| `lit_export` | 导出 BibTeX、RIS、CSL-JSON、CSV 或检索记录 |
| `lit_curate` | 待人审队列、补全元数据、预印本与正式版关联、合并重复、文件改路径后重绑 |
| `lit_wiki_pending` | 列出还没有 wiki 页的论文并生成一批编译请求 |
| `lit_wiki_apply` | 校验 agent 写的结果并落盘 |
| `lit_wiki_lint` | 体检 wiki |

## 运行要求

- 环境 `phi:python@1`；脚本只用 Python 标准库，没有第三方依赖。联网仅用于 `lit_search`、`lit_expand` 和元数据补全。
- 本地 PDF 读文字需要 `pdftotext`（poppler）；没有文字层时可选 `tesseract` 或 `LIT_PDF_OCR_CMD`，缺失时只给出警告。
- 可选环境变量：`OPENALEX_MAILTO`、`CROSSREF_MAILTO`、`CONTACT_EMAIL`（礼貌性联系邮箱）。

## 测试

```bash
python3 -m pytest resources/plugins/literature-library/tests -q
```

全部离线。`tests/test_plugin_contract.py` 检查清单、技能前置信息，以及每个工具参数都是对应脚本的真实参数。

## 边界

不做引文核查、证据评分、假设或组学分析；不接 Google Scholar、CNKI 和商业文献网关；不内嵌 MCP 服务，
也不打包任何密钥。设计取舍见 [DESIGN.md](DESIGN.md)，版本变化见 [CHANGELOG.md](CHANGELOG.md)。
