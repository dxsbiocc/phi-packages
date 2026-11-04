---
name: lit-ingest
description: >-
  Ingest literature into the project SQLite library from hits.jsonl, BibTeX,
  RIS, DOI/PMID lists, PDF text sidecars, or authorized local PDFs. Digital
  PDFs use pdftotext when available; optional OCR soft-fail (tesseract /
  LIT_PDF_OCR_CMD). Deduplicates before insert. Records pub_files (path +
  content_hash). Optional soft OA URL verify. Path-watch rebind via
  watch_rebind.py. TextIn is not a required parser / not the KB backbone.
phi:
  environment: phi:python@1
  attachTo:
    - LiteratureLibrary
  scripts:
    - name: ingest
      description: Parse input file and upsert into library.sqlite with dedup.
      run: [python, ./scripts/ingest.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          input_path:
            type: string
            format: input-path
          format:
            type: string
            description: hits|bib|ris|ids|pdf_text|pdf (optional auto-detect)
          run_id:
            type: string
          no_enrich:
            type: boolean
            description: Skip auto-enrich after local/ids ingest (default false)
          enrich_mock:
            type: boolean
            description: Use offline mock fixture when auto-enriching
          verify_oa:
            type: boolean
            description: Soft-verify candidate OA/PDF URLs; never store PDF bytes
          oa_mock_fixture:
            type: string
            format: input-path
            description: JSON map url→mock HEAD result for offline OA verify
        required:
          - project_root
          - collection_id
          - input_path
      approval: write
      output: ./schemas/ingest-result.json
      timeoutSeconds: 120
    - name: rebind_files
      description: >-
        Rebind pub_files.file_path by content_hash after files move (dry-run
        default; --confirm writes).
      run: [python, ./scripts/rebind_files.py]
      args:
        type: object
        additionalProperties: false
        properties:
          project_root:
            type: string
            format: project-path
          collection_id:
            type: string
          scan_dir:
            type: string
            format: input-path
          confirm:
            type: boolean
          sample_limit:
            type: integer
        required:
          - project_root
          - scan_dir
      approval: write
      timeoutSeconds: 300
---

# lit-ingest (0.5.3)

入库前去重（见 `lit-db`）。文件源：`.jsonl` hits、`.bib`、`.ris`、DOI/PMID 清单、PDF **文本** sidecar、授权本地 `.pdf`（只认文件名/sidecar，**不**解析二进制）。

**Local files (0.5.1):**
- 从 sidecar 正文 / 文件名提取 DOI、PMID、PMCID、arXiv
- 写入 `pub_files`（`file_path` + `content_hash` + `ingest_origin=local_file`）；不把未验证 OA 链接当已有 PDF
- 空壳（缺 title/abstract）默认自动 `lit-enrich`（`--no_enrich` 可关；`--enrich_mock` 离线）
- `--verify_oa`：对候选 URL 做软 HEAD/GET 分类 `verified_oa|unknown|paywall_or_blocked|not_pdf`，结果进 `publications.oa_status` / `oa_json`；**从不**把付费墙正文当库内全文

```bash
python scripts/ingest.py \
  --project_root ./demo-project \
  --collection_id col-... \
  --input_path tests/fixtures/pdf_sidecar_zheng.txt \
  --format pdf_text --verify_oa --oa_mock_fixture tests/fixtures/oa_mock_head.json

python scripts/ingest.py \
  --project_root ./demo-project \
  --collection_id col-... \
  --input_path ./papers/10.1016_j.cell.2017.05.035.pdf --format pdf
```

**Path drift (0.5.3):** `rebind_files.py` scans a directory, matches `content_hash`, updates `file_path` (dry-run default; `--confirm` writes).

```bash
python scripts/rebind_files.py --project_root P --scan_dir ./papers
python scripts/rebind_files.py --project_root P --collection_id C --scan_dir ./papers --confirm
```

工具名：`lit_ingest` / `lit_rebind_files`。
