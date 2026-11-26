---
name: nsfc-proposal
description: Draft, structure, and review National Natural Science Foundation of China (NSFC) grant proposals across Youth (青年C类), General (面上), and Key (重点) programs, supporting both Mechanism-Driven (机制探索型) and Method/Algorithm-Driven (方法与算法构建型) tracks. Covers 科学问题属性, <= 400-char Chinese abstract (targeting 350-385 chars with Python len() verification), English abstract (<= 4000 chars), 立项依据, 研究内容/目标/关键科学问题, 研究方案, 创新点, and 研究基础 (<= 5 representative works). Use for 国自然, NSFC标书, 基金申请书, 青年基金, 面上项目, 基金评审. Not for journal manuscripts (use manuscript-writing) or standalone literature reviews (use literature-review).
license: MIT
metadata: {"version": "1.1.0", "category": "文献与基金", "license": "MIT"}
---

# NSFC Grant Proposal Writing & Review (国自然标书写作与评审)

Design, draft, and review high-scoring NSFC proposals (青年科学基金 C 类 / 面上项目 / 重点项目) with sharp **scientific-question framing, strict logic chain, verified preliminary data, and zero fabricated citations or results**.

## 1. Scope & Boundary

- **Use for**: Framing core scientific hypotheses, **科学问题属性** (4 categories), `<= 400`-char Chinese abstract + `<= 4000`-char English abstract + 5 bilingual keywords, 立项依据, 研究内容/目标/关键科学问题, 研究方案, 可行性, 创新点, and 研究基础 (`<= 5` representative works).
- **Not for (Handoff)**:
  - Writing journal manuscripts -> `manuscript-writing`.
  - Multi-paper systematic literature review -> `literature-review`.
  - Polishing English prose only -> `manuscript-polishing`.

## 2. Epistemic Contract & Stop Gates

1. **Zero Fabricated Preliminary Data or Publications**: Never invent pilot `n`, fold-changes, `P`-values, or applicant papers. Tag missing evidence `[待补充前期数据: <具体指标/图表>]` and unverified literature `[待核验文献: <方向>]`.
2. **400-Char Chinese Abstract Verification Gate**: NSFC counts Chinese characters + punctuation `<= 400`. Always target **350–385 chars**, and verify exact length with Python `len()` when available (otherwise warn author to verify in the NSFC portal).
3. **Input Stop Gate**: If program tier (`青年 C 类` / `面上` / `重点`), research track (`机制探索型` vs `方法/算法构建型`), or core idea/preliminary data are missing, **STOP and ask for these 3 inputs** (or emit a labeled structural skeleton with `[待填]` slots).

## 3. Phi Runtime & Context Routing

- **Workspace Grounding**: Check workspace CVs, preliminary figures/tables, or `.phi-artifact.json` to ground 研究基础 in real evidence. Use `lit_search`/`pubmed` if enabled; otherwise mark `[待核验文献]`.
- **On-Demand Track References** (load ONLY the track matching the user's proposal):
  - **机制探索型**（湿实验 / 临床 / 生命医学） -> read [references/mechanism-track.md](references/mechanism-track.md).
  - **方法 / 算法构建型**（生信 / AI / 工程 / 交叉） -> read [references/method-track.md](references/method-track.md).

## 4. Core Workflow & Program-Tier Matrix

| Dimension | 青年科学基金（C 类） | 面上项目 | 重点项目 |
| --- | --- | --- | --- |
| **Aim Count** | **2–3 个聚焦子课题**（求深忌贪大） | **3–4 个递进子课题**（完整闭环） | **3–4 个大模块 + 协同验证** |
| **Preliminary Data** | 切入点苗头（1–2 组核心图） | 扎实直接前期数据（2–4 组图） | 系统性积累 + 平台/队列先发优势 |
| **Representative Works** | `<= 5` 篇（突出独立衔接） | `<= 5` 篇（突出连续深耕） | `<= 5` 篇（突出高影响引领） |

1. **Lock Tier, Track & 科学问题属性**: Choose 1 of the 4 attributes (*鼓励探索、突出原创* / *聚焦前沿、独辟蹊径* / *需求牵引、突破瓶颈* / *共性导向、交叉融通*) + 150–200-char rationale; load the matching track reference.
2. **Align Logic Spine**: `立项依据（痛点 + 前期证据 -> 假说） -> 关键科学问题（2–3 个本质难点，非操作步骤） -> 研究内容 -> 技术路线 -> 特色与创新 -> 研究基础`.
3. **Draft / Audit & Stress-Test**: Verify Chinese abstract `350–385` chars via `len()`, English abstract `<= 4000` chars, and eliminate slogan/AI prose (*"填补空白"*, *"开辟崭新道路"*).

## 5. Output Contract

1. **Proposal Alignment Card** (Program Tier, Research Track, 科学问题属性, Core Scientific Question).
2. **Drafted Section(s) or Scored 5-Dimension Review Report** (with exact character count for abstracts).
3. **Consistency & Gap Manifest** (`假说 <-> 内容 <-> 关键科学问题 <-> 创新点` alignment + `[待补充前期数据]` / `[待核验文献]` list).
