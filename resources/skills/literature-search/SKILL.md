---
name: literature-search
description: Search biomedical literature and keep a sourced note. Use when the user asks for papers, preprints, citations, PubMed, PMC, bioRxiv, medRxiv, or ClinicalTrials.gov, or wants an evidence summary with identifiers. Uses the installed pubmed, biorxiv, and clinical-trials connectors and does not invent citations. Also trigger on 文献, 论文, and 预印本.
allowed-tools: Read Bash
metadata: {"version": "1.0.0", "skill-author": "Phi Packages"}
---

# Literature search

## Overview

Answer a literature question from connector results. Every paper or trial in the answer needs an identifier that came back from a tool in this session. If a connector is not installed or returns nothing, say that and stop. Do not fill the gap from memory.

Source choice and the note shape are in [references/sources.md](references/sources.md). A blank note is in [assets/note-template.md](assets/note-template.md).

## When to use this skill

Use it for papers, preprints, citations, PubMed, PMC, bioRxiv, medRxiv, and clinical-trial records.

Use a different skill when the question is not a literature search:

- Target, disease, or drug associations: the `open-targets` connector.
- Cancer genomics studies already in cBioPortal: the `cbioportal` connector.
- Turning a downloaded PDF or Office file into Markdown: `markitdown`.
- Writing methods or results from the user's own analysis: `scientific-writing`.

## Procedure

1. Restate the question as a search: organism or disease, molecule or gene, method, and date bounds. Ask only when a missing bound would change the source.
2. Pick the source from `references/sources.md`. Published papers start at `pubmed`. Unpublished manuscripts start at `biorxiv`. Intervention studies start at `clinical-trials`.
3. Use the tools exposed by that installed connector. Read the tool schema before calling it. Do not invent tool names.
4. Keep the identifiers the tool returned (PMID, PMCID, DOI, NCT number, preprint DOI). Drop a hit you cannot tie to an identifier from the response.
5. Write the note from `assets/note-template.md`. Quote only short phrases that appear in the returned title or abstract. Paraphrase the rest.
6. If the user wants the full text, prefer a PMC full text the `pubmed` connector already returned. Otherwise pass a local PDF through `markitdown`. Do not claim a paper was read when only the abstract was returned.

## Report

Lead with the answer to the question, then the short list of sources with identifiers. State which connectors were called and which were unavailable. Separate "reported by the abstract" from "not in the retrieved text".
