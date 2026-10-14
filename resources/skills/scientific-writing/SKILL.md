---
name: scientific-writing
description: Draft methods, results, and figure captions from completed analysis outputs. Use when the user asks for a paper section, methods writeup, results paragraph, figure or table legend, or a plain-language summary of results. Uses only numbers and claims that appear in the project files. Also trigger on 方法学, 结果, and 图注.
allowed-tools: Read Write Edit Bash
metadata: {"version": "1.0.0", "skill-author": "Phi Packages"}
---

# Scientific writing

## Overview

Write the prose the user asked for from files in the project. Numbers, sample sizes, test names, and software versions come from those files or from the skill that produced them. A missing number stays missing.

Use a Phi wrapper and its packaged environment when the analysis has not been run yet. Use the agent's built-in search tools when the section needs citations. Use `seaborn`, `matplotlib`, or `omics-visualization` when the user wants a figure rather than a caption. Patterns for each section are in [references/section-patterns.md](references/section-patterns.md).

## When to use this skill

Use it for methods, results, figure or table legends, and plain-language summaries of a completed analysis.

Do not use it to invent a result, to run the analysis, or to search the literature unless the user also asked for citations.

## Procedure

1. Read the plan, the command that ran, and the output tables or figure files the user pointed at. If none are named, look for the newest analysis outputs in the project and say which files were used.
2. Extract the facts the section needs: sample size, groups, filters, software name and version, test, estimate, and the path of the figure or table. Keep a scratch list with the source path for each fact.
3. Draft only the section that was requested, following `references/section-patterns.md`.
4. Attach the source path in a short note under the draft, not in the prose, so the user can check the sentence against the file.
5. If a sentence needs a fact that is not in the files, write "[not in the outputs]" in place of the number. Do not estimate it.

## Style

- Past tense for methods that already ran. Present tense for a plan that has not run, and label it as a plan.
- Name the comparison (which groups, which direction) before the statistic.
- One claim per sentence. Keep abbreviations defined at first use in that section.
- Do not upgrade "associated" to "caused" unless the user supplied a design that supports that word and the outputs say so.

## Report

Give the draft first. Then list the files it is based on and any "[not in the outputs]" gaps. Do not describe the draft as a finished paper.
