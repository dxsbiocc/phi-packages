---
name: reproducible-analysis
description: Plan a scientific analysis that can be rerun. Use when the user wants an analysis, pipeline, notebook, or methods plan and the inputs, commands, environment, and outputs must be explicit. Chooses installed Phi skills such as scanpy, nextflow, scikit-learn, pysam, or rdkit instead of an improvised workflow. Also trigger on 可重复 and 分析方案.
allowed-tools: Read Write Edit Bash
metadata: {"version": "1.0.0", "skill-author": "Phi Packages"}
---

# Reproducible analysis

## Overview

Write a plan another person can run with the same files and the same Phi skills. Pick skills that already exist. Record inputs, parameters, environment, commands, and outputs before changing data.

Use `find-skills` to choose the skill. Use `scientific-writing` after the run, when the user wants methods or results prose. The plan outline is [references/plan-template.md](references/plan-template.md).

## When to use this skill

Use it when the user asks for an analysis, a pipeline, a notebook, or a methods plan, and the result needs to be rerunnable.

Do not use it to execute a Nextflow module design that `nextflow` already covers, or to draw a figure a plotting skill already covers. Those skills own the execution detail. This skill owns the plan that points at them.

## Procedure

1. State the question as a comparison or measurement, the sample unit, and the decision the result should support.
2. List the input files that exist in the project. Record format, path, and what each row or record is. Do not invent a path.
3. Call `find-skills` and name the skill that owns the analysis. Typical matches: `scanpy` for scRNA-seq, `scvelo` for velocity, `pysam` for alignment and variant files, `nextflow` for a multi-step pipeline, `scikit-learn` or `scikit-survival` for tabular models, `rdkit` for molecules, `seaborn` or `matplotlib` for figures, `omics-visualization` for template omics figures.
4. Write the plan from the template. Every step names the skill, the command or script entry point that skill documents, the parameters that change the result, and the output file.
5. Name the environment the skill declares (`phi:python@1`, `phi:r@1`, or the skill-local environment). Do not add a new `environment.yml` unless `create-skills` is in use and the checklist there is met.
6. Run only the steps the user asked to run. Leave the rest as a plan. After a run, fill in the outputs that were actually written.

## Rules

- One plan file, checked in only when the user wants it in the project. Do not scatter parameters across chat and scripts.
- Pin versions the skill already pins. Do not add a floating `pip install` when the skill says the library is in `phi:python@1`.
- Record random seeds, genome builds, and column names that the result depends on.
- If a required input is missing, stop and name the file. Do not synthesize a stand-in dataset and treat it as the user's data.

## Report

Lead with the question and the skill that will run it. Then the input list, the command sequence, the environment, and the expected outputs. Mark steps that were run and steps that are still a plan.
