# create-template

Agent skill for **creating and registering** reusable plot templates in the
sibling [`omics-visualization`](../omics-visualization) skill (same
visualization plugin, `skills/omics-visualization/`).

Use when a figure image should become a catalog template, or when the user asks
to 做成模板 / 复现为绘图模板 / 保存模板 / 扩展模板 / 新增模板. For plotting
user data with an existing template, use `omics-visualization`
(`skill://omics-visualization`) instead.

Validators such as `python3 scripts/validate_template_contracts.py` and
`python3 scripts/qa_single_plot.py` run from the installed omics-visualization
skill root (`../omics-visualization/`). This skill does not ship those scripts.

See [SKILL.md](SKILL.md).
