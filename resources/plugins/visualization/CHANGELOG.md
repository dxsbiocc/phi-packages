# Changelog

## 1.2.3 (2026-11-14)

- 模板画廊同时列出每个模板目录下的 `preview_<name>.png` CONFIG 变体图（沿用同一模板 id，标题带变体名），插件页现在能看到已有的 82 张变体；`scripts/layouts/` 的布局审核图仍排除在外，逃逸出模板目录的符号链接和 `preview_.png` 一类无名文件会被跳过。
- `create-template` 注册清单补充可选的变体图说明。

## 1.2.2 (2026-11-11)

- `tree-phylo-rings` 把扇形开口锁成精确 90°（原先 ggtree `open.angle` 会多出约半个 tip）。

## 1.2.1 (2026-11-11)

- 压缩全部模板预览图，插件下载体积更小。
- 新增原创内容的 MIT 许可，以及 loading.io 调色板的来源声明。
- 修复 `bar/upset` 和 `scatter/umap_atlas` 的分段标记，智能体现在能读到它们的配置和数据准备段。
- 修复 `scatter/variate_panel`、`radar/basic` 和 `line/ridge` 的示例说明，文件头里的命令可以照写直接运行。

## 1.2.0 (2026-11-06)

- 模板扩充到 211 个，智能体按目录精确选择模板。
- 提供 CONFIG 预览变体，并拆分容易混淆的图表类型。
