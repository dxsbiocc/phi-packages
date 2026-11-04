# 科研绘图插件

这个内置插件提供：

- `agents/Visualization.md`：专属的 Visualization 绘图智能体。
- `skills/omics-visualization/`：用已有模板绘制用户数据的科研绘图模板与脚本。
- `skills/create-template/`：把图做成可复用模板并登记进 omics-visualization 目录。不负责用已有模板画用户数据，本身没有脚本。
- `resources/runtime/environments/phi-r/`：智能体与 omics-visualization 技能共享的内置 R 环境
  `phi:r@1`（`phi-r`），同时提供 `Rscript` 和 Python。

插件不声明私有托管环境；绘图技能与其他内置技能共享 `phi:r@1`。
`create-template` 不声明环境。`python3 scripts/validate_template_contracts.py` 等校验命令在已安装的 `omics-visualization` 技能根目录运行。

## 绘图模板画廊

`skills/omics-visualization/references/template_gallery.json` 列出目录中的全部绘图模板及其 `preview.png`。查看已安装的科研绘图插件时，插件页读取这份清单，按图族展示模板画廊。预览图是透明 PNG，画廊用白底，避免把空白看成黑底。这些图是模板示例，不是用户数据的作图结果。`scripts/layouts/` 里的布局审核图不是绘图模板，不进入画廊。

新增或改名模板后，在 omics-visualization 技能根目录重新生成并校验：

```bash
python3 scripts/build_template_gallery.py
python3 scripts/validate_template_contracts.py --json
```
