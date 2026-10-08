# 科研绘图插件

这个内置插件提供：

- `agents/Visualization.md`：专属的 Visualization 绘图智能体。
- `skills/omics-visualization/`：科研绘图模板与脚本。
- `resources/runtime/environments/phi-r/`：智能体与技能共享的内置 R 环境
  `phi:r@1`（`phi-r`），同时提供 `Rscript` 和 Python。

插件不声明私有托管环境；它与其他内置技能共享 `phi:r@1`。
