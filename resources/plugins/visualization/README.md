# Scientific Plotting Plugin

This built-in plugin provides:

- `agents/Visualization.md`: the dedicated Visualization plotting agent.
- `skills/omics-visualization/`: the scientific plotting templates and scripts that draw user data with existing templates.
- `skills/create-template/`: turns a figure into a reusable template and registers it in the omics-visualization catalog. It does not draw user data with existing templates and has no scripts of its own.
- `resources/runtime/environments/phi-r/`: the built-in R environment `phi:r@1` (`phi-r`) shared by the agent and the omics-visualization skill, providing both `Rscript` and Python.

The plugin declares no private hosted environment; the plotting skills share `phi:r@1` with the other built-in skills.
`create-template` declares no environment. Validation commands such as `python3 scripts/validate_template_contracts.py` run from the installed `omics-visualization` skill root.

## Plot template gallery

`skills/omics-visualization/references/template_gallery.json` lists every plotting template in the catalog together with its `preview.png`. When the installed scientific plotting plugin is viewed, the plugin page reads this list and shows the template gallery by figure family. Previews are transparent PNGs, and the gallery uses a white background so that blank areas are not mistaken for black. These images are template examples, not plots of user data. The layout audit images under `scripts/layouts/` are not plotting templates and do not enter the gallery.

After adding or renaming a template, regenerate and validate from the omics-visualization skill root:

```bash
python3 scripts/build_template_gallery.py
python3 scripts/validate_template_contracts.py --json
```