# phi-r

Shared R 4.4 environment for R notebooks, scanpy's R interop, and the bundled
visualization plugin. Notebooks use IRkernel from this prefix as the `phi-r`
kernel; Phi writes its kernelspec with `<prefix>/bin/R` and this environment's
`environmentVariables`. Visualization script tools use both `python` and
`Rscript` from the same prefix.

Only `python=3.12` and `r-base=4.4` are pinned in `environment.yml`. Python is
present solely for the visualization skill's standard-library scripts. The lock
files pin every other package.

Channels are `conda-forge`, then `bioconda`. `bioconductor-*` packages come from bioconda. Every other package comes from conda-forge.

## Packages

| Use                | Packages                                                                                                                                                    |
| ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| script interpreter | `python=3.12` (standard library only)                                                                                                                       |
| notebooks          | `r-base=4.4`, `r-irkernel`                                                                                                                                  |
| general data work  | `r-tidyverse`                                                                                                                                               |
| scanpy R interop   | `bioconductor-singlecellexperiment`, `r-seurat` (with `r-seuratobject`)                                                                                     |
| visualization      | the explicit `r-*` and `bioconductor-complexheatmap` declarations in `environment.yml`                                                                      |
| source packages    | `GenomeInfoDbData`, then eight ordered visualization packages: `gground`, `ggideogram`, `ggcor`, `linkET`, `ggsankey`, `ggsvg`, `gridGeometry`, `ggmagnify` |

`r-ggplot2` is deliberately unpinned: `scripts/lib/common.R` patches
`ggideogram` 0.1.0 for ggplot2 4.x. `gridGeometry` immediately precedes
`ggmagnify` because the latter imports it. The other visualization source
packages do not depend on one another.

`r-waffle` is deliberately absent: its `r-rttf2pt1` dependency has no
`osx-arm64` build. The waffle template uses `geom_tile`. No host executables are
part of this environment.

### scanpy R interop without zellkonverter

`resources/skills/scanpy/references/r_interop.md` converts with `zellkonverter`. It is not installed: zellkonverter reads and writes `.h5ad` through basilisk, which creates its own Python environments at run time. An environment prefix is read-only after it is built, so that would fail (or write outside Phi's control). Use `SingleCellExperiment` and `Seurat` here, and do the `.h5ad` side in `phi-python`.

### GenomeInfoDbData

`SingleCellExperiment` loads `GenomeInfoDb`, which needs the data package `GenomeInfoDbData`. The bioconda package `bioconductor-genomeinfodbdata` contains no data: its post-link script downloads the tarball from Bioconductor at install time, which fails in a Phi build (the script cannot find `yq`, and it would fetch an unlocked file anyway). The package is therefore installed as a pinned `sourcePackages` entry from the Bioconductor GitHub mirror (version 1.2.15, commit `b5339e0`). The conda package stays in the lock as a dependency of `bioconductor-genomeinfodb`.

The builder defers that conda post-link failure only because the failed R package
has an exact, pinned `sourcePackages` replacement. Installing and verifying the
replacement still must succeed; unrelated link-script failures remain fatal.

## Locks

Solved on 2026-10-07 against the lock baselines: macOS arm64 11.0, macOS x64
10.15, and linux glibc 2.17 / linux 4.18.

| Platform       | Packages | Download bytes | lockSha256 (12) |
| -------------- | -------: | -------------: | --------------- |
| `darwin-arm64` |      481 |      666039033 | `20684ee37dc8`  |
| `darwin-x64`   |      482 |      601133421 | `ec53d72ec9df`  |
| `linux-x64`    |      478 |      731177933 | `ff7587a90f78`  |

Regenerate all three platform locks:

```bash
bun run runtime:lock -- --spec resources/runtime/environments/phi-r/environment.yml
```

That writes `locks/darwin-arm64.txt`, `locks/darwin-x64.txt`, and `locks/linux-x64.txt` next to `environment.yml`.

## Known platform differences

With the macOS 10.15 baseline, `darwin-x64` resolves `r-seurat` 5.4.0 and `r-seuratobject` 5.3.0; `darwin-arm64` and `linux-x64` get 5.5.1 and 5.4.0. `r-base` is 4.4.3 on all three.
