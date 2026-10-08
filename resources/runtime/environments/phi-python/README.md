# phi-python

Shared Python 3.12 environment for the official Python skills. Notebooks use `ipykernel` from this prefix as the default kernel. Office skills run entirely inside this managed environment and do not require host-installed document-conversion tools.

Only `python=3.12` is pinned in `environment.yml`. The lock files pin every other package.

Channels are `conda-forge`, then `bioconda`. `harmonypy`, `bbknn`, and `pysam` come from bioconda. Every other package comes from conda-forge.

`coreutils` remains shared runtime tooling for compatibility. The Office
skills do not depend on it.

## Skills

| Skill           | Packages                                                                    |
| --------------- | --------------------------------------------------------------------------- |
| pptx            | `python-pptx`, `nodejs`, `pptxgenjs`, `defusedxml`, `lxml`, `pillow`, `git` |
| xlsx            | `openpyxl`, `python-calamine`, `xlsxwriter`, `defusedxml`, `lxml`, `git`    |
| docx            | `python-docx`, `lxml`                                                       |
| pdf             | `reportlab`, `pdfplumber`, `pypdf`, `poppler`                               |
| markitdown      | `markitdown`, `requests`, `openai`, `python-dotenv`, `tesseract`            |
| matplotlib      | `numpy`, `scipy`, `pandas`, `matplotlib`                                    |
| scikit-learn    | `scikit-learn`                                                              |
| scanpy          | `scanpy`, `leidenalg`, `python-igraph`, `harmonypy`, `bbknn`                |
| scvelo          | `scvelo`                                                                    |
| rdkit           | `rdkit`                                                                     |
| anndata         | `anndata`                                                                   |
| seaborn         | `seaborn`                                                                   |
| networkx        | `networkx`                                                                  |
| shap            | `shap`                                                                      |
| pysam           | `pysam`                                                                     |
| scikit-survival | `scikit-survival`                                                           |
| notebooks       | `python=3.12`, `ipykernel`                                                  |

## Locks

Solved on 2026-10-07 against the lock baselines: macOS arm64 11.0, macOS x64 10.15, and linux glibc 2.17 / linux 4.18.

Regenerate all three platform locks:

```bash
bun run runtime:lock -- --spec resources/runtime/environments/phi-python/environment.yml
```

That writes `locks/darwin-arm64.txt`, `locks/darwin-x64.txt`, and `locks/linux-x64.txt` next to `environment.yml`.

## Known platform differences

Locks are solved against fixed baselines (linux glibc 2.17, macOS arm64 11.0, macOS x64 10.15;
see `scripts/runtime/lock-env.ts`). With the macOS 10.15 baseline, `darwin-x64` resolves
`harmonypy` 0.2.0 while `darwin-arm64` and `linux-x64` get 2.0.2; the two versions differ in
API. Skills that call Harmony directly should check `harmonypy.__version__`, or the x64
baseline can be raised when Intel Macs on 10.15 no longer need support.
