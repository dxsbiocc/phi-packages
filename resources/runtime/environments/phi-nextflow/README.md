# phi-nextflow

Managed Nextflow for the wrapper executor (`phi:nextflow@1`). Local wrapper runs spawn `nextflow` from this prefix by default; a host nextflow is used only when the user sets one explicitly in Settings → Environment, and only if it meets the wrappers' minimum version (it is then labelled "host (unmanaged)"). Remote runs keep the remote host's nextflow. Docker and Singularity / Apptainer are host dependencies, declared under `host:` and not installed here; the executor puts the one the run's profile needs on PATH.

Only `nextflow>=25.04.3` is constrained in `environment.yml`: that is the highest `manifest.nextflowVersion` among the bundled wrappers (`workflows/rna-seq`). The lock files pin every package.

Channels are `conda-forge`, then `bioconda`. `nextflow`, `nf-core`, and `nf-test` come from bioconda. Every other package comes from conda-forge.

## Contents

| Package     | Why                                                                                        |
| ----------- | ------------------------------------------------------------------------------------------ |
| `nextflow`  | the executor's launcher (26.04.6 in all three locks)                                       |
| `openjdk`   | Nextflow is a JVM program; nextflow 26.04.x requires `openjdk >=17,<26` (25.0.x is locked) |
| `coreutils` | declared by the nextflow package                                                           |
| `curl`      | declared by the nextflow package                                                           |
| `nf-core`   | authoring wrapper modules (nextflow skill); 4.1.0                                          |
| `nf-test`   | wrapper module tests (`tests/main.nf.test`); 0.9.5                                         |

`-profile conda` does not use this prefix for task software. Nextflow solves each process's `environment.yml` at run time with the bundled micromamba into `<runtime root>/nextflow-conda` (see `src/main/agent/wrappers/composition/conda-profile.ts`).

## Locks

Solved on 2026-10-01 against the lock baselines: macOS arm64 11.0, macOS x64 10.15, and linux glibc 2.17 / linux 4.18.

| Platform       | Packages | Download |
| -------------- | -------- | -------- |
| `darwin-arm64` | 166      | 361 MiB  |
| `darwin-x64`   | 160      | 348 MiB  |
| `linux-x64`    | 194      | 347 MiB  |

Download sizes are the `# download-bytes:` header of each lock.

Regenerate all three platform locks:

```bash
npm run runtime:lock -- --spec resources/runtime/environments/phi-nextflow/environment.yml
```

That writes `locks/darwin-arm64.txt`, `locks/darwin-x64.txt`, and `locks/linux-x64.txt` next to `environment.yml`.

## Known platform differences

With the macOS 10.15 baseline, `darwin-x64` resolves `openjdk` 25.0.1, `coreutils` 9.5, and `curl` 8.18.0; `darwin-arm64` and `linux-x64` get `openjdk` 25.0.2, `coreutils` 9.12, and `curl` 8.22.0. Nextflow, nf-core, and nf-test are the same everywhere.
