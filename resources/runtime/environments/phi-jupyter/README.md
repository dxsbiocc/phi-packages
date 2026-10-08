# phi-jupyter

The Jupyter Server that hosts notebook sessions. Phi runs `<prefix>/bin/jupyter server` with this environment's `environmentVariables`; it never falls back to a host `jupyter`. No JupyterLab UI is installed.

Kernels do not live here. They run in the analysis environments: `ipykernel` in `phi-python` (the default kernel) and IRkernel in `phi-r`. Phi writes their kernelspecs into `<runtime root>/jupyter/kernels/` and starts the server with `--KernelSpecManager.kernel_dirs` set to that directory, and with `JUPYTER_PATH`, `JUPYTER_DATA_DIR`, `JUPYTER_CONFIG_DIR`, and `JUPYTER_RUNTIME_DIR` under `<runtime root>/jupyter/`, so `~/Library/Jupyter`, `~/.local/share/jupyter`, and `~/.jupyter` are not read.

Only `python=3.12` is pinned in `environment.yml`. The lock files pin every other package. The only channel is `conda-forge`.

## Packages

| Use         | Packages                                                                              |
| ----------- | ------------------------------------------------------------------------------------- |
| server      | `jupyter_server` (brings `jupyter_core`, the `jupyter` command, and `jupyter_client`) |
| interpreter | `python=3.12`                                                                         |

## Locks

Solved on 2026-10-01 against the lock baselines: macOS arm64 11.0, macOS x64 10.15, and linux glibc 2.17 / linux 4.18.

Regenerate all three platform locks:

```bash
npm run runtime:lock -- --spec resources/runtime/environments/phi-jupyter/environment.yml
```

That writes `locks/darwin-arm64.txt`, `locks/darwin-x64.txt`, and `locks/linux-x64.txt` next to `environment.yml`.

## Known platform differences

`darwin-x64` resolves `python` 3.12.12 under the macOS 10.15 baseline; `darwin-arm64` and `linux-x64` get 3.12.14. `jupyter_server` is 2.21.1 on all three.
