# Building didgpu with CUDA support

didgpu is designed to compile with or without the NVIDIA CUDA Toolkit.
Without it, the package builds CPU-only and the CUDA backend reports
itself unavailable (`didgpu_backend_info()` shows `cuda = FALSE`,
`didgpu_has_cuda_support()` returns `FALSE`). With it, src/cuda_*.cu
files are compiled by `nvcc` and linked into the package DLL.

This document covers Windows; macOS/Linux is similar but uses the
Unix paths in `src/Makevars`.

---

## Requirements (Windows)

1. **NVIDIA GPU** with a recent driver. Confirm via `nvidia-smi`.
   This dev box runs an RTX 4000 Ada (compute capability 8.9, 12 GB).
2. **CUDA Toolkit >= 12.0** providing `nvcc.exe`, cuBLAS, and cuSOLVER.
   - Driver-only is NOT enough; the toolkit is a separate install.
3. **Rtools 4.4** (already required by any source-built R package on
   Windows).
4. **Visual Studio Build Tools** with the MSVC compiler. nvcc on
   Windows uses MSVC `cl.exe` to compile host-side code; even though
   linking the package DLL goes through Rtools g++, the `.cu` -> `.obj`
   path requires `cl.exe`. Install via the Visual Studio installer
   ("Desktop development with C++" workload).

---

## Install steps

### 1. Install the CUDA Toolkit

Download from <https://developer.nvidia.com/cuda-downloads> or via
winget:

```powershell
winget install --id Nvidia.CUDA --exact --accept-source-agreements
```

The installer typically lays the toolkit down at
`C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.x\`. It also
sets `CUDA_PATH` and `CUDA_PATH_V12_X` environment variables and
appends `%CUDA_PATH%\bin` to `PATH`.

Verify:

```powershell
nvcc --version
where nvcc
```

If `nvcc` is not on PATH, set `CUDA_HOME` to the toolkit root before
building didgpu. The Makevars.win falls back to `$CUDA_PATH` if
`CUDA_HOME` is unset.

### 2. Verify MSVC is available

```powershell
where cl
```

If not found, install via Visual Studio Build Tools and re-open the
shell so the environment variables apply. nvcc will pick up `cl.exe`
automatically if it's on PATH.

### 3. Rebuild didgpu

```r
remove.packages("didgpu")
install.packages(
  "C:/Users/ammonsj/DID GPU/didgpu",
  repos = NULL, type = "source",
  INSTALL_opts = "--no-multiarch"
)
```

Watch the build log for a line like:

```
  [nvcc] cuda_saxpy.cu
```

That's the nvcc invocation. If it does not appear, `nvcc` was not
detected (`HAVE_NVCC` was empty in Makevars.win) and the package
silently built CPU-only.

### 4. Confirm CUDA is wired in

```r
library(didgpu)
didgpu_has_cuda_support()
#> [1] TRUE

# SAXPY round-trip: y <- 2 * x + y on GPU
didgpu_run_saxpy(2.0, c(1, 2, 3), c(10, 20, 30))
#> [1] 12 24 36
```

If `didgpu_has_cuda_support()` returns FALSE after a fresh install
with nvcc on PATH, check:

- `cuda_saxpy.o` appears under
  `C:/Users/ammonsj/AppData/Local/R/win-library/4.4/didgpu/libs/`
  (it should, if the linker pulled it in).
- `R CMD INSTALL` was run from a shell where `nvcc` resolves.

---

## What's compiled vs. what's stub

Build artefact | Source | Built when
---|---|---
`didgpu_init.o` | `src/didgpu_init.cpp` | Always
`cpu_hello.o`   | `src/cpu_hello.cpp`   | Always (smoke target)
`RcppExports.o` | `src/RcppExports.cpp` | Always (Rcpp glue)
`cuda_saxpy.o`  | `src/cuda_saxpy.cu`   | Only when nvcc detected

The CUDA path is intentionally minimal at present: a SAXPY kernel that
proves nvcc -> g++ link works end-to-end. The real didgpu CUDA
implementation is Phase 3 of the project plan (see
`NOTES_did_gpu_checkpointed.Rmd`); when that lands it will live in
new src/cuda_*.cu files compiled by the same Makevars rule.

---

## Known issues on Windows

- **The winget CUDA package may be incomplete.** Installing via
  `winget install --id Nvidia.CUDA` on this dev box gave us a CUDA
  13.2 install whose `include/` directory is missing the `crt/`
  subdirectory. `cuda_runtime.h` `#include`s `"crt/host_config.h"`,
  which then fails. The "real" `host_config.h` at the include root is
  a deprecated forwarder that itself includes `crt/host_config.h`.
  Workaround: run the full CUDA installer from
  <https://developer.nvidia.com/cuda-downloads> selecting the
  **Toolkit** component (not just nvcc), which puts the crt/ headers
  in place. winget's minimal selection apparently omits them.
- **Rtools g++ vs. MSVC name mangling.** Any host function in a `.cu`
  file that you intend to call from `.cpp` must be `extern "C"`. The
  reference pattern is in `src/cuda_saxpy.cu`.
- **Compute capabilities baked into Makevars.** The default
  `NVCC_FLAGS` target sm_75, sm_80, sm_86, sm_89. If your GPU is older
  (sm_70 or below), edit `src/Makevars.win` before installing, or
  newer GPUs may take a JIT-compile hit at first kernel launch.
- **`CUDA_PATH_V12_X` is set, `CUDA_PATH` is not.** Some installers
  leave only the versioned variable. Add an explicit `CUDA_PATH` to
  your environment or set `CUDA_HOME` directly:

  ```powershell
  setx CUDA_HOME "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6"
  ```
