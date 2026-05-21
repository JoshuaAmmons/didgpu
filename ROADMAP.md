# didgpu Roadmap

> Living document. Update as decisions land. To resume after a context break,
> read top-to-bottom — the **Current state** and **Resume here** sections
> describe exactly what's done and what's next.
>
> **Companion:** [`ESTIMATES.md`](ESTIMATES.md) — per-phase time estimates
> (O / M / P / PERT-expected, in focused work hours) and calendar
> translation table.

## Goal

Ship a GPU-accelerated R package where a Windows + RStudio + NVIDIA-GPU user
gets a one-line install with zero toolchain setup on their end:

```r
install.packages("didgpu", repos = "https://jdammons.r-universe.dev")
library(didgpu)
fit <- didgpu_cs(panel, "Y", "unit", "period", "D", backend = "cuda")
```

End user needs only:

- a current NVIDIA driver (they have it because they have a CUDA GPU)
- R + RStudio

End user does NOT need: CUDA Toolkit, MSVC Build Tools, Rtools, WSL, admin
rights, environment-variable surgery.

## Architecture

### Two-DLL design on Windows

The package installs two shared libraries side-by-side in `libs/x64/`:

```
didgpu/
├── DESCRIPTION
├── R/                       # R source (data.table-based fallback + GPU dispatch)
├── inst/
│   └── include/
│       └── didgpu_cuda_api.h    # C-ABI contract (the boundary)
├── libs/x64/
│   ├── didgpu.dll               # Built by Rtools g++ (MinGW). R-loadable.
│   ├── didgpu_cuda.dll          # Built by nvcc + MSVC link.exe. Owns CUDA.
│   ├── cudart64_*.dll           # Bundled CUDA runtime DLLs (redist)
│   ├── cublas64_*.dll
│   ├── cusolver64_*.dll
│   └── curand64_*.dll
```

The two DLLs communicate via the **pure-C ABI** defined in
`inst/include/didgpu_cuda_api.h`. Strict rules:

1. `extern "C"` only — no C++ name mangling across the boundary.
2. POD types only (no `std::*`, no Rcpp types).
3. No C++ exceptions cross the boundary; all functions return `int` error
   codes; use `didgpu_cuda_last_error()` for human messages.
4. **All pointers are HOST pointers.** `didgpu_cuda.dll` handles `cudaMalloc`,
   H2D, D2H, `cudaFree` internally. This is what insulates the MinGW-built
   `didgpu.dll` from CUDA runtime symbols (`__security_cookie`, `_Init_thread_*`,
   etc.) that broke the single-DLL Windows build.
5. Matrices are row-major; sizes are `int`; RNG seeds are `uint64_t`.
6. ABI versioning via `DIDGPU_CUDA_ABI_VERSION` (bump on incompatible changes).

### Single-SO design on Linux/macOS

Linux/Mac use only g++ (no MSVC), so no toolchain bridge is needed. Everything
compiles into one `libdidgpu.so`. The C-API in `didgpu_cuda_api.h` is still the
*logical* boundary between layers — same header, same rules — but
implementation is one shared object.

## End-user UX matrix

| Platform | GPU? | UX |
|---|---|---|
| Windows + RStudio | yes | `install.packages(...)` — uses Windows binary with bundled CUDA runtime. |
| Windows + RStudio | no  | Same install, falls back to R/data.table backend silently. |
| Linux + RStudio | yes | Source install via r-universe; requires CUDA Toolkit + driver. |
| Linux + RStudio | no  | Source install; R/data.table backend only. |
| macOS | n/a | R/data.table backend only (NVIDIA dropped Mac). |

## Phase plan

Numbered. Pre-Windows-binary phases (0a, 0b, 1, 2) can all run on WSL2 Linux —
fastest path to working CUDA code. Windows-binary work (3, 4) comes after.

### Phase 0a — WSL2 Ubuntu dev env (in progress)
- [x] Verify GPU passthrough: `wsl.exe -d Ubuntu -- nvidia-smi` works,
      RTX 4000 Ada visible.
- [x] Write `tools/setup-wsl-env.sh` (one-shot installer).
- [ ] Run `tools/setup-wsl-env.sh` — installs gcc, R, CUDA Toolkit 12.6,
      r-base-dev + headers, Rcpp/RcppEigen/data.table/testthat/broom/quadprog.
      ~15-25 min, ~2-3 GB download.

### Phase 0b — Linux build with CUDA
- [x] `src/Makevars` (Linux) — rewritten with proper CUDA detection,
      `.DEFAULT_GOAL := all`, `-lcurand`, sm_90.
- [ ] Run `R CMD INSTALL .` in WSL2. Verify all 6 `.cu` files compile and
      `libdidgpu.so` links.
- [ ] Smoke test: `didgpu_run_saxpy()` and `didgpu_cuda_did()` produce
      results matching the R backend within tolerance.

### Phase 1 — Wire existing CUDA scaffolds through R (3 orphan scaffolds)
Tasks #79-81. Each scaffold has a host launcher in `src/cuda_*.cu` already.
- [ ] `cuda_cs_inner` — Rcpp wrapper + R dispatch in `R/cs.R`.
- [ ] `cuda_fect_svd` (truncated + softthreshold) — wire to `R/fect_ife.R`
      and `R/fect_mc.R`.
- [ ] `cuda_testmechs_bootstrap` — wire to `R/testmechs_bootstrap.R`.
- [ ] Parity tests: `backend = "cuda"` == `backend = "r"` within numerical
      tolerance for each.

### Phase 2 — Fill kernel gaps (the biggest user-facing wins)
Tasks #82-87. New kernels that don't exist yet.
- [ ] **Cluster bootstrap** kernel (`didgpu_cuda_cluster_bootstrap`).
      Used by both didgpu and didgpu_cs. Biggest single performance win
      because B = 1000+ replicates × n_clusters is embarrassingly parallel.
- [ ] **Multiplier (wild) bootstrap** kernel
      (`didgpu_cuda_multiplier_bootstrap`). Even faster than cluster bs.
- [ ] **CS IPW + DR augmentation** (per-(g,t) logistic regression for
      propensity, weighted outcome regression for augmentation).
- [ ] **CS aggregation + influence accumulation** — element-wise GPU
      reductions; keeps the full CS pipeline resident on GPU.
- [ ] **fect_ife fused alternation** — FE step + truncated SVD step in one
      host loop, data stays on GPU across iterations.
- [ ] **fect_mc CV inner loop** — K-fold parallel on GPU.

### Phase 3 — Batched LOO + benchmarks
Tasks #88-89.
- [ ] `didgpu_cuda_loo_batched` — fits K leave-one-cohort-out replicates in
      one GPU launch (just adds a leave-out axis to the existing per-(g,t)
      batched problem).
- [ ] Benchmark suite (CUDA vs R, all 5 families) on RTX 4000 Ada.

### Phase 4 — Verification + docs
Tasks #90-91.
- [ ] Full equivalence test grid: backend="cuda" == backend="r" within
      tolerance, across parameter grids for every estimator.
- [ ] README + NEWS updates with GPU coverage matrix + benchmark table.

### Phase 5 — Windows two-DLL build
Task #95. Tackled AFTER kernels are correct on Linux.
- [ ] `src/Makevars.win` — refactor to drive the two-DLL build.
- [ ] `tools/build-win-binary.ps1` — orchestrates:
    1. nvcc compiles `.cu` files into `didgpu_cuda.dll` via MSVC link.exe.
    2. `dlltool` creates `libdidgpu_cuda.dll.a` (MinGW import library).
    3. Rtools g++ builds `didgpu.dll`, links against the import library.
    4. Bundle CUDA runtime DLLs (cudart, cublas, cusolver, curand) from
       NVIDIA's redist archives into `inst/libs/x64/`.
    5. `R CMD INSTALL --build` produces `didgpu_<ver>.zip`.
- [ ] **Refactor `didgpu_init.cpp`** to use the host-pointer C-API
      (`didgpu_cuda_api.h`) instead of calling cudaMalloc/cudaMemcpy
      directly. This is what makes the MinGW side CUDA-free. The CUDA
      calls move into a new `src/cuda_host_wrappers.cu` file that
      compiles only into `didgpu_cuda.dll`.

### Phase 6 — Distribution
Task #96.
- [ ] Push to GitHub.
- [ ] r-universe entry — auto-builds Linux/Mac source binaries (no CUDA).
- [ ] Custom Windows binary URL — points at GitHub Release asset built
      by Phase 5.
- [ ] r-universe routing → end user `install.packages(...)` just works.

### Phase 7 — Fresh-Windows-VM UX validation
Task #97.
- [ ] Fresh Windows VM (or fresh R library path), only the NVIDIA driver
      installed. `install.packages("didgpu", repos="https://...")`. Verify
      `didgpu_cuda_available()` returns TRUE and a GPU pipeline runs
      end-to-end. This is the actual UX contract test.

## Distribution pipeline (Phase 6 detail)

```
GitHub (source repo)
   │
   ├── r-universe.dev/jdammons/builds source binaries automatically
   │     ├── Linux: source install, CUDA optional at user side
   │     └── macOS: source install, CPU only
   │
   └── tools/build-win-binary.ps1 (run on dev box with NVIDIA GPU + MSVC + Rtools)
         ├── nvcc → didgpu_cuda.dll
         ├── dlltool → libdidgpu_cuda.dll.a
         ├── Rtools g++ → didgpu.dll
         ├── bundle CUDA runtime DLLs into inst/libs/x64/
         └── publish didgpu_<ver>.zip to GitHub Release
              ↑
              r-universe pulls this URL for the Windows binary slot
```

End user types `install.packages("didgpu", repos="https://jdammons.r-universe.dev")`,
r-universe serves the right binary, RStudio installs it, done.

## Reality check: verified vs scaffolded

**Important context discovered 2026-05-19:** Several earlier tasks marked
"completed" tracked code-written, not code-verified-on-GPU. The task list
made it look like CUDA was working when in fact no CUDA path has ever been
end-to-end verified to produce correct output. Be careful with this
distinction going forward.

### Verified and working today (the real package)

These are tested, cross-validated, and shippable:

- **R/data.table backend** across all 5 estimator families:
  - `didgpu()` — DIDmultiplegtDYN-equivalent
  - `didgpu_cs()` — Callaway-Sant'Anna (OR / IPW / DR)
  - `didgpu_fect()` — fe / ife / mc
  - `didgpu_test_sharp_null()` — TestMechs (CS / ARP / FSST)
  - `didgpu_honest_did()` — HonestDiD
- **60× faster than `DIDmultiplegtDYN` reference** at 200K rows
  — *R-vs-R*, not GPU.
- **Cross-validation** against `did`/`DRDID` reference packages (max abs
  diff < 0.25 on event-study estimates).
- **300+ tests passing** including 21-scenario fuzz harness, parallel
  bootstrap, checkpoint round-trips, multivalued treatment.
- **Full S3 method coverage** — print, summary, coef, confint, vcov, plot,
  plus broom tidy/glance/augment.
- **Checkpoint/resume infrastructure** — `didgpu_resume()`,
  `didgpu_bootstrap_more()`, atomic writes via saveRDS tmp+rename.
- **Leave-one-out (`didgpu_loo`)** — all families, tornado plot, 8 tests
  passing.
- **R CMD check passes** with only the declared CUDA `.cu` WARN and
  `SystemRequirements: GNU make` NOTE.

### Scaffolded but not verified

These have `.cu` files and (sometimes) R-side dispatch branches, but the
CUDA path has **never been compiled, run, or output-verified** on this box
or any other:

- `src/cuda_didkernel.cu` — the DiD U-statistic kernel
- `src/cuda_fect_fe.cu` — fect_fe iterative demeaning
- `src/cuda_fect_svd.cu` — truncated + soft-threshold SVD
- `src/cuda_cs_inner.cu` — Callaway-Sant'Anna per-(g,t) inner regression
- `src/cuda_testmechs_bootstrap.cu` — partial-density bootstrap
- `src/cuda_saxpy.cu` — hello-world kernel

The `R/backend.R` dispatch has `if (backend == "cuda")` branches in
several places, but these branches **silently fall back to R** when the
CUDA build is broken (which it is on Windows right now). The "Overnight
100x stress tests" (task #75) ran entirely on the R backend; the GPU
paths were never exercised.

### Implication for shipping

If we wanted to publish *today*, we could ship the R-only package as
something like `didfast` or rename to be honest about scope. The R
backend alone is a meaningful contribution — a fast, cross-validated R
implementation of 5 estimator families with checkpoint/resume,
parallel bootstrap, and broom integration is genuinely useful even
without GPU.

The "GPU" branding is currently aspirational, not factual. Phases 0a-7
in this document are what makes it factual.

## Current state

**Today: 2026-05-19.** Active branch: `master` (no feature branch yet).

### What's done

| File | What's there |
|---|---|
| `R/loo.R` | LOO across all families. 8 tests passing. |
| `tests/testthat/test-loo.R` | LOO test suite. |
| `inst/include/didgpu_cuda_api.h` | C-API contract pinned at ABI v1. |
| `src/Makevars` | Linux build, ready for WSL2. |
| `tools/setup-wsl-env.sh` | One-shot WSL2 setup script. |
| `src/Makevars.win` | Partial Windows build — compiles `.cu` files but
                       BLOCKED on MinGW↔MSVC link. Will be refactored in Phase 5. |
| `src/nvcc_wrapper.bat` | cmd-side helper for nvcc + vcvars. Reusable in Phase 5. |
| `src/cuda_shim/crt/` | Real NVIDIA crt headers from cuda_crt redist
                         (patches around broken CUDA 13.2.1 install). Phase 5. |
| `~/cuda_local/` (user-local) | Shadow CUDA tree for the broken install.
                                  Junctions to Program Files + redist additions. |
| `~/cuda_crt_redist/` | NVIDIA redist downloads (cache; can be deleted). |

### What's not done / blocked

- Phase 0a — waiting for user to run `tools/setup-wsl-env.sh`.
- Everything after that is sequential on Phase 0a.

### Decisions log (so we don't relitigate)

- **Decision (2026-05-19):** Dev on WSL2 Linux first; ship Windows binary
  via a separate Phase 5 two-DLL build. Reason: MinGW (Rtools) cannot link
  MSVC-compiled CUDA objects (`__security_cookie`, `_Init_thread_*`,
  `__GSHandlerCheck` undefined). Documented in
  [this conversation's context]. The two-DLL split is the industry-standard
  workaround.
- **Decision (2026-05-19):** Small-K QP solvers (Cox-Shi, ARP, FSST,
  HonestDiD; K ≤ 10) stay on CPU. GPU kernel-launch latency exceeds the
  CPU runtime for these. User-confirmed.
- **Decision (2026-05-19):** Keep `backend = "r"` selectable even on
  CUDA-equipped machines. Useful for debugging; matches CRAN requirement
  for non-CUDA fallback.
- **Decision (2026-05-19):** No deadline pressure; do this right. The
  one-line `install.packages(...)` UX is the bar.
- **Decision (2026-05-19):** Target user is Windows + RStudio + GPU. Linux
  is dev-only.

## Resume here

After a context break, the natural next step is:

1. **Check whether `tools/setup-wsl-env.sh` has been run.** Either ask the
   user or test:
   ```
   wsl.exe -d Ubuntu -- bash -c "command -v nvcc && command -v R"
   ```
   If both print paths, Phase 0a is done.

2. **If Phase 0a done, run Phase 0b:**
   ```
   wsl.exe -d Ubuntu -- bash -c "cd '/mnt/c/Users/ammonsj/DID GPU/didgpu' && R CMD INSTALL ."
   ```

3. **If install succeeds, smoke test:** run
   `tests/testthat/test-cuda-saxpy.R` (or equivalent) under WSL.

4. **Then Phase 1:** wire the 3 orphan scaffolds. Start with
   `cuda_cs_inner` (most user-facing). Tasks #79, #80, #81 in TaskList.

5. **Don't touch Phase 5 (Windows two-DLL) until Phases 1-4 are real on
   Linux.** It's tempting to context-switch when a kernel works — resist.
   Get all the kernels right on one platform first.

## Reference: file locations

- This roadmap: `ROADMAP.md` (repo root)
- C-API contract: `inst/include/didgpu_cuda_api.h`
- Linux build: `src/Makevars`
- Windows build (partial): `src/Makevars.win`
- WSL setup: `tools/setup-wsl-env.sh`
- Task list: in TaskCreate/TaskList (#79-#97 are current; #77-78 deleted)
- User-local CUDA shadow tree: `C:\Users\ammonsj\cuda_local\`
- NVIDIA redist downloads: `C:\Users\ammonsj\cuda_crt_redist\`
