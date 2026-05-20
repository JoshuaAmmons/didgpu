# Windows GPU build — status & remaining work (#95)

## TL;DR

The Windows CUDA toolchain is **fully assembled with no admin rights**,
and **everything compiles** — all 7 `.cu` kernels (nvcc + MSVC) and all
the R-facing C++ (MinGW g++, including `cuda_runtime.h`). The single
remaining blocker is the **MinGW ↔ MSVC link barrier**, whose designed
solution is the two-DLL split. **The admin password is no longer
needed** — a system-wide CUDA install would not change the link error.

## What works now (no admin)

- **`tools/assemble-cuda-local.ps1`** downloads the CUDA 12.6 redist
  component ZIPs (nvcc, cudart, cccl, cuBLAS, cuSOLVER, cuRAND) and
  merges them into `C:\Users\ammonsj\cuda_local` — a complete toolkit
  in a user-writable folder. No installer, no UAC. This fully replaces
  the broken/partial Program Files CUDA 13.2 install.
- **`tools/build-didgpu-win.ps1`** points the build at `cuda_local`
  (`CUDA_PATH`) and runs `R CMD INSTALL` with Rtools44 + VS Build Tools.
- `src/Makevars.win` fixes:
  - pattern rule `%.o: %.cu` (was explicit per-file rules that silently
    dropped `cuda_bootstrap.cu`);
  - `-allow-unsupported-compiler` (CUDA 12.6 vs the newer MSVC 14.44).
- Result: every kernel compiles for sm_75/80/86/89/90; `didgpu_init.cpp`
  compiles against the CUDA headers under MinGW g++.

## The remaining barrier: MinGW can't link MSVC objects

The single-DLL approach bundles the nvcc/MSVC-compiled `.cu` objects
into the MinGW-linked `didgpu.dll`. The link fails because the MSVC
objects reference symbols MinGW's `ld` cannot supply:

1. `__security_cookie` / `__security_check_cookie` / `__GSHandlerCheck`
   — MSVC stack-protector. **Fixed** via no-op stubs in
   `src/win_msvc_stubs.cpp`.
2. MSVC C++ math overloads (`?sqrt@@YAMM@Z` = `float sqrt(float)`,
   `?sin@@YAMM@Z`, `fabsf`, `frexpf`, ...) — dozens of mangled symbols
   from MSVC's CRT. **Not stubbable** in practice.
3. Cross-references to our own `didgpu_cuda_*` functions don't resolve
   cleanly under the mixed-toolchain link.

These are fundamental: MinGW and MSVC have different C++ runtimes/ABIs.
No amount of stubbing makes a MinGW DLL safely host MSVC objects.

## The fix: two-DLL split (already designed)

`inst/include/didgpu_cuda_api.h` defines a pure-C ABI for exactly this:

- **`didgpu_cuda.dll`** — built entirely by `nvcc --shared` (which uses
  MSVC's `link.exe`, so MSVC resolves all its own math/runtime symbols).
  Exports the `extern "C"` `didgpu_cuda_*` functions (via a `.def` file
  or `__declspec(dllexport)`).
- **`didgpu.dll`** — the R-facing package DLL, built by MinGW as today,
  but linking against a MinGW import lib for `didgpu_cuda.dll`
  (`dlltool -d didgpu_cuda.def -l libdidgpu_cuda.dll.a`) and calling
  only the C-ABI functions. No MSVC objects in this DLL → no ABI clash.
- **Bundling**: `didgpu_cuda.dll` + the CUDA runtime DLLs (cudart64,
  cublas64, cusolver64, curand64) ship in the package's `libs/x64/` so
  they load alongside `didgpu.dll`.

### Concrete Makevars.win changes for the two-DLL build
1. Add a rule: `didgpu_cuda.dll: $(CUDA_OBJECTS)` →
   `nvcc --shared -o didgpu_cuda.dll $^ -L"$(CUDA_HOME)/lib/x64"
   -lcudart -lcublas -lcusolver -lcurand -Xlinker /DEF:didgpu_cuda.def`.
2. Generate the MinGW import lib via `dlltool` from `didgpu_cuda.def`.
3. Remove `$(CUDA_OBJECTS)` from `OBJECTS`; set
   `PKG_LIBS = -L. -ldidgpu_cuda`.
4. `$(SHLIB): didgpu_cuda.dll libdidgpu_cuda.dll.a`.
5. An install step (or `Makevars` `all:` extension) to copy
   `didgpu_cuda.dll` + the four CUDA runtime DLLs into the install
   `libs/x64/`.
6. Author `src/didgpu_cuda.def` listing the exported `didgpu_cuda_*`
   symbols (saxpy, run_one_event_time, fect_fe, cs_inner_or,
   cs_inner_logit, fect_svd_truncated, fect_svd_softthreshold,
   cluster_bootstrap, multiplier_bootstrap, testmechs_bootstrap).

This is the bulk of remaining Phase 5 (#95) work and overlaps Phase 6
(#96, bundling for distribution). It is well-defined but multi-step,
with several link/runtime-path failure points to work through
iteratively.

## Note for the curious

`src/win_msvc_stubs.cpp` is guarded `#if _WIN32 && HAS_CUDA`. In the
two-DLL design `didgpu.dll` carries no MSVC objects, so the stub TU is
empty and harmless there; it can be removed once two-DLL lands.
