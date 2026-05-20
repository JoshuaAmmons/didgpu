# Overnight session status — 2026-05-20

Autonomous session run while Joshua was away. **The headline goal is
done: the Windows + NVIDIA GPU build works, with no admin rights.** The
admin password you were going to fetch is **not needed** — the
user-local CUDA toolkit + the two-DLL split do it entirely in user space.

## TL;DR

- **#95 Windows two-DLL GPU build — DONE & verified.** Full multi-arch
  (Turing→Hopper) builds, links, bundles, loads, and computes
  **bit-exact vs R** on the RTX 4000 Ada. The prior session left this as
  "needs a Windows session, can't do from WSL"; it's now resolved.
- **#96 distribution — prepared & a binary built.** A 2.92 MB GPU-lite
  `didgpu_0.1.0.zip` is built + validated, ready to upload; `RELEASE_PLAN.md`
  has the exact publish commands. Only the publish step (your GitHub auth)
  remains.
- **#99 + #100 — DONE.** Two new pre-trends diagnostics, pure-R, verified.
- **Linux suite stays green: 976 pass / 0 fail** (50 files; 82 skips are
  reference-package-gated). Windows GPU verification all bit-exact.

## What landed this session (all pushed to `main`)

| Commit | What |
|--------|------|
| `9ccfc4b` | **#95 host-pointer C-ABI refactor** — moved all CUDA device-memory mgmt off the MinGW side into the `.cu` (nvcc/MSVC) side, so `didgpu.dll` has zero CUDA-runtime calls. Verified bit-identical on Linux first. |
| `1cbb057` | **#95 two-DLL build works end to end** — root-caused the runtime load failure to missing transitive deps (`cusolver→cusparse→nvJitLink`); added them to the toolkit + bundle. Loads + computes on Windows. |
| `5c4766d` | **#95 docs** — README/NEWS flipped to "CUDA live on Windows + Linux, no admin." |
| `d07c4bf` | **#96 `DIDGPU_LITE` build** — cudart-only, **2.32 MB vs 1.1 GB** (drops only the size-gated fect SVD → transparent R fallback). |
| `84b4ccd` | **#96 RELEASE_PLAN.md** — channel analysis (r-universe = CPU baseline; GitHub Releases = GPU binary). |
| `405fc73` | **#96 win.binary built + validated** — `didgpu_0.1.0.zip` 2.92 MB, confirmed it packages the CUDA DLLs. |
| `4cd5941` | **#99 `didgpu_equivalence()`** — pre-trends TOST equivalence test (per-horizon + joint + `breakdown_delta`). |
| `e420c2f` | **#100 `didgpu_joint_placebo()`** — windowed/subset joint placebo test (full-window == `p_jointplacebo` exactly). |
| `421e26f` | **#96 README install + release configs** — three install paths; copy-paste r-universe `packages.json` + `gh release` commands. |

## How it works (the two-DLL split)

MinGW (Rtools) can't link nvcc/MSVC objects (different C++ runtimes/ABIs).
So didgpu ships as two DLLs across a pure-C ABI:
`didgpu_cuda.dll` (built by `nvcc --shared`, MSVC resolves its own
symbols, owns all CUDA state) + `didgpu.dll` (MinGW, R-facing, zero CUDA
calls, links a `dlltool` import lib). The full dependency closure
(cudart + — in the full build — cublas/cublasLt/cusolver/cusparse/nvJitLink)
is bundled into `libs/x64`. Details in `WINDOWS_BUILD_STATUS.md`.

## Verification (RTX 4000 Ada, sm_89 and full multi-arch)

CUDA == R, identical on Windows native and Linux:
de Chaisemartin DID effects **bit-exact**; saxpy **bit-exact**; truncated
SVD recon 6e-11; soft-threshold SVD 2e-8; fect_fe 1e-6. Lite build
verified to load with only cudart and fall back correctly for fect SVD.

## Your machine right now

`didgpu` is installed as the **full multi-arch GPU build** (everything,
incl. fect GPU SVD). The GPU-lite release artifact is at
`C:\Users\ammonsj\didgpu_release\didgpu_0.1.0.zip`.

## What needs YOU (can't be done unsupervised)

1. **Publish the GPU binary** — `gh release create v0.1.0 …` then the repo
   is one-line installable (commands in `RELEASE_PLAN.md` appendix).
2. **r-universe** — create `jdammons.r-universe.dev` with the `packages.json`
   from `RELEASE_PLAN.md` for the CPU baseline.
3. **#97 fresh-Windows-VM UX validation** — needs a clean VM.

## Deliberately not done

- **#101 `didgpu_bacon()`** (Goodman-Bacon decomposition) — on the roadmap
  but lower-priority and substantial (2×2 enumeration + weights, ideally
  cross-validated against `bacondecomp`). Left for a supervised session so
  the output format / validation target can be chosen with you.

## Note

The Price DAOs project was not touched. The recurring autonomous cron is
being stopped now that the prompt's goals are complete (per its own
"stop scheduling further work" instruction).
