# Distribution plan (#96)

## The goal

A Windows + RStudio + NVIDIA-GPU user runs **one command** and is
immediately working with GPU acceleration — no CUDA Toolkit install, no
admin rights, no compiler.

## What we can ship

The package builds in three meaningful flavours (all from the same source
tree; see `WINDOWS_BUILD_STATUS.md` for the two-DLL mechanics):

| Flavour | How | `libs/x64` payload | GPU coverage |
|---|---|---|---|
| **CPU-only** | build with no `nvcc` on PATH | ~1 MB (`didgpu.dll`) | none — every `backend="cuda"` call falls back to R |
| **GPU lite** | `DIDGPU_LITE=1` + nvcc | **2.32 MB** (`didgpu.dll` + `didgpu_cuda.dll` + `cudart64`) | everything **except** fect GPU SVD (which falls back to R, transparently) |
| **GPU full** | nvcc, default | **~1.1 GB** (+ cuBLAS 100 MB, cuBLASLt 507 MB, cuSOLVER 132 MB, cuSPARSE 274 MB, nvJitLink 37 MB, cuRAND 60 MB) | all paths incl. fect GPU SVD on very large panels |

All three are numerically identical where they overlap (verified
bit-exact for DID; ≤1e-6 for fect/SVD). The **GPU lite** flavour is the
one to distribute: it keeps every high-value GPU path (the CS cluster
bootstrap is the 179–228× headline; CS inner, DID, TestMechs, fect_fe all
run on the GPU) and drops only the size-gated fect SVD, which already
prefers CPU LAPACK on the small matrices typical of fect and only ever
engaged GPU on `n_units ≥ 2000` balanced panels.

Why not ship full? 1.1 GB is impractical: r-universe rejects it on size,
and a 1.1 GB download to save CPU time on a niche, size-gated path is a
bad trade. Power users who genuinely want fect GPU SVD on huge panels can
build full from source with the documented toolchain.

## Channel analysis

### r-universe (https://jdammons.r-universe.dev)
- Builds **from source on its own CI** (Linux/Windows/macOS runners).
- Those runners have **no NVIDIA GPU and no CUDA toolkit**, so the Windows
  build there resolves to **CPU-only** (Makevars.win sees no `nvcc`,
  `HAVE_CUDA_BUILD` stays empty, package builds clean without GPU).
- **Verdict:** r-universe gives users a correct, always-installable
  **CPU baseline** via the standard `install.packages(..., repos=)`. It
  **cannot** produce the GPU build. That's fine — it's the safety net.

### GitHub Releases (prebuilt GPU win.binary)
- We build the **GPU lite** win.binary locally (with the no-admin
  toolchain) and upload the `.zip` as a release asset.
- Users install the GPU build with one line (no toolchain needed):
  ```r
  install.packages(
    "https://github.com/JoshuaAmmons/didgpu/releases/download/v0.1.0/didgpu_0.1.0.zip",
    repos = NULL, type = "win.binary")
  ```
- At 2.32 MB of libs the asset is small and downloads instantly.
- **Verdict:** this is the GPU distribution channel.

## Recommended strategy

1. **r-universe** publishes the source package → CPU baseline, auto-built,
   always works. (Set up `jdammons.r-universe.dev` by adding the repo to a
   `universe` packages list — a one-time GitHub step the maintainer does.)
2. **GitHub Releases** hosts the prebuilt **GPU lite** win.binary per
   supported R minor (4.4 today; add 4.5 when the user upgrades). README
   gives the one-line installer above.
3. README documents both: "CPU from r-universe; GPU from the Releases
   one-liner; build full-fat GPU from source if you need fect SVD on huge
   panels."

## Concrete build commands

GPU **lite** win.binary (the release artifact), full multi-arch:
```powershell
$env:CUDA_PATH   = "C:\Users\<you>\cuda_local"   # from assemble-cuda-local.ps1
$env:DIDGPU_LITE = "1"                            # cudart-only
# (leave DIDGPU_CUDA_ARCH unset for the full Turing..Hopper gencode list)
$env:PATH = "C:\rtools44\usr\bin;C:\rtools44\x86_64-w64-mingw32.static.posix\bin;$env:PATH"
& "C:\Program Files\R\R-4.4.1\bin\R.exe" CMD INSTALL --build `
    --no-multiarch "C:\Users\<you>\DID GPU\didgpu"
# -> produces didgpu_0.1.0.zip in the CWD (this is the upload asset)
```

GPU **full** (for a power-user source build): same, omit `DIDGPU_LITE`.

CPU baseline: nothing special — r-universe does it, or any machine without
`nvcc`.

## R-version matrix

A win.binary is tied to the R minor (4.4 ABI). Build + release one `.zip`
per supported R minor. Today: R 4.4. The CUDA arch span (sm_75..sm_90,
Turing→Hopper) is already baked into `didgpu_cuda.dll`, so one binary
covers essentially all current NVIDIA GPUs regardless of R version.

## Open decisions for the maintainer (you)

- **Confirm GPU lite as the distributed flavour** (recommended). If you
  want fect GPU SVD in the shipped binary, we'd need to solve the 1.1 GB
  problem differently (e.g. delay-load cuSOLVER so the heavy DLLs are an
  optional side-download) — more work, niche payoff.
- **Account/automation:** publishing needs your GitHub auth (create the
  Release, upload the asset) and the r-universe setup. Not done here (no
  external account actions taken).
- **CI:** a GitHub Actions Windows runner can't build the GPU variant (no
  GPU/CUDA), so the GPU `.zip` is produced on your machine (or any machine
  with the no-admin toolchain) and uploaded manually / via `gh release`.

## Status

- [x] Source builds CPU / GPU-lite / GPU-full from one tree.
- [x] GPU-lite verified: 2.32 MB libs, loads, all non-SVD GPU paths
      bit-exact, fect SVD falls back correctly.
- [x] GPU-lite win.binary **built + validated**: full multi-arch
      `didgpu_0.1.0.zip` (**2.92 MB**) produced via `R CMD INSTALL --build`,
      confirmed to contain `libs/x64/{didgpu.dll, didgpu_cuda.dll,
      cudart64_12.dll}` (so `--build` does package the bundled DLLs — no
      silent breakage). Sits at `C:\Users\ammonsj\didgpu_release\` (outside
      the repo), ready to upload.
- [ ] Upload that `.zip` to a GitHub Release (needs your auth), e.g.
      `gh release create v0.1.0 C:\Users\ammonsj\didgpu_release\didgpu_0.1.0.zip`
- [ ] Add the repo to r-universe for the CPU baseline.
- [ ] README install section for both channels.
