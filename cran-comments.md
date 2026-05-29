# cran-comments.md

## Resubmission of 0.1.1 (now 0.1.2)

This is a resubmission of the first version (the previous pretest for
0.1.1 returned 1 ERROR + 1 NOTE on both Windows and Debian r-devel; see
fixes below). Bumped to 0.1.2 per CRAN convention on resubmission.

### What the pretest flagged, and what we fixed

1. **ERROR — `test-parallel.R` exceeded the 2-core cap.**
   The test honoured `_R_CHECK_LIMIT_CORES_` explicitly but the env-var
   detection did not fire in CRAN's pretest environment, so 4 worker
   processes were spawned and `.check_ncores()` (correctly) rejected.
   *Fix:* both `test_that` blocks in `tests/testthat/test-parallel.R`
   now call `skip_on_cran()`. The parallel path is exercised in our
   GitHub Actions CI matrix on every push.

2. **NOTE — invalid file URIs in README.md.**
   Three README links pointed at files excluded from the source tarball
   (`WINDOWS_BUILD_STATUS.md`, `BENCHMARKS.md`) or out-of-tree
   (`../NOTES_did_gpu_checkpointed.Rmd`).
   *Fix:* the first two now use absolute GitHub URLs so they resolve
   from CRAN's rendered README; the out-of-tree link was removed.

3. **NOTE — "Possibly misspelled words in DESCRIPTION".**
   The flagged tokens are domain terminology and proper names
   (`Chaisemartin`, `D'Haultfoeuille`, `multivalued`, `checkpointing`,
   `reimplementation`, `backends`, `resumable`, `de`) plus software
   names (`Rcpp`, `CUDA`). Software names are now single-quoted per
   CRAN convention (`'Rcpp'`, `'CUDA'`). The remaining tokens are
   intentional and correct; we politely ask that the NOTE be tolerated.

4. **NOTE (Windows only) — "Overall checktime 21 min > 10 min".**
   The package compiles a non-trivial amount of CUDA kernel code via
   `nvcc` and runs an end-to-end test suite covering five estimator
   families; the long checktime is dominated by the install step
   itself, which we cannot meaningfully shrink without disabling the
   GPU backend. CPU-bound tests run quickly; expensive GPU and
   reference-parity tests already `skip_on_cran()`. We ask that the
   checktime NOTE be tolerated.

`didgpu` is a heterogeneity-robust dynamic difference-in-differences
estimator (de Chaisemartin and D'Haultfoeuille, 2024) with per-cell
checkpointing, cluster bootstrap, and *optional* GPU acceleration via
NVIDIA CUDA. It is a clean-room reimplementation of the inner algorithms
of the already-on-CRAN reference package `DIDmultiplegtDYN`; bit-for-bit
numerical equivalence is verified by an extensive randomized
differential test suite (tens of thousands of panels across the
weighted x flag matrix, all backends agree).

`didgpu` is a heterogeneity-robust dynamic difference-in-differences
estimator (de Chaisemartin and D'Haultfoeuille, 2024) with per-cell
checkpointing, cluster bootstrap, and *optional* GPU acceleration via
NVIDIA CUDA. It is a clean-room reimplementation of the inner algorithms
of the already-on-CRAN reference package `DIDmultiplegtDYN`; bit-for-bit
numerical equivalence is verified by an extensive randomized
differential test suite (tens of thousands of panels across the
weighted x flag matrix, all backends agree).

## Test environments

- Local Windows 11, R 4.6.0, Rtools45, optional CUDA Toolkit 13.2 +
  RTX PRO Blackwell (sm_120) GPU available.
- Working build also exercised on Linux via r-universe
  (https://joshuaammons.r-universe.dev/didgpu).

## R CMD check results

Local `R CMD check --as-cran` (Windows + Rtools45, R 4.6.0):

```
Status: 0 ERRORs, 0 WARNINGs, 1 NOTE
```

The single NOTE is the standard new-submission notice:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Joshua Ammons <jdammons89@gmail.com>'
New submission
```

(Two further NOTEs that fire only in my local environment are also
present and are explained below for transparency; they are not relevant
on CRAN's build farms.)

### Local-only notes (will NOT fire on CRAN)

* `top-level files`: 'pandoc' is not installed on my local Windows
  machine, so `R CMD check` cannot validate `README.md` / `NEWS.md`.
  CRAN's machines have pandoc.

* `compiled code`: when the optional NVIDIA CUDA Toolkit is present at
  install time, `src/Makevars.win` bundles a small set of CUDA runtime
  DLLs (`cudart64_*`, `cublas64_*`, etc.) into `inst/libs/x64` so the
  GPU backend works on the user's machine without additional `PATH`
  configuration. CRAN's build farm does not carry the CUDA toolkit, so
  the bundling step is a no-op there and this NOTE does not fire.

## CUDA is an *optional* system requirement

The package is fully functional **without** NVIDIA CUDA installed. The
pure-R and Rcpp (C++) backends provide all features and pass all tests.

The CUDA backend is detected at install time by a `configure` /
`configure.win` script and skipped cleanly when absent. The C++/R code
remains identical; only the `.cu` files are skipped from compilation
when no CUDA toolkit is found. The package installs and loads cleanly
on CUDA-less systems.

This optional system requirement is declared in DESCRIPTION:

```
SystemRequirements: GNU make; optional NVIDIA CUDA Toolkit (>= 12.0)
                    with cuBLAS and cuSOLVER for the GPU backend
```

Tests that require GPU code are guarded by both `skip_on_cran()` (a
standard testthat skip) and an internal `didgpu_has_cuda_support()`
runtime check, so they never run on CRAN's machines.

## Downstream dependencies

This is a new package; there are no reverse dependencies on CRAN.

## Repository / Maintainer

* URL:        https://github.com/JoshuaAmmons/didgpu
* BugReports: https://github.com/JoshuaAmmons/didgpu/issues
* Maintainer: Joshua Ammons <jdammons89@gmail.com>
