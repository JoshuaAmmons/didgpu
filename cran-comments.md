# cran-comments.md

This is a first submission of `didgpu` (version 0.1.0).

## Test environments

- Local Windows 10, R 4.4.1 (release)
- (Add R-hub / win-builder / GitHub Actions runs here before submission)

## R CMD check results

Local `R CMD check --as-cran` reports:

```
Status: 1 WARNING, 1 NOTE  (after addressing CRAN feasibility URL note)
```

### WARNING (intentional)

```
* checking if this is a source package ... WARNING
Subdirectory 'src' contains:
  cuda_didkernel.cu cuda_saxpy.cu
These are unlikely file names for src files.
```

The `.cu` extension is NVIDIA's standard, documented extension for
CUDA C++ source files (the same way `.f90` is the standard extension
for Fortran 90). These two files are CUDA C++ kernels that are
*conditionally* compiled by `nvcc` (the NVIDIA CUDA compiler) at
install time when a CUDA Toolkit is present on the user's machine.

The R CMD check heuristic that produced this WARNING accepts the
extensions `.c`, `.cc`, `.cpp`, `.h`, `.hpp`, `.f`, `.f90`, `.f95`,
`.m`, `.mm`, but not `.cu`. So the WARNING fires whenever a package
ships CUDA sources, even when their compilation is fully optional.

This is the canonical pattern previously used by the long-lived
`gputools` package on CRAN (which shipped CUDA-only kernels with
`.cu` files in `src/` for many years). The didgpu build is even more
defensive: `src/Makevars` and `src/Makevars.win` auto-detect `nvcc`
and skip the `.cu` files entirely when it is absent, so the package
installs successfully on machines without any CUDA toolchain (which
is the default install experience). The `.cu` files cannot affect
non-CUDA installs.

We considered the alternative of shipping `.cu` files in `inst/cuda/`
and copying them into `src/` via a `configure` script at install
time, but Rtools' shell environment on Windows makes that pattern
fragile, and it does not gain anything in functional terms — the
`.cu` files are still present in the tarball regardless.

### NOTE (intentional)

```
* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.
```

Declared in `DESCRIPTION` via `SystemRequirements: GNU make. ...`.

## Downstream dependencies

This is a new package, so there are no reverse dependencies.

## Maintainer / contact

Joshua Ammons <jdammons89@gmail.com>
