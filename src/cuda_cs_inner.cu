// ============================================================================
// CUDA scaffold for the Callaway-Sant'Anna per-(g, t) inner regressions.
//
// The CS framework computes ATT(g, t) for each (cohort g, calendar t)
// cell. With C cohorts and T calendar times, there are up to C * T
// cells, each one an independent regression on a subset of units
// (treated cohort + control group).
//
// All three inner estimators (OR / IPW / DR) share the same compute
// pattern at the cell level:
//   1. Build the per-cell design matrix X (n_cell x k) and outcome
//      change delta (n_cell)
//   2. Solve a normal-equations problem: beta = (X' X)^{-1} X' delta
//      (OR), or weighted variants (IPW / DR)
//   3. Predict, compute ATT, accumulate
//
// GPU acceleration approach:
//   - Stack all C * T cells as a "batched least-squares" problem,
//     one batch per cell.
//   - Use cuBLAS gemmStridedBatched or gemmBatched for X' X and X' y.
//   - Use cuSOLVER cusolverDnDpotrsBatched for the Cholesky solve.
//   - Per-cell cell-size varies; pad to a common max size or use
//     variable-length batched routines.
//
// For the typical scale (C in the tens, T in the tens, k a handful of
// covariates), the per-batch work is small; the win comes from amortising
// kernel-launch overhead across hundreds of cells in one launch.
//
// STATUS: scaffolded. The host launcher signature is defined; the
// cuBLAS / cuSOLVER calls themselves are TODOs to be filled in once
// the package is built against a real CUDA Toolkit. The R-side
// fallback (R/cs_methods.R) handles per-cell regression via base R
// `qr.solve()` and `stats::glm.fit()` and is already what runs locally.
// ============================================================================

#ifdef HAS_CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

// Host launcher for the batched-OR inner regression.
//
// Inputs:
//   d_X     : device buffer, n_total x (k + 1) (row-major), augmented
//             with intercept column.
//   d_y     : device buffer, n_total (outcome change delta per row).
//   d_offsets : device buffer, n_cells + 1 (prefix-sum of per-cell row
//               counts; cell c lives at [offsets[c], offsets[c+1])).
//   d_D     : device buffer, n_total (1 = treated row, 0 = control row).
//   d_att_out : device buffer, n_cells (written by the kernel).
//
// Returns 0 on success, nonzero CUDA error code otherwise.
//
// TODO: implement using cublasGemmStridedBatchedEx + cusolverDnDpotrsBatched.
extern "C" int didgpu_cuda_cs_batched_or(
    const double* d_X,
    const double* d_y,
    const int*    d_offsets,
    const int*    d_D,
    int n_total, int n_cells, int k_features,
    double* d_att_out) {
  (void)d_X; (void)d_y; (void)d_offsets; (void)d_D;
  (void)n_total; (void)n_cells; (void)k_features; (void)d_att_out;
  // Scaffolded; the production implementation should:
  //   1. For each cell c, form X_c' X_c (k+1 x k+1) and X_c' y_c (k+1)
  //      using cublasDgemmStridedBatched (or the variable-stride
  //      version since cell sizes can differ).
  //   2. Cholesky-decompose every X' X via cusolverDnDpotrfBatched.
  //   3. Solve via cusolverDnDpotrsBatched to get beta_c per cell.
  //   4. Compute fitted_c = X_c %*% beta_c via batched gemm.
  //   5. ATT_c = mean(y_c[D=1] - fitted_c[D=1]) — done via per-cell
  //      reduction (custom kernel or cub::DeviceReduce).
  //   6. Write d_att_out[c] = ATT_c.
  //
  // The kernel-launch-overhead savings come from doing all cells in
  // one call. For C * T = 100 cells with 5 covariates and 200 obs per
  // cell, the GPU should finish in well under 1ms; CPU is ~10ms.
  return -1;  // not implemented
}

// Host launcher for the batched-DR inner regression (same shape as OR
// but adds the propensity-score IPW correction).
extern "C" int didgpu_cuda_cs_batched_dr(
    const double* d_X,
    const double* d_y,
    const int*    d_offsets,
    const int*    d_D,
    int n_total, int n_cells, int k_features,
    double* d_att_out) {
  (void)d_X; (void)d_y; (void)d_offsets; (void)d_D;
  (void)n_total; (void)n_cells; (void)k_features; (void)d_att_out;
  // TODO: same primitives as the OR case, plus a per-cell logistic
  // regression for the propensity score (Newton-Raphson is the
  // standard; ~10 iterations, each a Cholesky solve).
  return -1;  // not implemented
}

#endif  // HAS_CUDA
