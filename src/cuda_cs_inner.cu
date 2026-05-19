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
//   1. Build the per-cell design matrix X (n_cell x p) and outcome
//      change delta (n_cell)
//   2. Solve a normal-equations problem: beta = (X' X)^{-1} X' delta
//      (OR), or weighted variants (IPW / DR)
//   3. Predict, compute ATT, accumulate
//
// GPU acceleration approach (Phase 2 target):
//   - Stack all C * T cells as a "batched least-squares" problem,
//     one batch per cell.
//   - Use cuBLAS gemmStridedBatched or gemmBatched for X' X and X' y.
//   - Use cuSOLVER cusolverDnDpotrsBatched for the Cholesky solve.
//   - Per-cell cell-size varies; offsets[] tells the kernel where
//     each cell lives in the concatenated buffers.
//
// For the typical scale (C in the tens, T in the tens, p a handful of
// covariates), the per-batch work is small; the win comes from
// amortising kernel-launch overhead across hundreds of cells in one
// launch.
//
// STATUS: scaffold. Signature matches inst/include/didgpu_cuda_api.h
// (canonical Phase-1 ABI). Body returns -1 ("not implemented") so
// the R-side dispatch falls back to the per-cell R loop. Phase 2
// (tasks #82-#85) fills in the cuBLAS/cuSOLVER calls.
// ============================================================================

#ifdef HAS_CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

// Host launcher for the CS batched inner regression.
//
// Signature matches `didgpu_cuda_cs_inner_batched` in
// inst/include/didgpu_cuda_api.h. All pointers are HOST pointers; the
// kernel allocates device memory internally per the C-ABI contract.
//
// Cell c lives at rows [X_offsets[c], X_offsets[c+1]) in X_concat /
// Y_concat / W_concat. Y_offsets and W_offsets are passed separately
// because the contract permits weights/outcomes to live in arrays
// indexed slightly differently than X — but in the canonical layout
// the three offset arrays are identical pointers.
//
// est_method: 0 = OR, 1 = IPW, 2 = DR
//
// out_att      : length n_cells (one ATT estimate per cell)
// out_influence: optional, length n_units * n_cells (row-major,
//                unit-major). Pass NULL to skip influence-function
//                computation. Phase 2 will populate this for
//                multiplier-bootstrap support.
//
// Returns: 0 on success; nonzero error code on failure. The R-side
// dispatch treats any nonzero return as "fall back to R impl".
extern "C" int didgpu_cuda_cs_inner_batched(
    const double* X_concat, const int* X_offsets,
    const double* Y_concat, const int* Y_offsets,
    const double* W_concat, const int* W_offsets,
    int n_cells, int p, int n_units,
    int est_method,
    double* out_att, double* out_influence) {
  (void)X_concat;    (void)X_offsets;
  (void)Y_concat;    (void)Y_offsets;
  (void)W_concat;    (void)W_offsets;
  (void)n_cells;     (void)p;            (void)n_units;
  (void)est_method;
  (void)out_att;     (void)out_influence;
  // -----------------------------------------------------------------
  // Phase 2 implementation plan (tasks #82-#85):
  //   1. cudaMalloc + H2D-copy the three concatenated buffers plus the
  //      offsets array (single array — Y_offsets and W_offsets equal
  //      X_offsets in the canonical layout).
  //   2. For each cell c, form X_c' X_c (p x p) and X_c' y_c (p) via
  //      cublasDgemmStridedBatched with strides derived from offsets.
  //      Per-cell sizes vary, so we use the variable-stride variant
  //      (or pad to max_cell_size with a mask).
  //   3. Cholesky-decompose every X' X via cusolverDnDpotrfBatched.
  //   4. Solve via cusolverDnDpotrsBatched to get beta_c per cell.
  //   5. For OR (method=0): ATT_c = mean(delta[D=1]) - mean(X_c[D=1]
  //      @ beta_c). For IPW (method=1): include the propensity-weight
  //      reweighting term. For DR (method=2): combine OR fitted-
  //      counterfactual with IPW-weighted residual correction.
  //   6. If out_influence != NULL, also write per-unit per-cell IF
  //      values (Phase 2 task #85 — required for multiplier-bootstrap
  //      SEs).
  //
  // Until that work lands, returning -1 keeps the R-side dispatch
  // honest: any backend = "cuda" call silently falls back to the
  // existing per-cell R loop in .cs_compute_att_gt.
  // -----------------------------------------------------------------
  return -1;  // Phase 2: not implemented
}

#endif  // HAS_CUDA
