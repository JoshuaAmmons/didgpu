/* ============================================================================
 * didgpu_cuda_api.h
 *
 * Pure-C ABI contract between the two halves of the didgpu binary on Windows
 * (and the single-DLL build on Linux/Mac, which uses the same header).
 *
 * On Windows:
 *   didgpu.dll        — built by Rtools g++ (MinGW). R-facing, no CUDA deps.
 *   didgpu_cuda.dll   — built by nvcc (MSVC under the hood). Owns all CUDA
 *                       state: device memory, streams, cuBLAS/cuSOLVER
 *                       handles. MinGW links against an import library.
 *
 * On Linux/Mac:
 *   didgpu.so         — single shared object containing both layers; the
 *                       C-API here is still the implementation boundary so
 *                       the two-tier separation is preserved logically.
 *
 * ABI rules — STRICTLY ENFORCED:
 *   1. extern "C" only. No C++ name mangling across the boundary.
 *   2. POD types only: int, double, float, pointers to POD, const char*.
 *      NO C++ types (no std::*), no Rcpp types, no STL containers.
 *   3. No exceptions cross the boundary. All functions return int error
 *      codes; 0 = success, non-zero = error. Use didgpu_cuda_last_error()
 *      to get a human-readable message for the most recent failure.
 *   4. All pointer arguments are HOST pointers. The CUDA DLL allocates and
 *      frees device memory internally. The caller never sees device
 *      pointers. (This is what makes the Rcpp side CUDA-free.)
 *   5. Matrix layouts are ROW-MAJOR (so they match the kernel-side
 *      conventions). The R side has to transpose if it has column-major
 *      data — usually a one-time cost per call.
 *   6. Sizes are int (signed 32-bit). No size_t (different width on Win64
 *      vs Linux LP64 — int is portable).
 *   7. RNG seeds are uint64_t for cross-platform consistency.
 *
 * Versioning:
 *   didgpu_cuda_abi_version() returns a monotonically increasing integer.
 *   The Rcpp side checks this at package load time and refuses to dispatch
 *   if didgpu_cuda.dll is from a different ABI generation. Bump this
 *   integer whenever any signature in this header changes incompatibly.
 *
 * Linkage:
 *   On Windows, didgpu_cuda.dll EXPORTS these symbols (via .def file).
 *   On Linux/Mac, these are ordinary extern "C" symbols visible across
 *   translation units within libdidgpu.so.
 * ============================================================================
 */

#ifndef DIDGPU_CUDA_API_H
#define DIDGPU_CUDA_API_H

#include <stdint.h>

/* DIDGPU_CUDA_API: visibility/export annotation.
 *
 *   - When building didgpu_cuda.dll, define DIDGPU_CUDA_BUILDING_DLL to
 *     expand to __declspec(dllexport) on Windows.
 *   - When using the header from didgpu.dll (the consumer), the default is
 *     __declspec(dllimport) on Windows.
 *   - On Linux/Mac it's just `extern`.
 */
#if defined(_WIN32) || defined(_WIN64)
#  if defined(DIDGPU_CUDA_BUILDING_DLL)
#    define DIDGPU_CUDA_API __declspec(dllexport)
#  elif defined(DIDGPU_CUDA_STATIC)
#    define DIDGPU_CUDA_API
#  else
#    define DIDGPU_CUDA_API __declspec(dllimport)
#  endif
#else
#  define DIDGPU_CUDA_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * ABI version. Bump when any signature changes incompatibly.
 * --------------------------------------------------------------------------- */
#define DIDGPU_CUDA_ABI_VERSION 1

DIDGPU_CUDA_API int didgpu_cuda_abi_version(void);

/* ---------------------------------------------------------------------------
 * Error reporting.
 * Standard codes used throughout:
 *   0  success
 *  -1  CUDA runtime error (cuda*  call failed; see last_error)
 *  -2  cuBLAS/cuSOLVER/cuRAND error
 *  -3  argument / shape mismatch detected by the DLL
 *  -4  out of GPU memory
 *  -5  no CUDA-capable device available
 *  -9  unspecified failure
 * --------------------------------------------------------------------------- */
DIDGPU_CUDA_API const char* didgpu_cuda_last_error(void);

/* ---------------------------------------------------------------------------
 * Probe and device info.
 * --------------------------------------------------------------------------- */

/* Returns 1 if a CUDA-capable device is present and initialization
 * succeeded; 0 otherwise. Never throws. Cheap to call repeatedly. */
DIDGPU_CUDA_API int didgpu_cuda_available(void);

/* Writes one-line device summary into buf ("RTX 4000 Ada, 12.3 GB, CC 8.9").
 * Returns 0 on success. buf must be at least 256 bytes. */
DIDGPU_CUDA_API int didgpu_cuda_device_info(char* buf, int buflen);

/* Optional: hint that we're done with all GPU work for a while. Frees the
 * internal workspace pool. Safe to call any time; subsequent CUDA calls
 * will re-allocate on demand. */
DIDGPU_CUDA_API int didgpu_cuda_release_workspace(void);

/* ---------------------------------------------------------------------------
 * didgpu (DIDmultiplegtDYN-equivalent) — per (event-time, direction) DiD.
 *
 * Runs ONE inner solve for a single (k, direction) and returns the scalar
 * DiD point estimate in *out_did. All inputs are host arrays as produced
 * by R/core_r.R::.core_one_event_time(). See R-side comments there for the
 * meaning of each variable.
 *
 * Lengths:
 *   outcome, N_gt, row_to_g, row_to_t, cohort_key : length n_rows
 *   F_g, S_g, T_g, L_g                            : length n_groups
 *   out_did                                       : length 1 (scalar)
 * --------------------------------------------------------------------------- */
DIDGPU_CUDA_API int didgpu_cuda_did_one_event_time(
    const double* outcome,
    const double* N_gt,
    const int* row_to_g,
    const int* row_to_t,
    const int* cohort_key,
    const int* F_g,
    const int* S_g,
    const int* T_g,
    const int* L_g,
    int n_rows, int n_groups, int n_cohorts,
    int k, int direction, double G_over_Ninc,
    double* out_did);

/* ---------------------------------------------------------------------------
 * fect family kernels.
 *
 * Y, M and any matrix output is ROW-MAJOR with shape (n_units, n_periods),
 * i.e. element (i, t) at offset i * n_periods + t. The R side has to
 * transpose column-major R matrices on the way in/out (cheap).
 * --------------------------------------------------------------------------- */

/* fect_fe: iterative two-way demeaning over controls-only cells.
 *   M[i, t] = 1 means cell is excluded (treated/missing); 0 means included.
 *   Outputs: alpha (length n_units), xi (length n_periods),
 *            *out_iter (iterations used), *out_delta (final max-change). */
DIDGPU_CUDA_API int didgpu_cuda_fect_fe(
    const double* Y_rm, const int* M_rm,
    int n_units, int n_periods,
    double tol, int max_iter,
    double* out_alpha, double* out_xi,
    int* out_iter, double* out_delta);

/* Truncated SVD: A = U * diag(S) * V^T, keeping only the top r components.
 * A is (m, n) row-major; out_U_rm is (m, r); out_S is (r,); out_V_rm is (n, r). */
DIDGPU_CUDA_API int didgpu_cuda_fect_svd_truncated(
    const double* A_rm,
    int m, int n, int r,
    double* out_U_rm, double* out_S, double* out_V_rm);

/* Soft-thresholded SVD reconstruction for matrix completion (Athey et al.):
 *   A_hat = U * diag(max(s - lambda, 0)) * V^T
 * out_A_rm is the (m, n) reconstructed matrix. */
DIDGPU_CUDA_API int didgpu_cuda_fect_svd_softthreshold(
    const double* A_rm,
    int m, int n, double lambda,
    double* out_A_rm);

/* Fused ife alternation: alternates fect_fe step + truncated SVD step until
 * convergence. Returns the final alpha, xi, and rank-r factor matrices. */
DIDGPU_CUDA_API int didgpu_cuda_fect_ife(
    const double* Y_rm, const int* M_rm,
    int n_units, int n_periods, int r,
    double tol, int max_iter,
    double* out_alpha, double* out_xi,
    double* out_U_rm, double* out_S, double* out_V_rm,
    int* out_iter, double* out_delta);

/* mc CV inner loop: K-fold CV over a grid of lambda values, returns the
 * lambda that minimizes held-out MSE. */
DIDGPU_CUDA_API int didgpu_cuda_fect_mc_cv(
    const double* Y_rm, const int* M_rm,
    int n_units, int n_periods,
    const double* lambda_grid, int n_lambdas,
    int n_folds, uint64_t seed,
    double* out_lambda_best, double* out_cv_mse_grid);

/* ---------------------------------------------------------------------------
 * Callaway-Sant'Anna inner regression (per (g, t) cell), batched.
 *
 * Cells are passed as concatenated arrays plus per-cell offset arrays.
 * For cell c, the rows are at positions [X_offsets[c], X_offsets[c+1]) in
 * X_concat/Y_concat/W_concat.
 *
 *   est_method: 0=OR, 1=IPW, 2=DR (Sant'Anna-Zhao 2020 doubly-robust)
 *   n_cells: total number of (g, t) cells
 *   p:       number of covariates (columns of X)
 *
 * Outputs:
 *   out_att (length n_cells)            — ATT(g, t)
 *   out_influence (n_units * n_cells)   — per-unit, per-cell influence
 *                                          function values, row-major
 *                                          (unit-major). NULL to skip. */
DIDGPU_CUDA_API int didgpu_cuda_cs_inner_batched(
    const double* X_concat, const int* X_offsets,
    const double* Y_concat, const int* Y_offsets,
    const double* W_concat, const int* W_offsets,
    int n_cells, int p, int n_units,
    int est_method,
    double* out_att, double* out_influence);

/* ---------------------------------------------------------------------------
 * Bootstrap primitives (used by didgpu AND didgpu_cs).
 * --------------------------------------------------------------------------- */

/* Cluster bootstrap: draws B resamples of clusters and recomputes the
 * point estimate from per-unit influence functions.
 *
 *   influence: (n_units, n_dims) row-major
 *   cluster_id: (n_units,) — integer cluster IDs in [0, n_clusters)
 *   B: number of bootstrap replicates
 *   seed: RNG seed (must be same across runs for determinism)
 *
 * Output:
 *   out_estimates: (B, n_dims) row-major; out_estimates[b * n_dims + d]
 *                  is replicate b's estimate for dim d.
 *
 * Determinism: results are bit-identical for fixed (seed, n_units, n_clusters,
 * B, n_dims) on the same GPU architecture; we use a fixed-stride per-replicate
 * cuRAND substream. */
DIDGPU_CUDA_API int didgpu_cuda_cluster_bootstrap(
    const double* influence, int n_units, int n_dims,
    const int* cluster_id, int n_clusters,
    int B, uint64_t seed,
    double* out_estimates);

/* Multiplier (wild) bootstrap: per-unit influence * Rademacher (or N(0,1)).
 *
 *   mult_kind: 0 = Rademacher (+1/-1), 1 = N(0, 1)
 *
 * Same output shape as cluster bootstrap. Much faster than cluster bs
 * for large B because it's pure SAXPY-style ops with no resampling. */
DIDGPU_CUDA_API int didgpu_cuda_multiplier_bootstrap(
    const double* influence, int n_units, int n_dims,
    int B, int mult_kind, uint64_t seed,
    double* out_estimates);

/* ---------------------------------------------------------------------------
 * TestMechs partial-density bootstrap.
 *
 * D, M, Y: integer-coded observations of length n. D in {0, 1};
 *          M in {0, ..., K-1}; Y in {0, ..., dy-1}.
 *   method: 0 = nonparametric (multinomial resample),
 *           1 = Bayesian Dirichlet posterior draw.
 *   out_beta_obs: (B, 2 * K * dy) row-major. Each row is the partial-density
 *                 vector for one bootstrap draw — stacked [D=0 block, D=1
 *                 block], each of length K*dy. Each per-D block must sum to 1.
 * --------------------------------------------------------------------------- */
DIDGPU_CUDA_API int didgpu_cuda_testmechs_bootstrap(
    const int* D, const int* M, const int* Y,
    int n, int K, int dy,
    int B, int method, uint64_t seed,
    double* out_beta_obs);

/* ---------------------------------------------------------------------------
 * Leave-one-out batched refit (Phase 3 target).
 *
 * Given the K leave-out subsets implied by `cohort_id` and `leave_outs`,
 * runs the per-(g,t) inner regression for ALL K leave-out replicates in
 * one GPU launch. Used by didgpu_loo when backend = "cuda".
 *
 *   cohort_id:  per-unit cohort labels (length n_units), -1 = never-treated
 *   leave_outs: distinct cohort values to leave out (length K)
 *
 * Output:
 *   out_headlines: (K,) — headline estimate with each cohort left out. */
DIDGPU_CUDA_API int didgpu_cuda_loo_batched(
    const double* X_concat, const int* X_offsets,
    const double* Y_concat, const int* Y_offsets,
    const double* W_concat, const int* W_offsets,
    int n_cells, int p, int n_units,
    const int* cohort_id,
    const int* leave_outs, int K,
    int est_method,
    double* out_headlines);

#ifdef __cplusplus
}
#endif

#endif /* DIDGPU_CUDA_API_H */
