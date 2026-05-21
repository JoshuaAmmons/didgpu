// ============================================================================
// CUDA SVD primitive for fect_ife and fect_mc.
//
// Uses cuSOLVER's cusolverDnDgesvdj (Jacobi SVD) for dense double-
// precision matrices. cuSOLVER is part of the standard CUDA Toolkit
// (linked via -lcusolver in Makevars).
//
// Two entry points:
//   1. didgpu_cuda_fect_svd_truncated: rank-r truncated SVD, used by ife.
//      Returns L = U_r * sqrt(D_r) and F = sqrt(D_r) * V_r^T directly,
//      avoiding the cost of building L and F separately on host.
//
//   2. didgpu_cuda_fect_svd_softthreshold: full SVD with soft-threshold
//      on singular values, then reconstructs Y_hat = U * D_st * V^T.
//      Used by mc.
//
// Both entry points handle the column-major <-> row-major conversion
// internally (cuSOLVER expects column-major; the R-side passes
// row-major).
//
// STATUS: scaffold. Compiles when nvcc + cuSOLVER are present. Untested
// locally without the full CUDA Toolkit install, but structured to
// mirror the R reference implementations in R/fect_ife.R and
// R/fect_mc.R line-by-line.
// ============================================================================

#ifdef HAS_CUDA

#if defined(DIDGPU_LITE)
// ===========================================================================
// LITE build (-DDIDGPU_LITE): compiled WITHOUT cuSOLVER/cuBLAS so the
// resulting didgpu_cuda.dll depends only on cudart64 (~0.5 MB) instead of
// dragging in cuBLAS + cuBLASLt + cuSOLVER + cuSPARSE + nvJitLink (~1.1 GB
// of redistributable DLLs). This is the build used for the distributable
// Windows binary — see RELEASE_PLAN.md.
//
// The fect GPU SVD path is therefore unavailable: both entry points return
// the sentinel -99, which the Rcpp wrapper turns into R_NilValue, so the R
// side falls back to LAPACK svd(). This costs nothing in practice — the SVD
// path is size-gated (only ever engaged on very large balanced panels) and
// the fallback is numerically identical, just on the CPU.
// ===========================================================================
extern "C" int didgpu_cuda_fect_svd_truncated(
    const double*, int, int, int, double*, double*) { return -99; }
extern "C" int didgpu_cuda_fect_svd_softthreshold(
    const double*, int, int, double, double*, int*) { return -99; }

#else  // ---- full build (with cuSOLVER/cuBLAS) --------------------------

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cmath>
#include <algorithm>

// Helper: column-major <-> row-major in place. cuSOLVER works in
// column-major; the rest of didgpu uses row-major. We transpose once
// on input and once on output.
//
// For an (m x n) matrix, layout swap means swapping the dimensions:
// transposed[(j, i)] = original[(i, j)].
__global__ void k_transpose_rm_to_cm(
    const double* __restrict__ src, double* __restrict__ dst,
    int m, int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = m * n;
  if (idx >= total) return;
  const int i = idx / n;            // row in row-major
  const int j = idx % n;            // col in row-major
  dst[j * m + i] = src[idx];        // column-major: col j, row i
}

__global__ void k_transpose_cm_to_rm(
    const double* __restrict__ src, double* __restrict__ dst,
    int m, int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = m * n;
  if (idx >= total) return;
  const int i = idx % m;            // row in column-major
  const int j = idx / m;            // col in column-major
  dst[i * n + j] = src[idx];        // row-major: row i, col j
}

// Element-wise sqrt of the first r entries of S (skips the remaining
// singular values; called once per SVD).
__global__ void k_sqrt_first_r(const double* __restrict__ S,
                                 double* __restrict__ S_sqrt, int r) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= r) return;
  const double v = S[i];
  S_sqrt[i] = (v > 0.0) ? sqrt(v) : 0.0;
}

// Scale columns of a column-major (m x r) matrix: out[i, j] =
// in[i, j] * scale[j]. Used to build L = U_r * diag(sqrt(D_r)) from
// the SVD's U output.
__global__ void k_scale_cols(const double* __restrict__ in,
                               const double* __restrict__ scale,
                               double* __restrict__ out,
                               int m, int r) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = m * r;
  if (idx >= total) return;
  const int j = idx / m;     // column in column-major
  out[idx] = in[idx] * scale[j];
}

// Scale ROWS of a row-major (r x n) matrix: out[i, j] = in[i, j] *
// scale[i]. Used to build F = diag(sqrt(D_r)) * V_r^T.
__global__ void k_scale_rows_rm(const double* __restrict__ in,
                                  const double* __restrict__ scale,
                                  double* __restrict__ out,
                                  int r, int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = r * n;
  if (idx >= total) return;
  const int i = idx / n;     // row in row-major
  out[idx] = in[idx] * scale[i];
}


// Fused transpose+scale: build F = sqrt(D_r) * V_r^T directly from
// the col-major V matrix output by cuSOLVER, without an intermediate
// buffer.
//
// Inputs:
//   V_cm   : (n x r) column-major (the first r columns of cuSOLVER's
//            V output). ld = n.
//   sqrt_S : length r vector of sqrt(singular values).
// Output:
//   F_rm   : (r x n) row-major. F_rm[a, b] = sqrt_S[a] * V_cm[b, a].
//
// Index identity: F_rm[a*n + b] in row-major is at the same linear
// index as V_cm[b + a*n] in col-major (since both equal a*n + b).
// So the kernel reads V_cm[idx] for idx in [0, r*n), interprets the
// "row" coordinate of F as idx / n, multiplies by sqrt_S of that row,
// and writes to F_rm[idx]. No actual data shuffle needed beyond the
// element-wise multiplication.
__global__ void k_build_F_from_Vcm(const double* __restrict__ V_cm,
                                     const double* __restrict__ sqrt_S,
                                     double*       __restrict__ F_rm,
                                     int r, int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = r * n;
  if (idx >= total) return;
  const int a = idx / n;       // row of F == column of V_cm
  F_rm[idx] = sqrt_S[a] * V_cm[idx];
}


// ---------------------------------------------------------------------------
// Truncated rank-r SVD: input M (m x n, row-major), output L (m x r)
// and F (r x n) such that L * F is a rank-r approximation of M.
// Specifically L = U_r * sqrt(D_r), F = sqrt(D_r) * V_r^T.
// ---------------------------------------------------------------------------
// INTERNAL device-side helper: every pointer here is a device pointer.
// The host-facing didgpu_cuda_fect_svd_truncated() wrapper below owns the
// device memory and exposes only host pointers across the C ABI.
static int fect_svd_truncated_dev(
    const double* d_M_rm,   // device, row-major (m x n)
    int m, int n, int r,
    double* d_L_out_rm,     // device, row-major (m x r), preallocated
    double* d_F_out_rm) {   // device, row-major (r x n), preallocated

  if (r <= 0 || r > std::min(m, n)) return -1;

  // ----------------------------------------------------------------
  // ALL local variables declared up front. Required because we use
  // `goto fail` for error handling, and C++ forbids jumping past a
  // variable that has a non-trivial initializer into a scope where
  // it's expected to be initialized. Hoisting every declaration here
  // (with explicit nullptr init) makes every goto target safe and
  // also lets `fail:` free every potentially-allocated pointer.
  // ----------------------------------------------------------------
  cudaError_t e;
  cusolverStatus_t st;
  cusolverDnHandle_t h = nullptr;
  gesvdjInfo_t params = nullptr;
  double* d_M_cm   = nullptr;
  double* d_U      = nullptr;
  double* d_V      = nullptr;
  double* d_S      = nullptr;
  double* d_work   = nullptr;
  int*    d_info   = nullptr;
  double* d_sqrt_S = nullptr;
  double* d_L_cm   = nullptr;
  // (formerly d_F_temp) — no longer needed; k_build_F_from_Vcm fuses
  // the transpose-from-col-major and the per-row sqrt_S scaling into
  // one element-wise pass over d_V.
  int lwork = 0;
  int total = 0;

  st = cusolverDnCreate(&h);
  if (st != CUSOLVER_STATUS_SUCCESS) return -2;

  // cuSOLVER wants column-major; transpose input.
  e = cudaMalloc((void**)&d_M_cm, sizeof(double) * m * n);
  if (e != cudaSuccess) { cusolverDnDestroy(h); return -3; }
  total = m * n;
  k_transpose_rm_to_cm<<<(total + 255) / 256, 256>>>(d_M_rm, d_M_cm, m, n);
  cudaDeviceSynchronize();

  // Allocate SVD outputs.
  e = cudaMalloc((void**)&d_U, sizeof(double) * m * m);
  if (e != cudaSuccess) goto fail;
  e = cudaMalloc((void**)&d_V, sizeof(double) * n * n);
  if (e != cudaSuccess) goto fail;
  e = cudaMalloc((void**)&d_S, sizeof(double) * std::min(m, n));
  if (e != cudaSuccess) goto fail;

  // Jacobi SVD parameters + workspace query.
  cusolverDnCreateGesvdjInfo(&params);
  cusolverDnXgesvdjSetTolerance(params, 1e-7);
  cusolverDnXgesvdjSetMaxSweeps(params, 100);

  st = cusolverDnDgesvdj_bufferSize(
      h, CUSOLVER_EIG_MODE_VECTOR, /*econ=*/1,
      m, n, d_M_cm, m, d_S, d_U, m, d_V, n,
      &lwork, params);
  if (st != CUSOLVER_STATUS_SUCCESS) goto fail;

  e = cudaMalloc((void**)&d_work, sizeof(double) * lwork);
  if (e != cudaSuccess) goto fail;
  e = cudaMalloc((void**)&d_info, sizeof(int));
  if (e != cudaSuccess) goto fail;

  st = cusolverDnDgesvdj(
      h, CUSOLVER_EIG_MODE_VECTOR, /*econ=*/1,
      m, n, d_M_cm, m, d_S, d_U, m, d_V, n,
      d_work, lwork, d_info, params);
  cudaFree(d_work); d_work = nullptr;
  cudaFree(d_info); d_info = nullptr;
  if (st != CUSOLVER_STATUS_SUCCESS) goto fail;

  // Compute sqrt(D_r) on device via a small kernel (no host roundtrip).
  cudaMalloc((void**)&d_sqrt_S, sizeof(double) * r);
  k_sqrt_first_r<<<(r + 31) / 32, 32>>>(d_S, d_sqrt_S, r);
  cudaDeviceSynchronize();

  // L = U[:, 1..r] * diag(sqrt(D_r)). U is m x m column-major; we want
  // the first r columns scaled by sqrt(D_r). Use a column-wise scaling
  // kernel, then transpose to row-major into d_L_out_rm.
  cudaMalloc((void**)&d_L_cm, sizeof(double) * m * r);
  k_scale_cols<<<(m * r + 255) / 256, 256>>>(d_U, d_sqrt_S, d_L_cm, m, r);
  cudaDeviceSynchronize();
  // Transpose L_cm (m x r col-major) -> L_rm (m x r row-major).
  k_transpose_cm_to_rm<<<(m * r + 255) / 256, 256>>>(d_L_cm, d_L_out_rm, m, r);

  // F = diag(sqrt(D_r)) * V[:, 1..r]^T. V is n x n column-major; we
  // want V[:, 1..r] (n x r) then transpose to (r x n), then scale
  // ROWS by sqrt(D_r). Fused into a single kernel via the index
  // identity F_rm[a*n + b] == V_cm[a*n + b] when V_cm is (n x r)
  // col-major (so its element (b, a) lives at b + a*n = a*n + b).
  k_build_F_from_Vcm<<<(r * n + 255) / 256, 256>>>(
      d_V, d_sqrt_S, d_F_out_rm, r, n);
  cudaDeviceSynchronize();

  cudaFree(d_L_cm);
  cudaFree(d_sqrt_S);
  cudaFree(d_U); cudaFree(d_V); cudaFree(d_S); cudaFree(d_M_cm);
  cusolverDnDestroyGesvdjInfo(params);
  cusolverDnDestroy(h);
  return 0;

fail:
  if (d_L_cm)   cudaFree(d_L_cm);
  if (d_sqrt_S) cudaFree(d_sqrt_S);
  if (d_work)   cudaFree(d_work);
  if (d_info)   cudaFree(d_info);
  if (d_U)      cudaFree(d_U);
  if (d_V)      cudaFree(d_V);
  if (d_S)      cudaFree(d_S);
  if (d_M_cm)   cudaFree(d_M_cm);
  if (params)   cusolverDnDestroyGesvdjInfo(params);
  cusolverDnDestroy(h);
  return -10;
}


// ---------------------------------------------------------------------------
// Soft-threshold scaling kernel: D_st[i] = max(D[i] - lambda, 0).
// Used by mc to soft-threshold the singular values in place.
// ---------------------------------------------------------------------------
__global__ void k_soft_threshold(double* d, int n, double lambda) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  const double v = d[i] - lambda;
  d[i] = (v > 0.0) ? v : 0.0;
}


// ---------------------------------------------------------------------------
// Full-SVD with soft-threshold reconstruction: Y_hat = U * D_st * V^T.
// Used by mc's per-iter update.
//
// Same scaffold structure as the truncated SVD above. Once the
// scaffolds are filled in, this is the kernel chain mc needs per iter.
// ---------------------------------------------------------------------------
// INTERNAL device-side helper (device pointers only); wrapped by the
// host-facing didgpu_cuda_fect_svd_softthreshold() below.
static int fect_svd_softthreshold_dev(
    const double* d_Y_complete_rm,
    int m, int n,
    double lambda,
    double* d_Y_hat_rm,
    int* out_n_nonzero) {
  // TODO(phase4): wire cudaError_t checks on every cudaMalloc/cudaMemcpy
  // call below. Currently this function silently ignores allocation
  // failures — fine for the scaffold, not OK for production.
  cusolverStatus_t st;
  cublasHandle_t hb = nullptr;
  cusolverDnHandle_t hs = nullptr;
  cublasCreate(&hb);
  st = cusolverDnCreate(&hs);
  if (st != CUSOLVER_STATUS_SUCCESS) { cublasDestroy(hb); return -1; }

  // Transpose input to column-major.
  double* d_Y_cm = nullptr;
  cudaMalloc((void**)&d_Y_cm, sizeof(double) * m * n);
  const int total = m * n;
  k_transpose_rm_to_cm<<<(total + 255) / 256, 256>>>(d_Y_complete_rm, d_Y_cm, m, n);
  cudaDeviceSynchronize();

  // Full SVD (econ = 1 returns U as m x min(m,n), V as n x min(m,n)).
  const int k = (m < n) ? m : n;
  double *d_U = nullptr, *d_V = nullptr, *d_S = nullptr;
  cudaMalloc((void**)&d_U, sizeof(double) * m * k);
  cudaMalloc((void**)&d_V, sizeof(double) * n * k);
  cudaMalloc((void**)&d_S, sizeof(double) * k);

  gesvdjInfo_t params = nullptr;
  cusolverDnCreateGesvdjInfo(&params);
  cusolverDnXgesvdjSetTolerance(params, 1e-7);
  cusolverDnXgesvdjSetMaxSweeps(params, 100);

  int lwork = 0;
  st = cusolverDnDgesvdj_bufferSize(
      hs, CUSOLVER_EIG_MODE_VECTOR, /*econ=*/1,
      m, n, d_Y_cm, m, d_S, d_U, m, d_V, n, &lwork, params);
  if (st != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_U); cudaFree(d_V); cudaFree(d_S); cudaFree(d_Y_cm);
    cusolverDnDestroyGesvdjInfo(params);
    cusolverDnDestroy(hs); cublasDestroy(hb);
    return -2;
  }
  double* d_work = nullptr;
  cudaMalloc((void**)&d_work, sizeof(double) * lwork);
  int* d_info = nullptr;
  cudaMalloc((void**)&d_info, sizeof(int));

  st = cusolverDnDgesvdj(
      hs, CUSOLVER_EIG_MODE_VECTOR, /*econ=*/1,
      m, n, d_Y_cm, m, d_S, d_U, m, d_V, n,
      d_work, lwork, d_info, params);
  cudaFree(d_work); cudaFree(d_info);
  if (st != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_U); cudaFree(d_V); cudaFree(d_S); cudaFree(d_Y_cm);
    cusolverDnDestroyGesvdjInfo(params);
    cusolverDnDestroy(hs); cublasDestroy(hb);
    return -3;
  }

  // Soft-threshold the singular values in place: D_st = max(D - lambda, 0).
  k_soft_threshold<<<(k + 255) / 256, 256>>>(d_S, k, lambda);
  cudaDeviceSynchronize();

  // Count the non-zero singular values on host (one-time small cost).
  if (out_n_nonzero) {
    double* h_S = new double[k];
    cudaMemcpy(h_S, d_S, sizeof(double) * k, cudaMemcpyDeviceToHost);
    int nz = 0;
    for (int i = 0; i < k; ++i) if (h_S[i] > 0) ++nz;
    *out_n_nonzero = nz;
    delete[] h_S;
  }

  // Y_hat = U * diag(D_st) * V^T. First scale columns of U by D_st:
  double* d_U_scaled = nullptr;
  cudaMalloc((void**)&d_U_scaled, sizeof(double) * m * k);
  k_scale_cols<<<(m * k + 255) / 256, 256>>>(d_U, d_S, d_U_scaled, m, k);
  cudaDeviceSynchronize();
  // Then Y_hat_cm = U_scaled %*% V^T via cuBLAS Dgemm.
  // U_scaled is m x k col-major, V is n x k col-major (we want V^T = k x n).
  double* d_Y_hat_cm = nullptr;
  cudaMalloc((void**)&d_Y_hat_cm, sizeof(double) * m * n);
  const double one = 1.0, zero = 0.0;
  cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_T,
              m, n, k,
              &one, d_U_scaled, m,
              d_V, n,
              &zero, d_Y_hat_cm, m);

  // Transpose Y_hat_cm (m x n col-major) to row-major output.
  k_transpose_cm_to_rm<<<(m * n + 255) / 256, 256>>>(d_Y_hat_cm, d_Y_hat_rm, m, n);

  cudaFree(d_U_scaled); cudaFree(d_Y_hat_cm);
  cudaFree(d_U); cudaFree(d_V); cudaFree(d_S); cudaFree(d_Y_cm);
  cusolverDnDestroyGesvdjInfo(params);
  cusolverDnDestroy(hs); cublasDestroy(hb);
  return 0;
}


// ===========================================================================
// Host-facing entry points (pure-C ABI). These own all device memory: they
// upload the HOST inputs, call the device-side helpers above, download the
// HOST outputs, and free everything. No device pointer crosses the boundary,
// so the MinGW-built didgpu.dll never links the CUDA runtime.
// ===========================================================================

// Truncated rank-r SVD. M_rm is HOST (m x n) row-major; out_L_rm (m x r)
// and out_F_rm (r x n) are HOST buffers the caller has preallocated, filled
// such that out_L_rm * out_F_rm is a rank-r approximation of M.
extern "C" int didgpu_cuda_fect_svd_truncated(
    const double* M_rm,
    int m, int n, int r,
    double* out_L_rm,
    double* out_F_rm) {
  if (r <= 0 || r > std::min(m, n)) return -1;

  cudaError_t e;
  double *d_M = nullptr, *d_L = nullptr, *d_F = nullptr;
  e = cudaMalloc((void**)&d_M, sizeof(double) * m * n);
  if (e != cudaSuccess) { return -3; }
  e = cudaMalloc((void**)&d_L, sizeof(double) * m * r);
  if (e != cudaSuccess) { cudaFree(d_M); return -3; }
  e = cudaMalloc((void**)&d_F, sizeof(double) * r * n);
  if (e != cudaSuccess) { cudaFree(d_M); cudaFree(d_L); return -3; }

  e = cudaMemcpy(d_M, M_rm, sizeof(double) * m * n, cudaMemcpyHostToDevice);
  if (e != cudaSuccess) { cudaFree(d_M); cudaFree(d_L); cudaFree(d_F); return -3; }

  int rc = fect_svd_truncated_dev(d_M, m, n, r, d_L, d_F);
  if (rc != 0) { cudaFree(d_M); cudaFree(d_L); cudaFree(d_F); return rc; }

  e = cudaMemcpy(out_L_rm, d_L, sizeof(double) * m * r, cudaMemcpyDeviceToHost);
  if (e == cudaSuccess)
    e = cudaMemcpy(out_F_rm, d_F, sizeof(double) * r * n, cudaMemcpyDeviceToHost);
  cudaFree(d_M); cudaFree(d_L); cudaFree(d_F);
  return (e == cudaSuccess) ? 0 : -3;
}

// Soft-thresholded SVD reconstruction. Y_complete_rm is HOST (m x n)
// row-major; out_Y_hat_rm (m x n) is a HOST buffer the caller preallocated.
extern "C" int didgpu_cuda_fect_svd_softthreshold(
    const double* Y_complete_rm,
    int m, int n,
    double lambda,
    double* out_Y_hat_rm,
    int* out_n_nonzero) {
  cudaError_t e;
  double *d_Y = nullptr, *d_Yhat = nullptr;
  e = cudaMalloc((void**)&d_Y, sizeof(double) * m * n);
  if (e != cudaSuccess) { return -3; }
  e = cudaMalloc((void**)&d_Yhat, sizeof(double) * m * n);
  if (e != cudaSuccess) { cudaFree(d_Y); return -3; }

  e = cudaMemcpy(d_Y, Y_complete_rm, sizeof(double) * m * n,
                 cudaMemcpyHostToDevice);
  if (e != cudaSuccess) { cudaFree(d_Y); cudaFree(d_Yhat); return -3; }

  int rc = fect_svd_softthreshold_dev(d_Y, m, n, lambda, d_Yhat, out_n_nonzero);
  if (rc != 0) { cudaFree(d_Y); cudaFree(d_Yhat); return rc; }

  e = cudaMemcpy(out_Y_hat_rm, d_Yhat, sizeof(double) * m * n,
                 cudaMemcpyDeviceToHost);
  cudaFree(d_Y); cudaFree(d_Yhat);
  return (e == cudaSuccess) ? 0 : -3;
}

#endif  // DIDGPU_LITE (full build with cuSOLVER/cuBLAS)
#endif  // HAS_CUDA
