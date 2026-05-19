// ============================================================================
// CUDA kernel for fect_fe: iterative two-way demeaning of the controls-only
// outcome matrix, until alpha + xi converge.
//
// Algorithm (matches R/fect_fe.R::.fect_fe_fit):
//   Inputs:
//     Y[n_units * n_periods]   - row-major outcome matrix (NaN = missing)
//     M[n_units * n_periods]   - row-major treatment mask (1 = drop from fit)
//   Outputs:
//     alpha[n_units]            - per-unit FE
//     xi[n_periods]             - per-time FE
//
// Per iteration:
//   1. R[i, t] = Y[i, t] - xi[t]   (control cells only; NaN for treated/missing)
//   2. alpha[i] = mean(R[i, :], skipping NaN)
//   3. R[i, t] = Y[i, t] - alpha[i] - xi[t]   (same masking)
//   4. xi[t] = colmean of R, skipping NaN
//
// On the GPU, each kernel launches one block per row (for alpha) or
// per column (for xi), with threads cooperating on the reduction.
//
// Convergence check is computed on host after each iter (small cost
// compared to the per-iter passes through the n_units * n_periods grid).
//
// STATUS: scaffold. Compiles when nvcc is on PATH and CUDA headers are
// present (we've already shown that build works for the existing
// cuda_didkernel.cu). Untested on this dev machine (no admin to install
// the full CUDA Toolkit), but written to mirror the R reference impl
// line-by-line so the first run on a CUDA-capable box should be a
// straightforward smoke test.
// ============================================================================

#ifdef HAS_CUDA
#include <cuda_runtime.h>
#include <cmath>

// Block size for per-row / per-col reductions. 256 is a safe default
// for any reasonable panel; rows / cols longer than this loop within
// the block.
#define DIDGPU_FECT_BLOCK 256

// ---------------------------------------------------------------------------
// Row-mean reduction: alpha[i] = mean(Y[i, t] - xi[t]) over t where
// M[i, t] == 0 and Y[i, t] is finite. One block per row.
// ---------------------------------------------------------------------------
__global__ void k_fect_row_mean(
    const double* __restrict__ Y,
    const int*    __restrict__ M,
    const double* __restrict__ xi,
    double*       __restrict__ alpha,
    int n_units, int n_periods) {

  __shared__ double sh_sum[DIDGPU_FECT_BLOCK];
  __shared__ int    sh_cnt[DIDGPU_FECT_BLOCK];

  const int i = blockIdx.x;
  if (i >= n_units) return;

  double local_sum = 0.0;
  int    local_cnt = 0;
  const double* yrow = Y + i * n_periods;
  const int*    mrow = M + i * n_periods;

  for (int t = threadIdx.x; t < n_periods; t += blockDim.x) {
    if (mrow[t] == 0) {
      const double y = yrow[t];
      if (!isnan(y)) {
        local_sum += y - xi[t];
        local_cnt += 1;
      }
    }
  }
  sh_sum[threadIdx.x] = local_sum;
  sh_cnt[threadIdx.x] = local_cnt;
  __syncthreads();

  // Reduce within block.
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      sh_sum[threadIdx.x] += sh_sum[threadIdx.x + s];
      sh_cnt[threadIdx.x] += sh_cnt[threadIdx.x + s];
    }
    __syncthreads();
  }

  if (threadIdx.x == 0) {
    alpha[i] = (sh_cnt[0] > 0) ? (sh_sum[0] / static_cast<double>(sh_cnt[0]))
                                : 0.0;
  }
}


// ---------------------------------------------------------------------------
// Col-mean reduction: xi[t] = mean(Y[i, t] - alpha[i]) over i where
// M[i, t] == 0 and Y[i, t] is finite. One block per column.
// ---------------------------------------------------------------------------
__global__ void k_fect_col_mean(
    const double* __restrict__ Y,
    const int*    __restrict__ M,
    const double* __restrict__ alpha,
    double*       __restrict__ xi,
    int n_units, int n_periods) {

  __shared__ double sh_sum[DIDGPU_FECT_BLOCK];
  __shared__ int    sh_cnt[DIDGPU_FECT_BLOCK];

  const int t = blockIdx.x;
  if (t >= n_periods) return;

  double local_sum = 0.0;
  int    local_cnt = 0;

  for (int i = threadIdx.x; i < n_units; i += blockDim.x) {
    const int idx = i * n_periods + t;
    if (M[idx] == 0) {
      const double y = Y[idx];
      if (!isnan(y)) {
        local_sum += y - alpha[i];
        local_cnt += 1;
      }
    }
  }
  sh_sum[threadIdx.x] = local_sum;
  sh_cnt[threadIdx.x] = local_cnt;
  __syncthreads();

  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      sh_sum[threadIdx.x] += sh_sum[threadIdx.x + s];
      sh_cnt[threadIdx.x] += sh_cnt[threadIdx.x + s];
    }
    __syncthreads();
  }

  if (threadIdx.x == 0) {
    xi[t] = (sh_cnt[0] > 0) ? (sh_sum[0] / static_cast<double>(sh_cnt[0]))
                              : 0.0;
  }
}


// ---------------------------------------------------------------------------
// Host launcher: run the demeaning iteration to convergence on device-
// resident Y, M. Returns alpha, xi on device. Caller is responsible for
// H2D transfer of inputs and D2H of outputs.
//
// Convergence check (max |alpha_new - alpha_old| and max |xi_new -
// xi_old|) is done on host between launches; this keeps the kernel
// itself stateless and easy to compose.
//
// Returns: 0 on success, non-zero on CUDA error. On success, *out_iter
// holds the iteration count and *out_delta holds the final max-abs
// change.
// ---------------------------------------------------------------------------
extern "C" int didgpu_cuda_fect_fe(
    const double* d_Y,
    const int*    d_M,
    double*       d_alpha,
    double*       d_xi,
    int n_units, int n_periods,
    double tol, int max_iter,
    int* out_iter, double* out_delta) {

  // Zero-initialise alpha and xi.
  cudaError_t e;
  e = cudaMemset(d_alpha, 0, n_units   * sizeof(double));
  if (e != cudaSuccess) return static_cast<int>(e);
  e = cudaMemset(d_xi,    0, n_periods * sizeof(double));
  if (e != cudaSuccess) return static_cast<int>(e);

  // Scratch host buffers for convergence check.
  double* h_alpha     = new double[n_units];
  double* h_xi        = new double[n_periods];
  double* h_alpha_old = new double[n_units];
  double* h_xi_old    = new double[n_periods];
  for (int i = 0; i < n_units;   ++i) { h_alpha_old[i] = 0.0; }
  for (int t = 0; t < n_periods; ++t) { h_xi_old[t]    = 0.0; }

  int iter = 0;
  double delta = 1e30;
  for (iter = 1; iter <= max_iter; ++iter) {
    // Update alpha.
    k_fect_row_mean<<<n_units, DIDGPU_FECT_BLOCK>>>(
      d_Y, d_M, d_xi, d_alpha, n_units, n_periods);
    e = cudaGetLastError();
    if (e != cudaSuccess) { delete[] h_alpha; delete[] h_xi;
                            delete[] h_alpha_old; delete[] h_xi_old;
                            return static_cast<int>(e); }
    // Update xi.
    k_fect_col_mean<<<n_periods, DIDGPU_FECT_BLOCK>>>(
      d_Y, d_M, d_alpha, d_xi, n_units, n_periods);
    e = cudaGetLastError();
    if (e != cudaSuccess) { delete[] h_alpha; delete[] h_xi;
                            delete[] h_alpha_old; delete[] h_xi_old;
                            return static_cast<int>(e); }
    // Pull alpha, xi for convergence check.
    cudaMemcpy(h_alpha, d_alpha, n_units   * sizeof(double),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(h_xi,    d_xi,    n_periods * sizeof(double),
               cudaMemcpyDeviceToHost);
    double d_a = 0.0, d_x = 0.0;
    for (int i = 0; i < n_units; ++i) {
      const double dd = std::fabs(h_alpha[i] - h_alpha_old[i]);
      if (dd > d_a) d_a = dd;
      h_alpha_old[i] = h_alpha[i];
    }
    for (int t = 0; t < n_periods; ++t) {
      const double dd = std::fabs(h_xi[t] - h_xi_old[t]);
      if (dd > d_x) d_x = dd;
      h_xi_old[t] = h_xi[t];
    }
    delta = (d_a > d_x) ? d_a : d_x;
    if (delta < tol) break;
  }

  delete[] h_alpha; delete[] h_xi;
  delete[] h_alpha_old; delete[] h_xi_old;
  if (out_iter)  *out_iter  = iter;
  if (out_delta) *out_delta = delta;
  return 0;
}

#endif  // HAS_CUDA
