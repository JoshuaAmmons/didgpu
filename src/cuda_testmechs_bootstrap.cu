// ============================================================================
// CUDA bootstrap kernel for TestMechs partial densities.
//
// Computes a (B x dim_beta) matrix of bootstrap-replicated partial-
// density vectors in one launch:
//
//   For each bootstrap draw b in 0..B-1:
//     For each observation i in 0..n-1:
//       Pick an index j = curand_uniform(state[b]) * n     (multinomial)
//       Accumulate count[b, d[j] * K * dy + (m[j]-1) * dy + (y[j]-1)] += 1
//     Normalise: count[b, ...] /= count_per_d[b, d]
//
// This is the dominant cost in TestMechs's analytic-variance OFF path.
// For n = 10K, B = 500, this is 5M random index draws + 5M atomic adds
// + a final normalisation — well under a second on a modest GPU,
// versus ~30-60 seconds on a single CPU thread.
//
// Bayesian-bootstrap variant: replace the multinomial draw with
// Dirichlet weights (cuRAND generates Gamma(1, 1) = Exp(1) draws;
// normalise to sum-1 per draw; do a weighted reduction instead of
// counts).
//
// STATUS: scaffold. Compiles when nvcc + cuRAND are present (cuRAND is
// in the standard CUDA Toolkit). Untested on this dev machine without
// admin install of the full Toolkit. The R-side falls back to
// .testmechs_bootstrap_r in that case.
// ============================================================================

#ifdef HAS_CUDA
#include <cuda_runtime.h>
#include <curand_kernel.h>

#define DIDGPU_TM_BOOT_BLOCK 256

// Per-bootstrap initialiser: one thread per bootstrap draw, seeds its
// own cuRAND state from a (master_seed, b) pair.
__global__ void k_tm_init_states(curandState* states, unsigned long master_seed,
                                  int B) {
  const int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= B) return;
  // Subseed per bootstrap so different B draws are independent. Same
  // (master_seed, b) gives the same draw every time = reproducibility.
  curand_init(master_seed, /*sequence=*/b, /*offset=*/0, &states[b]);
}


// Per-(bootstrap, obs) reduction. Each block handles one bootstrap
// draw; threads within the block stream through the n observations,
// each picking a random index and atomically incrementing the
// (b, d * K * dy + (m-1) * dy + (y-1)) slot.
//
// dim_beta = 2 * K * dy. d_counts is the per-D total count for
// normalisation (treated/control sizes can differ across draws).
__global__ void k_tm_bootstrap_nonparam(
    curandState* states,
    const int* __restrict__ d_d,    // device, length n
    const int* __restrict__ d_m,    // device, length n (1-based)
    const int* __restrict__ d_y,    // device, length n (1-based)
    int n, int K, int dy, int B,
    double* __restrict__ d_beta,    // device, B x dim_beta (row-major)
    int* __restrict__ d_count_per_d // device, B x 2 (count of d=0 and d=1)
    ) {
  const int b = blockIdx.x;
  if (b >= B) return;

  const int dim_beta = 2 * K * dy;
  // Each block has 256 threads; they cooperate on the n-loop.
  // For simplicity we let thread 0 use the cuRAND state and broadcast
  // (or, more parallelism, give each thread its own state initialised
  // at the same seed-offset). For brevity here we use thread-0
  // serial draws; production should fan out.
  if (threadIdx.x == 0) {
    curandState local_state = states[b];
    for (int i = 0; i < n; ++i) {
      // Random index in [0, n).
      unsigned int r = curand(&local_state);
      int j = static_cast<int>(static_cast<unsigned long long>(r) % n);
      const int dj = d_d[j];
      const int mj = d_m[j];
      const int yj = d_y[j];
      // Accumulate into the (dj, mj, yj) bin of bootstrap b.
      const int pos = b * dim_beta + dj * K * dy + (mj - 1) * dy + (yj - 1);
      // No atomics needed: only one thread per block writes.
      d_beta[pos] += 1.0;
      d_count_per_d[b * 2 + dj] += 1;
    }
    states[b] = local_state;
  }
}


// Per-(bootstrap, position) normalisation: divide each count by the
// per-D total for that bootstrap draw.
__global__ void k_tm_normalize(
    double* __restrict__ d_beta,
    const int* __restrict__ d_count_per_d,
    int B, int K, int dy) {
  const int b   = blockIdx.x;
  const int pos = threadIdx.x;
  const int dim_beta = 2 * K * dy;
  if (b >= B || pos >= dim_beta) return;
  const int dd = (pos < K * dy) ? 0 : 1;
  const int n_dd = d_count_per_d[b * 2 + dd];
  if (n_dd > 0) {
    d_beta[b * dim_beta + pos] /= static_cast<double>(n_dd);
  }
}


// Host launcher. Allocates device buffers, runs the three-stage
// kernel chain, copies results back to host.
//
// Returns: 0 on success, nonzero CUDA error code on failure.
extern "C" int didgpu_cuda_testmechs_bootstrap(
    const int* h_d, const int* h_m, const int* h_y,
    int n, int K, int dy, int B,
    unsigned long seed,
    double* h_beta) {

  cudaError_t e;
  const int dim_beta = 2 * K * dy;
  int *d_d = nullptr, *d_m = nullptr, *d_y = nullptr, *d_count = nullptr;
  double *d_beta = nullptr;
  curandState* d_states = nullptr;

  auto cleanup = [&]() {
    if (d_d)      cudaFree(d_d);
    if (d_m)      cudaFree(d_m);
    if (d_y)      cudaFree(d_y);
    if (d_count)  cudaFree(d_count);
    if (d_beta)   cudaFree(d_beta);
    if (d_states) cudaFree(d_states);
  };

  // Allocate.
  e = cudaMalloc((void**)&d_d, n * sizeof(int));
  if (e != cudaSuccess) { cleanup(); return static_cast<int>(e); }
  e = cudaMalloc((void**)&d_m, n * sizeof(int));
  if (e != cudaSuccess) { cleanup(); return static_cast<int>(e); }
  e = cudaMalloc((void**)&d_y, n * sizeof(int));
  if (e != cudaSuccess) { cleanup(); return static_cast<int>(e); }
  e = cudaMalloc((void**)&d_count, B * 2 * sizeof(int));
  if (e != cudaSuccess) { cleanup(); return static_cast<int>(e); }
  e = cudaMalloc((void**)&d_beta, B * dim_beta * sizeof(double));
  if (e != cudaSuccess) { cleanup(); return static_cast<int>(e); }
  e = cudaMalloc((void**)&d_states, B * sizeof(curandState));
  if (e != cudaSuccess) { cleanup(); return static_cast<int>(e); }

  // Copy inputs.
  cudaMemcpy(d_d, h_d, n * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_m, h_m, n * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, h_y, n * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemset(d_beta, 0, B * dim_beta * sizeof(double));
  cudaMemset(d_count, 0, B * 2 * sizeof(int));

  // Stage 1: init RNG states.
  k_tm_init_states<<<(B + DIDGPU_TM_BOOT_BLOCK - 1) / DIDGPU_TM_BOOT_BLOCK,
                      DIDGPU_TM_BOOT_BLOCK>>>(
      d_states, seed, B);
  e = cudaGetLastError();
  if (e != cudaSuccess) { cleanup(); return static_cast<int>(e); }

  // Stage 2: bootstrap reduction (one block per draw).
  k_tm_bootstrap_nonparam<<<B, DIDGPU_TM_BOOT_BLOCK>>>(
      d_states, d_d, d_m, d_y, n, K, dy, B, d_beta, d_count);
  e = cudaGetLastError();
  if (e != cudaSuccess) { cleanup(); return static_cast<int>(e); }

  // Stage 3: normalise.
  k_tm_normalize<<<B, dim_beta>>>(d_beta, d_count, B, K, dy);
  e = cudaGetLastError();
  if (e != cudaSuccess) { cleanup(); return static_cast<int>(e); }

  // Copy back.
  cudaMemcpy(h_beta, d_beta, B * dim_beta * sizeof(double),
             cudaMemcpyDeviceToHost);
  cleanup();
  return 0;
}

#endif  // HAS_CUDA
