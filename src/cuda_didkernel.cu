// ============================================================================
// CUDA kernels for the binary, no-controls DiD U-statistic.
//
// One launch processes a SINGLE bootstrap iteration. To run a 100-rep
// bootstrap in parallel, launch 100 streams (see launch_did_iter_batched
// for the scaffold).
//
// Algorithm mirrors R/core_r.R::.core_one_event_time. The data layout
// is column-major arrays of length n_rows (= G * T after balancing).
// Per-group quantities (F_g, S_g, T_g, L_g, d_sq) are stored once per
// group and looked up via row_to_group[row_idx].
//
// Inputs (all on device):
//   outcome     [n_rows]   double — Y
//   N_gt        [n_rows]   double — 0 or 1 in the binary case
//   row_to_g    [n_rows]   int    — group index for each row
//   row_to_t    [n_rows]   int    — time index for each row (1-based)
//   row_to_dsq  [n_rows]   int    — d_sq (baseline treatment) for each row
//   F_g         [n_groups] int    — first-switch period per group
//   S_g         [n_groups] int    — 1 (in), 0 (out), -1 (never)
//   T_g         [n_groups] int    — last usable time per group
//   L_g         [n_groups] int    — post-switch horizon per group
// Outputs:
//   did         [effects]  double — per-event-time DID estimates
//
// Compile only when nvcc is available (Makevars sets HAS_CUDA=1).
// ============================================================================

#include <cuda_runtime.h>

#define CUDA_CHECK(call) do { cudaError_t e = call; if (e != cudaSuccess) return (int)e; } while (0)

// --------------------------------------------------------------------------
// Kernel 1: per-row "lag-by-k" gather.
//   diff_y_k[r] = outcome[r] - outcome[lag_idx], where lag_idx is the row
//   in this group at time (row_to_t[r] - k). If no such row exists (early
//   in the panel), diff_y_k[r] = NAN.
//
// Assumes the panel is sorted by (group, time) so that row r and row r-k
// in the same group are k rows apart in storage. This is what .prep_panel
// guarantees via setkeyv(d, c("group_XX", "time_XX")).
// --------------------------------------------------------------------------
__global__ void k_lag_diff(
    const double* __restrict__ outcome,
    const int*    __restrict__ row_to_g,
    int n_rows, int k,
    double* __restrict__ diff_y_k)
{
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r >= n_rows) return;
  int lag_r = r - k;
  if (lag_r < 0 || row_to_g[lag_r] != row_to_g[r]) {
    diff_y_k[r] = NAN;
  } else {
    diff_y_k[r] = outcome[r] - outcome[lag_r];
  }
}

// --------------------------------------------------------------------------
// Kernel 2: build the per-row never_change_k mask and the per-row
// "candidate switcher cell" mask. Final dist_k mask depends on whether
// the (time, d_sq) cohort has any controls, so we defer the final
// gating until after the cohort-sum kernel.
// --------------------------------------------------------------------------
__global__ void k_make_masks(
    const int* __restrict__ row_to_g,
    const int* __restrict__ row_to_t,
    const double* __restrict__ N_gt,
    const double* __restrict__ diff_y_k,
    const int* __restrict__ F_g,
    const int* __restrict__ S_g,
    const int* __restrict__ L_g,
    int n_rows, int k, int direction,
    int* __restrict__ never_change_k,
    int* __restrict__ candidate_dist_k)
{
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r >= n_rows) return;
  int g  = row_to_g[r];
  int t  = row_to_t[r];
  int fg = F_g[g];
  bool diff_ok = !isnan(diff_y_k[r]);
  bool gt_ok   = N_gt[r] > 0.0;

  never_change_k[r] = (t < fg && gt_ok && diff_ok) ? 1 : 0;

  int sg = S_g[g];
  int lg = L_g[g];
  bool is_cand = (t == fg + k - 1)
              && (k <= lg)
              && (sg == direction)
              && diff_ok
              && gt_ok;
  candidate_dist_k[r] = is_cand ? 1 : 0;
}

// --------------------------------------------------------------------------
// Kernel 3: per-(time, d_sq) cohort sums. We compute both:
//   N_t_control[(t, d_sq)] = sum_r N_gt[r] * never_change_k[r] in that cohort
//   N_t_switch_cand[(t, d_sq)] = sum_r N_gt[r] * candidate_dist_k[r] in cohort
//
// For up-to-O(few) d_sq levels (typically 1 or 2 in binary) and up to
// O(T) times, we flatten the cohort key to a single index and use
// atomicAdd into a scratch array of size n_cohorts. n_cohorts is small
// enough (T * |d_sq|) that this is fast.
//
// Caller computes cohort_key[r] = (row_to_t[r] - 1) * n_dsq_levels +
// dsq_index_of(row_to_dsq[r]) on the host before launch.
// --------------------------------------------------------------------------
__global__ void k_cohort_sums(
    const double* __restrict__ N_gt,
    const int*    __restrict__ never_change_k,
    const int*    __restrict__ candidate_dist_k,
    const int*    __restrict__ cohort_key,
    int n_rows,
    double* __restrict__ N_t_control,
    double* __restrict__ N_t_switch_cand)
{
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r >= n_rows) return;
  int c = cohort_key[r];
  double n  = N_gt[r];
  if (never_change_k[r]) atomicAdd(&N_t_control[c],     n);
  if (candidate_dist_k[r]) atomicAdd(&N_t_switch_cand[c], n);
}

// --------------------------------------------------------------------------
// Kernel 4: finalise dist_k (gate by N_t_control > 0) and compute the
// per-row kernel value:
//   kernel[r] = (G / N_inc) * 1[t in (k+1)..T_g] * N_gt
//             * (dist_k - (N_t_switch / N_t_control) * never_change_k)
//             * diff_y_k
// where dist_k is the gated version and N_t_switch is the *gated* per-
// cohort switcher mass (sum of dist_k after gating, which equals
// N_t_switch_cand when N_t_control > 0 for that cohort, else 0).
// --------------------------------------------------------------------------
__global__ void k_finalize_dist_and_kernel(
    const double* __restrict__ N_gt,
    const int*    __restrict__ candidate_dist_k,
    const int*    __restrict__ never_change_k,
    const double* __restrict__ diff_y_k,
    const int*    __restrict__ row_to_g,
    const int*    __restrict__ row_to_t,
    const int*    __restrict__ cohort_key,
    const int*    __restrict__ T_g,
    const double* __restrict__ N_t_control,
    const double* __restrict__ N_t_switch_cand,
    int n_rows, int k, double G_over_Ninc,
    int*    __restrict__ dist_k_final,
    double* __restrict__ kernel_val)
{
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r >= n_rows) return;
  int c = cohort_key[r];
  double n_ctrl = N_t_control[c];
  int dist = (candidate_dist_k[r] && n_ctrl > 0.0) ? 1 : 0;
  dist_k_final[r] = dist;

  int g = row_to_g[r];
  int t = row_to_t[r];
  bool in_window = (t >= k + 1) && (t <= T_g[g]);

  double n_switch = (n_ctrl > 0.0) ? N_t_switch_cand[c] : 0.0;
  double ratio = (n_ctrl > 0.0) ? (n_switch / n_ctrl) : 0.0;

  double inner = (double)dist - ratio * (double)never_change_k[r];
  double v = G_over_Ninc * (in_window ? 1.0 : 0.0) * N_gt[r] * inner * diff_y_k[r];
  if (isnan(v)) v = 0.0;
  kernel_val[r] = v;
}

// --------------------------------------------------------------------------
// Kernel 5: per-group reduction of kernel_val -> U_g, then sum U_g -> DID.
// For simplicity we do a single-pass atomicAdd into a per-group accumulator
// (n_groups is small) and then a single thread reduces to scalar DID.
// --------------------------------------------------------------------------
__global__ void k_pergroup_sum(
    const double* __restrict__ kernel_val,
    const int*    __restrict__ row_to_g,
    int n_rows,
    double* __restrict__ U_g)
{
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r >= n_rows) return;
  atomicAdd(&U_g[row_to_g[r]], kernel_val[r]);
}

__global__ void k_final_did(
    const double* __restrict__ U_g,
    int n_groups,
    double* __restrict__ did_out)
{
  // One block, n_threads = next-pow-2(n_groups) or 256 with stride.
  __shared__ double s[256];
  int tid = threadIdx.x;
  double acc = 0.0;
  for (int g = tid; g < n_groups; g += blockDim.x) acc += U_g[g];
  s[tid] = acc;
  __syncthreads();
  for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
    if (tid < stride) s[tid] += s[tid + stride];
    __syncthreads();
  }
  if (tid == 0) *did_out = s[0] / (double)n_groups;
}

// --------------------------------------------------------------------------
// Host-side launcher: run all kernels for a single (k, direction) pair on
// one bootstrap iteration. Returns the DID value via did_out (host-resident).
//
// Caller is responsible for:
//   - allocating + populating all device arrays (outcome, N_gt, row_to_*,
//     F_g, S_g, T_g, L_g, cohort_key)
//   - providing a workspace big enough for diff_y_k, never_change_k,
//     candidate_dist_k, dist_k_final, kernel_val, U_g, N_t_control,
//     N_t_switch_cand (sizes documented in the param list)
//   - choosing direction (1 = in, 0 = out)
//   - knowing N_inc out-of-band so it can pass G_over_Ninc
//     (= G / N_inc, scalar). N_inc itself is derivable from the
//     candidate_dist_k * N_gt sum after gating; we leave it to the
//     caller for now so the kernel chain stays linear.
// --------------------------------------------------------------------------
extern "C" int didgpu_cuda_run_one_event_time(
    const double* outcome,           // [n_rows]
    const double* N_gt,              // [n_rows]
    const int*    row_to_g,          // [n_rows]
    const int*    row_to_t,          // [n_rows]
    const int*    cohort_key,        // [n_rows]
    const int*    F_g,               // [n_groups]
    const int*    S_g,               // [n_groups]
    const int*    T_g,               // [n_groups]
    const int*    L_g,               // [n_groups]
    int n_rows, int n_groups, int n_cohorts,
    int k, int direction, double G_over_Ninc,
    // Workspace (caller-allocated):
    double* diff_y_k,                // [n_rows]
    int*    never_change_k,          // [n_rows]
    int*    candidate_dist_k,        // [n_rows]
    int*    dist_k_final,            // [n_rows]
    double* kernel_val,              // [n_rows]
    double* U_g,                     // [n_groups]
    double* N_t_control,             // [n_cohorts]
    double* N_t_switch_cand,         // [n_cohorts]
    double* did_out_device           // [1]
){
  const int block = 256;
  const int grid_rows = (n_rows + block - 1) / block;

  // Zero accumulators.
  CUDA_CHECK(cudaMemsetAsync(N_t_control,     0, n_cohorts * sizeof(double), 0));
  CUDA_CHECK(cudaMemsetAsync(N_t_switch_cand, 0, n_cohorts * sizeof(double), 0));
  CUDA_CHECK(cudaMemsetAsync(U_g,             0, n_groups * sizeof(double), 0));

  k_lag_diff<<<grid_rows, block>>>(outcome, row_to_g, n_rows, k, diff_y_k);

  k_make_masks<<<grid_rows, block>>>(
      row_to_g, row_to_t, N_gt, diff_y_k,
      F_g, S_g, L_g, n_rows, k, direction,
      never_change_k, candidate_dist_k);

  k_cohort_sums<<<grid_rows, block>>>(
      N_gt, never_change_k, candidate_dist_k, cohort_key,
      n_rows, N_t_control, N_t_switch_cand);

  k_finalize_dist_and_kernel<<<grid_rows, block>>>(
      N_gt, candidate_dist_k, never_change_k, diff_y_k,
      row_to_g, row_to_t, cohort_key, T_g,
      N_t_control, N_t_switch_cand,
      n_rows, k, G_over_Ninc,
      dist_k_final, kernel_val);

  k_pergroup_sum<<<grid_rows, block>>>(kernel_val, row_to_g, n_rows, U_g);

  k_final_did<<<1, 256>>>(U_g, n_groups, did_out_device);

  CUDA_CHECK(cudaDeviceSynchronize());
  return 0;
}
