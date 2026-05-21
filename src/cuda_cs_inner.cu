// ============================================================================
// CUDA kernel for the Callaway-Sant'Anna per-(g, t) inner regression.
//
// Implements the OR (outcome-regression) estimator end-to-end. IPW
// (est_method = 1) and DR (est_method = 2) still return -3 ("not yet
// implemented"); the R fallback handles them. Phase 2 follow-ups
// (still under task #84) will add the propensity-score logistic
// kernel and the DR augmentation step.
//
// OR algorithm (matches R/cs_methods.R::.cs_inner_or):
//   For each cell c with rows [X_offsets[c], X_offsets[c+1]):
//     1. Split rows into treated (W = 1) and control (W = 0).
//     2. If p == 1 (intercept only): ATT = mean(Y_t) - mean(Y_c).
//     3. Else fit OLS on controls: beta = (X_c' X_c)^{-1} X_c' Y_c
//        via in-thread Cholesky decomposition + forward/backward
//        substitution.
//     4. ATT = mean(Y_t - X_t @ beta).
//   If Cholesky fails (singular Gram matrix) or n_t == 0 or n_c == 0,
//   the cell's ATT is NaN.
//
// Design choice — one thread per cell:
//   CS cells are typically small (n_c ~ 10-100, p < 16). The work per
//   cell is dominated by the X_c' X_c reduction which is O(n_c * p^2)
//   = a few thousand FLOPs. With n_cells ~ 50-200, total work is well
//   under 1 ms even with one thread per cell. Multi-thread cooperation
//   inside a cell would add atomic-add overhead and shared-memory
//   pressure for marginal gain.
//   The hard upper bound is p <= 16 (the per-thread XtX[16*16]
//   buffer); the launcher rejects larger p.
//
// Influence functions (out_influence) are not yet populated in this
// kernel — that's task #85. The Rcpp wrapper passes
// want_influence = TRUE today and gets back NULL (kernel returns -3
// for that path); the R fallback path remains the only source of IF
// data for now.
// ============================================================================

#ifdef HAS_CUDA
#include <cuda_runtime.h>


// Hard upper bound on per-cell covariate count. Each thread keeps
// XtX (p*p doubles) and a few smaller buffers in local memory; the
// 16 limit is generous (CS typically uses p <= 5) and keeps the
// per-thread footprint under 4 KB.
#define DIDGPU_CS_MAX_P 16


__global__ void k_cs_inner_or(
    const double* __restrict__ X_concat,   // (n_total, p) row-major
    const double* __restrict__ Y_concat,   // (n_total)
    const double* __restrict__ W_concat,   // (n_total) -- D in {0, 1}
    const int*    __restrict__ X_offsets,  // (n_cells + 1)
    int n_cells, int p,
    double* __restrict__ out_att,          // (n_cells)
    double* __restrict__ out_IF_per_row) { // (n_total) — per-row IF, or NULL
  const int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= n_cells) return;

  const int row_start = X_offsets[c];
  const int row_end   = X_offsets[c + 1];

  // Per-thread scratch. With DIDGPU_CS_MAX_P = 16 these are 16*16*8 =
  // 2 KB plus a couple of vectors. The compiler may demote XtX to
  // local (global) memory if p is fully dynamic; access patterns are
  // cache-friendly so the latency hit is small.
  double XtX[DIDGPU_CS_MAX_P * DIDGPU_CS_MAX_P];
  double XtY[DIDGPU_CS_MAX_P];
  double L  [DIDGPU_CS_MAX_P * DIDGPU_CS_MAX_P];
  double z  [DIDGPU_CS_MAX_P];
  double beta[DIDGPU_CS_MAX_P];

  // Initialize.
  #pragma unroll
  for (int i = 0; i < DIDGPU_CS_MAX_P * DIDGPU_CS_MAX_P; ++i) XtX[i] = 0.0;
  #pragma unroll
  for (int i = 0; i < DIDGPU_CS_MAX_P; ++i) { XtY[i] = 0.0; beta[i] = 0.0; }

  int n_t = 0, n_c = 0;
  double sum_Yt = 0.0;

  // --- Pass 1: accumulate Gram and cross-products on controls;
  //              count treated/controls and sum_Y_t. ---
  for (int r = row_start; r < row_end; ++r) {
    const double w = W_concat[r];
    if (w > 0.5) {
      ++n_t;
      sum_Yt += Y_concat[r];
    } else {
      ++n_c;
      const double y = Y_concat[r];
      const int base = r * p;
      for (int a = 0; a < p; ++a) {
        const double xa = X_concat[base + a];
        XtY[a] += xa * y;
        for (int b = 0; b < p; ++b) {
          XtX[a * p + b] += xa * X_concat[base + b];
        }
      }
    }
  }

  if (n_t == 0 || n_c == 0) {
    out_att[c] = nan("");
    if (out_IF_per_row) {
      for (int r = row_start; r < row_end; ++r) out_IF_per_row[r] = 0.0;
    }
    return;
  }

  if (p == 1) {
    // No covariates: simple difference of means. (X is the intercept
    // column = all 1.0, so XtX[0] = n_c and XtY[0] = sum(Y_c).)
    const double mean_Yc = XtY[0] / static_cast<double>(n_c);
    const double mean_Yt = sum_Yt / static_cast<double>(n_t);
    const double att = mean_Yt - mean_Yc;
    out_att[c] = att;
    if (out_IF_per_row) {
      // IF[treated row r] = (Y[r] - mean_Yc) - att = Y[r] - mean_Yt
      // IF[control row r] = 0  (OR puts no per-control IF mass at the
      // ATT level — control units enter only through the projection)
      for (int r = row_start; r < row_end; ++r) {
        out_IF_per_row[r] = (W_concat[r] > 0.5) ? (Y_concat[r] - mean_Yt) : 0.0;
      }
    }
    return;
  }

  // --- Pass 2: Cholesky decompose XtX = L * L^T, then solve. ---
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = XtX[i * p + j];
      for (int k = 0; k < j; ++k) s -= L[i * p + k] * L[j * p + k];
      if (i == j) {
        if (s <= 0.0) {  // Rank-deficient or numerically singular.
          out_att[c] = nan("");
          if (out_IF_per_row) {
            for (int r = row_start; r < row_end; ++r) out_IF_per_row[r] = 0.0;
          }
          return;
        }
        L[i * p + i] = sqrt(s);
      } else {
        L[i * p + j] = s / L[j * p + j];
      }
    }
  }
  // Forward sub: L * z = XtY.
  for (int i = 0; i < p; ++i) {
    double s = XtY[i];
    for (int j = 0; j < i; ++j) s -= L[i * p + j] * z[j];
    z[i] = s / L[i * p + i];
  }
  // Backward sub: L^T * beta = z.
  for (int i = p - 1; i >= 0; --i) {
    double s = z[i];
    for (int j = i + 1; j < p; ++j) s -= L[j * p + i] * beta[j];
    beta[i] = s / L[i * p + i];
  }

  // --- Pass 3: ATT = mean(Y_t - X_t @ beta). ---
  double sum_fittedt = 0.0;
  for (int r = row_start; r < row_end; ++r) {
    if (W_concat[r] > 0.5) {
      double fitted = 0.0;
      const int base = r * p;
      for (int a = 0; a < p; ++a) fitted += X_concat[base + a] * beta[a];
      sum_fittedt += fitted;
    }
  }
  const double att = (sum_Yt - sum_fittedt) / static_cast<double>(n_t);
  out_att[c] = att;

  // --- Pass 4 (optional): per-row influence function. ---
  // Matches R/cs_methods.R::.cs_inner_or:
  //   IF[treated row r] = (Y[r] - fitted[r]) - att
  //   IF[control row r] = 0
  if (out_IF_per_row) {
    for (int r = row_start; r < row_end; ++r) {
      if (W_concat[r] > 0.5) {
        double fitted = 0.0;
        const int base = r * p;
        for (int a = 0; a < p; ++a) fitted += X_concat[base + a] * beta[a];
        out_IF_per_row[r] = (Y_concat[r] - fitted) - att;
      } else {
        out_IF_per_row[r] = 0.0;
      }
    }
  }
}


// Host launcher (Linux internal form, used by the Rcpp wrapper).
//
// Differs from the public ABI in inst/include/didgpu_cuda_api.h in
// that influence is returned as a PER-ROW vector (length n_total)
// rather than the (n_units x n_cells) public layout. The Rcpp
// wrapper does the unit-major expansion using a row-to-unit map it
// receives from the R side.
//
// Per-row IF is the natural output of the kernel because the kernel
// only knows about row indices within each cell, not the unit
// identity of those rows. Pushing the row→unit scatter to the host
// keeps the kernel simple and avoids inflating the kernel signature
// with a per-row unit_id buffer.
extern "C" int didgpu_cuda_cs_inner_or(
    const double* h_X_concat, const int* h_X_offsets,
    const double* h_Y_concat,
    const double* h_W_concat,
    int n_cells, int p,
    double* h_out_att,
    double* h_out_IF_per_row /* length n_total, NULL to skip */) {

  if (n_cells <= 0 || p <= 0 || p > DIDGPU_CS_MAX_P) return -3;
  const int n_total = h_X_offsets[n_cells];
  if (n_total <= 0) return -3;

  cudaError_t e;
  double* d_X = nullptr;
  double* d_Y = nullptr;
  double* d_W = nullptr;
  int*    d_off = nullptr;
  double* d_att = nullptr;
  double* d_IF  = nullptr;
  auto cleanup = [&]() {
    if (d_X)   cudaFree(d_X);
    if (d_Y)   cudaFree(d_Y);
    if (d_W)   cudaFree(d_W);
    if (d_off) cudaFree(d_off);
    if (d_att) cudaFree(d_att);
    if (d_IF)  cudaFree(d_IF);
  };

  e = cudaMalloc((void**)&d_X,   sizeof(double) * n_total * p);
  if (e != cudaSuccess) { cleanup(); return -4; }
  e = cudaMalloc((void**)&d_Y,   sizeof(double) * n_total);
  if (e != cudaSuccess) { cleanup(); return -4; }
  e = cudaMalloc((void**)&d_W,   sizeof(double) * n_total);
  if (e != cudaSuccess) { cleanup(); return -4; }
  e = cudaMalloc((void**)&d_off, sizeof(int)    * (n_cells + 1));
  if (e != cudaSuccess) { cleanup(); return -4; }
  e = cudaMalloc((void**)&d_att, sizeof(double) * n_cells);
  if (e != cudaSuccess) { cleanup(); return -4; }
  if (h_out_IF_per_row) {
    e = cudaMalloc((void**)&d_IF, sizeof(double) * n_total);
    if (e != cudaSuccess) { cleanup(); return -4; }
  }

  e = cudaMemcpy(d_X, h_X_concat, sizeof(double) * n_total * p,
                  cudaMemcpyHostToDevice);
  if (e != cudaSuccess) { cleanup(); return -1; }
  e = cudaMemcpy(d_Y, h_Y_concat, sizeof(double) * n_total,
                  cudaMemcpyHostToDevice);
  if (e != cudaSuccess) { cleanup(); return -1; }
  e = cudaMemcpy(d_W, h_W_concat, sizeof(double) * n_total,
                  cudaMemcpyHostToDevice);
  if (e != cudaSuccess) { cleanup(); return -1; }
  e = cudaMemcpy(d_off, h_X_offsets, sizeof(int) * (n_cells + 1),
                  cudaMemcpyHostToDevice);
  if (e != cudaSuccess) { cleanup(); return -1; }

  const int block = 64;
  const int grid  = (n_cells + block - 1) / block;
  k_cs_inner_or<<<grid, block>>>(
      d_X, d_Y, d_W, d_off, n_cells, p, d_att, d_IF);
  e = cudaGetLastError();
  if (e != cudaSuccess) { cleanup(); return -1; }

  e = cudaMemcpy(h_out_att, d_att, sizeof(double) * n_cells,
                  cudaMemcpyDeviceToHost);
  if (e != cudaSuccess) { cleanup(); return -1; }
  if (h_out_IF_per_row) {
    e = cudaMemcpy(h_out_IF_per_row, d_IF, sizeof(double) * n_total,
                    cudaMemcpyDeviceToHost);
  }
  cleanup();
  return (e != cudaSuccess) ? -1 : 0;
}


// ============================================================================
// IPW (est_method = 1) and DR (est_method = 2) inner regression.
//
// Mirrors R/cs_methods.R::.cs_inner_ipw and .cs_inner_dr branch-for-
// branch, including:
//   * the no-covariate (p == 1) closed form with its SPECIAL influence
//     function (different from the general w1*δ - w0*δ - att form);
//   * a per-cell logistic regression (IRLS) for the propensity score,
//     replicating stats::glm.fit (mustart = (y+0.5)/2, maxit = 25,
//     deviance convergence eps = 1e-8, WLS weight mu*(1-mu));
//   * propensity trimming to [0.01, 0.99];
//   * the fallbacks: IPW on a rank-deficient design -> constant
//     propensity (= no-cov form); DR on a rank-deficient propensity
//     design -> OR (.cs_inner_or).
// ============================================================================

// Solve an SPD system A x = b (p x p, row-major) via Cholesky. A is
// read (not overwritten); b holds the RHS on entry and the solution
// on exit. Returns false if A is not positive-definite (rank-deficient),
// which the callers treat as glm.fit's "NA coefficients" signal.
__device__ inline bool cs_spd_solve(const double* A, double* b, int p) {
  double L[DIDGPU_CS_MAX_P * DIDGPU_CS_MAX_P];
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = A[i * p + j];
      for (int k = 0; k < j; ++k) s -= L[i * p + k] * L[j * p + k];
      if (i == j) {
        if (s <= 0.0) return false;
        L[i * p + i] = sqrt(s);
      } else {
        L[i * p + j] = s / L[j * p + j];
      }
    }
  }
  double y[DIDGPU_CS_MAX_P];
  for (int i = 0; i < p; ++i) {
    double s = b[i];
    for (int j = 0; j < i; ++j) s -= L[i * p + j] * y[j];
    y[i] = s / L[i * p + i];
  }
  for (int i = p - 1; i >= 0; --i) {
    double s = y[i];
    for (int j = i + 1; j < p; ++j) s -= L[j * p + i] * b[j];
    b[i] = s / L[i * p + i];
  }
  return true;
}

// IRLS logistic fit of D ~ X over rows [r0, r1). X is row-major with p
// columns (intercept included as column 0). Writes beta[0..p). Returns
// false if any WLS solve is singular (rank-deficient design).
__device__ bool cs_logit_fit(
    const double* X, const double* D,
    int r0, int r1, int p, double* beta) {
  for (int j = 0; j < p; ++j) beta[j] = 0.0;
  double dev_old = 0.0;
  bool have_dev = false;
  for (int iter = 0; iter < 25; ++iter) {
    double A[DIDGPU_CS_MAX_P * DIDGPU_CS_MAX_P];
    double b[DIDGPU_CS_MAX_P];
    for (int t = 0; t < p * p; ++t) A[t] = 0.0;
    for (int t = 0; t < p; ++t) b[t] = 0.0;
    double dev = 0.0;
    for (int r = r0; r < r1; ++r) {
      const double* xr = X + r * p;
      const double y = D[r];
      double eta;
      if (iter == 0) {
        const double mus = (y + 0.5) * 0.5;        // glm.fit mustart
        eta = log(mus / (1.0 - mus));              // logit link
      } else {
        eta = 0.0;
        for (int j = 0; j < p; ++j) eta += xr[j] * beta[j];
      }
      double mu = 1.0 / (1.0 + exp(-eta));
      if (mu < 1e-10)        mu = 1e-10;
      if (mu > 1.0 - 1e-10)  mu = 1.0 - 1e-10;
      const double W = mu * (1.0 - mu);            // IRLS weight
      const double bcoef = W * eta + (y - mu);     // X_i row's b-contribution
      for (int a = 0; a < p; ++a) {
        b[a] += xr[a] * bcoef;
        for (int c2 = 0; c2 < p; ++c2) A[a * p + c2] += W * xr[a] * xr[c2];
      }
      dev += (y > 0.5) ? (-2.0 * log(mu)) : (-2.0 * log(1.0 - mu));
    }
    if (!cs_spd_solve(A, b, p)) return false;      // singular -> fallback
    for (int j = 0; j < p; ++j) beta[j] = b[j];
    if (have_dev && fabs(dev - dev_old) / (fabs(dev) + 0.1) < 1e-8) break;
    dev_old = dev;
    have_dev = true;
  }
  return true;
}

// No-covariate IPW/DR closed form (the p == 1 branch and the IPW
// rank-deficient fallback). att = mean(δ_t) - mean(δ_c); special IF.
__device__ void cs_nocov_diff(
    const double* Y, const double* W, int r0, int r1,
    int n_t, int n_c, int c,
    double* out_att, double* out_IF_per_row) {
  double sum_t = 0.0, sum_c = 0.0;
  for (int r = r0; r < r1; ++r) {
    if (W[r] > 0.5) sum_t += Y[r]; else sum_c += Y[r];
  }
  const double mean_t = sum_t / n_t;
  const double mean_c = sum_c / n_c;
  const double att = mean_t - mean_c;
  out_att[c] = att;
  if (out_IF_per_row) {
    for (int r = r0; r < r1; ++r) {
      out_IF_per_row[r] = (W[r] > 0.5)
          ? (Y[r] - mean_t - 0.5 * att)
          : (-(Y[r] - mean_c) - 0.5 * att);
    }
  }
}

__global__ void k_cs_inner_logit(
    const double* __restrict__ X_concat,   // (n_total, p) row-major
    const double* __restrict__ Y_concat,   // delta
    const double* __restrict__ W_concat,   // D in {0, 1}
    const int*    __restrict__ X_offsets,  // (n_cells + 1)
    int n_cells, int p, int est_method,    // 1 = IPW, 2 = DR
    double* __restrict__ out_att,
    double* __restrict__ out_IF_per_row) {
  const int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= n_cells) return;

  const int r0 = X_offsets[c];
  const int r1 = X_offsets[c + 1];
  int n_t = 0, n_c = 0;
  for (int r = r0; r < r1; ++r) { if (W_concat[r] > 0.5) ++n_t; else ++n_c; }

  if (n_t == 0 || n_c == 0) {
    out_att[c] = nan("");
    if (out_IF_per_row) for (int r = r0; r < r1; ++r) out_IF_per_row[r] = 0.0;
    return;
  }

  const double E_D = static_cast<double>(n_t) / static_cast<double>(n_t + n_c);
  const double n   = static_cast<double>(n_t + n_c);

  // ---- no-covariate branch: IPW == DR == simple difference ----
  if (p == 1) {
    cs_nocov_diff(Y_concat, W_concat, r0, r1, n_t, n_c, c,
                  out_att, out_IF_per_row);
    return;
  }

  // ---- DR outcome model: OR on controls (also used by DR fallback) ----
  double beta_or[DIDGPU_CS_MAX_P];
  bool or_ok = false;
  double mean_dc = 0.0;
  if (est_method == 2) {
    double Aor[DIDGPU_CS_MAX_P * DIDGPU_CS_MAX_P];
    double bor[DIDGPU_CS_MAX_P];
    for (int t = 0; t < p * p; ++t) Aor[t] = 0.0;
    for (int t = 0; t < p; ++t) bor[t] = 0.0;
    double sum_dc = 0.0;
    for (int r = r0; r < r1; ++r) {
      if (W_concat[r] > 0.5) continue;             // controls only
      const double* xr = X_concat + r * p;
      const double yv = Y_concat[r];
      sum_dc += yv;
      for (int a = 0; a < p; ++a) {
        bor[a] += xr[a] * yv;
        for (int c2 = 0; c2 < p; ++c2) Aor[a * p + c2] += xr[a] * xr[c2];
      }
    }
    mean_dc = sum_dc / n_c;
    or_ok = cs_spd_solve(Aor, bor, p);
    if (or_ok) for (int j = 0; j < p; ++j) beta_or[j] = bor[j];
  }

  // ---- propensity model: logistic D ~ X ----
  double beta_ps[DIDGPU_CS_MAX_P];
  const bool ps_ok = cs_logit_fit(X_concat, W_concat, r0, r1, p, beta_ps);

  // ---- fallbacks when the propensity design is rank-deficient ----
  if (!ps_ok) {
    if (est_method == 1) {
      // IPW -> constant propensity (= no-cov simple difference).
      cs_nocov_diff(Y_concat, W_concat, r0, r1, n_t, n_c, c,
                    out_att, out_IF_per_row);
    } else {
      // DR -> OR (.cs_inner_or): att = mean over treated of
      // (δ - m_hat); IF[treated] = (δ - m_hat) - att, IF[control] = 0.
      double sum_resid_t = 0.0;
      for (int r = r0; r < r1; ++r) {
        if (W_concat[r] <= 0.5) continue;
        const double* xr = X_concat + r * p;
        double m = mean_dc;
        if (or_ok) { m = 0.0; for (int j = 0; j < p; ++j) m += xr[j] * beta_or[j]; }
        sum_resid_t += (Y_concat[r] - m);
      }
      const double att = sum_resid_t / n_t;
      out_att[c] = att;
      if (out_IF_per_row) {
        for (int r = r0; r < r1; ++r) {
          if (W_concat[r] > 0.5) {
            const double* xr = X_concat + r * p;
            double m = mean_dc;
            if (or_ok) { m = 0.0; for (int j = 0; j < p; ++j) m += xr[j] * beta_or[j]; }
            out_IF_per_row[r] = (Y_concat[r] - m) - att;
          } else {
            out_IF_per_row[r] = 0.0;
          }
        }
      }
    }
    return;
  }

  // ---- main IPW / DR path (propensity ok) ----
  // att = (1/n) sum_i (w1_i - w0_i) * v_i, where v = δ (IPW) or
  // v = δ - m_hat (DR); w1_i = D_i/E_D, w0_i = (1-D_i)(p/(1-p))/E_D.
  double sum_att = 0.0;
  for (int r = r0; r < r1; ++r) {
    const double* xr = X_concat + r * p;
    double eta = 0.0;
    for (int j = 0; j < p; ++j) eta += xr[j] * beta_ps[j];
    double ph = 1.0 / (1.0 + exp(-eta));
    if (ph < 0.01) ph = 0.01;
    if (ph > 0.99) ph = 0.99;
    double v = Y_concat[r];
    if (est_method == 2) {
      double m = mean_dc;
      if (or_ok) { m = 0.0; for (int j = 0; j < p; ++j) m += xr[j] * beta_or[j]; }
      v -= m;
    }
    const bool treated = (W_concat[r] > 0.5);
    const double w1 = treated ? (1.0 / E_D) : 0.0;
    const double w0 = treated ? 0.0 : ((ph / (1.0 - ph)) / E_D);
    sum_att += (w1 - w0) * v;
  }
  const double att = sum_att / n;
  out_att[c] = att;

  if (out_IF_per_row) {
    for (int r = r0; r < r1; ++r) {
      const double* xr = X_concat + r * p;
      double eta = 0.0;
      for (int j = 0; j < p; ++j) eta += xr[j] * beta_ps[j];
      double ph = 1.0 / (1.0 + exp(-eta));
      if (ph < 0.01) ph = 0.01;
      if (ph > 0.99) ph = 0.99;
      double v = Y_concat[r];
      if (est_method == 2) {
        double m = mean_dc;
        if (or_ok) { m = 0.0; for (int j = 0; j < p; ++j) m += xr[j] * beta_or[j]; }
        v -= m;
      }
      const bool treated = (W_concat[r] > 0.5);
      const double w1 = treated ? (1.0 / E_D) : 0.0;
      const double w0 = treated ? 0.0 : ((ph / (1.0 - ph)) / E_D);
      out_IF_per_row[r] = (w1 - w0) * v - att;
    }
  }
}

// Host launcher for IPW / DR (est_method 1 / 2). Same buffer plumbing
// as didgpu_cuda_cs_inner_or; dispatches to k_cs_inner_logit.
extern "C" int didgpu_cuda_cs_inner_logit(
    const double* h_X_concat, const int* h_X_offsets,
    const double* h_Y_concat,
    const double* h_W_concat,
    int n_cells, int p, int est_method,
    double* h_out_att,
    double* h_out_IF_per_row) {

  if (est_method != 1 && est_method != 2) return -3;
  if (n_cells <= 0 || p <= 0 || p > DIDGPU_CS_MAX_P) return -3;
  const int n_total = h_X_offsets[n_cells];
  if (n_total <= 0) return -3;

  cudaError_t e;
  double* d_X = nullptr; double* d_Y = nullptr; double* d_W = nullptr;
  int* d_off = nullptr; double* d_att = nullptr; double* d_IF = nullptr;
  auto cleanup = [&]() {
    if (d_X) cudaFree(d_X); if (d_Y) cudaFree(d_Y); if (d_W) cudaFree(d_W);
    if (d_off) cudaFree(d_off); if (d_att) cudaFree(d_att);
    if (d_IF) cudaFree(d_IF);
  };

  e = cudaMalloc((void**)&d_X,   sizeof(double) * n_total * p); if (e) { cleanup(); return -4; }
  e = cudaMalloc((void**)&d_Y,   sizeof(double) * n_total);     if (e) { cleanup(); return -4; }
  e = cudaMalloc((void**)&d_W,   sizeof(double) * n_total);     if (e) { cleanup(); return -4; }
  e = cudaMalloc((void**)&d_off, sizeof(int) * (n_cells + 1));  if (e) { cleanup(); return -4; }
  e = cudaMalloc((void**)&d_att, sizeof(double) * n_cells);     if (e) { cleanup(); return -4; }
  if (h_out_IF_per_row) {
    e = cudaMalloc((void**)&d_IF, sizeof(double) * n_total);    if (e) { cleanup(); return -4; }
  }

  e = cudaMemcpy(d_X, h_X_concat, sizeof(double) * n_total * p, cudaMemcpyHostToDevice); if (e) { cleanup(); return -1; }
  e = cudaMemcpy(d_Y, h_Y_concat, sizeof(double) * n_total, cudaMemcpyHostToDevice);     if (e) { cleanup(); return -1; }
  e = cudaMemcpy(d_W, h_W_concat, sizeof(double) * n_total, cudaMemcpyHostToDevice);     if (e) { cleanup(); return -1; }
  e = cudaMemcpy(d_off, h_X_offsets, sizeof(int) * (n_cells + 1), cudaMemcpyHostToDevice); if (e) { cleanup(); return -1; }

  const int block = 64;
  const int grid  = (n_cells + block - 1) / block;
  k_cs_inner_logit<<<grid, block>>>(
      d_X, d_Y, d_W, d_off, n_cells, p, est_method, d_att, d_IF);
  e = cudaGetLastError();
  if (e != cudaSuccess) { cleanup(); return -1; }

  e = cudaMemcpy(h_out_att, d_att, sizeof(double) * n_cells, cudaMemcpyDeviceToHost);
  if (e != cudaSuccess) { cleanup(); return -1; }
  if (h_out_IF_per_row) {
    e = cudaMemcpy(h_out_IF_per_row, d_IF, sizeof(double) * n_total, cudaMemcpyDeviceToHost);
  }
  cleanup();
  return (e != cudaSuccess) ? -1 : 0;
}


// Public-ABI form (matches inst/include/didgpu_cuda_api.h). For
// est_method = 0 (OR) it forwards to didgpu_cuda_cs_inner_or with
// no IF output (since the public ABI doesn't carry the row→unit
// scatter map). For est_method = 1, 2 returns -3 (not implemented).
//
// The Rcpp wrapper calls didgpu_cuda_cs_inner_or directly rather
// than this entry point, so this exists primarily for ABI
// compatibility with the future Windows two-DLL build.
extern "C" int didgpu_cuda_cs_inner_batched(
    const double* h_X_concat, const int* h_X_offsets,
    const double* h_Y_concat, const int* /*Y_offsets*/,
    const double* h_W_concat, const int* /*W_offsets*/,
    int n_cells, int p, int /*n_units*/,
    int est_method,
    double* h_out_att, double* /*out_influence*/) {
  if (est_method != 0) return -3;  // IPW + DR: follow-up.
  return didgpu_cuda_cs_inner_or(
      h_X_concat, h_X_offsets, h_Y_concat, h_W_concat,
      n_cells, p, h_out_att, /*IF=*/nullptr);
}

#endif  // HAS_CUDA
