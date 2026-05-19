// Native routine registration and the R-facing entry points that
// dispatch to CPU or CUDA implementations.
//
// The CUDA path is only compiled when HAS_CUDA is defined (set by
// Makevars/Makevars.win when nvcc is detected). Otherwise the CUDA
// stubs return NA so callers can detect-and-fall-back at runtime.

#include <Rcpp.h>

#ifdef HAS_CUDA
#include <cuda_runtime.h>
extern "C" int didgpu_cuda_saxpy(int n, float a, const float* x, float* y);
extern "C" int didgpu_cuda_run_one_event_time(
    const double* outcome, const double* N_gt,
    const int* row_to_g, const int* row_to_t, const int* cohort_key,
    const int* F_g, const int* S_g, const int* T_g, const int* L_g,
    int n_rows, int n_groups, int n_cohorts,
    int k, int direction, double G_over_Ninc,
    double* diff_y_k, int* never_change_k, int* candidate_dist_k,
    int* dist_k_final, double* kernel_val, double* U_g,
    double* N_t_control, double* N_t_switch_cand,
    double* did_out_device);
extern "C" int didgpu_cuda_fect_fe(
    const double* d_Y, const int* d_M,
    double* d_alpha, double* d_xi,
    int n_units, int n_periods,
    double tol, int max_iter,
    int* out_iter, double* out_delta);
extern "C" int didgpu_cuda_cs_inner_batched(
    const double* X_concat, const int* X_offsets,
    const double* Y_concat, const int* Y_offsets,
    const double* W_concat, const int* W_offsets,
    int n_cells, int p, int n_units,
    int est_method,
    double* out_att, double* out_influence);
extern "C" int didgpu_cuda_cs_inner_or(
    const double* X_concat, const int* X_offsets,
    const double* Y_concat,
    const double* W_concat,
    int n_cells, int p,
    double* out_att,
    double* out_IF_per_row);
extern "C" int didgpu_cuda_fect_svd_truncated(
    const double* d_M_rm,
    int m, int n, int r,
    double* d_L_out_rm,
    double* d_F_out_rm);
extern "C" int didgpu_cuda_fect_svd_softthreshold(
    const double* d_Y_complete_rm,
    int m, int n,
    double lambda,
    double* d_Y_hat_rm,
    int* out_n_nonzero);
extern "C" int didgpu_cuda_testmechs_bootstrap(
    const int* h_d, const int* h_m, const int* h_y,
    int n, int K, int dy, int B,
    unsigned long seed,
    double* h_beta);
extern "C" int didgpu_cuda_cluster_bootstrap(
    const double* h_IF, int n_units, int n_dims,
    const int* h_cluster_id, int n_clusters,
    int B, unsigned long long seed,
    double* h_out_estimates);
extern "C" int didgpu_cuda_multiplier_bootstrap(
    const double* h_IF, int n_units, int n_dims,
    int B, int mult_kind,
    unsigned long long seed,
    double* h_out_estimates);
#endif

// [[Rcpp::export]]
bool didgpu_has_cuda_support() {
#ifdef HAS_CUDA
  return true;
#else
  return false;
#endif
}

// Run one (event-time, direction) DID kernel on the GPU and return the
// scalar DID estimate. All inputs are R vectors of the right length and
// type; this function does the H2D/D2H transfers and workspace alloc.
//
// Arguments mirror the .core_one_event_time() interface for the binary
// no-controls case; see R/core_r.R.
//
// Returns the per-event-time DID estimate (sum_g U_g / G).
//
// [[Rcpp::export]]
double didgpu_cuda_did(
    Rcpp::NumericVector outcome,
    Rcpp::NumericVector N_gt,
    Rcpp::IntegerVector row_to_g,
    Rcpp::IntegerVector row_to_t,
    Rcpp::IntegerVector cohort_key,
    Rcpp::IntegerVector F_g,
    Rcpp::IntegerVector S_g,
    Rcpp::IntegerVector T_g,
    Rcpp::IntegerVector L_g,
    int n_cohorts,
    int k, int direction, double G_over_Ninc) {
#ifdef HAS_CUDA
  const int n_rows = outcome.size();
  const int n_groups = F_g.size();
  if (N_gt.size() != n_rows || row_to_g.size() != n_rows
      || row_to_t.size() != n_rows || cohort_key.size() != n_rows)
    Rcpp::stop("vector length mismatch");
  if (S_g.size() != n_groups || T_g.size() != n_groups || L_g.size() != n_groups)
    Rcpp::stop("group-vector length mismatch");

  // Device allocations.
  double *d_outcome=nullptr, *d_Ngt=nullptr, *d_diff=nullptr, *d_kernel=nullptr;
  double *d_Ug=nullptr, *d_Nctrl=nullptr, *d_Nswitch=nullptr, *d_did=nullptr;
  int *d_rtog=nullptr, *d_rtot=nullptr, *d_ckey=nullptr;
  int *d_Fg=nullptr, *d_Sg=nullptr, *d_Tg=nullptr, *d_Lg=nullptr;
  int *d_nck=nullptr, *d_cdist=nullptr, *d_dist=nullptr;

  auto fail = [&](const char* msg) -> double {
    if (d_outcome) cudaFree(d_outcome);
    if (d_Ngt)     cudaFree(d_Ngt);
    if (d_diff)    cudaFree(d_diff);
    if (d_kernel)  cudaFree(d_kernel);
    if (d_Ug)      cudaFree(d_Ug);
    if (d_Nctrl)   cudaFree(d_Nctrl);
    if (d_Nswitch) cudaFree(d_Nswitch);
    if (d_did)     cudaFree(d_did);
    if (d_rtog)    cudaFree(d_rtog);
    if (d_rtot)    cudaFree(d_rtot);
    if (d_ckey)    cudaFree(d_ckey);
    if (d_Fg)      cudaFree(d_Fg);
    if (d_Sg)      cudaFree(d_Sg);
    if (d_Tg)      cudaFree(d_Tg);
    if (d_Lg)      cudaFree(d_Lg);
    if (d_nck)     cudaFree(d_nck);
    if (d_cdist)   cudaFree(d_cdist);
    if (d_dist)    cudaFree(d_dist);
    Rcpp::stop("CUDA error: %s", msg);
    return 0.0;
  };

  cudaError_t e;
  #define ALLOC(p, n, T) do { e = cudaMalloc((void**)&p, (n) * sizeof(T)); \
                              if (e != cudaSuccess) return fail(cudaGetErrorString(e)); } while(0)
  #define H2D(dst, src, n, T) do { e = cudaMemcpy(dst, src, (n) * sizeof(T), \
                                                    cudaMemcpyHostToDevice); \
                                    if (e != cudaSuccess) return fail(cudaGetErrorString(e)); } while(0)

  ALLOC(d_outcome, n_rows, double);
  ALLOC(d_Ngt,     n_rows, double);
  ALLOC(d_rtog,    n_rows, int);
  ALLOC(d_rtot,    n_rows, int);
  ALLOC(d_ckey,    n_rows, int);
  ALLOC(d_Fg,      n_groups, int);
  ALLOC(d_Sg,      n_groups, int);
  ALLOC(d_Tg,      n_groups, int);
  ALLOC(d_Lg,      n_groups, int);
  ALLOC(d_diff,    n_rows, double);
  ALLOC(d_nck,     n_rows, int);
  ALLOC(d_cdist,   n_rows, int);
  ALLOC(d_dist,    n_rows, int);
  ALLOC(d_kernel,  n_rows, double);
  ALLOC(d_Ug,      n_groups, double);
  ALLOC(d_Nctrl,   n_cohorts, double);
  ALLOC(d_Nswitch, n_cohorts, double);
  ALLOC(d_did,     1, double);

  H2D(d_outcome, &outcome[0],     n_rows,   double);
  H2D(d_Ngt,     &N_gt[0],        n_rows,   double);
  H2D(d_rtog,    &row_to_g[0],    n_rows,   int);
  H2D(d_rtot,    &row_to_t[0],    n_rows,   int);
  H2D(d_ckey,    &cohort_key[0],  n_rows,   int);
  H2D(d_Fg,      &F_g[0],         n_groups, int);
  H2D(d_Sg,      &S_g[0],         n_groups, int);
  H2D(d_Tg,      &T_g[0],         n_groups, int);
  H2D(d_Lg,      &L_g[0],         n_groups, int);

  int ec = didgpu_cuda_run_one_event_time(
      d_outcome, d_Ngt, d_rtog, d_rtot, d_ckey,
      d_Fg, d_Sg, d_Tg, d_Lg,
      n_rows, n_groups, n_cohorts,
      k, direction, G_over_Ninc,
      d_diff, d_nck, d_cdist, d_dist, d_kernel, d_Ug,
      d_Nctrl, d_Nswitch, d_did);
  if (ec != 0) return fail("kernel chain failed");

  double did = 0.0;
  e = cudaMemcpy(&did, d_did, sizeof(double), cudaMemcpyDeviceToHost);
  if (e != cudaSuccess) return fail(cudaGetErrorString(e));

  cudaFree(d_outcome); cudaFree(d_Ngt); cudaFree(d_diff); cudaFree(d_kernel);
  cudaFree(d_Ug); cudaFree(d_Nctrl); cudaFree(d_Nswitch); cudaFree(d_did);
  cudaFree(d_rtog); cudaFree(d_rtot); cudaFree(d_ckey);
  cudaFree(d_Fg); cudaFree(d_Sg); cudaFree(d_Tg); cudaFree(d_Lg);
  cudaFree(d_nck); cudaFree(d_cdist); cudaFree(d_dist);

  return did;
  #undef ALLOC
  #undef H2D
#else
  (void)outcome; (void)N_gt; (void)row_to_g; (void)row_to_t; (void)cohort_key;
  (void)F_g; (void)S_g; (void)T_g; (void)L_g;
  (void)n_cohorts; (void)k; (void)direction; (void)G_over_Ninc;
  Rcpp::stop("didgpu was built without CUDA support. "
             "Install the NVIDIA CUDA Toolkit so nvcc is on PATH, then reinstall.");
#endif
}

// Run the fect_fe iterative two-way demeaning on the GPU.
// Inputs: Y (n_units x n_periods, row-major NumericMatrix; NaN = missing)
//         M (same shape; 1 = treated/excluded, 0 = control/included)
// Returns a list with alpha (length n_units), xi (length n_periods),
// iter (int), delta (double).
//
// [[Rcpp::export]]
Rcpp::List didgpu_cuda_fect_fe_r(
    Rcpp::NumericMatrix Y,
    Rcpp::IntegerMatrix M,
    double tol,
    int max_iter) {
#ifdef HAS_CUDA
  const int n_units   = Y.nrow();
  const int n_periods = Y.ncol();
  if (M.nrow() != n_units || M.ncol() != n_periods)
    Rcpp::stop("Y and M must have the same dimensions");

  // R stores matrices column-major, but our kernel wants row-major
  // (Y[i, t] at index i * n_periods + t). Transpose on the host before
  // sending to the GPU.
  std::vector<double> Y_rm(static_cast<size_t>(n_units) * n_periods);
  std::vector<int>    M_rm(static_cast<size_t>(n_units) * n_periods);
  for (int i = 0; i < n_units; ++i) {
    for (int t = 0; t < n_periods; ++t) {
      Y_rm[i * n_periods + t] = Y(i, t);
      M_rm[i * n_periods + t] = M(i, t);
    }
  }

  double *d_Y = nullptr, *d_alpha = nullptr, *d_xi = nullptr;
  int    *d_M = nullptr;
  auto fail = [&](const char* msg) -> Rcpp::List {
    if (d_Y)     cudaFree(d_Y);
    if (d_M)     cudaFree(d_M);
    if (d_alpha) cudaFree(d_alpha);
    if (d_xi)    cudaFree(d_xi);
    Rcpp::stop("CUDA error: %s", msg);
    return Rcpp::List::create();
  };

  cudaError_t e;
  e = cudaMalloc((void**)&d_Y,     n_units * n_periods * sizeof(double));
  if (e != cudaSuccess) return fail(cudaGetErrorString(e));
  e = cudaMalloc((void**)&d_M,     n_units * n_periods * sizeof(int));
  if (e != cudaSuccess) return fail(cudaGetErrorString(e));
  e = cudaMalloc((void**)&d_alpha, n_units * sizeof(double));
  if (e != cudaSuccess) return fail(cudaGetErrorString(e));
  e = cudaMalloc((void**)&d_xi,    n_periods * sizeof(double));
  if (e != cudaSuccess) return fail(cudaGetErrorString(e));

  e = cudaMemcpy(d_Y, Y_rm.data(),
                  n_units * n_periods * sizeof(double),
                  cudaMemcpyHostToDevice);
  if (e != cudaSuccess) return fail(cudaGetErrorString(e));
  e = cudaMemcpy(d_M, M_rm.data(),
                  n_units * n_periods * sizeof(int),
                  cudaMemcpyHostToDevice);
  if (e != cudaSuccess) return fail(cudaGetErrorString(e));

  int iter = 0;
  double delta = 0.0;
  int ec = didgpu_cuda_fect_fe(d_Y, d_M, d_alpha, d_xi,
                                  n_units, n_periods,
                                  tol, max_iter, &iter, &delta);
  if (ec != 0) return fail("kernel chain failed");

  Rcpp::NumericVector alpha(n_units), xi(n_periods);
  cudaMemcpy(&alpha[0], d_alpha, n_units * sizeof(double),
              cudaMemcpyDeviceToHost);
  cudaMemcpy(&xi[0],    d_xi,    n_periods * sizeof(double),
              cudaMemcpyDeviceToHost);
  cudaFree(d_Y); cudaFree(d_M); cudaFree(d_alpha); cudaFree(d_xi);

  return Rcpp::List::create(
    Rcpp::_["alpha"] = alpha,
    Rcpp::_["xi"]    = xi,
    Rcpp::_["iter"]  = iter,
    Rcpp::_["delta"] = delta);
#else
  (void)Y; (void)M; (void)tol; (void)max_iter;
  Rcpp::stop("didgpu was built without CUDA support. "
             "Install the NVIDIA CUDA Toolkit and reinstall.");
#endif
}


// [[Rcpp::export]]
Rcpp::NumericVector didgpu_run_saxpy(double a, Rcpp::NumericVector x, Rcpp::NumericVector y) {
#ifdef HAS_CUDA
  if (x.size() != y.size()) Rcpp::stop("x and y must be the same length");
  int n = x.size();
  std::vector<float> xf(n), yf(n);
  for (int i = 0; i < n; ++i) { xf[i] = static_cast<float>(x[i]); yf[i] = static_cast<float>(y[i]); }
  int err = didgpu_cuda_saxpy(n, static_cast<float>(a), xf.data(), yf.data());
  if (err) Rcpp::stop("CUDA SAXPY failed with code %d", err);
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = static_cast<double>(yf[i]);
  return out;
#else
  (void)a; (void)x; (void)y;
  Rcpp::stop("didgpu was built without CUDA support. Reinstall after installing the NVIDIA CUDA Toolkit so nvcc is on PATH.");
#endif
}


// Multiplier (wild) bootstrap on per-unit influence functions.
//
// Computes B bootstrap replicates of
//
//   out[b, d] = sum_i xi[i, b] * IF[i, d]
//
// where xi[i, b] is i.i.d. Rademacher (mult_kind = 0; +1/-1 with
// equal probability) or N(0, 1) (mult_kind = 1).
//
// Inputs:
//   IF        : R numeric matrix (n_units, n_dims).
//   B         : number of bootstrap replicates.
//   mult_kind : 0 = Rademacher, 1 = N(0, 1).
//   seed      : cuRAND seed.
//
// Returns NULL on any CUDA failure; a (B, n_dims) numeric matrix
// otherwise. As with the other bootstrap entry points, cuRAND and
// R's MT19937 produce different per-replicate draws — column SDs
// converge to the same population SE.
//
// [[Rcpp::export]]
SEXP didgpu_cuda_multiplier_bootstrap_r(Rcpp::NumericMatrix IF,
                                          int B,
                                          int mult_kind,
                                          int seed) {
#ifdef HAS_CUDA
  const int n_units = IF.nrow();
  const int n_dims  = IF.ncol();
  if (B <= 0)
    Rcpp::stop("B must be positive; got %d", B);
  if (mult_kind != 0 && mult_kind != 1)
    Rcpp::stop("mult_kind must be 0 (Rademacher) or 1 (N(0,1)); got %d",
               mult_kind);

  std::vector<double> IF_rm(static_cast<size_t>(n_units) * n_dims);
  for (int i = 0; i < n_units; ++i)
    for (int j = 0; j < n_dims; ++j)
      IF_rm[static_cast<size_t>(i) * n_dims + j] = IF(i, j);

  std::vector<double> out(static_cast<size_t>(B) * n_dims);
  int rc = didgpu_cuda_multiplier_bootstrap(
      IF_rm.data(), n_units, n_dims,
      B, mult_kind,
      static_cast<unsigned long long>(seed),
      out.data());
  if (rc != 0) return R_NilValue;

  Rcpp::NumericMatrix out_mat(B, n_dims);
  for (int b = 0; b < B; ++b)
    for (int d = 0; d < n_dims; ++d)
      out_mat(b, d) = out[static_cast<size_t>(b) * n_dims + d];
  return out_mat;
#else
  (void)IF; (void)B; (void)mult_kind; (void)seed;
  return R_NilValue;
#endif
}


// Cluster bootstrap on per-unit influence functions.
//
// Computes, on the GPU, B bootstrap replicates of the influence-
// function-weighted point estimate
//
//   out[b, d] = sum_i  w_{i, b} * IF[i, d]
//
// where w_{i, b} is the count of times unit i's cluster was sampled
// in replicate b's draw of n_clusters clusters with replacement.
//
// This is the "delta-method shortcut" — orders of magnitude faster
// than re-running the full estimator per replicate, with the same
// asymptotic distribution (Hansen 2022, Ch.10).
//
// Inputs:
//   IF         : R numeric matrix (n_units, n_dims), row-major from
//                R's POV (R matrices are column-major, but we convert
//                row-major on the way in to match the kernel ABI).
//   cluster_id : R integer vector (n_units), 0-based cluster IDs in
//                [0, n_clusters).
//   n_clusters : number of distinct clusters (must equal
//                max(cluster_id) + 1; caller enforces).
//   B          : number of bootstrap replicates.
//   seed       : RNG seed for cuRAND. Like .testmechs_bootstrap_cuda,
//                the cuRAND stream differs from R's MT19937 so per-
//                replicate output differs from a R-side cluster
//                bootstrap, but the columnwise SDs converge to the
//                same population SE.
//
// Returns NULL on any CUDA error. On success returns a B x n_dims
// numeric matrix of bootstrap estimates.
//
// [[Rcpp::export]]
SEXP didgpu_cuda_cluster_bootstrap_r(Rcpp::NumericMatrix IF,
                                      Rcpp::IntegerVector cluster_id,
                                      int n_clusters,
                                      int B,
                                      int seed) {
#ifdef HAS_CUDA
  const int n_units = IF.nrow();
  const int n_dims  = IF.ncol();
  if (cluster_id.size() != n_units)
    Rcpp::stop("cluster_id length (%d) must equal nrow(IF) (%d)",
               (int)cluster_id.size(), n_units);
  if (n_clusters <= 0)
    Rcpp::stop("n_clusters must be positive; got %d", n_clusters);
  if (B <= 0) Rcpp::stop("B must be positive; got %d", B);

  // Convert R column-major IF to row-major host buffer.
  std::vector<double> IF_rm(static_cast<size_t>(n_units) * n_dims);
  for (int i = 0; i < n_units; ++i)
    for (int j = 0; j < n_dims; ++j)
      IF_rm[static_cast<size_t>(i) * n_dims + j] = IF(i, j);

  std::vector<double> out(static_cast<size_t>(B) * n_dims);

  int rc = didgpu_cuda_cluster_bootstrap(
      IF_rm.data(), n_units, n_dims,
      &cluster_id[0], n_clusters,
      B, static_cast<unsigned long long>(seed),
      out.data());
  if (rc != 0) return R_NilValue;

  // Pack row-major host buffer back to a column-major R matrix.
  Rcpp::NumericMatrix out_mat(B, n_dims);
  for (int b = 0; b < B; ++b)
    for (int d = 0; d < n_dims; ++d)
      out_mat(b, d) = out[static_cast<size_t>(b) * n_dims + d];
  return out_mat;
#else
  (void)IF; (void)cluster_id; (void)n_clusters; (void)B; (void)seed;
  return R_NilValue;
#endif
}


// TestMechs nonparametric bootstrap on the GPU.
//
// Inputs:
//   d, m, y : integer-coded observation triples (each length n).
//             d in {0, 1}, m in {1, .., K}, y in {1, .., dy}.
//             1-based for m, y to match the R-side conventions.
//   B       : number of bootstrap replicates.
//   K, dy   : alphabet sizes (caller is responsible for setting these
//             correctly; the kernel assumes they're right).
//   seed    : RNG seed (used by cuRAND, NOT the same stream as R's
//             Mersenne-Twister — outputs differ from .testmechs_
//             bootstrap_r at the same seed, but both are valid bootstraps).
//
// Returns NULL on any CUDA error. On success returns a B x dim_beta
// numeric matrix where dim_beta = 2 * K * dy (per-D block of K*dy
// (m, y) cells, stacked [D=0, D=1]). Each per-D block sums to 1 within
// a row.
//
// The "bayes" method has no GPU path yet (the kernel only does
// nonparametric multinomial resampling); the R-side helper falls back
// to .testmechs_bootstrap_r in that case.
//
// [[Rcpp::export]]
SEXP didgpu_cuda_testmechs_bootstrap_r(
    Rcpp::IntegerVector d,
    Rcpp::IntegerVector m,
    Rcpp::IntegerVector y,
    int B, int K, int dy,
    int seed) {
#ifdef HAS_CUDA
  const int n = d.size();
  if (m.size() != n || y.size() != n)
    Rcpp::stop("d, m, y must be the same length; got %d, %d, %d",
               n, (int)m.size(), (int)y.size());
  if (B  <= 0) Rcpp::stop("B must be positive; got %d", B);
  if (K  <= 0) Rcpp::stop("K must be positive; got %d", K);
  if (dy <= 0) Rcpp::stop("dy must be positive; got %d", dy);

  const int dim_beta = 2 * K * dy;
  std::vector<double> beta(static_cast<size_t>(B) * dim_beta, 0.0);

  int rc = didgpu_cuda_testmechs_bootstrap(
      &d[0], &m[0], &y[0],
      n, K, dy, B,
      static_cast<unsigned long>(seed),
      beta.data());
  if (rc != 0) return R_NilValue;

  // Pack row-major host buffer (B x dim_beta) into a column-major R
  // matrix of the same shape.
  Rcpp::NumericMatrix out(B, dim_beta);
  for (int b = 0; b < B; ++b)
    for (int j = 0; j < dim_beta; ++j)
      out(b, j) = beta[static_cast<size_t>(b) * dim_beta + j];
  return out;
#else
  (void)d; (void)m; (void)y; (void)B; (void)K; (void)dy; (void)seed;
  return R_NilValue;
#endif
}


// Rank-r truncated SVD on the GPU. Returns L (m x r) and F (r x n)
// such that L * F is a rank-r approximation of M (m x n). This is the
// LR decomposition the alternating fect_ife loop consumes directly.
//
// Inputs:
//   M : R numeric matrix (m x n). NaN cells are passed through to the
//       kernel; the caller (fect_ife) is responsible for any masking.
//   r : truncation rank, 1 <= r <= min(m, n).
//
// Returns NULL if CUDA reports any error (caller falls back to
// .fect_svd_r). On success returns list(L, F, status = 0).
//
// [[Rcpp::export]]
SEXP didgpu_cuda_fect_svd_truncated_r(Rcpp::NumericMatrix M, int r) {
#ifdef HAS_CUDA
  const int m = M.nrow();
  const int n = M.ncol();
  if (r <= 0 || r > std::min(m, n)) {
    Rcpp::stop("r must satisfy 1 <= r <= min(m, n); got r=%d, m=%d, n=%d",
               r, m, n);
  }

  // R matrices are column-major; the kernel wants row-major.
  std::vector<double> M_rm(static_cast<size_t>(m) * n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      M_rm[static_cast<size_t>(i) * n + j] = M(i, j);

  double *d_M = nullptr, *d_L = nullptr, *d_F = nullptr;
  cudaError_t e;
  auto cleanup_and_return_null = [&]() -> SEXP {
    if (d_M) cudaFree(d_M);
    if (d_L) cudaFree(d_L);
    if (d_F) cudaFree(d_F);
    return R_NilValue;
  };

  e = cudaMalloc((void**)&d_M, sizeof(double) * m * n);
  if (e != cudaSuccess) return cleanup_and_return_null();
  e = cudaMalloc((void**)&d_L, sizeof(double) * m * r);
  if (e != cudaSuccess) return cleanup_and_return_null();
  e = cudaMalloc((void**)&d_F, sizeof(double) * r * n);
  if (e != cudaSuccess) return cleanup_and_return_null();

  e = cudaMemcpy(d_M, M_rm.data(), sizeof(double) * m * n,
                  cudaMemcpyHostToDevice);
  if (e != cudaSuccess) return cleanup_and_return_null();

  int rc = didgpu_cuda_fect_svd_truncated(d_M, m, n, r, d_L, d_F);
  if (rc != 0) return cleanup_and_return_null();

  std::vector<double> L_rm(static_cast<size_t>(m) * r);
  std::vector<double> F_rm(static_cast<size_t>(r) * n);
  cudaMemcpy(L_rm.data(), d_L, sizeof(double) * m * r,
              cudaMemcpyDeviceToHost);
  cudaMemcpy(F_rm.data(), d_F, sizeof(double) * r * n,
              cudaMemcpyDeviceToHost);
  cudaFree(d_M); cudaFree(d_L); cudaFree(d_F);

  // Convert row-major host -> column-major R matrices.
  Rcpp::NumericMatrix L(m, r), F(r, n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < r; ++j)
      L(i, j) = L_rm[static_cast<size_t>(i) * r + j];
  for (int i = 0; i < r; ++i)
    for (int j = 0; j < n; ++j)
      F(i, j) = F_rm[static_cast<size_t>(i) * n + j];

  Rcpp::List result;
  result["L"]      = L;
  result["F"]      = F;
  result["status"] = 0;
  return result;
#else
  (void)M; (void)r;
  return R_NilValue;
#endif
}


// Soft-thresholded SVD reconstruction on the GPU. Computes
//   Y_hat = U %*% diag(max(s - lambda, 0)) %*% V^T
// where (U, s, V) come from svd(Y_complete). Used by fect_mc's per-
// iteration update.
//
// Inputs:
//   Y_complete : R numeric matrix (m x n), the matrix being
//                completed. Caller is responsible for filling treated
//                cells with their current estimate before the call.
//   lambda     : soft-threshold parameter (>= 0).
//
// Returns NULL on CUDA error. On success returns list(Y_hat, n_nonzero,
// status = 0).
//
// [[Rcpp::export]]
SEXP didgpu_cuda_fect_svd_softthreshold_r(Rcpp::NumericMatrix Y_complete,
                                            double lambda) {
#ifdef HAS_CUDA
  const int m = Y_complete.nrow();
  const int n = Y_complete.ncol();
  if (lambda < 0.0) Rcpp::stop("lambda must be non-negative; got %f", lambda);

  std::vector<double> Y_rm(static_cast<size_t>(m) * n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      Y_rm[static_cast<size_t>(i) * n + j] = Y_complete(i, j);

  double *d_Y = nullptr, *d_Yhat = nullptr;
  cudaError_t e;
  auto cleanup_and_return_null = [&]() -> SEXP {
    if (d_Y)    cudaFree(d_Y);
    if (d_Yhat) cudaFree(d_Yhat);
    return R_NilValue;
  };

  e = cudaMalloc((void**)&d_Y,    sizeof(double) * m * n);
  if (e != cudaSuccess) return cleanup_and_return_null();
  e = cudaMalloc((void**)&d_Yhat, sizeof(double) * m * n);
  if (e != cudaSuccess) return cleanup_and_return_null();
  e = cudaMemcpy(d_Y, Y_rm.data(), sizeof(double) * m * n,
                  cudaMemcpyHostToDevice);
  if (e != cudaSuccess) return cleanup_and_return_null();

  int n_nonzero = 0;
  int rc = didgpu_cuda_fect_svd_softthreshold(d_Y, m, n, lambda,
                                                d_Yhat, &n_nonzero);
  if (rc != 0) return cleanup_and_return_null();

  std::vector<double> Yhat_rm(static_cast<size_t>(m) * n);
  cudaMemcpy(Yhat_rm.data(), d_Yhat, sizeof(double) * m * n,
              cudaMemcpyDeviceToHost);
  cudaFree(d_Y); cudaFree(d_Yhat);

  Rcpp::NumericMatrix Y_hat(m, n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      Y_hat(i, j) = Yhat_rm[static_cast<size_t>(i) * n + j];

  Rcpp::List result;
  result["Y_hat"]     = Y_hat;
  result["n_nonzero"] = n_nonzero;
  result["status"]    = 0;
  return result;
#else
  (void)Y_complete; (void)lambda;
  return R_NilValue;
#endif
}


// Batched per-(g, t) CS inner regression on the GPU.
//
// Inputs (constructed by R-side .cs_inner_batched_cuda):
//   X_concat   : concatenated row-major design matrices for all cells.
//                Length = sum_c (n_c * p). For cells without
//                covariates pass a length-(sum_c * 1) intercept column
//                and set p = 1.
//   X_offsets  : length n_cells + 1; offsets[c] = starting row index
//                of cell c in the row-stacked layout. offsets[n_cells]
//                = total rows.
//   Y_concat   : concatenated delta values (length sum_c n_c).
//   W_concat   : concatenated treatment-indicator weights (D in
//                {0, 1}). Length sum_c n_c.
//   est_method : 0 = OR, 1 = IPW, 2 = DR
//   n_units    : panel-level unique-unit count (sets the row dimension
//                of out_influence when influence functions are requested).
//   want_influence : if TRUE, allocate and fill an n_units x n_cells IF
//                    matrix. If FALSE, that work is skipped (faster).
//
// Returns:
//   NULL if the kernel reports a non-zero status code (caller is
//   expected to fall back to the R per-cell loop). Otherwise a list
//   with components:
//     att        — numeric vector of length n_cells
//     influence  — n_units x n_cells matrix, or NULL if !want_influence
//     status     — 0 (success)
//
// This wrapper is the Phase-1 plumbing for task #79. The underlying
// kernel currently returns -1 (Phase 2), so in practice this function
// always returns NULL today and the R side runs the R fallback. Once
// Phase 2 (#82-#85) fills the kernel, no R-side or Rcpp changes are
// required — this seam stays stable.
//
// [[Rcpp::export]]
SEXP didgpu_cuda_cs_inner_batched_r(
    Rcpp::NumericVector X_concat,
    Rcpp::IntegerVector X_offsets,
    Rcpp::NumericVector Y_concat,
    Rcpp::NumericVector W_concat,
    Rcpp::IntegerVector unit_id_per_row,
    int p, int n_units, int est_method,
    bool want_influence) {
#ifdef HAS_CUDA
  const int n_cells = X_offsets.size() - 1;
  if (n_cells <= 0) Rcpp::stop("X_offsets must have at least 2 elements");
  if (p <= 0)       Rcpp::stop("p must be positive");

  const int n_total = X_offsets[n_cells];
  if (Y_concat.size() != n_total)
    Rcpp::stop("Y_concat length (%d) does not match X_offsets last element (%d)",
               (int)Y_concat.size(), n_total);
  if (W_concat.size() != n_total)
    Rcpp::stop("W_concat length (%d) does not match X_offsets last element (%d)",
               (int)W_concat.size(), n_total);
  if (X_concat.size() != n_total * p)
    Rcpp::stop("X_concat length (%d) does not match n_total * p (%d * %d = %d)",
               (int)X_concat.size(), n_total, p, n_total * p);
  if (unit_id_per_row.size() != n_total)
    Rcpp::stop("unit_id_per_row length (%d) does not match n_total (%d)",
               (int)unit_id_per_row.size(), n_total);

  // Phase 2 #84: only OR is implemented in the kernel. IPW / DR fall
  // back via the -3 return.
  if (est_method != 0) return R_NilValue;

  std::vector<double> att(n_cells, NA_REAL);
  std::vector<double> IF_per_row;
  double* IF_ptr = nullptr;
  if (want_influence) {
    IF_per_row.assign(n_total, 0.0);
    IF_ptr = IF_per_row.data();
  }

  int rc = didgpu_cuda_cs_inner_or(
      &X_concat[0], &X_offsets[0],
      &Y_concat[0],
      &W_concat[0],
      n_cells, p,
      att.data(), IF_ptr);

  if (rc != 0) return R_NilValue;

  Rcpp::NumericVector att_out(att.begin(), att.end());
  Rcpp::List result;
  result["att"]    = att_out;
  result["status"] = 0;
  if (want_influence) {
    // Scatter per-row IF into the (n_units x n_cells) layout. For
    // each cell c with rows [X_offsets[c], X_offsets[c+1]), each row
    // r maps to unit unit_id_per_row[r] (0-based). IF[unit, cell]
    // is set to IF_per_row[r] for that row's unit.
    Rcpp::NumericMatrix IF(n_units, n_cells);
    for (int c = 0; c < n_cells; ++c) {
      const int row_start = X_offsets[c];
      const int row_end   = X_offsets[c + 1];
      for (int r = row_start; r < row_end; ++r) {
        const int u = unit_id_per_row[r];
        if (u >= 0 && u < n_units) IF(u, c) = IF_per_row[r];
      }
    }
    result["influence"] = IF;
  } else {
    result["influence"] = R_NilValue;
  }
  return result;
#else
  (void)X_concat; (void)X_offsets; (void)Y_concat; (void)W_concat;
  (void)unit_id_per_row;
  (void)p; (void)n_units; (void)est_method; (void)want_influence;
  return R_NilValue;
#endif
}
