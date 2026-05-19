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
