// ============================================================================
// Rcpp port of .core_one_event_time() for the binary, no-controls case.
//
// Input shape: the panel must be passed AS PREPPED columns (output of
// .prep_panel() in R) sorted by (group, time). The C++ kernel does NOT
// do any panel preparation — it expects:
//
//   outcome[n_rows]        : Y at each (group, time) cell (NA -> NaN)
//   N_gt[n_rows]           : weight (0 = row excluded; user weight or 1)
//   group_id[n_rows]       : 0-based group index 0..n_groups-1
//   time_id[n_rows]        : 0-based time index 0..n_times-1
//   cohort_id[n_rows]      : 0-based (time, d_sq) cohort index 0..n_cohorts-1
//   F_g[n_groups]          : per-group first-switch period (0-based)
//   S_g[n_groups]          : per-group switcher direction (1=in, 0=out, -1=NA)
//   T_g[n_groups]          : per-group last-usable period (0-based)
//   L_g[n_groups]          : per-group post-switch horizon length
//   group_offset[n_groups+1]: prefix-sum of rows per group; the rows of
//                             group g live at [group_offset[g], group_offset[g+1])
//
// Output: a list with att (double), N_inc (int), N_eff (int),
//   U_g (numeric of length n_groups).
//
// Matches the R-side .core_one_event_time() bit-for-bit when
// - controls=NULL, normalized=FALSE, same_switchers=FALSE,
//   only_never_switchers=FALSE (the simple-binary path).
//
// The kernel intentionally does NOT support the optional features above
// — the R orchestrator detects unsupported feature combinations and
// falls back to .core_one_event_time() (the R version).
// ============================================================================

#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

// Helper: NaN-safe check.
static inline bool is_na(double x) { return ISNAN(x); }

// [[Rcpp::export]]
List didgpu_cpp_core_one_event_time(
    NumericVector outcome,
    NumericVector N_gt,
    IntegerVector group_id,
    IntegerVector time_id,
    IntegerVector cohort_id,
    IntegerVector F_g,
    IntegerVector S_g,
    IntegerVector T_g,
    IntegerVector L_g,
    IntegerVector group_offset,
    int n_cohorts,
    int k,
    int direction) {

  const int n_rows   = outcome.size();
  const int n_groups = F_g.size();

  if (N_gt.size()       != n_rows)   stop("N_gt length mismatch");
  if (group_id.size()   != n_rows)   stop("group_id length mismatch");
  if (time_id.size()    != n_rows)   stop("time_id length mismatch");
  if (cohort_id.size()  != n_rows)   stop("cohort_id length mismatch");
  if (S_g.size()        != n_groups) stop("S_g length mismatch");
  if (T_g.size()        != n_groups) stop("T_g length mismatch");
  if (L_g.size()        != n_groups) stop("L_g length mismatch");
  if (group_offset.size() != n_groups + 1) stop("group_offset length mismatch");
  if (k < 1) stop("k must be >= 1");

  // -------- pass 1: per-row diff_y_k and never_change_k masks --------
  std::vector<double> diff_y_k(n_rows, NA_REAL);
  std::vector<int>    never_change_k(n_rows, 0);

  for (int g = 0; g < n_groups; ++g) {
    const int start = group_offset[g];
    const int end   = group_offset[g + 1];
    const int fg    = F_g[g];
    // Long difference: outcome[i] - outcome[i - k] (within group, row-wise).
    // Panel is sorted by (group, time); we use raw row offset = relative time.
    for (int i = start; i < end; ++i) {
      const int lag_i = i - k;
      if (lag_i >= start) {
        const double yi = outcome[i];
        const double yk = outcome[lag_i];
        if (!is_na(yi) && !is_na(yk)) {
          diff_y_k[i] = yi - yk;
        }
      }
      // never_change_k: time < F_g, N_gt > 0, valid diff.
      if (time_id[i] < fg && N_gt[i] > 0 && !is_na(diff_y_k[i])) {
        never_change_k[i] = 1;
      }
    }
  }

  // -------- pass 2: per-cohort control mass N_t_control --------
  std::vector<double> N_t_control(n_cohorts, 0.0);
  for (int i = 0; i < n_rows; ++i) {
    if (never_change_k[i] == 1) {
      N_t_control[cohort_id[i]] += N_gt[i];
    }
  }

  // -------- pass 3: per-row dist_k mask --------
  // dist_k[i] = 1 iff:
  //   time_id[i] == F_g - 1 + k        (the switcher's i-th post-period)
  //   k <= L_g                          (horizon reachable)
  //   S_g != -1, S_g == direction
  //   !is.na(diff_y_k[i])
  //   N_gt[i] > 0
  //   N_t_control[cohort] > 0
  std::vector<int> dist_k(n_rows, 0);
  for (int g = 0; g < n_groups; ++g) {
    const int start = group_offset[g];
    const int end   = group_offset[g + 1];
    const int fg = F_g[g];
    const int sg = S_g[g];
    const int Lg = L_g[g];
    if (sg < 0 || sg != direction || k > Lg) continue;
    const int target_t = fg - 1 + k;
    for (int i = start; i < end; ++i) {
      if (time_id[i] != target_t) continue;
      if (is_na(diff_y_k[i]))     continue;
      if (N_gt[i] <= 0)           continue;
      if (N_t_control[cohort_id[i]] <= 0) continue;
      dist_k[i] = 1;
    }
  }

  // -------- pass 4: N_inc, N_t_switch (per cohort), per-row contribution --------
  double N_inc = 0.0;
  std::vector<double> N_t_switch(n_cohorts, 0.0);
  for (int i = 0; i < n_rows; ++i) {
    if (dist_k[i] == 1) {
      N_inc += N_gt[i];
      N_t_switch[cohort_id[i]] += N_gt[i];
    }
  }

  if (N_inc == 0.0) {
    return List::create(
      _["att"]   = NA_REAL,
      _["N_inc"] = 0,
      _["N_eff"] = 0,
      _["U_g"]   = NumericVector(n_groups, 0.0));
  }

  // -------- pass 5: kernel + U_g --------
  // kernel[i] = (G / N_inc) * [time in (k+1)..T_g] * N_gt * (dist - ratio * never) * diff_y_k
  const double G_over_Ninc = static_cast<double>(n_groups) / N_inc;
  NumericVector U_g(n_groups, 0.0);
  int N_eff = 0;
  for (int g = 0; g < n_groups; ++g) {
    const int start = group_offset[g];
    const int end   = group_offset[g + 1];
    const int Tg    = T_g[g];
    double ug = 0.0;
    for (int i = start; i < end; ++i) {
      const int t = time_id[i];
      // Contribution window: time in (k+1)..T_g, inclusive.
      // R uses 1-based time; here time_id is 0-based but the mapping
      // preserved t = F_g + k - 1 boundary above. The contribution
      // window in 0-based indexing is t >= k (i.e., t in [k, Tg]).
      if (t < k || t > Tg) continue;
      // ratio at this row's cohort.
      const int c = cohort_id[i];
      const double nt_c = N_t_control[c];
      const double ratio = (nt_c > 0.0) ? (N_t_switch[c] / nt_c) : 0.0;
      const double dy = diff_y_k[i];
      if (is_na(dy)) continue;
      const double kernel =
        G_over_Ninc * N_gt[i] *
        (static_cast<double>(dist_k[i]) - ratio * static_cast<double>(never_change_k[i])) *
        dy;
      ug += kernel;
      // N_eff: contributing rows. A row contributes if it's a switcher
      // OR a never-change row with at least one cohort-mate switcher.
      const bool is_switcher = (dist_k[i] == 1);
      const bool is_active_control =
        (never_change_k[i] == 1) && (N_t_switch[c] > 0.0);
      if (N_gt[i] > 0.0 && (is_switcher || is_active_control)) {
        N_eff += 1;
      }
    }
    U_g[g] = ug;
  }

  // -------- final: att = sum(U_g) / G --------
  double sum_U = 0.0;
  for (int g = 0; g < n_groups; ++g) sum_U += U_g[g];
  const double att = sum_U / static_cast<double>(n_groups);

  return List::create(
    _["att"]   = att,
    _["N_inc"] = static_cast<int>(N_inc),
    _["N_eff"] = static_cast<int>(N_eff),
    _["U_g"]   = U_g);
}
