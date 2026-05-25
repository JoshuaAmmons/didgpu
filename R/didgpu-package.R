#' didgpu: GPU-Accelerated, Checkpointed Difference-in-Differences
#'
#' A clean-room reimplementation of the de Chaisemartin and D'Haultfoeuille
#' dynamic difference-in-differences estimator, designed for long-running
#' econometric work. Three architectural commitments distinguish it from
#' the reference 'DIDmultiplegtDYN':
#'
#' 1. **Per-cell checkpointing.** Each bootstrap iteration is saved to disk
#'    as it completes. A crash at hour N preserves N-worth of work.
#' 2. **Resumable runs.** Re-invoking with the same `checkpoint_dir` skips
#'    completed cells. Identical seed plus identical config produces
#'    identical aggregated output.
#' 3. **Pluggable backend.** Pure R for correctness, Rcpp+Eigen for CPU
#'    throughput, optional CUDA (cuBLAS + cuSOLVER) for the bootstrap
#'    inner loop. Same numerics across backends.
#'
#' @section Status: This package is in early development. Currently
#'   supported: binary, non-absorbing treatment; the core estimator
#'   in pure R; checkpointing layer; CPU backend. Not yet supported:
#'   multivalued/continuous treatment, controls, weights, the
#'   `by_path` framework, the GPU backend (scaffolded only).
#'
#' @docType package
#' @name didgpu-package
#' @aliases didgpu-package
#' @useDynLib didgpu, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom data.table ":=" .SD .N as.data.table CJ copy data.table
#'   setDT setkeyv setnames shift
#' @importFrom stats cov pchisq pnorm qnorm rnorm sd
#' @importFrom graphics abline segments
#' @importFrom utils packageVersion read.csv
#' @importFrom tools md5sum
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom MASS ginv
"_PACKAGE"

# Suppress R CMD check NOTE about data.table NSE columns we use internally.
# Extend as more columns are introduced.
utils::globalVariables(c(
  # data.table NSE generics
  ".SD", ".N", ".I", "..keep", "..mycontrols_XX",
  # high-level columns used in simulators / utilities
  "G", "T_", "D", "Y", "F_g", "k", "b", "cell_id",
  "unit", "period", "k_evt", "tau_k",
  # prep-panel columns (the "_XX" suffix convention from the reference)
  "outcome_XX", "treatment_XX", "group_XX", "time_XX", "N_gt_XX",
  "d_sq_XX", "F_g_XX", "S_g_XX", "T_g_XX", "L_g_XX", "L_g_placebo_XX",
  "avg_post", "never_change_d_XX", "diff_y_XX", "weight_XX_input",
  "trends_np_XX",
  # per-event-time / per-placebo scratch columns built inside the core
  "diff_y_k_XX", "dist_k_XX", "never_change_k_XX",
  "N_t_control", "N_t_switch", "ratio_XX", "kernel_XX", "contrib_mask_XX",
  "diff_y_pl_k_XX", "dist_k_pl_XX", "never_change_k_pl_XX",
  "N_t_control_pl", "N_t_switch_pl", "ratio_pl_XX", "kernel_pl_XX",
  "contrib_pl_mask_XX",
  # normalized = TRUE scratch (cumulative treatment-change magnitude)
  "sum_temp_XX", "sum_treat_until_XX", "delta_cum_XX",
  "sum_temp_pl_XX", "sum_treat_until_pl_XX", "delta_cum_pl_XX",
  # predict_het scratch
  "Yg_Fg_min1_XX", "feasible_het_XX", "Yg_Fg_i_XX",
  "prod_het_XX", "gr_id_XX", "weight_XX",
  # controls.R scratch
  "ever_change_d_XX", "diff_y_XX_orig", "fd_X_all_non_missing_XX",
  "._num_XX", "._den_XX", "._sumnum_XX", "._sumden_XX",
  ".__diff_X_k_XX", ".__diff_X_pl_k_XX",
  # dont_drop_larger_lower scratch in .prep_panel
  "d_sq_tmp", "diff_from_sq_tmp",
  "ever_strict_increase_tmp", "ever_strict_decrease_tmp",
  # continuous scratch in .prep_panel
  "d_sq_XX_orig", "treatment_XX_orig", "S_g_het_XX",
  # same_switchers pre-pass scratch
  "still_switcher_XX", "N_g_control_check_XX",
  "diff_y_last_XX", "never_change_d_last_XX",
  "N_gt_control_last_XX", "N_g_control_last_m_XX", "diff_y_relev_XX",
  # same_switchers_pl pre-pass scratch
  "still_switcher_pl_XX", "N_g_control_check_pl_XX",
  "diff_y_last_pl_XX", "never_change_d_last_pl_XX",
  "N_gt_control_last_pl_XX", "N_g_control_last_m_pl_XX",
  "diff_y_relev_pl_XX",
  # by_path scratch
  "group_path_XX", "time_path_XX", "treatment_path_XX",
  "group_path_int_XX", "time_path_int_XX",
  "d_sq_path_XX", "F_g_path_XX", "path",
  # cpu-backend scratch
  ".g0", ".t0", ".c0", "N",
  # Callaway-Sant'Anna scratch
  "Y_XX", "G_XX", "T_XX", "D_XX", "Y_pre", "Y_t", "delta_XX",
  # summarize.R scratch
  "Y_", "G_", "T_", "D_", "d_sq", "F_g", "d_sq2",
  "avg_post", "direction",
  # did_static.R (.didm_point) scratch
  "D_lag", "Y_lag", "dY", "trans",
  # cs_continuous.R scratch
  "dy", "dose",
  # twfe.R scratch (data.table update-join)
  "i.vv"
))
