# ============================================================================
# CUDA glue layer: convert a prepped panel into the flat-array layout
# the CUDA kernel needs, then call the native function and return the
# DID estimate.
#
# This is the bridge between `.prep_panel()` / `.core_one_event_time()`
# (R-level) and `didgpu_cuda_did()` (Rcpp/CUDA in src/didgpu_init.cpp).
#
# Once CUDA is built (didgpu_has_cuda_support() is TRUE), the backend
# dispatch can swap `.core_one_event_time` for `.cuda_one_event_time`
# in the per-event-time loop. The result must match the r-backend
# bit-for-bit (within FP tolerance).
# ============================================================================


#' Run one (k, direction) DID on the GPU
#'
#' Drop-in alternative to `.core_one_event_time` that dispatches to the
#' CUDA kernel. Only callable when `didgpu_has_cuda_support()` is TRUE.
#'
#' @param prepped Output of `.prep_panel()`.
#' @param k Integer event-time.
#' @param direction 1 (in) or 0 (out).
#' @return A list with `att`, `N_inc`, `U_g`, `N_eff` (same shape as
#'   `.core_one_event_time` returns).
#' @keywords internal
#' @noRd
.cuda_one_event_time <- function(prepped, k, direction = 1L) {
  if (!isTRUE(didgpu_has_cuda_support())) {
    stop("CUDA backend not built. Install the NVIDIA CUDA Toolkit and ",
         "reinstall didgpu so nvcc compiles src/cuda_*.cu.")
  }
  k <- as.integer(k); direction <- as.integer(direction)

  d <- prepped
  data.table::setkeyv(d, c("group_XX", "time_XX"))
  n_rows <- nrow(d)

  # Build per-group quantities (one row per group).
  g <- unique(d[, list(group_XX, F_g_XX, S_g_XX, T_g_XX, L_g_XX)])
  data.table::setkeyv(g, "group_XX")
  n_groups <- nrow(g)

  # S_g_XX is NA for never-switchers; the kernel expects an int sentinel.
  # -1 = never-switcher (won't match direction 0 or 1).
  S_g <- as.integer(ifelse(is.na(g$S_g_XX), -1L, g$S_g_XX))

  # row_to_g maps 0-based group index per row. group_XX is already a
  # consecutive 1-based integer (prep_panel ensures this), so subtract 1.
  row_to_g <- as.integer(d$group_XX - 1L)
  row_to_t <- as.integer(d$time_XX)

  # Cohort key: flatten (time, d_sq) -> single 0-based int.
  dsq_levels <- sort(unique(d$d_sq_XX))   # sort() drops NA: NA is not a level
  n_dsq <- length(dsq_levels)
  dsq_idx <- match(d$d_sq_XX, dsq_levels) - 1L  # 0-based
  t_idx <- as.integer(d$time_XX - min(d$time_XX))  # 0-based
  n_cohorts <- as.integer((max(t_idx) + 1L) * n_dsq)
  cohort_key <- as.integer(t_idx * n_dsq + dsq_idx)

  # Groups unobserved at the global first period have d_sq_XX == NA (see
  # core_r.R, "Baseline (period-1) treatment"), so match() yields NA and
  # cohort_key would be NA_integer_ -- which reaches the CUDA kernels as
  # INT_MIN and causes an illegal memory access (CUDA error 700) the
  # moment k_finalize_dist_and_kernel reads N_t_control[cohort_key[r]].
  # On the CPU path these rows fall out of every cohort mask
  # (`d_sq_XX == l` is NA, never TRUE) and contribute exactly zero.
  # Reproduce that here: park them in a dedicated padding cohort. Their
  # kernel contribution is identically zero -- candidate_dist is 0 for
  # them (S_g is the -1 never-switcher sentinel, never == direction) and
  # the padding cohort's switcher mass is therefore 0, so ratio == 0 and
  # kernel_val == 0 -- but the reads are now in bounds.
  na_key <- is.na(cohort_key)
  if (any(na_key)) {
    cohort_key[na_key] <- n_cohorts
    n_cohorts <- n_cohorts + 1L
  }

  # First compute N_inc on the host (we need it for G_over_Ninc; the
  # CUDA kernel takes it as a scalar). N_inc = sum(N_gt * gated dist).
  # The "gated" version requires knowing N_t_control > 0, which we can
  # check up front with the same recipe.
  d2 <- data.table::copy(d)
  if (k == 1L) {
    d2[, diff_y_k_XX := diff_y_XX]
  } else {
    d2[, diff_y_k_XX := outcome_XX -
        data.table::shift(outcome_XX, k, type = "lag"),
       by = group_XX]
  }
  d2[, never_change_k_XX := as.integer(time_XX < F_g_XX &
                                        N_gt_XX > 0 &
                                        !is.na(diff_y_k_XX))]
  d2[is.na(never_change_k_XX), never_change_k_XX := 0L]
  d2[, N_t_control := sum(N_gt_XX * never_change_k_XX),
     by = c("time_XX", "d_sq_XX")]
  d2[, dist_k_XX := as.integer(
       time_XX == (F_g_XX + k - 1L) &
       k <= L_g_XX &
       !is.na(S_g_XX) & S_g_XX == direction &
       !is.na(diff_y_k_XX) &
       N_gt_XX > 0 &
       !is.na(N_t_control) & N_t_control > 0
     )]
  d2[is.na(dist_k_XX), dist_k_XX := 0L]
  N_inc <- sum(d2$N_gt_XX * d2$dist_k_XX, na.rm = TRUE)
  if (N_inc == 0) {
    return(list(att = NA_real_, N_inc = 0L, N_eff = 0L,
                U_g = numeric(n_groups)))
  }
  G_over_Ninc <- n_groups / N_inc

  did <- didgpu_cuda_did(
    outcome     = as.numeric(d$outcome_XX),
    N_gt        = as.numeric(d$N_gt_XX),
    row_to_g    = row_to_g,
    row_to_t    = row_to_t,
    cohort_key  = cohort_key,
    F_g         = as.integer(g$F_g_XX),
    S_g         = S_g,
    T_g         = as.integer(g$T_g_XX),
    L_g         = as.integer(g$L_g_XX),
    n_cohorts   = n_cohorts,
    k           = k,
    direction   = direction,
    G_over_Ninc = G_over_Ninc
  )

  # N_eff: we cheap-out and return N_inc as a placeholder. The CUDA
  # kernel doesn't yet expose the "contributing rows" count.
  list(att = did, N_inc = as.integer(N_inc),
       N_eff = as.integer(N_inc), U_g = numeric(n_groups))
}
