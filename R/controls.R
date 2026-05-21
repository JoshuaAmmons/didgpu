# ============================================================================
# Controls support — WORK IN PROGRESS, not wired into the backend.
#
# This file holds the pre-fit logic that computes the per-baseline-level
# regression coefficients (theta_d) from the FWL machinery the reference
# uses. To complete:
#
# 1. (Done in this file) `.prefit_controls(prepped, controls)` returns
#    theta_d[l, k] from the OLS of sqrt(N_gt)*diff_y on
#    sqrt(N_gt)*(diff_X - avg_diff_X[time, d_sq]) restricted to
#    never-switcher rows of baseline level l.
#
# 2. (Not done) `.core_one_event_time` must accept the prefit and, for
#    each event-time k, build `diff_X_k = X - shift(X, k)` per control,
#    then adjust `diff_y_k := diff_y_k - sum_c theta_d[l, c] * diff_X_c_k`
#    for rows with d_sq == l. Reference: core.R:267-269. The current
#    `.apply_controls_adjustment` only handles k = 1 (it modifies
#    diff_y_XX in place, which only flows through to the k=1 case).
#
# 3. (Not done) Add `controls` arg to `didgpu()`, plumb through
#    `args$controls` to `.backend_r_impl`, call `.prefit_controls` after
#    `.prep_panel`, pass the prefit to the per-k kernel calls.
#
# 4. (Not done) Variance correction `part2_switch` (core.R:504-533).
#    Affects SE only; bootstrap SEs work without it.
#
# Algorithm (verified against did_multiplegt_dyn_core.R:225-273 and
# did_multiplegt_main.R:348-491):
#
#   For each control variable c in 1..K:
#     diff_X_c    = X_c - lag(X_c)         by group
#     diff_X_c_N  = sqrt(N_gt) * diff_X_c
#
#   For each (time, d_sq):
#     avg_diff_X_c[time, d_sq] = mean(diff_X_c restricted to
#                                       ever_change_d == 0 & non-missing)
#     resid_X_c = sqrt(N_gt) * (diff_X_c - avg_diff_X_c[time, d_sq])
#
#   For each baseline level l:
#     Restrict to ever_change_d == 0 & d_sq == l & non-missing rows.
#     X_mat = [resid_X_1, ..., resid_X_K]
#     y     = sqrt(N_gt) * diff_y
#     theta_d[l, ] = ginv(X_mat' X_mat) X_mat' y
#
#   Inside the per-event-time core, BEFORE building diff_y_k:
#     For each row with d_sq == l and each control c, compute
#       diff_X_c_k = X_c - shift(X_c, k) by group
#     and subtract sum_c theta_d[l, c] * diff_X_c_k from diff_y_k.
#
# Variance correction (part2_switch) is deferred — bootstrap SEs work fine
# without it.
# ============================================================================


#' Compute pre-fit controls coefficients per baseline level
#'
#' @param prepped Output of `.prep_panel()`.
#' @param controls Character vector of column names in prepped.
#' @return A list with `theta` (named list keyed by baseline level,
#'   each a numeric vector of length K) and `diff_X` (a list of
#'   the per-row first-differenced control columns, length K, each
#'   length n_rows; the kernel uses these to adjust diff_y).
#' @keywords internal
#' @noRd
.prefit_controls <- function(prepped, controls) {
  d <- prepped
  K <- length(controls)
  if (K == 0L) return(list(theta = list(), diff_X = list(), controls = controls))

  # ever_change_d_XX is PER-(group, time): 1 if this row is at or after
  # the group's first switch, 0 before. So a switcher unit has
  # ever_change=0 for its pre-switch rows AND those rows count as
  # "controls" for the regression below. Reference: main.R:193-196.
  d[, ever_change_d_XX := as.integer(time_XX >= F_g_XX)]
  d[is.na(F_g_XX), ever_change_d_XX := 0L]

  # First-difference each control by group.
  diff_X <- vector("list", K)
  for (k in seq_len(K)) {
    cname <- controls[k]
    if (!cname %in% names(d))
      stop("control variable not in prepped panel: ", cname)
    d[, paste0("diff_X_", k, "_XX") := get(cname) -
        data.table::shift(get(cname), 1L, type = "lag"),
      by = group_XX]
    diff_X[[k]] <- d[[paste0("diff_X_", k, "_XX")]]
  }

  # `fd_X_all_non_missing_XX`: 1 iff every control's diff_X is non-NA.
  # The reference uses this as part of the row mask for the control
  # regression and the per-(time, d_sq) means (main.R:350-355, 367-380).
  d[, fd_X_all_non_missing_XX := 1L]
  for (k in seq_len(K)) {
    d[is.na(get(paste0("diff_X_", k, "_XX"))), fd_X_all_non_missing_XX := 0L]
  }

  # Per (time, d_sq[, trends_nonparam]) mean of diff_X_k among control
  # units. Mask matches the reference's mask: ever_change == 0 AND
  # non-NA diff_y AND all control diffs present.
  cohort_cols <- c("time_XX", "d_sq_XX")
  if ("trends_np_XX" %in% names(d)) cohort_cols <- c(cohort_cols, "trends_np_XX")
  resid_X <- vector("list", K)
  mask_ctrl <- d$ever_change_d_XX == 0L &
               !is.na(d$diff_y_XX) &
               d$fd_X_all_non_missing_XX == 1L &
               d$N_gt_XX > 0
  for (k in seq_len(K)) {
    diff_col   <- paste0("diff_X_", k, "_XX")
    avg_col    <- paste0("avg_diff_X_", k, "_XX")
    resid_col  <- paste0("resid_X_", k, "_XX")

    d[, "._num_XX" := ifelse(mask_ctrl, N_gt_XX * get(diff_col), 0)]
    d[, "._den_XX" := ifelse(mask_ctrl, N_gt_XX, 0)]
    d[, "._sumnum_XX" := sum(._num_XX, na.rm = TRUE),
      by = cohort_cols]
    d[, "._sumden_XX" := sum(._den_XX, na.rm = TRUE),
      by = cohort_cols]
    d[, (avg_col) := ifelse(._sumden_XX > 0, ._sumnum_XX / ._sumden_XX, NA_real_)]
    d[, (resid_col) := sqrt(N_gt_XX) * (get(diff_col) - get(avg_col))]
    d[is.na(get(resid_col)), (resid_col) := 0]
    d[, c("._num_XX", "._den_XX", "._sumnum_XX", "._sumden_XX") := NULL]

    resid_X[[k]] <- d[[resid_col]]
  }

  # Per baseline level: OLS of diff_y_w on resid_X among controls.
  diff_y_w <- sqrt(d$N_gt_XX) * d$diff_y_XX
  diff_y_w[is.na(diff_y_w)] <- 0

  dsq_levels <- sort(unique(d$d_sq_XX))
  theta <- list()
  for (l in dsq_levels) {
    mask <- d$ever_change_d_XX == 0L &
            d$d_sq_XX == l &
            !is.na(d$diff_y_XX) &
            d$fd_X_all_non_missing_XX == 1L &
            d$N_gt_XX > 0
    n_obs <- sum(mask)
    if (n_obs <= K + 1L) {
      # Not enough rows; theta_l = 0 (no adjustment).
      theta[[as.character(l)]] <- rep(0, K)
      next
    }
    X_mat <- do.call(cbind, lapply(seq_len(K),
                                    function(k) d[[paste0("resid_X_", k, "_XX")]][mask]))
    y_vec <- diff_y_w[mask]
    XtX <- crossprod(X_mat)
    XtY <- crossprod(X_mat, y_vec)
    theta_l <- tryCatch(
      as.numeric(MASS::ginv(XtX) %*% XtY),
      error = function(e) rep(0, K)
    )
    theta[[as.character(l)]] <- theta_l
  }

  list(theta = theta, diff_X = diff_X, controls = controls,
       dsq_levels = dsq_levels)
}


#' Apply the controls FWL adjustment to diff_y_XX in-place
#'
#' Modifies `d$diff_y_XX` by subtracting `sum_k theta_d[l, k] * diff_X_k`
#' for each row with `d_sq == l`. The kernel then uses the adjusted
#' diff_y as if no controls existed.
#'
#' @keywords internal
#' @noRd
.apply_controls_adjustment <- function(d, prefit) {
  K <- length(prefit$controls)
  if (K == 0L) return(d)
  d[, diff_y_XX_orig := diff_y_XX]
  for (l in prefit$dsq_levels) {
    key <- as.character(l)
    theta_l <- prefit$theta[[key]]
    if (is.null(theta_l) || all(theta_l == 0)) next
    mask <- d$d_sq_XX == l
    adj <- rep(0, nrow(d))
    for (k in seq_len(K)) {
      dx <- d[[paste0("diff_X_", k, "_XX")]]
      adj <- adj + theta_l[k] * ifelse(is.na(dx), 0, dx)
    }
    d[mask, diff_y_XX := diff_y_XX - adj[mask]]
  }
  d
}
