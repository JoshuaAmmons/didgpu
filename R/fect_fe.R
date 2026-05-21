# ============================================================================
# fect_fe: counterfactual ATT via two-way fixed effects.
#
# Algorithm (Liu, Wang & Xu 2024, simplified):
#   1. Build the (n_units x n_periods) outcome matrix Y and treatment
#      mask M (1 = treated, 0 = control).
#   2. Fit alpha_i (unit FE) + xi_t (time FE) on CONTROL cells only,
#      by iterative demeaning:
#         repeat:
#           Y_c = Y with M = 1 cells masked out
#           alpha_i = mean(Y_c[i, ] - xi_t,  na.rm = TRUE)
#           xi_t    = mean(Y_c[, t] - alpha_i, na.rm = TRUE)
#         until alpha + xi changes < tol.
#   3. Predict Y_hat[i, t] = alpha_i + xi_t for treated cells.
#   4. ATT = mean of (Y[i, t] - Y_hat[i, t]) over treated cells.
#   5. Per-event-time effects: ATT split by (t - F_g).
#
# This is the simplest fect estimator. The two harder ones — ife
# (interactive fixed effects, Bai 2009) and mc (matrix completion,
# Athey et al. 2021) — both reduce to the same primitive (iterative
# soft-thresholded SVD on the controls matrix), which is what the GPU
# path is for: dense SVD on the GPU via cuSOLVER's batched SVDJ.
#
# The CUDA kernel for the demeaning step is in src/cuda_fect_fe.cu
# (scaffolded; compiles when nvcc + headers are present).
# ============================================================================


# Reshape a long panel (one row per (group, time)) into wide matrices
# Y and M used by the fect kernel.
#' @keywords internal
#' @noRd
.fect_build_matrices <- function(df, outcome, group, time, treatment) {
  d <- data.table::as.data.table(df)
  d <- d[!is.na(get(group)) & !is.na(get(time)), ]
  data.table::setorderv(d, c(group, time))
  units   <- sort(unique(d[[group]]))
  periods <- sort(unique(d[[time]]))
  n_units   <- length(units)
  n_periods <- length(periods)
  Y <- matrix(NA_real_, n_units, n_periods,
              dimnames = list(as.character(units), as.character(periods)))
  M <- matrix(0L, n_units, n_periods,
              dimnames = list(as.character(units), as.character(periods)))
  # Vectorised fill: build row / col indices once.
  ri <- match(d[[group]], units)
  ci <- match(d[[time]],  periods)
  idx <- ri + (ci - 1L) * n_units
  Y[idx] <- as.numeric(d[[outcome]])
  M[idx] <- as.integer(d[[treatment]] > 0)
  M[is.na(M)] <- 0L
  list(Y = Y, M = M, units = units, periods = periods)
}


# Iterative two-way demeaning on the controls-only matrix.
# Returns (alpha, xi) such that Y_hat[i, t] = alpha[i] + xi[t].
#' @keywords internal
#' @noRd
.fect_fe_fit <- function(Y, M, tol = 1e-5, max_iter = 500L) {
  n_units   <- nrow(Y)
  n_periods <- ncol(Y)
  # Mask treated cells with NA so they're ignored.
  Y_c <- Y
  Y_c[M == 1L] <- NA_real_

  alpha <- numeric(n_units)
  xi    <- numeric(n_periods)
  prev_loss <- Inf

  for (iter in seq_len(max_iter)) {
    # Update alpha: row mean of (Y_c - xi) over non-NA cells.
    # Vectorised: subtract xi column-wise, then rowMeans with na.rm.
    R <- Y_c - matrix(xi, n_units, n_periods, byrow = TRUE)
    alpha_new <- rowMeans(R, na.rm = TRUE)
    # If a unit has 0 controls (all post-treatment), rowMeans returns NaN.
    alpha_new[is.nan(alpha_new)] <- 0

    # Update xi: column mean of (Y_c - alpha).
    R <- Y_c - matrix(alpha_new, n_units, n_periods, byrow = FALSE)
    xi_new <- colMeans(R, na.rm = TRUE)
    xi_new[is.nan(xi_new)] <- 0

    # Convergence: max abs change in alpha + xi.
    delta <- max(max(abs(alpha_new - alpha)), max(abs(xi_new - xi)))
    alpha <- alpha_new
    xi    <- xi_new

    # Numerical loss for diagnostic: sum of squared control residuals.
    resid <- Y_c - alpha - matrix(xi, n_units, n_periods, byrow = TRUE)
    loss <- sum(resid^2, na.rm = TRUE)

    if (delta < tol) break
    if (iter > 1L && abs(loss - prev_loss) < tol) break
    prev_loss <- loss
  }

  attr(alpha, "iter")  <- iter
  attr(alpha, "delta") <- delta
  list(alpha = alpha, xi = xi, iter = iter, delta = delta)
}


# Compute per-event-time ATT from a fitted (alpha, xi) and the panel.
# Returns a list with overall ATT, per-event-time ATT, sample sizes.
#' @keywords internal
#' @noRd
.fect_compute_att <- function(Y, M, fit, group, time, units, periods,
                                effects = NULL) {
  alpha <- fit$alpha
  xi    <- fit$xi
  Y_hat <- alpha + matrix(xi, nrow(Y), ncol(Y), byrow = TRUE)
  residual <- Y - Y_hat   # = Y - Y(0)_hat, this IS the treatment effect estimate
  # Overall ATT: mean of residual over treated cells.
  treated_mask <- (M == 1L) & !is.na(Y) & !is.na(Y_hat)
  ate <- mean(residual[treated_mask], na.rm = TRUE)
  n_treated <- sum(treated_mask)

  # Per-event-time decomposition. For each treated cell, the event-time
  # is (period - F_g) where F_g is the first treated period for that unit.
  # F_g_i = min t such that M[i, t] = 1; NA if never treated.
  F_g <- apply(M, 1L, function(row) {
    idx <- which(row == 1L)
    if (length(idx) == 0L) NA_integer_ else min(idx)
  })

  # Map (i, t) to event-time. For each treated cell, event_time = t - F_g_i.
  ev_results <- NULL
  if (!is.null(effects) && effects > 0L) {
    eff <- numeric(effects)
    n_eff <- integer(effects)
    n_switchers <- integer(effects)
    for (k in seq_len(effects)) {
      # Cells at event-time k - 1 (i.e., k periods after F_g, zero-indexed).
      # We use 1-indexed convention: k = 1 means the first post-treatment period,
      # which is t = F_g.
      cells <- treated_mask & FALSE
      switchers_this_k <- 0L
      for (i in seq_len(nrow(Y))) {
        if (is.na(F_g[i])) next
        t_idx <- F_g[i] + k - 1L
        if (t_idx >= 1L && t_idx <= ncol(Y) &&
            treated_mask[i, t_idx]) {
          cells[i, t_idx] <- TRUE
          switchers_this_k <- switchers_this_k + 1L
        }
      }
      n <- sum(cells)
      n_eff[k] <- as.integer(n)
      n_switchers[k] <- switchers_this_k
      if (n > 0L) {
        eff[k] <- mean(residual[cells], na.rm = TRUE)
      } else {
        eff[k] <- NA_real_
      }
    }
    ev_results <- list(effects = eff,
                       n_eff = n_eff,
                       n_switchers = n_switchers)
  }

  list(ate = ate,
       n_treated = n_treated,
       Y_hat = Y_hat,
       residual = residual,
       per_event = ev_results,
       F_g = F_g)
}


# R-side implementation of one fect_fe fit (one bootstrap iter).
# Returns a list with the same shape as a didgpu cell: ate, effects,
# n_eff_effects, etc.
#' @keywords internal
#' @noRd
.fect_fe_one_iter <- function(df, args, iter_seed) {
  df_use <- if (iter_seed == 0L) df else .cluster_resample(df, args, iter_seed)
  mats <- .fect_build_matrices(df_use, args$outcome, args$group,
                                 args$time, args$treatment)
  t0 <- Sys.time()
  # Backend dispatch: CUDA if available + requested AND the matrix is
  # large enough to amortise GPU overhead. The fect_fe kernel does a
  # per-iteration D2H copy for its convergence check, so it has the
  # same small-matrix pathology as the SVD path (see BENCHMARKS.md);
  # the size gate keeps backend = "cuda" from ever being slower than
  # the R demeaning loop.
  use_cuda <- identical(args$backend, "cuda") &&
              isTRUE(tryCatch(didgpu_has_cuda_support(),
                               error = function(e) FALSE)) &&
              .fect_cuda_svd_worthwhile(nrow(mats$Y), ncol(mats$Y))
  fit <- if (use_cuda) {
    cuda_res <- didgpu_cuda_fect_fe_r(mats$Y, mats$M,
                                        tol = args$tol %||% 1e-5,
                                        max_iter = args$max_iter %||% 500L)
    list(alpha = cuda_res$alpha,
         xi = cuda_res$xi,
         iter = cuda_res$iter,
         delta = cuda_res$delta)
  } else {
    .fect_fe_fit(mats$Y, mats$M,
                  tol = args$tol %||% 1e-5,
                  max_iter = args$max_iter %||% 500L)
  }
  res <- .fect_compute_att(mats$Y, mats$M, fit,
                             args$group, args$time, mats$units, mats$periods,
                             effects = args$effects)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  n_eff_effects <- if (is.null(res$per_event)) integer(0)
                    else res$per_event$n_eff
  n_inc_effects <- if (is.null(res$per_event)) integer(0)
                    else res$per_event$n_switchers
  effects_vec <- if (is.null(res$per_event)) numeric(0)
                  else res$per_event$effects

  list(
    effects        = effects_vec,
    ate            = res$ate,
    placebos       = numeric(0),
    n_effects      = length(effects_vec),
    n_placebos     = 0L,
    n_inc_effects  = n_inc_effects,
    n_inc_placebos = integer(0),
    n_eff_effects  = n_eff_effects,
    n_eff_placebos = integer(0),
    iter_seed      = as.integer(iter_seed),
    wall_seconds   = wall,
    backend        = if (use_cuda) "cuda" else "r",
    fect_method    = "fe",
    fect_iter      = fit$iter,
    fect_delta     = fit$delta
  )
}
