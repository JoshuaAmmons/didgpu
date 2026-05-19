# ============================================================================
# fect_ife: interactive fixed effects (Bai 2009).
#
# Algorithm:
#   Y[i, t] = alpha[i] + xi[t] + sum_{k=1..r} lambda[i, k] * F[k, t] + eps
#
#   Fit on control cells only (M[i, t] = 0); predict Y(0) for treated.
#
# Iteration (alternating):
#   1. fe step: alpha_i = rowmean(Y_c - xi - L * F)
#               xi_t    = colmean(Y_c - alpha - L * F)
#   2. SVD step: R = Y_c - alpha - xi
#                R = U D V^T   (truncated SVD, keep top r)
#                L_new = U_r * sqrt(D_r)         (n_units x r)
#                F_new = sqrt(D_r) * V_r^T       (r x n_periods)
#   3. Convergence: max |delta(alpha) U delta(xi) U delta(L*F)| < tol
#
# Predict: Y_hat[i, t] = alpha[i] + xi[t] + sum_k L[i, k] * F[k, t]
# ATT: mean over treated cells of (Y - Y_hat)
#
# The CUDA path lives in src/cuda_fect_svd.cu (cuSOLVER cusolverDnDgesvdj
# for the SVD step). The fe step reuses src/cuda_fect_fe.cu's
# row/col-mean kernels.
# ============================================================================


# Rank-r truncated SVD of a matrix M, treating NaN as 0 (for the
# residual matrix on control cells, we mask out treated cells with NaN
# beforehand and then convert to 0 for the SVD — same as fect's R
# implementation).
#' @keywords internal
#' @noRd
.fect_svd_r <- function(M, r) {
  # Replace NaN with 0 (control-cell positions where Y wasn't observed).
  M0 <- M
  M0[is.na(M0)] <- 0
  s <- svd(M0, nu = r, nv = r)
  list(L = s$u %*% diag(sqrt(s$d[seq_len(r)]), r, r),
       F = diag(sqrt(s$d[seq_len(r)]), r, r) %*% t(s$v),
       d = s$d[seq_len(r)])
}


# CUDA wrapper for the truncated-SVD step. Phase-1 wiring for task
# #80: dispatches to the cuSOLVER kernel and returns the same shape
# as .fect_svd_r (L, F). If CUDA is unavailable, the kernel returns
# nonzero, or anything else goes wrong, returns NULL so the caller
# falls back to .fect_svd_r.
#' @keywords internal
#' @noRd
.fect_svd_truncated_cuda <- function(M, r) {
  if (!isTRUE(tryCatch(didgpu_has_cuda_support(),
                       error = function(e) FALSE))) return(NULL)
  M0 <- M
  M0[is.na(M0)] <- 0
  result <- tryCatch(
    didgpu_cuda_fect_svd_truncated_r(M = M0, r = as.integer(r)),
    error = function(e) NULL)
  if (is.null(result)) return(NULL)
  if (is.null(result$L) || is.null(result$F)) return(NULL)
  # Pad with NA for d (cuSOLVER doesn't return S separately in this
  # kernel — the R fallback computes it from svd()). Phase 2 #86 will
  # extend the kernel to also return s.
  list(L = result$L, F = result$F, d = rep(NA_real_, r))
}


# Size gate for the fect CUDA SVD path.
#
# CRITICAL PERFORMANCE FINDING (see BENCHMARKS.md): for the small,
# tall-skinny matrices typical of fect panels (n_units in the
# hundreds, n_periods in the tens), the per-iteration cuSOLVER SVD is
# 100-300x SLOWER than base R's LAPACK svd(). cuSOLVER handle
# creation plus H2D/D2H transfer (~2-5 ms per call) dwarfs the SVD
# itself, and fect_ife / fect_mc call the SVD once per alternation
# iteration (dozens to hundreds of times).
#
# GPU SVD only starts to win when the matrix is large enough to
# amortise that fixed overhead — empirically n_units in the
# thousands AND a six-figure element count. Below that threshold we
# silently use the R svd() even when backend = "cuda", so that
# requesting the GPU backend NEVER makes fect slower than the CPU
# path. This mirrors how a good BLAS auto-dispatches small GEMMs to
# the CPU.
#
# The threshold is deliberately conservative: false negatives (using
# the CPU when the GPU would have marginally won) cost little, but
# false positives (using the GPU on a small matrix) cost 100x.
#' @keywords internal
#' @noRd
.fect_cuda_svd_worthwhile <- function(n_units, n_periods) {
  nu <- as.numeric(n_units)
  np <- as.numeric(n_periods)
  (nu >= 2000) && (nu * np >= 2e5)
}


# IFE fit: alternating fe-step + svd-step until convergence.
# Returns alpha, xi, L (n_units x r), F (r x n_periods), iter, delta.
#
# `use_cuda_svd`: when TRUE, every iteration's truncated SVD is
# computed on the GPU via .fect_svd_truncated_cuda. If a single CUDA
# call fails, the function silently falls back to .fect_svd_r for the
# remainder of the fit (the alternation may have started — switching
# back mid-loop is safe because both kernels produce the same shape).
#' @keywords internal
#' @noRd
.fect_ife_fit <- function(Y, M, r = 2L, tol = 1e-5, max_iter = 500L,
                          use_cuda_svd = FALSE) {
  n_units   <- nrow(Y)
  n_periods <- ncol(Y)
  Y_c <- Y
  Y_c[M == 1L] <- NA_real_

  alpha <- numeric(n_units)
  xi    <- numeric(n_periods)
  L     <- matrix(0, n_units,   r)
  F     <- matrix(0, r,         n_periods)
  prev_Y_hat <- matrix(0, n_units, n_periods)

  for (iter in seq_len(max_iter)) {
    # fe step on the residual (Y_c minus current factor structure).
    fac <- L %*% F
    R_for_fe <- Y_c - fac
    fe_fit <- .fect_fe_fit(R_for_fe, M, tol = tol, max_iter = max_iter)
    alpha <- fe_fit$alpha
    xi    <- fe_fit$xi

    # svd step on the residual after fe.
    R <- Y_c - alpha - matrix(xi, n_units, n_periods, byrow = TRUE)
    svd_res <- if (use_cuda_svd) {
      cuda_res <- .fect_svd_truncated_cuda(R, r)
      if (is.null(cuda_res)) {
        # CUDA failed once; stop trying for the rest of this fit.
        use_cuda_svd <- FALSE
        .fect_svd_r(R, r)
      } else cuda_res
    } else .fect_svd_r(R, r)
    L <- svd_res$L
    F <- svd_res$F

    # Convergence: max change in predicted Y over all cells.
    Y_hat <- alpha + matrix(xi, n_units, n_periods, byrow = TRUE) + L %*% F
    delta <- max(abs(Y_hat - prev_Y_hat))
    prev_Y_hat <- Y_hat

    if (delta < tol) break
  }

  list(alpha = alpha, xi = xi, L = L, F = F,
       iter = iter, delta = delta)
}


# Compute ATT + per-event-time effects from a fitted IFE model.
#' @keywords internal
#' @noRd
.fect_ife_compute_att <- function(Y, M, fit, effects = NULL) {
  n_units   <- nrow(Y)
  n_periods <- ncol(Y)
  Y_hat <- fit$alpha +
           matrix(fit$xi, n_units, n_periods, byrow = TRUE) +
           fit$L %*% fit$F
  residual <- Y - Y_hat
  treated_mask <- (M == 1L) & !is.na(Y) & !is.na(Y_hat)
  ate <- mean(residual[treated_mask], na.rm = TRUE)
  n_treated <- sum(treated_mask)

  F_g <- apply(M, 1L, function(row) {
    idx <- which(row == 1L)
    if (length(idx) == 0L) NA_integer_ else min(idx)
  })

  ev_results <- NULL
  if (!is.null(effects) && effects > 0L) {
    eff <- numeric(effects)
    n_eff <- integer(effects)
    n_switchers <- integer(effects)
    for (k in seq_len(effects)) {
      cells <- treated_mask & FALSE
      sw_k <- 0L
      for (i in seq_len(n_units)) {
        if (is.na(F_g[i])) next
        t_idx <- F_g[i] + k - 1L
        if (t_idx >= 1L && t_idx <= n_periods &&
            treated_mask[i, t_idx]) {
          cells[i, t_idx] <- TRUE
          sw_k <- sw_k + 1L
        }
      }
      n_eff[k] <- as.integer(sum(cells))
      n_switchers[k] <- sw_k
      eff[k] <- if (sum(cells) > 0L) mean(residual[cells], na.rm = TRUE)
                else NA_real_
    }
    ev_results <- list(effects = eff,
                       n_eff = n_eff,
                       n_switchers = n_switchers)
  }

  list(ate = ate, n_treated = n_treated,
       Y_hat = Y_hat, residual = residual,
       per_event = ev_results, F_g = F_g)
}


# Per-iter ife fit wrapping build/fit/compute_att.
#' @keywords internal
#' @noRd
.fect_ife_one_iter <- function(df, args, iter_seed) {
  df_use <- if (iter_seed == 0L) df else .cluster_resample(df, args, iter_seed)
  mats <- .fect_build_matrices(df_use, args$outcome, args$group,
                                 args$time, args$treatment)
  t0 <- Sys.time()
  # fect_ife ALWAYS uses the R svd() for its rank-r truncated SVD step.
  # Benchmarking (BENCHMARKS.md) showed the CUDA truncated-SVD path
  # loses to R at EVERY size tested, including large panels above the
  # mc size gate (2000x100: 0.09x; 8000x50: 0.26x). The reason: ife
  # needs only the top r singular triplets, and R's svd(nu=r, nv=r)
  # computes a cheap thin SVD, whereas the cuSOLVER kernel runs a full
  # gesvdj then truncates -- wasteful, and dominated by per-iteration
  # handle creation + H2D/D2H. Unlike mc (full SVD, GPU wins big at
  # scale), ife has no GPU-favourable regime, so the CUDA path is left
  # unwired here. The kernel + .fect_svd_truncated_cuda helper remain
  # for direct use / testing, just not in the production ife loop.
  use_cuda_svd <- FALSE
  fit <- .fect_ife_fit(mats$Y, mats$M, r = args$r %||% 2L,
                       tol = args$tol %||% 1e-5,
                       max_iter = args$max_iter %||% 500L,
                       use_cuda_svd = use_cuda_svd)
  res <- .fect_ife_compute_att(mats$Y, mats$M, fit, effects = args$effects)
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
    # ife always runs the SVD step on the host (see .fect_ife_one_iter
    # comment): the CUDA truncated-SVD path loses at every size.
    backend        = if (use_cuda_svd) "cuda" else "r",
    fect_method    = "ife",
    fect_iter      = fit$iter,
    fect_delta     = fit$delta,
    fect_r         = args$r %||% 2L
  )
}
