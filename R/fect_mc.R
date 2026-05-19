# ============================================================================
# fect_mc: matrix completion (Athey et al. 2021).
#
# Algorithm: iterative soft-thresholded SVD.
#
#   Y_complete = Y on control cells; Y_hat_old on treated cells.
#   repeat:
#     Y_complete = U D V^T   (full SVD of the n_units x n_periods matrix)
#     D_st = pmax(D - lambda, 0)         (soft-threshold)
#     Y_hat = U D_st V^T
#     update Y_complete: treated cells <- Y_hat[treated], control unchanged
#   until convergence (max |Y_hat - Y_hat_old| < tol).
#
# Counterfactual: Y_hat at treated cells. ATT = mean(Y - Y_hat) over treated.
#
# Picking lambda: in the reference fect package this is done by k-fold
# cross-validation over a grid. Here we accept a user-supplied lambda
# for now (default = a heuristic based on the SVD of the controls
# matrix). CV-based selection is a follow-up.
#
# The CUDA path reuses the same SVD primitive as fect_ife.
# ============================================================================


# Default lambda heuristic: a fraction of the largest singular value of
# the controls-only matrix. Conservative shrinkage; users should
# typically grid-search.
#' @keywords internal
#' @noRd
.fect_mc_default_lambda <- function(Y, M, frac = 0.1) {
  Y_c <- Y
  Y_c[M == 1L] <- 0  # zero treated cells for the initial scan
  Y_c[is.na(Y_c)] <- 0
  s <- svd(Y_c, nu = 0L, nv = 0L)
  frac * max(s$d)
}


# Cross-validated lambda selection for fect_mc.
#
# Splits the CONTROL cells (M == 0) into K folds at random. For each
# lambda in the grid:
#   - For each fold, mask the fold's cells (treat them as "missing"),
#     fit MC on the rest, and predict the held-out fold.
#   - Record out-of-sample MSE on the fold.
# Average across folds gives a CV-MSE curve. Pick the lambda minimising it.
#
# Default grid: geometric from 0.01 * sigma_max to 1.0 * sigma_max, 10 points.
# K defaults to 5.
#
# Reference: fect package's MC option does this same procedure.
#' @keywords internal
#' @noRd
.fect_mc_cv_lambda <- function(Y, M, K = 5L, n_grid = 10L,
                                  tol = 1e-5, max_iter = 200L,
                                  seed = 1L) {
  n_units   <- nrow(Y)
  n_periods <- ncol(Y)
  # Lambda grid: geometric from 0.01 * sigma_max to 1.0 * sigma_max.
  Y_c <- Y; Y_c[M == 1L] <- 0; Y_c[is.na(Y_c)] <- 0
  s <- svd(Y_c, nu = 0L, nv = 0L)
  sigma_max <- max(s$d)
  if (!is.finite(sigma_max) || sigma_max <= 0) return(0.1)
  lambdas <- exp(seq(log(0.01 * sigma_max), log(sigma_max),
                     length.out = n_grid))

  # Build fold assignment over control cells only.
  set.seed(seed)
  ctrl_idx <- which(M == 0L & !is.na(Y))
  if (length(ctrl_idx) < K * 5L) {
    # Too few control cells to do CV reliably; fall back to heuristic.
    return(.fect_mc_default_lambda(Y, M))
  }
  fold_id <- sample(rep(seq_len(K), length.out = length(ctrl_idx)))

  cv_mse <- numeric(length(lambdas))
  for (li in seq_along(lambdas)) {
    lambda <- lambdas[li]
    fold_mse <- numeric(K)
    for (k in seq_len(K)) {
      hold_out_idx <- ctrl_idx[fold_id == k]
      # Build a panel where the held-out cells are TREATED in M_cv
      # (so MC's fit ignores them).
      M_cv <- M
      M_cv[hold_out_idx] <- 1L
      fit <- .fect_mc_fit(Y, M_cv, lambda = lambda,
                            tol = tol, max_iter = max_iter)
      # Out-of-sample MSE: actual Y vs Y_hat on the held-out cells.
      held_Y    <- Y[hold_out_idx]
      held_Yhat <- fit$Y_hat[hold_out_idx]
      ok <- !is.na(held_Y) & !is.na(held_Yhat)
      fold_mse[k] <- if (any(ok)) mean((held_Y[ok] - held_Yhat[ok])^2)
                     else NA_real_
    }
    cv_mse[li] <- mean(fold_mse, na.rm = TRUE)
  }
  best <- which.min(cv_mse)
  if (length(best) == 0L) return(.fect_mc_default_lambda(Y, M))
  lambdas[best]
}


# Soft-threshold a vector: max(x - lambda, 0).
#' @keywords internal
#' @noRd
.soft_threshold <- function(x, lambda) {
  pmax(x - lambda, 0)
}


# CUDA wrapper for one soft-thresholded-SVD step. Phase-1 wiring for
# task #80: returns Y_hat (m x n) and the count of non-zero singular
# values, or NULL if anything goes wrong.
#' @keywords internal
#' @noRd
.fect_svd_softthreshold_cuda <- function(Y_complete, lambda) {
  if (!isTRUE(tryCatch(didgpu_has_cuda_support(),
                       error = function(e) FALSE))) return(NULL)
  result <- tryCatch(
    didgpu_cuda_fect_svd_softthreshold_r(Y_complete = Y_complete,
                                           lambda = lambda),
    error = function(e) NULL)
  if (is.null(result)) return(NULL)
  if (is.null(result$Y_hat)) return(NULL)
  list(Y_hat = result$Y_hat,
       n_nonzero = as.integer(result$n_nonzero %||% NA_integer_))
}


# MC fit: iterative soft-thresholded SVD.
#
# `use_cuda_svd`: when TRUE, every iteration's soft-thresholded SVD
# step runs on the GPU via .fect_svd_softthreshold_cuda. If any CUDA
# call fails, the function silently switches back to host svd() for
# the rest of the fit.
#' @keywords internal
#' @noRd
.fect_mc_fit <- function(Y, M, lambda = NULL, tol = 1e-5, max_iter = 500L,
                          use_cuda_svd = FALSE) {
  n_units   <- nrow(Y)
  n_periods <- ncol(Y)
  if (is.null(lambda)) lambda <- .fect_mc_default_lambda(Y, M)

  # Initialise Y_complete: control cells = Y, treated cells = 0.
  Y_complete <- Y
  Y_complete[M == 1L] <- 0
  Y_complete[is.na(Y_complete)] <- 0

  Y_hat <- matrix(0, n_units, n_periods)
  prev_Y_hat <- Y_hat
  n_nz <- 0L

  for (iter in seq_len(max_iter)) {
    Y_hat_new <- NULL
    if (use_cuda_svd) {
      cuda_res <- .fect_svd_softthreshold_cuda(Y_complete, lambda)
      if (!is.null(cuda_res)) {
        Y_hat_new <- cuda_res$Y_hat
        n_nz <- cuda_res$n_nonzero
      } else {
        use_cuda_svd <- FALSE   # disable for remainder of fit
      }
    }
    if (is.null(Y_hat_new)) {
      s <- svd(Y_complete)
      D_st <- .soft_threshold(s$d, lambda)
      nz <- D_st > 0
      n_nz <- sum(nz)
      if (any(nz)) {
        Y_hat_new <- s$u[, nz, drop = FALSE] %*%
                     diag(D_st[nz], sum(nz), sum(nz)) %*%
                     t(s$v[, nz, drop = FALSE])
      } else {
        Y_hat_new <- matrix(0, n_units, n_periods)
      }
    }
    Y_hat <- Y_hat_new
    # Update Y_complete: replace treated cells with current estimate.
    Y_complete[M == 1L] <- Y_hat[M == 1L]
    delta <- max(abs(Y_hat - prev_Y_hat))
    prev_Y_hat <- Y_hat
    if (delta < tol) break
  }

  list(Y_hat = Y_hat, lambda = lambda,
       iter = iter, delta = delta,
       n_nonzero_singular = as.integer(n_nz))
}


# Compute ATT + per-event-time effects from a fitted MC model.
#' @keywords internal
#' @noRd
.fect_mc_compute_att <- function(Y, M, fit, effects = NULL) {
  n_units   <- nrow(Y)
  n_periods <- ncol(Y)
  Y_hat <- fit$Y_hat
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


# Per-iter mc fit.
#' @keywords internal
#' @noRd
.fect_mc_one_iter <- function(df, args, iter_seed) {
  df_use <- if (iter_seed == 0L) df else .cluster_resample(df, args, iter_seed)
  mats <- .fect_build_matrices(df_use, args$outcome, args$group,
                                 args$time, args$treatment)
  t0 <- Sys.time()
  use_cuda <- identical(args$backend, "cuda") &&
              isTRUE(tryCatch(didgpu_has_cuda_support(),
                               error = function(e) FALSE))
  # CUDA path uses src/cuda_fect_mc.cu (cuSOLVER for SVD); untested
  # locally without nvcc.
  # Lambda selection: if user passed one, use it directly. Otherwise,
  # cross-validate on iter 0 (the point estimate) and reuse that lambda
  # for all bootstrap iters — full CV per bootstrap iter would be
  # prohibitively expensive.
  lambda <- if (!is.null(args$lambda)) {
    args$lambda
  } else {
    cached <- attr(df, "didgpu_mc_lambda_cv")
    if (!is.null(cached)) {
      cached
    } else if (iter_seed == 0L) {
      # CV on the point-estimate panel.
      l <- .fect_mc_cv_lambda(mats$Y, mats$M,
                                K = args$mc_cv_K %||% 5L,
                                n_grid = args$mc_cv_grid %||% 10L,
                                tol = args$tol %||% 1e-5,
                                max_iter = (args$max_iter %||% 500L) / 2L,
                                seed = args$seed %||% 1L)
      data.table::setattr(df, "didgpu_mc_lambda_cv", l)
      l
    } else {
      .fect_mc_default_lambda(mats$Y, mats$M)
    }
  }
  fit <- .fect_mc_fit(mats$Y, mats$M,
                        lambda = lambda,
                        tol = args$tol %||% 1e-5,
                        max_iter = args$max_iter %||% 500L,
                        use_cuda_svd = use_cuda)
  res <- .fect_mc_compute_att(mats$Y, mats$M, fit, effects = args$effects)
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
    fect_method    = "mc",
    fect_iter      = fit$iter,
    fect_delta     = fit$delta,
    fect_lambda    = fit$lambda,
    fect_n_nonzero = fit$n_nonzero_singular
  )
}
