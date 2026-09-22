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
  frac * .fect_mc_sigma_max(Y, M)
}


# Largest singular value of the matrix ACTUALLY being penalised, i.e. the
# residual after two-way fixed effects.
#
# The penalty must be scaled to the low-rank part, not to the level. The
# old version ran svd() on the raw outcome matrix with treated cells
# zeroed, so sigma_max -- and therefore lambda -- grew with the location
# of Y. That left the estimator level-dependent even once the fit itself
# removed fixed effects: on a known-zero DGP the ATT still drifted from
# -0.0321 to -0.0244 as a constant was added to the outcome.
#' @keywords internal
#' @noRd
.fect_mc_sigma_max <- function(Y, M) {
  Z <- Y
  obs <- !is.na(Z) & (M == 0L)
  fill <- if (any(obs)) mean(Z[obs]) else 0
  Z[M == 1L] <- NA_real_
  Z[is.na(Z)] <- fill
  R <- Z - .fect_mc_twoway(Z)
  s <- svd(R, nu = 0L, nv = 0L)
  mx <- max(s$d)
  if (!is.finite(mx) || mx <= 0) 0 else mx
}
# Cross-validated lambda selection for fect_mc.
#
# The held-out cells must look like the cells we actually have to
# predict. The treated block is a CONTIGUOUS RUN at the end of a unit's
# untreated history, so validation uses a rolling origin: for each unit,
# hold out the last `cv_nobs` untreated observations, then shift the
# origin back one period per fold.
#
# The previous version held out RANDOM SCATTERED control cells. A
# low-rank model interpolates isolated holes far more easily than it
# extrapolates a block, so scattered-hole MSE is minimised at too little
# shrinkage and the selected lambda left spurious factors in. On a DGP
# with NO factor structure, where full shrinkage is correct:
#
#   scattered CV -> lambda 3.15, keeps 3 singular values, ATT -0.031349
#   rolling CV   -> lambda 8.77, keeps 0 singular values, ATT -0.024425
#   fect                                                  ATT -0.024427
#
# Note fect works in lambda / (T * N) units, so its reported
# `lambda.cv` corresponds to `lambda * T * N` here.
#' @keywords internal
#' @noRd
.fect_mc_cv_lambda <- function(Y, M, n_grid = 10L, cv_nobs = 3L,
                                folds = 3L, tol = 1e-6, max_iter = 200L,
                                seed = 1L) {
  n_units   <- nrow(Y)
  n_periods <- ncol(Y)
  sigma_max <- .fect_mc_sigma_max(Y, M)
  if (!is.finite(sigma_max) || sigma_max <= 0) return(0.1)
  lambdas <- exp(seq(log(0.01 * sigma_max), log(sigma_max),
                     length.out = n_grid))

  # Rolling-origin hold-out sets, built once and reused for every lambda
  # so the grid is compared on identical cells.
  #
  # Only EVER-TREATED units are validated on. They are the units whose
  # counterfactuals have to be predicted, and their prediction task is
  # extrapolation past the end of a short untreated run. Never-treated
  # units have a full history, so holding out their last few periods is
  # a much easier, interpolation-like problem -- and because they are
  # usually the majority, including them dominates the MSE and selects
  # too little shrinkage. Measured against fect, with validation on all
  # units vs ever-treated only:
  #     long histories   all units lambda 0.244 -> ATT -0.027469
  #                      treated   lambda 5.257 -> ATT -0.024425  (fect -0.024427)
  #     short histories  treated   lambda 4.161 -> ATT +0.014847  (fect +0.014847)
  ever_treated <- rowSums(M == 1L, na.rm = TRUE) > 0L
  holds <- vector("list", folds)
  for (fd in seq_len(folds)) {
    idx <- integer(0)
    for (i in seq_len(n_units)) {
      if (!ever_treated[i]) next
      un <- which(M[i, ] == 0L & !is.na(Y[i, ]))
      if (length(un) < cv_nobs + fd) next
      last <- un[length(un) - (fd - 1L)]
      take <- un[un <= last & un > last - cv_nobs]
      idx <- c(idx, (take - 1L) * n_units + i)
    }
    holds[[fd]] <- idx
  }
  if (!length(unlist(holds))) return(.fect_mc_default_lambda(Y, M))

  cv_mse <- vapply(lambdas, function(lambda) {
    err <- numeric(0)
    for (fd in seq_len(folds)) {
      hold <- holds[[fd]]
      if (!length(hold)) next
      M_cv <- M
      M_cv[hold] <- 1L
      fit <- .fect_mc_fit(Y, M_cv, lambda = lambda,
                          tol = tol, max_iter = max_iter)
      e <- Y[hold] - fit$Y_hat[hold]
      err <- c(err, e[is.finite(e)])
    }
    if (length(err)) mean(err^2) else NA_real_
  }, numeric(1))

  if (all(is.na(cv_mse))) return(.fect_mc_default_lambda(Y, M))
  lambdas[which.min(cv_mse)]
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
# Exact additive two-way fixed effects of a COMPLETE matrix.
#
#   FE[i, t] = rowMean_i + colMean_t - grandMean
#
# Closed form (no iteration needed once the matrix has no missing
# cells). Adding a constant c to Z raises rowMean, colMean and
# grandMean each by c, so FE rises by c + c - c = c: the residual
# Z - FE is invariant to the level, which is what makes the MC
# estimator location-invariant.
#' @keywords internal
#' @noRd
.fect_mc_twoway <- function(Z) {
  gm <- mean(Z)
  a  <- rowMeans(Z) - gm
  x  <- colMeans(Z) - gm
  outer(a, rep(1, ncol(Z))) + outer(rep(1, nrow(Z)), x) + gm
}


# MC fit: two-way fixed effects PLUS an iterative soft-thresholded SVD
# of the residual.
#
# Athey et al. (2021) estimate Y = L + unit FE + time FE, penalising the
# nuclear norm of L ALONE. Penalising the raw outcome matrix instead
# shrinks the level, so the imputed counterfactual is biased toward zero
# and ATT = mean(Y - Y_hat) over treated cells inherits the units' level.
# That made the estimator depend on the location of Y: on a known-zero
# DGP, adding 100 to the outcome moved the reported ATT from +0.27 to
# +2.55, while fe, ife and fect::fect all returned -0.024 at every level.
# On a positive, trending outcome it manufactured large, monotonically
# rising, significant effects where fe and ife both found a null.
#
# Each iteration therefore:
#   1. takes two-way fixed effects off the current complete matrix,
#   2. soft-thresholds the SVD of the RESIDUAL,
#   3. adds the fixed effects back to form Y_hat,
#   4. re-imputes the treated cells from Y_hat.
#
# `use_cuda_svd`: when TRUE the SVD + soft-threshold step runs on the GPU
# via .fect_svd_softthreshold_cuda, applied to the residual so the GPU and
# host paths stay numerically identical. If any CUDA call fails, the
# function silently switches back to host svd() for the rest of the fit.
#' @keywords internal
#' @noRd
.fect_mc_fit <- function(Y, M, lambda = NULL, tol = 1e-5, max_iter = 500L,
                          use_cuda_svd = FALSE) {
  n_units   <- nrow(Y)
  n_periods <- ncol(Y)
  if (is.null(lambda)) lambda <- .fect_mc_default_lambda(Y, M)

  # Initialise Y_complete: control cells = Y, treated cells seeded from
  # the two-way fit on the controls so the first residual is sensible.
  Y_complete <- Y
  obs <- !is.na(Y_complete) & (M == 0L)
  seed_fill <- if (any(obs)) mean(Y_complete[obs]) else 0
  Y_complete[M == 1L] <- NA_real_
  Y_complete[is.na(Y_complete)] <- seed_fill

  Y_hat <- matrix(seed_fill, n_units, n_periods)
  prev_Y_hat <- Y_hat
  n_nz <- 0L
  delta <- NA_real_

  for (iter in seq_len(max_iter)) {
    FE <- .fect_mc_twoway(Y_complete)
    R  <- Y_complete - FE

    L_new <- NULL
    if (use_cuda_svd) {
      cuda_res <- .fect_svd_softthreshold_cuda(R, lambda)
      if (!is.null(cuda_res)) {
        L_new <- cuda_res$Y_hat
        n_nz  <- cuda_res$n_nonzero
      } else {
        use_cuda_svd <- FALSE   # disable for remainder of fit
      }
    }
    if (is.null(L_new)) {
      s <- svd(R)
      D_st <- .soft_threshold(s$d, lambda)
      nz <- D_st > 0
      n_nz <- sum(nz)
      if (any(nz)) {
        L_new <- s$u[, nz, drop = FALSE] %*%
                 diag(D_st[nz], sum(nz), sum(nz)) %*%
                 t(s$v[, nz, drop = FALSE])
      } else {
        L_new <- matrix(0, n_units, n_periods)
      }
    }

    Y_hat <- FE + L_new
    # Re-impute the treated cells; control cells keep their observed
    # values, and genuinely missing control cells follow Y_hat too.
    Y_complete[M == 1L] <- Y_hat[M == 1L]
    miss <- is.na(Y) & (M == 0L)
    if (any(miss)) Y_complete[miss] <- Y_hat[miss]

    delta <- max(abs(Y_hat - prev_Y_hat))
    prev_Y_hat <- Y_hat
    if (is.finite(delta) && delta < tol) break
  }

  list(Y_hat = Y_hat, lambda = lambda,
       iter = iter, delta = delta,
       converged = isTRUE(is.finite(delta) && delta < tol),
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
                               error = function(e) FALSE)) &&
              .fect_cuda_svd_worthwhile(nrow(mats$Y), ncol(mats$Y))
  # CUDA path uses cuSOLVER for the soft-thresholded SVD. Gated behind
  # .fect_cuda_svd_worthwhile: for small matrices the GPU SVD is far
  # slower than R's LAPACK (see BENCHMARKS.md), so backend = "cuda"
  # transparently falls back to svd() below the size threshold.
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
                                n_grid  = args$mc_cv_grid %||% 10L,
                                cv_nobs = args$mc_cv_nobs %||% 3L,
                                folds   = args$mc_cv_folds %||% 3L,
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
