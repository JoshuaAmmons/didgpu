# ============================================================================
# Robustness diagnostics for fect estimators.
#
# Two diagnostics, both standard in the fect / HonestDiD literature:
#
# 1. Placebo tests (`didgpu_fect_placebo`):
#    Refit the fect estimator pretending the M pre-treatment periods
#    immediately before F_g are "treated." If the estimated effect at
#    those PLACEBO cells is significantly different from zero, the
#    parallel-trends / no-anticipation assumption is suspect.
#
#    Output: per-placebo-horizon effect, SE, joint test p-value across
#    all placebos. Standard interpretation: small p-values are bad
#    (reject parallel trends); large p-values are reassuring but NOT
#    a proof.
#
# 2. Equivalence tests (`didgpu_fect_equivalence`):
#    Reverse the testing logic. Instead of H0: effect = 0, test
#    H0: |effect| > delta for a user-supplied tolerance delta.
#    REJECTING this null lets us conclude the placebo deviation is
#    below delta — i.e., the panel passes the user's bar for "close
#    enough to parallel trends."
#
#    Output: per-placebo-horizon equivalence p-value at the chosen
#    delta; an inverted-CI version reports the smallest delta at which
#    the panel passes.
#
# Both diagnostics share the placebo-fitting machinery defined here.
# ============================================================================


#' Run placebo (pre-treatment) tests on a fitted fect model
#'
#' For each unit g with first-switch period F_g, "treat" the M
#' pre-treatment periods immediately before F_g (i.e., periods
#' F_g - 1, F_g - 2, ..., F_g - M) and refit the fect estimator on
#' the remaining truly-pre-treatment cells. Effects at those placebo
#' cells should be statistically indistinguishable from zero if the
#' identifying assumption holds.
#'
#' @param df A panel data.frame (the same one used to fit).
#' @param outcome,group,time,treatment Column names.
#' @param method One of `"fe"`, `"ife"`, `"mc"`.
#' @param n_placebos Integer. Number of pre-treatment placebos to test
#'   per unit. Default `3L`. Equivalent to the `placebo` option in the
#'   fect package.
#' @param r,lambda,tol,max_iter Method-specific args forwarded to the
#'   underlying fit (see [didgpu_fect()]).
#' @param bootstrap_reps Integer. Number of bootstrap replicates for
#'   SE estimation. Default `0L` = no bootstrap, SEs are NA.
#' @param seed Integer. Bootstrap seed.
#' @param backend Backend for the underlying fits.
#' @return A data.frame with one row per placebo horizon: `horizon`
#'   (negative, -1 = period immediately before F_g), `estimate`,
#'   `se`, `p_value` (test of H0: estimate = 0), plus a `joint_p`
#'   attribute giving the joint p-value across all placebos.
#'
#' @examples
#' p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
#'                             tau_profile = c(0.5, 1.0), seed = 17L)
#' pl <- didgpu_fect_placebo(p, "Y", "unit", "period", "D",
#'                            method = "fe", n_placebos = 2L)
#' pl
#'
#' @export
didgpu_fect_placebo <- function(
    df, outcome, group, time, treatment,
    method = c("fe", "ife", "mc"),
    n_placebos = 3L,
    r = 2L, lambda = NULL,
    tol = 1e-5, max_iter = 500L,
    bootstrap_reps = 0L, seed = 1L,
    backend = "auto") {

  method <- match.arg(method)
  n_placebos <- as.integer(n_placebos)
  stopifnot(n_placebos >= 1L)

  mats <- .fect_build_matrices(df, outcome, group, time, treatment)
  Y <- mats$Y; M <- mats$M
  n_units   <- nrow(Y)
  n_periods <- ncol(Y)

  # Compute F_g per unit (first treated period) and the placebo cells.
  F_g <- apply(M, 1L, function(row) {
    idx <- which(row == 1L)
    if (length(idx) == 0L) NA_integer_ else min(idx)
  })

  fit_one_placebo_panel <- function(seed_offset) {
    # For each placebo horizon h in 1..n_placebos:
    # Mark cell (i, F_g[i] - h) as "treated" in M_pl, leaving all
    # subsequent cells out of the fit. Refit and pull the average
    # effect at those placebo cells.
    M_pl <- matrix(0L, n_units, n_periods)
    # Mask out the "real" post-treatment cells so they don't pollute
    # the placebo fit. We drop EVERYTHING from F_g onward.
    drop_mask <- matrix(0L, n_units, n_periods)
    for (i in seq_len(n_units)) {
      if (is.na(F_g[i])) next
      fg <- F_g[i]
      drop_mask[i, seq_len(n_periods) >= fg] <- 1L
      # Mark the placebo cells.
      for (h in seq_len(n_placebos)) {
        t_idx <- fg - h
        if (t_idx >= 1L) M_pl[i, t_idx] <- 1L
      }
    }
    Y_pl <- Y
    Y_pl[drop_mask == 1L] <- NA_real_

    fit <- switch(method,
      fe  = .fect_fe_fit(Y_pl, M_pl, tol = tol, max_iter = max_iter),
      ife = .fect_ife_fit(Y_pl, M_pl, r = r, tol = tol, max_iter = max_iter),
      mc  = .fect_mc_fit(Y_pl, M_pl, lambda = lambda,
                          tol = tol, max_iter = max_iter)
    )

    # Predicted Y for placebo cells, then residual.
    Y_hat <- switch(method,
      fe  = fit$alpha + matrix(fit$xi, n_units, n_periods, byrow = TRUE),
      ife = fit$alpha + matrix(fit$xi, n_units, n_periods, byrow = TRUE) +
            fit$L %*% fit$F,
      mc  = fit$Y_hat
    )
    residual <- Y_pl - Y_hat
    # Per-horizon average residual at the placebo cells.
    horizons <- -seq_len(n_placebos)
    eff <- numeric(n_placebos)
    n_cells <- integer(n_placebos)
    for (h in seq_len(n_placebos)) {
      cells_h <- matrix(FALSE, n_units, n_periods)
      for (i in seq_len(n_units)) {
        if (is.na(F_g[i])) next
        t_idx <- F_g[i] - h
        if (t_idx >= 1L && M_pl[i, t_idx] == 1L &&
            !is.na(residual[i, t_idx])) {
          cells_h[i, t_idx] <- TRUE
        }
      }
      n_cells[h] <- sum(cells_h)
      eff[h] <- if (n_cells[h] > 0L) mean(residual[cells_h], na.rm = TRUE)
                else NA_real_
    }
    list(horizons = horizons, estimate = eff, n_cells = n_cells)
  }

  # Point estimate (no bootstrap).
  point <- fit_one_placebo_panel(0L)

  # Bootstrap SEs (cluster on group). Build per-iter args bundle for
  # .cluster_resample compatibility.
  if (bootstrap_reps > 0L) {
    args_for_resample <- list(group = group, time = time, cluster = NULL)
    boot_mat <- matrix(NA_real_, nrow = bootstrap_reps,
                        ncol = n_placebos)
    for (b in seq_len(bootstrap_reps)) {
      df_b <- .cluster_resample(df, args_for_resample, iter_seed = seed + b)
      mats_b <- .fect_build_matrices(df_b, outcome, group, time, treatment)
      Y <- mats_b$Y; M <- mats_b$M
      n_units <- nrow(Y); n_periods <- ncol(Y)
      F_g <- apply(M, 1L, function(row) {
        idx <- which(row == 1L)
        if (length(idx) == 0L) NA_integer_ else min(idx)
      })
      res_b <- fit_one_placebo_panel(seed + b)
      boot_mat[b, ] <- res_b$estimate
    }
    se <- apply(boot_mat, 2L, function(col) {
      if (sum(!is.na(col)) >= 2L) stats::sd(col, na.rm = TRUE) else NA_real_
    })
  } else {
    se <- rep(NA_real_, n_placebos)
  }

  z <- ifelse(is.na(se) | se == 0, NA_real_, point$estimate / se)
  p_value <- ifelse(is.na(z), NA_real_, 2 * stats::pnorm(-abs(z)))

  out <- data.frame(
    horizon  = point$horizons,
    estimate = point$estimate,
    se       = se,
    p_value  = p_value,
    n_cells  = point$n_cells,
    stringsAsFactors = FALSE
  )
  # Joint p-value: chi-square of the placebo vector against the
  # bootstrap covariance (if available) at the standard 5% level.
  joint_p <- NA_real_
  if (bootstrap_reps > 0L && all(!is.na(point$estimate))) {
    V <- tryCatch(stats::cov(boot_mat[, , drop = FALSE]),
                  error = function(e) NULL)
    if (!is.null(V) && all(is.finite(diag(V)))) {
      inv_V <- tryCatch(solve(V), error = function(e) NULL)
      if (!is.null(inv_V)) {
        wald <- as.numeric(point$estimate %*% inv_V %*% point$estimate)
        if (is.finite(wald) && wald >= 0) {
          joint_p <- 1 - stats::pchisq(wald, df = n_placebos)
        }
      }
    }
  }
  attr(out, "joint_p") <- joint_p
  attr(out, "method")  <- method
  class(out) <- c("didgpu_fect_placebo", class(out))
  out
}


#' Equivalence test on fect placebos
#'
#' Tests H0: |effect at placebo horizon h| > delta for a user-supplied
#' tolerance delta. REJECTING this null lets us conclude the placebo
#' deviation is below delta. Implemented as a two-one-sided-tests
#' (TOST) procedure on the bootstrap distribution.
#'
#' @param placebo_result A data.frame returned by [didgpu_fect_placebo()]
#'   with non-NA `se` (i.e., from a run with `bootstrap_reps > 0`).
#' @param delta Numeric. The equivalence margin.
#' @return The input data.frame augmented with `equivalence_p` (the
#'   TOST p-value: small means the placebo deviation IS below delta)
#'   and `passes_at_delta` (logical, TRUE if equivalence_p < 0.05).
#'
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
#'                             tau_profile = c(0.5, 1.0), seed = 17L)
#' pl <- didgpu_fect_placebo(p, "Y", "unit", "period", "D",
#'                            method = "fe", n_placebos = 2L,
#'                            bootstrap_reps = 30L, seed = 1L)
#' eq <- didgpu_fect_equivalence(pl, delta = 0.5)
#' eq
#' }
#' @export
didgpu_fect_equivalence <- function(placebo_result, delta) {
  stopifnot(inherits(placebo_result, "didgpu_fect_placebo"))
  stopifnot(is.numeric(delta), length(delta) == 1L, delta > 0)
  if (any(is.na(placebo_result$se))) {
    warning("Equivalence test requires bootstrap SEs. Some placebo SEs ",
            "are NA; corresponding equivalence p-values will be NA.")
  }
  # TOST: p1 = P(estimate <= -delta), p2 = P(estimate >= +delta).
  # Combined: equivalence_p = max(P(Z > (delta - estimate)/se),
  #                               P(Z > (delta + estimate)/se))
  est <- placebo_result$estimate
  se  <- placebo_result$se
  p1 <- stats::pnorm((delta - est) / se, lower.tail = FALSE)
  p2 <- stats::pnorm((delta + est) / se, lower.tail = FALSE)
  eq_p <- pmax(p1, p2)
  out <- placebo_result
  out$equivalence_p   <- eq_p
  out$passes_at_delta <- eq_p < 0.05
  attr(out, "delta") <- delta
  out
}


#' Print method for didgpu_fect_placebo
#' @param x A placebo-test result.
#' @param ... Unused.
#' @return The input invisibly.
#' @export
print.didgpu_fect_placebo <- function(x, ...) {
  cat(sprintf("didgpu_fect placebo test (method = '%s')\n",
              attr(x, "method")))
  print(as.data.frame(x), row.names = FALSE)
  jp <- attr(x, "joint_p")
  if (!is.null(jp) && !is.na(jp)) {
    cat(sprintf("\nJoint test of placebos: chi-sq p = %.4g\n", jp))
  }
  cat(strrep("-", 50), "\n", sep = "")
  cat("Interpretation:\n")
  cat("  Small p_value at a horizon = reject the null that the\n")
  cat("  pre-treatment effect was zero (parallel-trends concern).\n")
  cat("  For equivalence testing at a user-supplied delta, see\n")
  cat("  didgpu_fect_equivalence().\n")
  invisible(x)
}
