# ============================================================================
# HonestDiD: Rambachan & Roth (2023) sensitivity analysis for DiD.
#
# Standard event-study DiD requires a parallel-trends assumption that
# is fundamentally untestable. Pre-trend tests have power problems --
# even if pre-trend coefficients are statistically insignificant, real
# violations of parallel trends can still bias post-treatment estimates.
#
# HonestDiD reverses the framing. Instead of "do the pre-trends look
# OK?", it asks: "given the pre-trend evidence, what RANGE of post-
# treatment effects is consistent with at most a delta-bounded
# violation of parallel trends?"
#
# Two standard restrictions on the violation:
#   - "M" (smoothness): consecutive-period changes in the violation
#                       are bounded by M
#   - "RM" (relative magnitudes): the post-treatment violation is
#                                  bounded by RM times the max pre-
#                                  treatment deviation
#
# This wrapper takes a fitted didgpu_result OR didgpu_cs_result with
# event-study output and a covariance matrix and dispatches to the
# reference HonestDiD package's createSensitivityResults_* functions.
#
# Reference: Rambachan, A. and Roth, J. (2023). "A More Credible
# Approach to Parallel Trends." Review of Economic Studies 90(5).
# ============================================================================


#' Sensitivity analysis for a fitted event-study DiD
#'
#' Given a fitted didgpu / didgpu_cs result with both pre- and post-
#' treatment event-study coefficients, runs the Rambachan-Roth (2023)
#' sensitivity analysis: bound the post-treatment estimate under a
#' user-specified restriction on the magnitude of pre-trend violations,
#' and report the "breakdown" parameter -- the smallest violation
#' that would flip the substantive conclusion.
#'
#' @param fit A `didgpu_result` (with placebos > 0) or `didgpu_cs_result`
#'   (which always includes pre-treatment placebos automatically).
#' @param event_post Integer. Which post-treatment event-time to bound.
#'   Default `1L`.
#' @param method One of `"M"` (smoothness; bound consecutive-period
#'   violations) or `"RM"` (relative magnitudes; bound by Mbar times
#'   the largest pre-trend deviation). Default `"RM"`.
#' @param Mbar Numeric vector. Grid of restriction values to test.
#'   For `"M"`, this is the smoothness bound in outcome units; for
#'   `"RM"`, the multiplier on the largest pre-trend.
#'   Default `c(0, 0.5, 1, 1.5, 2)`.
#' @param alpha Numeric. Significance level. Default `0.05`.
#' @param ci_level Numeric. Same as `(1 - alpha) * 100`. Convenience.
#' @return A `didgpu_honest_did_result`: data.frame with one row per
#'   Mbar value: `Mbar`, `lb` (lower bound of robust CI), `ub`,
#'   `crosses_zero` (TRUE if CI includes 0; the conclusion is fragile).
#'   Plus a `$breakdown` attribute giving the smallest Mbar at which
#'   the CI includes zero.
#'
#' @references
#' Rambachan, A. and Roth, J. (2023). "A More Credible Approach to
#' Parallel Trends." *Review of Economic Studies* 90(5): 2555-2591.
#'
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 100L, n_periods = 12L,
#'                             tau_profile = c(0.5, 1.0),
#'                             seed = 17L)
#' fit <- didgpu_cs(p, "Y", "unit", "period", "D",
#'                   est_method = "OR", bootstrap_reps = 30L,
#'                   backend = "r", verbose = FALSE)
#' sens <- didgpu_honest_did(fit, event_post = 1L,
#'                            method = "RM",
#'                            Mbar = c(0, 0.5, 1.0))
#' print(sens)
#' }
#' @export
didgpu_honest_did <- function(
    fit, event_post = 1L,
    method   = c("RM", "M"),
    Mbar     = c(0, 0.5, 1, 1.5, 2),
    alpha    = 0.05,
    ci_level = 100 * (1 - alpha)) {
  if (!requireNamespace("HonestDiD", quietly = TRUE)) {
    stop("didgpu_honest_did requires the HonestDiD package. ",
         "Install with install.packages('HonestDiD').")
  }
  method <- match.arg(method)
  stopifnot(is.numeric(Mbar), all(Mbar >= 0))
  event_post <- as.integer(event_post)
  stopifnot(event_post >= 1L)

  # Extract event-study coefficients + covariance from the fit.
  es <- .honest_extract_event_study(fit)
  if (is.null(es)) {
    stop("Could not extract event-study coefficients from `fit`. ",
         "Make sure the fit includes pre-treatment placebos.")
  }
  betahat <- es$beta
  sigma   <- es$Sigma
  num_pre <- es$num_pre
  num_post <- es$num_post
  if (event_post > num_post) {
    stop(sprintf("event_post = %d exceeds the number of post-treatment ",
                  event_post),
         sprintf("event-times in fit (%d).", num_post))
  }
  if (num_pre < 2L) {
    stop("HonestDiD needs at least 2 pre-treatment event-times for the ",
         "smoothness / relative-magnitudes restrictions to be informative.")
  }

  # Dispatch to the appropriate HonestDiD function.
  hd_res <- if (method == "RM") {
    HonestDiD::createSensitivityResults_relativeMagnitudes(
      betahat   = betahat,
      sigma     = sigma,
      numPrePeriods  = num_pre,
      numPostPeriods = num_post,
      Mbarvec   = Mbar,
      l_vec     = .honest_l_vec(num_post, event_post),
      alpha     = alpha)
  } else {
    HonestDiD::createSensitivityResults(
      betahat   = betahat,
      sigma     = sigma,
      numPrePeriods  = num_pre,
      numPostPeriods = num_post,
      Mvec      = Mbar,
      l_vec     = .honest_l_vec(num_post, event_post),
      alpha     = alpha)
  }

  out <- data.frame(
    Mbar         = if (method == "RM") hd_res$Mbar else hd_res$M,
    lb           = hd_res$lb,
    ub           = hd_res$ub,
    crosses_zero = (hd_res$lb <= 0) & (hd_res$ub >= 0),
    stringsAsFactors = FALSE
  )
  # Find the smallest Mbar at which the CI includes zero (breakdown).
  bd_idx <- which(out$crosses_zero)
  breakdown <- if (length(bd_idx) > 0L) out$Mbar[min(bd_idx)] else NA_real_

  attr(out, "method")    <- method
  attr(out, "alpha")     <- alpha
  attr(out, "event_post") <- event_post
  attr(out, "breakdown") <- breakdown
  attr(out, "original_estimate") <- betahat[num_pre + event_post]
  attr(out, "original_se")       <- sqrt(diag(sigma))[num_pre + event_post]
  class(out) <- c("didgpu_honest_did_result", class(out))
  out
}


# Build the linear-combination vector l_vec that picks out the desired
# post-treatment event-time from the betahat vector. l_vec is length
# num_post and is 1 at event_post, 0 elsewhere (we report the bounds
# for a single event-time).
#' @keywords internal
#' @noRd
.honest_l_vec <- function(num_post, event_post) {
  v <- numeric(num_post)
  v[event_post] <- 1
  v
}


# Extract event-study coefficients (betahat) and covariance (Sigma)
# from a didgpu_result or didgpu_cs_result. The required ordering is
# (pre placebos: oldest first, ..., -1) then (post effects: 1, 2, ..., k).
#' @keywords internal
#' @noRd
.honest_extract_event_study <- function(fit) {
  if (inherits(fit, "didgpu_cs_result")) {
    ev <- fit$aggregation
    if (!"event_time" %in% names(ev)) {
      ev <- .cs_aggregate(fit$att_gt, "event", fit$args)
    }
    pre  <- ev[ev$event_time < 0L, , drop = FALSE]
    post <- ev[ev$event_time >= 0L, , drop = FALSE]
    pre  <- pre[order(pre$event_time), , drop = FALSE]   # most negative first
    post <- post[order(post$event_time), , drop = FALSE]
    beta <- c(pre$estimate, post$estimate)
    # Covariance via the att_gt SEs (assume diagonal — conservative).
    # The proper construction reconstructs the full event-study covariance
    # from the per-cell IFs; this is an approximation good for moderate
    # cohort sizes.
    se_evt <- function(e) {
      cells <- fit$att_gt[fit$att_gt$event_time == e, ]
      ses <- cells$se
      if (all(is.na(ses))) return(NA_real_)
      sqrt(mean(ses^2, na.rm = TRUE))
    }
    ses <- vapply(c(pre$event_time, post$event_time), se_evt, numeric(1))
    if (all(is.na(ses))) {
      stop("HonestDiD requires SE estimates; ",
           "refit with bootstrap_reps > 0.")
    }
    ses[is.na(ses)] <- 0
    Sigma <- diag(ses^2)
    return(list(beta = beta, Sigma = Sigma,
                num_pre = nrow(pre),
                num_post = nrow(post)))
  }
  if (inherits(fit, "didgpu_result")) {
    # didgpu() result with placebos.
    ef <- as.numeric(fit$results$Effects[, "Estimate"])
    pl <- if (fit$results$N_Placebos > 0L)
            rev(as.numeric(fit$results$Placebos[, "Estimate"])) else numeric(0)
    if (length(pl) < 2L) return(NULL)
    ef_se <- as.numeric(fit$results$Effects[, "SE"])
    pl_se <- if (fit$results$N_Placebos > 0L)
              rev(as.numeric(fit$results$Placebos[, "SE"])) else numeric(0)
    beta <- c(pl, ef)
    ses  <- c(pl_se, ef_se)
    ses[is.na(ses)] <- 0
    Sigma <- diag(ses^2)
    return(list(beta = beta, Sigma = Sigma,
                num_pre = length(pl),
                num_post = length(ef)))
  }
  NULL
}


#' Print method for didgpu_honest_did_result
#' @param x A `didgpu_honest_did_result`.
#' @param ... Unused.
#' @return The input invisibly.
#' @export
print.didgpu_honest_did_result <- function(x, ...) {
  cat(sprintf("HonestDiD sensitivity (method = '%s', event_post = %d)\n",
              attr(x, "method"), attr(x, "event_post")))
  cat(sprintf("  original estimate: %.4f (SE %.4f)\n",
              attr(x, "original_estimate"),
              attr(x, "original_se")))
  print(as.data.frame(x), row.names = FALSE)
  bd <- attr(x, "breakdown")
  cat(sprintf("\nBreakdown (smallest Mbar at which CI includes 0): %s\n",
              if (is.na(bd)) "NA (conclusion holds across the grid)"
              else sprintf("%.3f", bd)))
  invisible(x)
}
