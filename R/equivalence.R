# ============================================================================
# Pre-trends equivalence (TOST) test for the de Chaisemartin didgpu()
# event study.
#
# The conventional pre-trend test checks H0: pre-trend = 0 and treats a
# non-significant result as support for parallel trends. That logic is
# backwards — "absence of evidence is not evidence of absence", and an
# underpowered design fails to reject precisely when it should worry us
# most. The equivalence (TOST) framing flips it: we pick a margin `delta`
# (the largest pre-trend we'd call economically negligible) and test, at
# each placebo horizon,
#
#     H0: |placebo effect| >= delta    vs.    H1: |placebo effect| < delta
#
# Rejecting H0 is *positive* evidence that the pre-treatment deviation is
# smaller than delta. This mirrors didgpu_fect_equivalence() for the fect
# family (see R/fect_robustness.R) and the equivalence-testing literature
# on pre-trends (e.g. Roth 2022; Hartman & Hidalgo 2018 for TOST in causal
# settings).
# ============================================================================


#' Pre-trends equivalence (TOST) test for a didgpu event study
#'
#' Runs a two-one-sided-tests (TOST) equivalence test on the placebo
#' (pre-treatment) estimates of a [didgpu()] fit. For a user-supplied
#' margin `delta`, each placebo horizon tests
#' \eqn{H_0: |\theta| \ge \delta} against \eqn{H_1: |\theta| < \delta};
#' rejecting \eqn{H_0} is positive evidence the pre-trend is within
#' \eqn{\pm}`delta`. The joint claim "all placebos lie within `delta`"
#' follows by the intersection-union principle: it holds iff every horizon
#' individually rejects.
#'
#' Unlike a conventional placebo test (where a *large* p-value is the
#' hoped-for result but only weakly informative), here a *small*
#' `equivalence_p` is the hoped-for result and is a genuine rejection.
#'
#' @param x A `didgpu_result` from [didgpu()] run with `placebo > 0` and
#'   `bootstrap_reps > 0` (equivalence testing needs standard errors).
#' @param delta Positive numeric equivalence margin on the outcome scale —
#'   the largest pre-trend you would consider economically negligible.
#' @param alpha One-sided significance level (default 0.05).
#' @return A `didgpu_equivalence` object (a data.frame with one row per
#'   placebo horizon: `event_time`, `estimate`, `std.error`,
#'   `equivalence_p`, `passes_at_delta`) carrying attributes `delta`,
#'   `alpha`, `joint_pass` (TRUE iff every horizon passes; NA if any SE is
#'   missing) and `breakdown_delta` (the smallest margin at which the joint
#'   equivalence would hold at level `alpha`).
#' @seealso [didgpu()], [didgpu_fect_equivalence()] for the fect-family
#'   analogue, [didgpu_honest_did()] for HonestDiD sensitivity bounds.
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
#'                            tau_profile = c(0.5, 1.0), seed = 7L)
#' p$D <- as.integer(p$D >= 0.5)
#' fit <- didgpu(p, "Y", "unit", "period", "D", effects = 3L, placebo = 3L,
#'               bootstrap_reps = 200L, verbose = FALSE)
#' didgpu_equivalence(fit, delta = 0.5)
#' }
#' @export
didgpu_equivalence <- function(x, delta, alpha = 0.05) {
  stopifnot(inherits(x, "didgpu_result"))
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta <= 0)
    stop("`delta` must be a single positive number (the equivalence margin).")
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single number in (0, 1).")

  pl <- x$results$Placebos
  if (is.null(pl) || nrow(pl) == 0L)
    stop("`x` has no placebo estimates. Re-run didgpu() with placebo > 0.")
  if (!all(c("Estimate", "SE") %in% colnames(pl)))
    stop("`x$results$Placebos` is missing the 'Estimate'/'SE' columns.")

  est <- as.numeric(pl[, "Estimate"])
  se  <- as.numeric(pl[, "SE"])
  if (all(is.na(se)))
    stop("All placebo SEs are NA. Re-run didgpu() with bootstrap_reps > 0 ",
         "so the equivalence test has standard errors.")
  if (any(is.na(se)))
    warning("Some placebo SEs are NA; their equivalence p-values will be NA.")

  # TOST: equivalence_p = max( P(Z > (delta - est)/se),
  #                            P(Z > (delta + est)/se) ).
  # The binding side reduces to P(Z > (delta - |est|)/se), so a small
  # value means we reject |theta| >= delta at both edges.
  p1   <- stats::pnorm((delta - est) / se, lower.tail = FALSE)
  p2   <- stats::pnorm((delta + est) / se, lower.tail = FALSE)
  eq_p <- pmax(p1, p2)
  passes <- eq_p < alpha

  # Smallest margin at which horizon j would pass at level alpha solves
  # P(Z > (d - |est|)/se) = alpha  =>  d = |est| + se * z_{1-alpha}.
  # The joint (IU) breakdown margin is the max over horizons.
  z         <- stats::qnorm(1 - alpha)
  min_delta <- abs(est) + se * z
  breakdown <- if (all(is.na(min_delta))) NA_real_
               else max(min_delta, na.rm = TRUE)

  out <- data.frame(
    event_time      = -seq_len(nrow(pl)),
    estimate        = est,
    std.error       = se,
    equivalence_p   = eq_p,
    passes_at_delta = passes,
    stringsAsFactors = FALSE
  )
  attr(out, "delta")           <- delta
  attr(out, "alpha")           <- alpha
  attr(out, "joint_pass")      <- all(passes)   # NA if any SE missing
  attr(out, "breakdown_delta") <- breakdown
  class(out) <- c("didgpu_equivalence", "data.frame")
  out
}


#' Print method for didgpu_equivalence
#' @param x A `didgpu_equivalence` object from [didgpu_equivalence()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.didgpu_equivalence <- function(x, ...) {
  delta <- attr(x, "delta"); alpha <- attr(x, "alpha")
  jp    <- attr(x, "joint_pass"); bd <- attr(x, "breakdown_delta")
  cat(sprintf("didgpu pre-trends equivalence test (delta = %.4g, alpha = %.3g)\n",
              delta, alpha))
  print(as.data.frame(x), row.names = FALSE)
  cat(strrep("-", 60), "\n", sep = "")
  joint_msg <- if (isTRUE(jp))
      sprintf("PASS - every placebo is within +/- %.4g at level %.3g", delta, alpha)
    else if (is.na(jp))
      "INDETERMINATE - some placebo SEs are missing"
    else
      "FAIL - at least one placebo is not shown to be within delta"
  cat(sprintf("Joint (intersection-union): %s\n", joint_msg))
  if (!is.na(bd))
    cat(sprintf("Smallest defensible margin at alpha = %.3g: delta >= %.4g\n",
                alpha, bd))
  cat("Interpretation: a SMALL equivalence_p rejects |pre-trend| >= delta,\n")
  cat("i.e. it is positive evidence the pre-trend lies within +/- delta.\n")
  invisible(x)
}
