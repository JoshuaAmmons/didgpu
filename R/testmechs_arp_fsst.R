# ============================================================================
# TestMechs ARP and FSST moment-inequality tests.
#
# Both methods are alternatives to Cox-Shi (CS) for testing whether
# the observed moments are consistent with the polytope implied by the
# null. They differ in how they compute the critical value:
#
#   - CS  (already implemented): chi-square approximation via projection
#         QP. Conservative; fast.
#   - ARP (Andrews, Roth, Pakes 2023): least-favorable critical value
#         via Monte Carlo on N(0, Sigma) draws. More precise for many
#         active constraints.
#   - FSST (Fang, Santos, Shaikh, Torgovitsky 2023): cone-based test
#         using a secondary bootstrap of cone statistics.
#
# Both are implemented here in their simplified "test only" form
# (return reject / pval); the production-grade versions in the
# reference HonestDiD / lpinfer packages add tuning parameters (kappa
# for hybrid ARP, lambda for FSST) that we expose via `...`.
# ============================================================================


# Least-favorable critical value for ARP. Computes by Monte Carlo:
# draw M samples from N(0, Sigma), evaluate the test statistic on each,
# return the alpha-quantile.
#' @keywords internal
#' @noRd
.testmechs_arp_cv <- function(Sigma, A, alpha = 0.05, M_mc = 5000L,
                                seed = 1L) {
  p <- ncol(A)
  if (nrow(Sigma) != p || ncol(Sigma) != p) {
    stop("Sigma must be p x p.")
  }
  # Cholesky factor of Sigma for fast N(0, Sigma) draws.
  chol_S <- tryCatch(chol(Sigma + diag(1e-8, p, p)),
                      error = function(e) NULL)
  if (is.null(chol_S)) return(NA_real_)

  set.seed(seed, kind = "Mersenne-Twister")
  Z <- matrix(stats::rnorm(M_mc * p), nrow = M_mc, ncol = p)
  draws <- Z %*% chol_S         # M_mc x p
  # Test statistic: max over rows of A of (negative slack). I.e.,
  # T(theta) = max(0, -A %*% theta).
  At <- t(A)
  stats_vec <- apply(draws, 1L, function(z) {
    slack <- as.numeric(z %*% At)
    max(c(0, -slack))
  })
  stats::quantile(stats_vec, probs = 1 - alpha, names = FALSE)
}


# ARP test for moment inequalities.
# Implements the conditional / least-favorable critical-value approach:
# T = max(0, -min(A %*% theta_hat)), then compare against the LF cv.
#' @keywords internal
#' @noRd
.testmechs_arp_test <- function(theta_hat, Sigma, A, alpha = 0.05,
                                  M_mc = 5000L, seed = 1L) {
  slack <- as.numeric(A %*% theta_hat)
  T_stat <- max(c(0, -slack))
  cv <- .testmechs_arp_cv(Sigma, A, alpha = alpha,
                           M_mc = M_mc, seed = seed)
  if (is.na(cv)) {
    return(list(method = "ARP", reject = NA, test_stat = T_stat,
                 cv = NA_real_, pval = NA_real_,
                 message = "Sigma Cholesky failed"))
  }
  reject <- T_stat > cv
  # Bootstrap-based p-value: fraction of draws whose T exceeds T_stat.
  pval <- if (is.na(cv)) NA_real_
          else {
            set.seed(seed, kind = "Mersenne-Twister")
            p <- ncol(A)
            chol_S <- chol(Sigma + diag(1e-8, p, p))
            Z <- matrix(stats::rnorm(M_mc * p), nrow = M_mc, ncol = p)
            draws <- Z %*% chol_S
            At <- t(A)
            ts <- apply(draws, 1L, function(z) {
              slack_z <- as.numeric(z %*% At)
              max(c(0, -slack_z))
            })
            mean(ts >= T_stat)
          }
  list(method = "ARP", reject = reject, test_stat = T_stat,
       cv = cv, pval = pval, M_mc = M_mc)
}


# FSST cone-based test. Implementation strategy:
#   1. Compute T_n = N_n^{-1/2} * min_{theta in polytope} (theta_hat - theta)' Sigma^{-1} (theta_hat - theta)
#      where N_n is the sample size.
#   2. Bootstrap-resample to compute the cone statistic on each replicate.
#   3. Reject if T_n exceeds the (1 - alpha)-quantile of the bootstrap distribution.
#
# Simplified version here: the inner projection uses quadprog (same as
# CS), the bootstrap reuses the supplied Sigma to draw N(0, Sigma)
# replicates as a Gaussian approximation to the full bootstrap.
# Production-grade FSST (lpinfer::lp_inference) requires a secondary
# bootstrap over the data — left as future work.
#' @keywords internal
#' @noRd
.testmechs_fsst_test <- function(theta_hat, Sigma, A, alpha = 0.05,
                                   B = 500L, seed = 1L) {
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("FSST requires quadprog.")
  }
  # Use the CS engine's projection routine as the cone-distance helper.
  base_cs <- .testmechs_cs_test(theta_hat, Sigma, A, alpha = alpha)
  T_n <- base_cs$test_stat

  # Gaussian-approx bootstrap distribution of T_n under the null.
  p <- ncol(A)
  chol_S <- tryCatch(chol(Sigma + diag(1e-8, p, p)),
                      error = function(e) NULL)
  if (is.null(chol_S)) {
    return(list(method = "FSST", reject = NA, test_stat = T_n,
                 cv = NA_real_, pval = NA_real_,
                 message = "Sigma Cholesky failed"))
  }
  set.seed(seed, kind = "Mersenne-Twister")
  ts <- numeric(B)
  for (b in seq_len(B)) {
    z <- as.numeric(stats::rnorm(p) %*% chol_S)
    # Cone distance: shift theta_hat by z, re-project.
    cs_b <- tryCatch(
      .testmechs_cs_test(theta_hat + z, Sigma, A, alpha = alpha),
      error = function(e) list(test_stat = NA_real_))
    ts[b] <- cs_b$test_stat
  }
  ts <- ts[!is.na(ts)]
  if (length(ts) < 2L) {
    return(list(method = "FSST", reject = NA, test_stat = T_n,
                 cv = NA_real_, pval = NA_real_,
                 message = "no valid bootstrap draws"))
  }
  cv <- stats::quantile(ts, probs = 1 - alpha, names = FALSE)
  pval <- mean(ts >= T_n)
  reject <- T_n > cv
  list(method = "FSST", reject = reject, test_stat = T_n,
       cv = cv, pval = pval, B = B)
}
