# ============================================================================
# Generic Cox-Shi (2023) moment-inequality test.
#
# Given:
#   theta_hat  -- p-vector of point estimates
#   Sigma      -- p x p covariance (typically from a bootstrap)
#   A          -- m x p inequality matrix; under H0, A %*% theta >= 0
#                 ("theta lies in the polytope")
#   A_eq       -- optional q x p equality matrix; under H0, A_eq %*% theta = b_eq
#   b_eq       -- q-vector of equality RHS
#
# The test:
#   1. Solve QP: minimise (theta - theta_hat)' Sigma^{-1} (theta - theta_hat)
#      subject to A %*% theta >= 0 and A_eq %*% theta = b_eq.
#      Solution: theta_proj.
#   2. Compute T_stat = (theta_hat - theta_proj)' Sigma^{-1} (theta_hat - theta_proj).
#   3. Approximate the asymptotic distribution under H0 as
#      chi-squared with dof equal to the rank-deficiency of A_active
#      (the rows of A whose corresponding constraint is binding at
#      theta_proj). For an empty active set, T_stat = 0.
#   4. p-value = 1 - pchisq(T_stat, dof).
#
# Reference: Cox, G. and Shi, X. (2023). "A simple uniformly valid test
# for inequalities." Quantitative Economics.
#
# For the mediation-specific polytope construction (which depends on
# the (D, M, Y) structure), see testmechs_sharp_null.R (work in
# progress; the generic CS test here is the engine that polytope-
# specific code will call into).
# ============================================================================


#' Run a generic Cox-Shi moment-inequality test
#'
#' Given a point-estimate vector, its covariance, and an inequality
#' polytope, computes the CS test statistic, its asymptotic chi-square
#' approximation, and a p-value for the null `H0: A %*% theta >= 0`.
#'
#' This is the engine that the TestMechs sharp-mediation test plugs
#' into; the polytope construction is specific to the mediation
#' problem and lives elsewhere.
#'
#' @param theta_hat Numeric vector of length p.
#' @param Sigma     p x p covariance matrix (positive-definite).
#' @param A         m x p inequality matrix.
#' @param A_eq      Optional q x p equality matrix.
#' @param b_eq      Optional q-vector of equality RHS.
#' @param alpha     Test level. Default `0.05`.
#' @return A list with `test_stat`, `cv` (chi-square critical value at
#'   `alpha`), `pval`, `dof`, `theta_proj`, `reject` (logical).
#'
#' @keywords internal
#' @noRd
.testmechs_cs_test <- function(theta_hat, Sigma, A,
                                  A_eq = NULL, b_eq = NULL,
                                  alpha = 0.05) {
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("The CS test requires the 'quadprog' package. ",
         "Install with install.packages('quadprog').")
  }
  p <- length(theta_hat)
  if (!is.matrix(Sigma) || nrow(Sigma) != p || ncol(Sigma) != p) {
    stop("Sigma must be a p x p matrix.")
  }
  if (!is.matrix(A) || ncol(A) != p) {
    stop("A must have p columns matching length(theta_hat).")
  }

  # Solve the QP: min 0.5 * (theta - theta_hat)' Sigma^{-1} (theta - theta_hat)
  # subject to A %*% theta >= 0 and A_eq %*% theta = b_eq.
  # quadprog::solve.QP takes Dmat (positive definite), dvec, Amat (transposed),
  # bvec, meq (number of equality constraints, which must be the FIRST meq).
  Sigma_inv <- tryCatch(solve(Sigma), error = function(e) MASS::ginv(Sigma))
  # quadprog requires Dmat positive-definite. If Sigma is singular,
  # use a tiny ridge.
  Dmat <- Sigma_inv
  Dmat <- (Dmat + t(Dmat)) / 2          # symmetrise
  eigs <- eigen(Dmat, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigs <= 0)) {
    ridge <- 1e-8 * max(abs(eigs))
    Dmat <- Dmat + diag(ridge, p, p)
  }
  dvec <- as.numeric(Sigma_inv %*% theta_hat)
  if (!is.null(A_eq) && !is.null(b_eq)) {
    Amat <- t(rbind(A_eq, A))
    bvec <- c(b_eq, rep(0, nrow(A)))
    meq  <- nrow(A_eq)
  } else {
    Amat <- t(A)
    bvec <- rep(0, nrow(A))
    meq  <- 0L
  }
  qp <- tryCatch(
    quadprog::solve.QP(Dmat = Dmat, dvec = dvec,
                       Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL)
  if (is.null(qp)) {
    # Infeasible or degenerate; treat T = 0 (the null trivially holds).
    return(list(test_stat = 0, cv = NA_real_, pval = 1,
                dof = 0L, theta_proj = theta_hat,
                reject = FALSE,
                infeasible = TRUE))
  }
  theta_proj <- qp$solution
  diff <- theta_hat - theta_proj
  T_stat <- as.numeric(t(diff) %*% Sigma_inv %*% diff)
  T_stat <- max(0, T_stat)

  # Active constraint set: rows of A where A %*% theta_proj is "close to" 0
  # (binding). The dof of the chi-square approximation equals the number
  # of binding constraints (intuitively: each binding moment is one
  # dimension along which the projection was forced).
  slack <- as.numeric(A %*% theta_proj)
  active <- abs(slack) < 1e-6
  dof <- max(1L, sum(active))   # avoid pchisq(., 0) = 0

  cv   <- stats::qchisq(1 - alpha, df = dof)
  pval <- stats::pchisq(T_stat, df = dof, lower.tail = FALSE)
  reject <- T_stat > cv

  list(test_stat = T_stat,
       cv        = cv,
       pval      = pval,
       dof       = as.integer(dof),
       theta_proj = theta_proj,
       reject    = reject,
       n_active_constraints = sum(active),
       infeasible = FALSE)
}
