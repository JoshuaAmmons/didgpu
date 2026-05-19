# ============================================================================
# Sharp null of full mediation: implementation for the binary-mediator
# case (M in {0, 1}), using the Cox-Shi (2023) moment-inequality test.
#
# Under the sharp null Y(1, m) = Y(0, m) for all m, with binary M and
# binned Y (d_y bins), the implied moment restrictions can be written
# as a polytope on a vector of compliance-type probabilities:
#
#   theta = (theta_at_y, theta_nt_y, theta_co_y_d : y in 1..d_y, d in {0,1})
#   where
#     theta_at_y    = P(type = always-taker, Y = y)
#     theta_nt_y    = P(type = never-taker,  Y = y)
#     theta_co_y_d  = P(type = complier, Y = y, D = d)
#
# Under no-defiers (a standard assumption; full mediation is even
# tighter), the observed cell probabilities P(Y = y, M = m, D = d) /
# P(D = d) factor through these compliance-type quantities:
#
#   beta[y, m=1, d=0] = pi_at * f_at(y)  (always-takers under D=0 have M=1)
#   beta[y, m=0, d=0] = pi_nt * f_nt(y) + pi_co * f_co_0(y)
#   beta[y, m=1, d=1] = pi_at * f_at(y) + pi_co * f_co_1(y)
#   beta[y, m=0, d=1] = pi_nt * f_nt(y)
#
# where f_*(y) = theta_*_y / pi_* are conditional Y densities within
# each type, and pi_* are type proportions. Under sharp null full
# mediation: f_co_0(y) = f_co_1(y), so the "complier delta" vanishes.
#
# This implementation builds A (the inequality matrix saying every
# theta_* >= 0) and A_eq (the equality matrix mapping theta to beta)
# and hands them to .testmechs_cs_test for the actual test.
#
# Limitations of this initial cut:
#   - Binary M only (K = 2). General K is on the roadmap.
#   - No regression adjustment (reg_formula path).
#   - No cluster bootstrap (uses i.i.d. nonparametric bootstrap).
#   - Default bootstrap; the analytic-variance path for CS requires
#     more delicate moment-derivative computations.
# ============================================================================


#' Build the polytope (A, A_eq, b_eq) for the binary-M full-mediation test
#'
#' theta layout (length 2 * d_y + 2 + 2 * d_y = 4 * d_y + 2):
#'   theta[1:d_y]                    = theta_at_y  (always-takers, Y = y)
#'   theta[(d_y+1):(2*d_y)]          = theta_nt_y  (never-takers, Y = y)
#'   theta[(2*d_y+1):(3*d_y)]        = theta_co_y  (compliers, Y = y) under null
#'                                                  (no D-dependence)
#'   theta[(3*d_y+1):(4*d_y)]        = nuisance: total per-D mass
#'                                    (collapses to identity under no
#'                                    regression adjustment)
#'   theta[4*d_y + 1] = pi_at, theta[4*d_y + 2] = pi_nt (type masses)
#'
#' Actually the cleanest formulation: parametrise by (theta_at_y, theta_nt_y,
#' theta_co_y_0, theta_co_y_1) — 4 * d_y free vars. Inequality: every
#' theta >= 0. Equalities: 4 * d_y mappings to observed beta cells +
#' two normalisation equalities (theta_at sums to pi_at; etc).
#'
#' Under sharp null: theta_co_y_0 = theta_co_y_1 for each y, which is
#' added as d_y additional equality constraints.
#'
#' @keywords internal
#' @noRd
.testmechs_build_polytope_binary_m <- function(beta_hat, d_y) {
  # beta layout (as built by .testmechs_partial_density):
  #   pos = dd * K * d_y + (mm - 1) * d_y + yy
  # with K = 2, mm in {1, 2}, dd in {0, 1}, yy in 1..d_y.
  # So beta has length 2 * 2 * d_y = 4 * d_y.
  K <- 2L
  stopifnot(length(beta_hat) == 2L * K * d_y)

  # theta_layout: (theta_at_y, theta_nt_y, theta_co_y_0, theta_co_y_1)
  # Total length p = 4 * d_y. Indices:
  #   at      [y]:  1..d_y
  #   nt      [y]:  d_y + 1 .. 2*d_y
  #   co_d=0  [y]:  2*d_y + 1 .. 3*d_y
  #   co_d=1  [y]:  3*d_y + 1 .. 4*d_y
  p <- 4L * d_y
  idx_at      <- 1L:d_y
  idx_nt      <- (d_y + 1L):(2L * d_y)
  idx_co_0    <- (2L * d_y + 1L):(3L * d_y)
  idx_co_1    <- (3L * d_y + 1L):(4L * d_y)

  # Equality matrix mapping theta to beta cells, factoring through
  # P(d). Let pi_d_0 = sum_{y, m} beta_hat[y, m, 0] = 1; pi_d_1 = 1.
  # We work per-d:
  #   beta[y, m=1, d=0] = theta_at[y]   (always-takers contribute under D=0 to M=1)
  #   beta[y, m=0, d=0] = theta_nt[y] + theta_co_0[y]
  #   beta[y, m=1, d=1] = theta_at[y] + theta_co_1[y]
  #   beta[y, m=0, d=1] = theta_nt[y]
  # That gives 4 * d_y equalities.
  A_eq <- matrix(0, nrow = 4L * d_y, ncol = p)
  b_eq <- numeric(4L * d_y)
  row <- 0L
  # beta indexing helper:
  bidx <- function(yy, mm, dd) dd * K * d_y + (mm - 1L) * d_y + yy
  for (y in 1L:d_y) {
    # eqn 1: beta[y, m=1, d=0] = theta_at[y]
    row <- row + 1L
    A_eq[row, idx_at[y]] <- 1
    b_eq[row] <- beta_hat[bidx(y, 2L, 0L)]
    # eqn 2: beta[y, m=0, d=0] = theta_nt[y] + theta_co_0[y]
    row <- row + 1L
    A_eq[row, idx_nt[y]]   <- 1
    A_eq[row, idx_co_0[y]] <- 1
    b_eq[row] <- beta_hat[bidx(y, 1L, 0L)]
    # eqn 3: beta[y, m=1, d=1] = theta_at[y] + theta_co_1[y]
    row <- row + 1L
    A_eq[row, idx_at[y]]   <- 1
    A_eq[row, idx_co_1[y]] <- 1
    b_eq[row] <- beta_hat[bidx(y, 2L, 1L)]
    # eqn 4: beta[y, m=0, d=1] = theta_nt[y]
    row <- row + 1L
    A_eq[row, idx_nt[y]] <- 1
    b_eq[row] <- beta_hat[bidx(y, 1L, 1L)]
  }

  # Sharp null: theta_co_y_0 == theta_co_y_1 for each y.
  A_eq_null <- matrix(0, nrow = d_y, ncol = p)
  b_eq_null <- numeric(d_y)
  for (y in 1L:d_y) {
    A_eq_null[y, idx_co_0[y]] <-  1
    A_eq_null[y, idx_co_1[y]] <- -1
  }

  # Inequality matrix: each theta >= 0 (p inequalities).
  A <- diag(p)
  # A %*% theta >= 0 (the .testmechs_cs_test convention).

  list(A    = A,
       A_eq = rbind(A_eq, A_eq_null),
       b_eq = c(b_eq, b_eq_null),
       p    = p,
       d_y  = d_y,
       theta_indices = list(at = idx_at, nt = idx_nt,
                             co_0 = idx_co_0, co_1 = idx_co_1))
}


# Wire the binary-M case end-to-end: bootstrap beta, build polytope,
# call .testmechs_cs_test.
#' @keywords internal
#' @noRd
.testmechs_test_sharp_null_binary_m_cs <- function(
    d, m, y, B, num_Ybins, seed, backend, alpha) {
  # Discretise Y.
  y_bin <- .testmechs_bin_y(y, num_Ybins)
  d_int <- as.integer(d)
  m_int <- as.integer(m)
  if (max(m_int) > 2L) {
    stop("Binary-M CS test requires M with exactly 2 levels; got ",
         max(m_int), ".")
  }
  # Make M 1-based (sample mediator might be 0/1).
  if (min(m_int) == 0L) m_int <- m_int + 1L

  # Point-estimate beta on the original data.
  point <- .testmechs_partial_density(d_int, m_int, y_bin)
  beta_hat <- point$beta
  d_y <- point$d_y

  # Bootstrap covariance.
  boot_mat <- .testmechs_bootstrap(d_int, m_int, y_bin, B = B,
                                    method = "nonparametric",
                                    seed = seed, backend = backend)
  Sigma <- .testmechs_sigma(boot_mat)

  # Build the polytope projection.
  poly <- .testmechs_build_polytope_binary_m(beta_hat, d_y)

  # Run the CS test on the IMPLIED theta. The theta_hat we plug in is
  # the unconstrained "data-implied" theta — i.e., we solve the
  # *unconstrained* system A_eq %*% theta = b_eq to get theta_hat
  # (via least squares; the system is over-determined under the null
  # because of the sharp-null equalities). Then the QP projects this
  # onto the non-negative orthant.
  theta_hat <- tryCatch(
    as.numeric(MASS::ginv(poly$A_eq) %*% poly$b_eq),
    error = function(e) rep(0, poly$p))
  # The Sigma we need is for theta_hat, not for beta_hat. Approximate
  # via delta method: Sigma_theta = J %*% Sigma %*% t(J), where
  # J = (A_eq^+) is the Moore-Penrose pseudo-inverse of A_eq restricted
  # to the beta-indexing rows.
  # For now, use J = pinv(A_eq[1:length(beta_hat), ]) — the equality
  # rows that map theta to beta. The sharp-null rows have b_eq = 0 so
  # contribute no variance.
  J <- tryCatch(MASS::ginv(poly$A_eq[1L:length(beta_hat), , drop = FALSE]),
                error = function(e) NULL)
  if (is.null(J)) {
    return(list(method = "CS", reject = NA, test_stat = NA_real_,
                 cv = NA_real_, pval = NA_real_, dof = NA_integer_,
                 message = "delta-method failed"))
  }
  Sigma_theta <- J %*% Sigma %*% t(J)

  cs <- .testmechs_cs_test(theta_hat = theta_hat,
                            Sigma = Sigma_theta,
                            A = poly$A,
                            alpha = alpha)
  list(
    method      = "CS",
    reject      = cs$reject,
    test_stat   = cs$test_stat,
    cv          = cs$cv,
    pval        = cs$pval,
    dof         = cs$dof,
    n_active    = cs$n_active_constraints,
    B           = B,
    K           = 2L,
    d_y         = d_y,
    backend     = backend,
    note        = "binary-M sharp-null test via CS; bootstrap-based covariance"
  )
}
