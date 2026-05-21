# ============================================================================
# Multi-level mediator (K > 2) extension of the TestMechs sharp-null test.
#
# With K levels of M and d_y bins of Y, the observed partial-density
# vector beta has length 2 * K * d_y.
#
# Under no-defiers, compliance types are pairs (M(0) = a, M(1) = b)
# with a <= b in some ordering (defiers excluded). For an ordered M
# (e.g. dose 0/1/2/3) we have K * (K + 1) / 2 types under no-defiers;
# we parametrise type probabilities pi_ab and the within-type
# Y-density f_ab(y).
#
# Under sharp full mediation: f_ab is the same for a == b (always-
# stayers see the same Y distribution under both D values). For a < b
# (compliers in some direction), the same restriction: f_ab(y | D = 0)
# = f_ab(y | D = 1). This is the "sharp" part — no direct effect of D
# on Y conditional on M.
#
# Implementation: build the polytope (A_eq, b_eq) for general K
# parallel to the binary-M case. Plug into .testmechs_cs_test.
#
# Status: v1 supports general K with the same "no-defiers" assumption
# as the binary case. For ORDERED mediators with an explicit ordering,
# additional shape constraints can be added later.
# ============================================================================


# Build the polytope for the multi-level case (K >= 2, ordered no-defiers).
#' @keywords internal
#' @noRd
.testmechs_build_polytope_multi_m <- function(beta_hat, K, d_y) {
  # Compliance types under no-defiers (assuming an ordered M):
  #   M(0) = a, M(1) = b with a <= b, a, b in 1..K.
  # Number of types T = K * (K + 1) / 2.
  type_pairs <- expand.grid(a = 1L:K, b = 1L:K)
  type_pairs <- type_pairs[type_pairs$a <= type_pairs$b, ]
  T_n <- nrow(type_pairs)

  # Parameters per type:
  #   theta_ab_y_d : P(type = (a, b), Y = y, D = d)
  # Under sharp null: theta_ab_y_0 = theta_ab_y_1 for every (a, b, y).
  # Total params: 2 * T_n * d_y.
  # Indexing: position(t, y, d) = d * T_n * d_y + (t - 1) * d_y + y.
  p <- 2L * T_n * d_y
  pos <- function(t, y, d) d * T_n * d_y + (t - 1L) * d_y + y

  # Observed cell beta[y, m, d] is the sum over types that produce M = m
  # under D = d. For type (a, b): if D = 0 then M = a, if D = 1 then M = b.
  # So beta[y, m, d=0] = sum over types with a == m of theta_ab_y_0
  #    beta[y, m, d=1] = sum over types with b == m of theta_ab_y_1
  A_eq <- matrix(0, nrow = 2L * K * d_y, ncol = p)
  b_eq <- numeric(2L * K * d_y)
  bidx <- function(yy, mm, dd) dd * K * d_y + (mm - 1L) * d_y + yy
  row <- 0L
  for (dd in c(0L, 1L)) {
    for (mm in 1L:K) {
      for (yy in 1L:d_y) {
        row <- row + 1L
        # Types contributing: if dd == 0, types with a == mm;
        #                     if dd == 1, types with b == mm.
        contrib <- if (dd == 0L) which(type_pairs$a == mm)
                   else which(type_pairs$b == mm)
        for (t in contrib) {
          A_eq[row, pos(t, yy, dd)] <- 1
        }
        b_eq[row] <- beta_hat[bidx(yy, mm, dd)]
      }
    }
  }

  # Sharp-null equalities: theta_ab_y_0 == theta_ab_y_1 for all (a, b, y).
  A_eq_null <- matrix(0, nrow = T_n * d_y, ncol = p)
  b_eq_null <- numeric(T_n * d_y)
  row <- 0L
  for (t in 1L:T_n) {
    for (yy in 1L:d_y) {
      row <- row + 1L
      A_eq_null[row, pos(t, yy, 0L)] <-  1
      A_eq_null[row, pos(t, yy, 1L)] <- -1
    }
  }

  # Inequality: all theta >= 0.
  A <- diag(p)

  list(A    = A,
       A_eq = rbind(A_eq, A_eq_null),
       b_eq = c(b_eq, b_eq_null),
       p    = p, K = K, d_y = d_y, T_n = T_n,
       type_pairs = type_pairs)
}


# Multi-level CS-test wrapper.
#' @keywords internal
#' @noRd
.testmechs_test_sharp_null_multi_m_cs <- function(
    d, m, y, B, num_Ybins, seed, backend, alpha) {
  y_bin <- .testmechs_bin_y(y, num_Ybins)
  d_int <- as.integer(d)
  m_int <- as.integer(m)
  if (min(m_int) == 0L) m_int <- m_int + 1L
  K <- max(m_int)
  if (K < 2L) stop("Need K >= 2 levels of mediator.")

  point <- .testmechs_partial_density(d_int, m_int, y_bin)
  beta_hat <- point$beta
  d_y <- point$d_y

  boot_mat <- .testmechs_bootstrap(d_int, m_int, y_bin, B = B,
                                    method = "nonparametric",
                                    seed = seed, backend = backend)
  Sigma <- .testmechs_sigma(boot_mat)

  poly <- .testmechs_build_polytope_multi_m(beta_hat, K, d_y)

  theta_hat <- tryCatch(
    as.numeric(MASS::ginv(poly$A_eq) %*% poly$b_eq),
    error = function(e) rep(0, poly$p))
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
    method = "CS", reject = cs$reject,
    test_stat = cs$test_stat, cv = cs$cv, pval = cs$pval,
    dof = cs$dof, n_active = cs$n_active_constraints,
    B = B, K = K, d_y = d_y, backend = backend,
    note = sprintf("multi-level M (K=%d) sharp-null test via CS", K)
  )
}
