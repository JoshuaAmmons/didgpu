# ============================================================================
# TestMechs: sharp test of full mediation (Kwon & Roth 2026, ReStud).
#
# Methodology: cross-sectional moment-inequality testing of the "sharp
# null of full mediation" -- given a binary treatment D, a discrete
# mediator M (K levels), and an outcome Y (binned to d_y levels), test
# whether D affects Y ONLY through M (i.e., Y(1, m) = Y(0, m) for all m).
#
# This is fundamentally DIFFERENT from the panel-DiD work in
# `didgpu()` / `didgpu_fect()`. There's no time dimension, no F_g, no
# event-time concept. The data is i.i.d. observations of (D, M, Y).
#
# Public API (mirrors the reference TestMechs package):
#   didgpu_test_sharp_null(df, d, m, y, method, B, ...)
#   didgpu_lb_frac_affected(df, d, m, y, B, ...)
#
# Three test methods:
#   "CS"   -- Cox-Shi (2023). QP + chi-sq cv. The simplest path.
#   "ARP"  -- Andrews-Roth-Pakes conditional/hybrid LP test.
#   "FSST" -- Fang-Santos-Shaikh-Torgovitsky. Cone test; needs 2x bootstrap.
#
# GPU acceleration story:
#   The bottleneck is the bootstrap loop computing the empirical
#   partial densities beta.obs = P(Y = y, M = m | D = d) on each
#   resample. With B = 500-2000 and n = O(10K), this dominates wall
#   time. CUDA kernel: one block per bootstrap draw, atomic-add into a
#   (B x dim_beta) matrix; or cuRAND multinomial draws + segmented
#   reduction. cuBLAS syrk forms Sigma in one call.
#
#   The LP/QP solvers for the moment-inequality step (Cox-Shi, ARP,
#   FSST) stay on CPU -- the matrices are tiny (K = 2..10, d_y = 2..10
#   means beta.obs has ~50-400 entries) and the existing
#   osqp/Rglpk/lpinfer solvers are well-tuned.
#
# STATUS: scaffolded. The public API is exposed; the actual estimator
# logic is stubbed. Real implementation work for the next session.
# ============================================================================


#' Sharp test of full mediation
#'
#' Tests whether a binary treatment `D` affects an outcome `Y` only
#' through a discrete mediator `M`. The null hypothesis is sharp full
#' mediation: `Y(1, m) = Y(0, m)` for every level `m`. Implemented via
#' three moment-inequality testing procedures.
#'
#' @param df A data.frame. Cross-sectional; one row per observation.
#' @param d Character. Column name of the binary treatment (0/1).
#' @param m Character. Column name of the discrete mediator.
#' @param y Character. Column name of the outcome (continuous outcomes
#'   are binned to `num_Ybins` quantile bins).
#' @param method One of `"CS"` (Cox-Shi), `"ARP"` (Andrews-Roth-Pakes),
#'   or `"FSST"` (Fang-Santos-Shaikh-Torgovitsky). Default `"CS"`.
#' @param B Integer. Number of bootstrap resamples for variance
#'   estimation. Default `500L`.
#' @param num_Ybins Integer or NULL. For continuous `y`, the number of
#'   quantile bins to use. Default `5L`. Ignored when `y` is already
#'   discrete.
#' @param cluster Character or NULL. Column name for cluster bootstrap.
#' @param reg_formula Optional one-sided formula for regression-adjusted
#'   partial densities.
#' @param alpha Numeric. Test level. Default `0.05`.
#' @param seed Integer. RNG seed.
#' @param backend One of `"auto"`, `"r"`, `"cuda"`. `"auto"` picks
#'   `"cuda"` when available.
#'
#' @return A list with `reject` (logical), `test_stat` (numeric),
#'   `cv` (critical value), `pval` (numeric), `method` (character).
#'
#' @references
#' Kwon, S. and Roth, J. (2026). "(Empirical) Bayes approaches to
#' parallel trends." *Review of Economic Studies*, forthcoming.
#'
#' @examples
#' \dontrun{
#' # NOT YET IMPLEMENTED. The example below shows the planned API.
#' df <- data.frame(D = sample(0:1, 1000, TRUE),
#'                  M = sample(1:3, 1000, TRUE),
#'                  Y = rnorm(1000))
#' res <- didgpu_test_sharp_null(df, "D", "M", "Y", method = "CS",
#'                                 B = 500L, seed = 1L)
#' res$pval
#' }
#' @export
didgpu_test_sharp_null <- function(
    df, d, m, y,
    method      = c("CS", "ARP", "FSST"),
    B           = 500L,
    num_Ybins   = 5L,
    cluster     = NULL,
    reg_formula = NULL,
    alpha       = 0.05,
    seed        = 1L,
    backend     = "auto") {

  method <- match.arg(method)
  stopifnot(is.data.frame(df) || data.table::is.data.table(df))
  for (nm in c("d", "m", "y")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v)) {
      stop("`", nm, "` must be a single non-empty character column name.")
    }
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  stopifnot(all(df[[d]] %in% c(0, 1)) || all(is.na(df[[d]])))
  B <- as.integer(B)
  stopifnot(B >= 1L)

  # Dispatch: all three test methods (CS / ARP / FSST) x all K >= 2.
  # The polytope construction is shared between binary and multi-M
  # via .testmechs_build_polytope_{binary,multi}_m; the test method
  # plugs the shared (theta_hat, Sigma_theta, A) into its own engine.
  K <- length(unique(stats::na.omit(df[[m]])))
  if (K < 2L) {
    stop("Mediator must have at least 2 levels.")
  }
  y_bin <- .testmechs_bin_y(df[[y]], num_Ybins %||% 5L)
  d_int <- as.integer(df[[d]])
  m_int <- as.integer(df[[m]])
  if (min(m_int) == 0L) m_int <- m_int + 1L
  point <- .testmechs_partial_density(d_int, m_int, y_bin)
  boot_mat <- .testmechs_bootstrap(d_int, m_int, y_bin, B = B,
                                    method = "nonparametric",
                                    seed = seed, backend = backend)
  Sigma_beta <- .testmechs_sigma(boot_mat)
  poly <- if (K == 2L) {
    .testmechs_build_polytope_binary_m(point$beta, point$d_y)
  } else {
    .testmechs_build_polytope_multi_m(point$beta, K, point$d_y)
  }
  theta_hat <- tryCatch(
    as.numeric(MASS::ginv(poly$A_eq) %*% poly$b_eq),
    error = function(e) rep(0, poly$p))
  J <- tryCatch(MASS::ginv(poly$A_eq[1L:length(point$beta), , drop = FALSE]),
                error = function(e) NULL)
  Sigma_theta <- if (is.null(J)) diag(poly$p) else J %*% Sigma_beta %*% t(J)
  Sigma_theta <- (Sigma_theta + t(Sigma_theta)) / 2

  test_res <- switch(method,
    "CS"   = .testmechs_cs_test(theta_hat, Sigma_theta, poly$A, alpha = alpha),
    "ARP"  = .testmechs_arp_test(theta_hat, Sigma_theta, poly$A, alpha = alpha,
                                   seed = seed),
    "FSST" = .testmechs_fsst_test(theta_hat, Sigma_theta, poly$A, alpha = alpha,
                                    B = B, seed = seed)
  )
  # The CS engine's return list doesn't carry a `method` field; tag it
  # consistently so downstream code can introspect.
  test_res$method <- method
  c(test_res, list(B = B, K = K, d_y = point$d_y, backend = backend,
                    note = sprintf("K = %d sharp-null via %s", K, method)))
}


#' Sharp lower bound on the fraction of always-takers affected
#'
#' Given a binary treatment `D`, a discrete mediator `M`, and an
#' outcome `Y`, computes a sharp lower bound on the fraction of
#' always-takers (units with `M(1) = M(0)`) whose outcome is moved by
#' the treatment. The bound is identified from the total-variation
#' distance between treated and control conditional distributions of
#' `Y` given `M`.
#'
#' @inheritParams didgpu_test_sharp_null
#' @param B Number of bootstrap resamples for the CI. Default `500L`.
#' @param at_group Optional character: name of an always-taker subgroup
#'   column. NULL pools over all groups.
#' @return A list with `lb` (lower bound), `ci_low`, `ci_high`,
#'   `method = "TV"`.
#'
#' @examples
#' \dontrun{
#' # NOT YET IMPLEMENTED. Planned API:
#' df <- data.frame(D = sample(0:1, 1000, TRUE),
#'                  M = sample(1:3, 1000, TRUE),
#'                  Y = rnorm(1000))
#' lb <- didgpu_lb_frac_affected(df, "D", "M", "Y", B = 500L)
#' lb$lb; lb$ci_low; lb$ci_high
#' }
#' @export
didgpu_lb_frac_affected <- function(
    df, d, m, y,
    B           = 500L,
    at_group    = NULL,
    num_Ybins   = 5L,
    cluster     = NULL,
    reg_formula = NULL,
    alpha       = 0.05,
    seed        = 1L,
    backend     = "auto") {

  stopifnot(is.data.frame(df) || data.table::is.data.table(df))
  for (nm in c("d", "m", "y")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v)) {
      stop("`", nm, "` must be a single non-empty character column name.")
    }
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  .testmechs_not_implemented("lb_frac_affected", NULL)
}


# Sentinel used until the real implementation lands.
.testmechs_not_implemented <- function(fn, method) {
  msg <- sprintf(
    paste0("didgpu_%s is scaffolded but not yet implemented.\n",
            "Planned implementation roadmap:\n",
            "  1. GPU bootstrap kernel: cuRAND multinomial draws +\n",
            "     atomic-add reduction to build the (B x dim_beta) matrix\n",
            "     of partial-density estimates in one launch.\n",
            "  2. cuBLAS syrk to form Sigma.obs (covariance of beta).\n",
            "  3. CPU LP/QP step%s using HonestDiD / osqp / lpinfer\n",
            "     (the matrices are tiny -- GPU here would be wasted).\n",
            "  4. Cross-validation against the reference TestMechs\n",
            "     package on the baranov_data example.\n",
            "Reference: TestMechs (Kwon & Roth 2026, ReStud)\n",
            "  https://github.com/jonathandroth/TestMechs"),
    fn,
    if (!is.null(method)) sprintf(" (method = '%s')", method) else ""
  )
  stop(msg, call. = FALSE)
}


# ----------------------------------------------------------------------------
# Helper: bin a continuous outcome y into d_y quantile bins. Used by
# both test_sharp_null and lb_frac_affected before the partial-density
# computation.
# ----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.testmechs_bin_y <- function(y, n_bins) {
  if (length(unique(stats::na.omit(y))) <= n_bins) {
    # Already discrete enough: just integer-encode.
    return(as.integer(factor(y, levels = sort(unique(y)))))
  }
  quants <- stats::quantile(y, probs = seq(0, 1, length.out = n_bins + 1L),
                              na.rm = TRUE, type = 7L)
  # Ensure unique breakpoints (ties get merged into fewer bins).
  quants <- unique(quants)
  if (length(quants) <= 1L) {
    return(rep(1L, length(y)))
  }
  out <- as.integer(cut(y, breaks = quants, include.lowest = TRUE,
                         labels = FALSE))
  out
}


# Compute the empirical partial-density vector beta.obs from a data
# frame, given d, m, y columns. The convention is the stacked vector:
#   beta = [P(Y = y, M = m | D = 0) for all (y, m),
#           P(Y = y, M = m | D = 1) for all (y, m)]
# of length 2 * K * d_y, where K = #distinct values of m and d_y =
# #distinct values of y.
#
# Returns the beta vector and the cell-count breakdown (useful for
# the bootstrap reductions later).
#' @keywords internal
#' @noRd
.testmechs_partial_density <- function(d_vec, m_vec, y_vec) {
  stopifnot(length(d_vec) == length(m_vec),
            length(m_vec) == length(y_vec))
  d_vec <- as.integer(d_vec)
  m_vec <- as.integer(m_vec)
  y_vec <- as.integer(y_vec)
  K   <- max(m_vec, na.rm = TRUE)
  d_y <- max(y_vec, na.rm = TRUE)
  beta <- numeric(2 * K * d_y)
  n_per_d <- table(d_vec)
  for (dd in c(0L, 1L)) {
    if (is.na(n_per_d[as.character(dd)])) next
    n_dd <- as.integer(n_per_d[as.character(dd)])
    for (mm in seq_len(K)) {
      for (yy in seq_len(d_y)) {
        count <- sum(d_vec == dd & m_vec == mm & y_vec == yy, na.rm = TRUE)
        # Index into beta: block dd of K * d_y, then position (m-1)*d_y + y.
        pos <- dd * K * d_y + (mm - 1L) * d_y + yy
        beta[pos] <- count / n_dd
      }
    }
  }
  list(beta = beta, K = K, d_y = d_y)
}
