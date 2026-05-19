# ============================================================================
# Bootstrap kernel for TestMechs partial densities.
#
# Computes a (B x dim_beta) matrix of bootstrap-replicated beta.obs
# vectors, where dim_beta = 2 * K * d_y and beta.obs[i] is the empirical
# partial density vector on the i-th resampled dataset.
#
# Two implementations:
#   .testmechs_bootstrap_r     -- pure-R reference (always works)
#   .testmechs_bootstrap_cuda  -- CUDA kernel (when nvcc is available)
#
# Both produce IDENTICAL outputs at the same RNG seed (Mersenne-Twister
# in R; cuRAND XORWOW or Mersenne-Twister with matching seeds on GPU).
# This makes parity testing trivial.
#
# Bootstrap variants:
#   - "nonparametric": resample n indices with replacement
#   - "bayes":         draw Dirichlet(1, ..., 1) weights, weighted reduction
# ============================================================================


#' Bootstrap-replicate the partial density vector
#'
#' For each bootstrap draw b = 1..B:
#'   1. Resample n observations (with replacement, or with Dirichlet weights).
#'   2. Compute beta_b = (P(Y = y, M = m | D = d) for all (y, m, d)).
#'
#' Returns a B x (2 * K * d_y) matrix.
#'
#' @param d Integer vector. Binary treatment (0/1).
#' @param m Integer vector. Mediator levels (1..K).
#' @param y Integer vector. Outcome bins (1..d_y).
#' @param B Integer. Number of bootstrap draws.
#' @param method One of `"nonparametric"` (default) or `"bayes"`.
#' @param seed Integer. RNG seed.
#' @param backend One of `"auto"`, `"r"`, `"cuda"`. `"auto"` picks
#'   `"cuda"` when available.
#' @return A `B x (2 * K * d_y)` numeric matrix.
#'
#' @keywords internal
#' @noRd
.testmechs_bootstrap <- function(d, m, y, B,
                                    method  = c("nonparametric", "bayes"),
                                    seed    = 1L,
                                    backend = "auto") {
  method <- match.arg(method)
  use_cuda <- backend %in% c("auto", "cuda") &&
              isTRUE(tryCatch(didgpu_has_cuda_support(),
                               error = function(e) FALSE))
  if (use_cuda && backend %in% c("auto", "cuda")) {
    return(.testmechs_bootstrap_cuda(d, m, y, B, method, seed))
  }
  .testmechs_bootstrap_r(d, m, y, B, method, seed)
}


# Pure-R reference implementation.
#' @keywords internal
#' @noRd
.testmechs_bootstrap_r <- function(d, m, y, B, method, seed) {
  n  <- length(d)
  K  <- max(m, na.rm = TRUE)
  dy <- max(y, na.rm = TRUE)
  dim_beta <- 2L * K * dy
  out <- matrix(0, nrow = B, ncol = dim_beta)
  # Force Mersenne-Twister and seed deterministically.
  set.seed(as.integer(seed), kind = "Mersenne-Twister")
  for (b in seq_len(B)) {
    if (method == "nonparametric") {
      idx <- sample.int(n, n, replace = TRUE)
      d_b <- d[idx]; m_b <- m[idx]; y_b <- y[idx]
    } else {
      # Bayesian bootstrap: draw Dirichlet(1,...,1) weights = normalized
      # Exp(1) draws, then compute weighted partial densities.
      w <- stats::rexp(n, 1)
      w <- w / sum(w) * n   # rescale so sum = n (same total as nonparam)
      d_b <- d; m_b <- m; y_b <- y
    }
    if (method == "nonparametric") {
      pd <- .testmechs_partial_density(d_b, m_b, y_b)
      out[b, ] <- pd$beta
    } else {
      # Weighted partial density: P_hat(Y=y, M=m | D=d) =
      #   sum_i w_i * 1[d_i=d, m_i=m, y_i=y] / sum_i w_i * 1[d_i=d]
      n_per_d <- c(`0` = sum(w[d == 0L]), `1` = sum(w[d == 1L]))
      for (dd in c(0L, 1L)) {
        n_dd <- as.numeric(n_per_d[as.character(dd)])
        if (is.na(n_dd) || n_dd == 0) next
        for (mm in seq_len(K)) {
          for (yy in seq_len(dy)) {
            mask <- d == dd & m == mm & y == yy
            sumw <- sum(w[mask])
            pos  <- dd * K * dy + (mm - 1L) * dy + yy
            out[b, pos] <- sumw / n_dd
          }
        }
      }
    }
  }
  out
}


# CUDA bootstrap (stub for now; the .cu kernel lives in
# src/cuda_testmechs_bootstrap.cu and is wired through an Rcpp helper).
# For now this falls back to the R implementation; the wiring will be
# completed in a follow-up session once the .cu file is fleshed out.
#' @keywords internal
#' @noRd
.testmechs_bootstrap_cuda <- function(d, m, y, B, method, seed) {
  # TODO: call didgpu_cuda_testmechs_bootstrap_r when wired in.
  # Until then, fall back to the R implementation.
  message("[testmechs] CUDA bootstrap not yet wired; falling back to R.")
  .testmechs_bootstrap_r(d, m, y, B, method, seed)
}


# Bootstrap covariance of beta.obs. Given the B x dim_beta matrix from
# .testmechs_bootstrap, returns a (dim_beta x dim_beta) covariance.
# This is the input to the CS / ARP / FSST moment-inequality tests.
#' @keywords internal
#' @noRd
.testmechs_sigma <- function(boot_matrix) {
  if (nrow(boot_matrix) < 2L) {
    stop("Cannot compute Sigma from fewer than 2 bootstrap draws.")
  }
  stats::cov(boot_matrix)
}
