# ============================================================================
# Callaway-Sant'Anna inner estimators: OR, IPW, DR.
#
# All three share the same per-(g, t) cell structure:
#   1. Identify treated units (F_g_unit == g) and control units (depends on
#      control_group: never-treated or not-yet-treated as of time t).
#   2. Compute the outcome change Δy_i = Y_i(t) - Y_i(g - 1) for each unit.
#   3. Apply the method-specific weighting / adjustment.
#   4. Return ATT(g, t), influence function (per unit), and sample sizes.
#
# Method-specific machinery:
#   OR  (outcome regression, Heckman-Ichimura-Todd):
#       Fit lm(Δy ~ X) on controls; predict for treated; ATT = mean(Δy - Δy_hat).
#       Without covariates this reduces to mean(Δy_t) - mean(Δy_c).
#
#   IPW (Abadie 2005):
#       Fit propensity glm(D ~ X, family = binomial); compute weights
#       w_c = p_hat / (1 - p_hat) for controls; ATT = mean(Δy_t) -
#       sum(w_c * Δy_c) / sum(w_c).
#
#   DR  (Sant'Anna & Zhao 2020, "doubly-robust DiD"):
#       Combines OR + IPW; the canonical DR formula:
#         ATT = mean(w1 * (Δy - m(X))) - mean(w0 * (Δy - m(X)))
#       where w1 = D / E[D], w0 = (1 - D) * (p(X) / (1 - p(X))) / E[D],
#       m(X) = E[Δy | X, D = 0] is the OR predictor.
#       Consistent if EITHER the propensity model OR the outcome model
#       is correctly specified — hence "doubly robust."
#
# All three return the influence function so the multiplier bootstrap
# can use it instead of refitting per resample.
# ============================================================================


# Identify control units at time t given the control_group choice.
#' @keywords internal
#' @noRd
.cs_control_units <- function(d, units, F_g_per_unit, g, t, control_group) {
  if (control_group == "never") {
    units[!is.finite(F_g_per_unit)]
  } else if (control_group == "notyet") {
    # Not yet treated as of time t: F_g > t (includes never-treated as
    # a limiting case). Exclude units in the SAME cohort g (those are
    # the treated).
    units[(F_g_per_unit > t) & (F_g_per_unit != g)]
  } else {
    stop("Unknown control_group: ", control_group)
  }
}


# Per-cell OR estimator. Without covariates, it's mean(Δy_t) - mean(Δy_c).
# With covariates, fit lm(Δy ~ X) on controls and predict for treated.
# Returns att, IF (length n_unit), n_treated, n_control.
#' @keywords internal
#' @noRd
.cs_inner_or <- function(delta, D_mask, X_treated, X_control, n_total) {
  n_t <- sum(D_mask)
  n_c <- length(delta) - n_t
  if (n_t == 0L || n_c == 0L) {
    return(list(att = NA_real_, IF = rep(0, n_total),
                 n_treated = n_t, n_control = n_c))
  }
  delta_t <- delta[D_mask]
  delta_c <- delta[!D_mask]

  if (is.null(X_treated) || is.null(X_control) || ncol(X_control) == 0L) {
    # No covariates: simple difference of means.
    m_hat_t <- mean(delta_c, na.rm = TRUE)
    fitted_c <- rep(m_hat_t, n_c)
    fitted_t <- rep(m_hat_t, n_t)
  } else {
    # Fit lm(Δy ~ X) on controls; predict.
    X_c <- cbind(1, X_control)
    qr_c <- tryCatch(qr(X_c), error = function(e) NULL)
    if (is.null(qr_c) || qr_c$rank < ncol(X_c)) {
      # Rank-deficient; fall back to intercept-only.
      m_hat <- mean(delta_c, na.rm = TRUE)
      fitted_c <- rep(m_hat, n_c)
      fitted_t <- rep(m_hat, n_t)
    } else {
      beta <- qr.solve(qr_c, delta_c)
      fitted_c <- as.numeric(X_c %*% beta)
      X_t <- cbind(1, X_treated)
      fitted_t <- as.numeric(X_t %*% beta)
    }
  }
  att <- mean(delta_t - fitted_t)
  # Influence function (unit-level contribution to ATT).
  IF <- numeric(n_total)
  IF[D_mask] <- (delta_t - fitted_t) - att
  # Control units' IF is zero in OR (they only affect the projection,
  # not the ATT directly).
  list(att = att, IF = IF,
       n_treated = as.integer(n_t),
       n_control = as.integer(n_c))
}


# Per-cell IPW estimator (Abadie 2005).
# Without covariates, propensity = n_t / (n_t + n_c) is constant, so
# weights cancel and IPW reduces to OR-without-X (= simple difference).
# With covariates, propensity glm gives unit-specific weights.
#' @keywords internal
#' @noRd
.cs_inner_ipw <- function(delta, D_mask, X, n_total) {
  n_t <- sum(D_mask)
  n_c <- length(delta) - n_t
  if (n_t == 0L || n_c == 0L) {
    return(list(att = NA_real_, IF = rep(0, n_total),
                 n_treated = n_t, n_control = n_c))
  }
  if (is.null(X) || ncol(X) == 0L) {
    # Constant propensity: IPW = simple difference.
    delta_t <- delta[D_mask]
    delta_c <- delta[!D_mask]
    att <- mean(delta_t) - mean(delta_c)
    IF <- numeric(n_total)
    IF[D_mask]  <- delta_t - mean(delta_t) - att / 2
    IF[!D_mask] <- -(delta_c - mean(delta_c)) - att / 2
    return(list(att = att, IF = IF,
                 n_treated = as.integer(n_t),
                 n_control = as.integer(n_c)))
  }
  # Fit propensity score by logistic regression.
  glm_fit <- tryCatch(
    suppressWarnings(stats::glm.fit(x = cbind(1, X),
                                    y = as.integer(D_mask),
                                    family = stats::binomial())),
    error = function(e) NULL)
  if (is.null(glm_fit) || any(is.na(glm_fit$coefficients))) {
    # Logistic failed; fall back to constant propensity.
    return(.cs_inner_ipw(delta, D_mask, X = NULL, n_total = n_total))
  }
  beta <- glm_fit$coefficients
  eta <- as.numeric(cbind(1, X) %*% beta)
  p_hat <- 1 / (1 + exp(-eta))
  # Trim extreme propensities for stability.
  p_hat <- pmin(pmax(p_hat, 0.01), 0.99)
  E_D <- mean(as.integer(D_mask))
  w1 <- as.integer(D_mask) / E_D
  w0 <- (1 - as.integer(D_mask)) * (p_hat / (1 - p_hat)) / E_D
  att <- mean(w1 * delta) - mean(w0 * delta)
  IF <- (w1 * delta) - (w0 * delta) - att
  list(att = att, IF = IF,
       n_treated = as.integer(n_t),
       n_control = as.integer(n_c))
}


# Per-cell DR estimator (Sant'Anna & Zhao 2020).
# Combines OR (outcome model m(X)) + IPW (propensity score p(X)).
#' @keywords internal
#' @noRd
.cs_inner_dr <- function(delta, D_mask, X, n_total) {
  n_t <- sum(D_mask)
  n_c <- length(delta) - n_t
  if (n_t == 0L || n_c == 0L) {
    return(list(att = NA_real_, IF = rep(0, n_total),
                 n_treated = n_t, n_control = n_c))
  }
  if (is.null(X) || ncol(X) == 0L) {
    # No covariates: DR reduces to the simple difference (and IPW).
    return(.cs_inner_ipw(delta, D_mask, X = NULL, n_total = n_total))
  }
  # Outcome model on controls.
  X_c <- cbind(1, X[!D_mask, , drop = FALSE])
  delta_c <- delta[!D_mask]
  qr_c <- tryCatch(qr(X_c), error = function(e) NULL)
  if (is.null(qr_c) || qr_c$rank < ncol(X_c)) {
    m_hat <- rep(mean(delta_c, na.rm = TRUE), n_total)
  } else {
    beta_or <- qr.solve(qr_c, delta_c)
    m_hat <- as.numeric(cbind(1, X) %*% beta_or)
  }
  # Propensity model on all units.
  glm_fit <- tryCatch(
    suppressWarnings(stats::glm.fit(x = cbind(1, X),
                                    y = as.integer(D_mask),
                                    family = stats::binomial())),
    error = function(e) NULL)
  if (is.null(glm_fit) || any(is.na(glm_fit$coefficients))) {
    # Propensity failed; fall back to OR.
    return(.cs_inner_or(delta, D_mask,
                         X_treated = X[D_mask, , drop = FALSE],
                         X_control = X[!D_mask, , drop = FALSE],
                         n_total = n_total))
  }
  beta_ps <- glm_fit$coefficients
  eta <- as.numeric(cbind(1, X) %*% beta_ps)
  p_hat <- pmin(pmax(1 / (1 + exp(-eta)), 0.01), 0.99)
  E_D <- mean(as.integer(D_mask))
  w1 <- as.integer(D_mask) / E_D
  w0 <- (1 - as.integer(D_mask)) * (p_hat / (1 - p_hat)) / E_D
  # DR formula: ATT = E[w1 * (Y - m)] - E[w0 * (Y - m)]
  resid <- delta - m_hat
  att <- mean(w1 * resid) - mean(w0 * resid)
  IF <- (w1 * resid) - (w0 * resid) - att
  list(att = att, IF = IF,
       n_treated = as.integer(n_t),
       n_control = as.integer(n_c))
}


# Method dispatch.
#' @keywords internal
#' @noRd
.cs_inner_dispatch <- function(method, delta, D_mask, X, n_total) {
  switch(method,
    "OR"  = .cs_inner_or(delta, D_mask,
                          X_treated = if (!is.null(X)) X[D_mask, , drop = FALSE] else NULL,
                          X_control = if (!is.null(X)) X[!D_mask, , drop = FALSE] else NULL,
                          n_total = n_total),
    "IPW" = .cs_inner_ipw(delta, D_mask, X, n_total),
    "DR"  = .cs_inner_dr(delta, D_mask, X, n_total),
    stop("Unknown method: ", method)
  )
}


# ============================================================================
# CUDA batched-inner: try to compute ALL cells in one GPU call.
#
# Phase-1 plumbing for task #79. The underlying kernel is scaffolded
# (returns -1), so this function returns NULL today and the caller
# falls back to the per-cell R loop. Once Phase 2 (#82-#85) fills the
# kernel, this function returns the same shape as the per-cell loop
# would have produced (an att vector + per-cell influence-function
# vectors) and the orchestrator skips the loop.
#
# Inputs:
#   cells      list of per-cell data; each element has $delta (numeric),
#              $D_mask (logical), $X (matrix or NULL), $units (integer
#              unit IDs in the order they appear in delta).
#   method     "OR" / "IPW" / "DR"
#   all_units  integer vector of all unit IDs in canonical order (the
#              influence-function rows are indexed by position in this
#              vector).
#
# Returns NULL if CUDA unavailable, the kernel returns nonzero, or any
# error is thrown — the caller is responsible for falling back.
#' @keywords internal
#' @noRd
.cs_inner_batched_cuda <- function(cells, method, all_units) {
  if (!isTRUE(tryCatch(didgpu_has_cuda_support(),
                       error = function(e) FALSE))) return(NULL)
  if (length(cells) == 0L) return(NULL)

  method_int <- switch(method, "OR" = 0L, "IPW" = 1L, "DR" = 2L, NA_integer_)
  if (is.na(method_int)) return(NULL)

  # Marshal cells -> concatenated buffers. If any cell has no
  # covariate matrix, fall back to a single intercept column. Mixed
  # (some-with-X some-without) is rejected — the canonical CS layout
  # uses a single design across all cells.
  has_X <- vapply(cells, function(c) !is.null(c$X), logical(1))
  if (any(has_X) && !all(has_X)) return(NULL)

  p <- if (any(has_X)) ncol(cells[[1]]$X) + 1L else 1L  # +1 for intercept
  n_cells <- length(cells)
  n_per_cell <- vapply(cells, function(c) length(c$delta), integer(1))
  offsets <- as.integer(c(0L, cumsum(n_per_cell)))
  n_total <- offsets[n_cells + 1L]
  if (n_total == 0L) return(NULL)

  X_concat <- numeric(n_total * p)
  Y_concat <- numeric(n_total)
  W_concat <- numeric(n_total)
  unit_id_per_row <- integer(n_total)
  unit_to_idx <- stats::setNames(seq_along(all_units) - 1L,
                                  as.character(all_units))
  for (c_idx in seq_len(n_cells)) {
    ce <- cells[[c_idx]]
    rows <- (offsets[c_idx] + 1L):offsets[c_idx + 1L]
    Y_concat[rows] <- ce$delta
    W_concat[rows] <- as.numeric(ce$D_mask)
    unit_id_per_row[rows] <- unit_to_idx[as.character(ce$units)]
    base <- offsets[c_idx] * p
    if (any(has_X)) {
      Xm <- cbind(1.0, ce$X)         # intercept first
      for (r in seq_len(nrow(Xm))) {
        X_concat[(base + (r - 1L) * p + 1L):(base + r * p)] <- Xm[r, ]
      }
    } else {
      X_concat[(base + 1L):(base + n_per_cell[c_idx])] <- 1.0
    }
  }

  result <- tryCatch(
    didgpu_cuda_cs_inner_batched_r(
      X_concat        = X_concat,
      X_offsets       = offsets,
      Y_concat        = Y_concat,
      W_concat        = W_concat,
      unit_id_per_row = unit_id_per_row,
      p               = as.integer(p),
      n_units         = length(all_units),
      est_method      = method_int,
      want_influence  = TRUE),
    error = function(e) NULL)

  if (is.null(result)) return(NULL)
  # Sanity-check return shape.
  if (is.null(result$att) || length(result$att) != n_cells) return(NULL)
  result
}
