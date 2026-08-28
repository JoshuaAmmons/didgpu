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
# ---------------------------------------------------------------------------
# Per-cell 2x2 estimators and their INFLUENCE FUNCTIONS.
#
# These mirror DRDID (the package that `did` itself calls) function for
# function, because bit-for-bit agreement with did::att_gt() is the
# contract:
#     est_method = "reg" -> DRDID::reg_did_panel
#     est_method = "ipw" -> DRDID::std_ipw_did_panel
#     est_method = "dr"  -> DRDID::drdid_panel
#
# The influence function of an ATT is NOT just the treated units'
# demeaned residual. It has three parts:
#   (a) the treated arm, normalised by mean(w.treat) = E[D];
#   (b) the comparison arm, normalised by mean(w.cont) -- which for IPW/DR
#       is E[p(X)(1-D)/(1-p(X))], NOT E[D];
#   (c) estimation-effect terms for the nuisance parameters: the OLS
#       outcome regression (asy.lin.rep.wols) and/or the propensity score
#       (asy.lin.rep.ps). These load onto CONTROL units.
#
# Previous versions dropped (c) entirely and used the wrong normaliser in
# (b): OR set every control unit's influence to zero, and IPW/DR omitted
# the 1/E[D] and 1/E[p(1-D)/(1-p)] scaling. Measured against did::att_gt
# on a 200-unit panel, the resulting multiplier-bootstrap SEs were ~0.13x
# (OR) and ~0.50x (IPW/DR) of the correct width -- confidence intervals
# two to eight times too narrow. The cluster bootstrap never touches
# these, which is why it was unaffected and correct throughout.
#
# Reference: Callaway and Sant'Anna (2021, J. Econometrics 225, Theorem
# 2); Sant'Anna and Zhao (2020). SE convention: sd(psi) * sqrt(n - 1) / n.
# ---------------------------------------------------------------------------

# Design matrix with intercept, matching DRDID's `int.cov`.
#' @keywords internal
#' @noRd
.cs_int_cov <- function(X, n) {
  if (is.null(X) || !is.matrix(X) || ncol(X) == 0L) {
    matrix(1, nrow = n, ncol = 1L)
  } else {
    cbind(1, X)
  }
}

# Weighted least squares of `y` on `Xm` over rows `keep`. NULL if singular.
#' @keywords internal
#' @noRd
.cs_wls <- function(Xm, y, w, keep) {
  Xk <- Xm[keep, , drop = FALSE]
  XpX <- crossprod(Xk * w[keep], Xk)
  rc <- tryCatch(rcond(XpX), error = function(e) 0)
  if (!is.finite(rc) || rc < .Machine$double.eps) return(NULL)
  as.numeric(solve(XpX, crossprod(Xk * w[keep], y[keep])))
}

# Propensity score by logistic regression on the full cell sample.
#' @keywords internal
#' @noRd
.cs_pscore <- function(Xm, D, w) {
  fit <- tryCatch(
    suppressWarnings(stats::glm.fit(x = Xm, y = D, weights = w,
                                    family = stats::binomial())),
    error = function(e) NULL)
  if (is.null(fit) || anyNA(fit$coefficients)) return(NULL)
  pmin(fit$fitted.values, 1 - 1e-6)
}

#' @keywords internal
#' @noRd
.cs_na_cell <- function(n_t, n_c, n_total) {
  list(att = NA_real_, IF = rep(0, n_total),
       n_treated = as.integer(n_t), n_control = as.integer(n_c))
}

#' @keywords internal
#' @noRd
.cs_safe_inv <- function(M) {
  rc <- tryCatch(rcond(M), error = function(e) 0)
  if (!is.finite(rc) || rc < .Machine$double.eps) return(NULL)
  solve(M)
}


# OR / "reg": mirrors DRDID::reg_did_panel.
#' @keywords internal
#' @noRd
.cs_inner_or <- function(delta, D_mask, X_treated, X_control, n_total,
                          X = NULL) {
  n_t <- sum(D_mask); n_c <- length(delta) - n_t
  if (n_t == 0L || n_c == 0L) return(.cs_na_cell(n_t, n_c, n_total))
  n  <- length(delta)
  D  <- as.numeric(D_mask)
  iw <- rep(1, n)
  if (is.null(X) && !is.null(X_control) && ncol(X_control) > 0L) {
    X <- matrix(NA_real_, nrow = n, ncol = ncol(X_control))
    X[D_mask, ] <- X_treated
    X[!D_mask, ] <- X_control
  }
  int.cov <- .cs_int_cov(X, n)

  reg.coeff <- .cs_wls(int.cov, delta, iw, keep = (D == 0))
  if (is.null(reg.coeff)) return(.cs_na_cell(n_t, n_c, n_total))
  out.delta <- as.numeric(int.cov %*% reg.coeff)

  # DRDID uses w.cont = w.treat = D: the comparison arm is the regression
  # prediction evaluated on the TREATED units.
  w.treat <- iw * D
  w.cont  <- iw * D
  reg.att.treat <- w.treat * delta
  reg.att.cont  <- w.cont * out.delta
  eta.treat <- mean(reg.att.treat) / mean(w.treat)
  eta.cont  <- mean(reg.att.cont)  / mean(w.cont)
  att <- eta.treat - eta.cont

  weights.ols <- iw * (1 - D)
  wols.x  <- weights.ols * int.cov
  wols.eX <- weights.ols * (delta - out.delta) * int.cov
  XpXinv <- .cs_safe_inv(crossprod(wols.x, int.cov) / n)
  if (is.null(XpXinv)) return(.cs_na_cell(n_t, n_c, n_total))
  asy.lin.rep.ols <- wols.eX %*% XpXinv

  inf.treat  <- (reg.att.treat - w.treat * eta.treat) / mean(w.treat)
  inf.cont.1 <- (reg.att.cont - w.cont * eta.cont)
  M1 <- colMeans(w.cont * int.cov)
  inf.cont.2 <- asy.lin.rep.ols %*% M1
  inf.control <- (inf.cont.1 + inf.cont.2) / mean(w.cont)

  list(att = att, IF = as.numeric(inf.treat - inf.control),
       n_treated = as.integer(n_t), n_control = as.integer(n_c))
}


# IPW: mirrors DRDID::std_ipw_did_panel (Hajek / standardised weights).
#' @keywords internal
#' @noRd
.cs_inner_ipw <- function(delta, D_mask, X, n_total, trim.level = 0.995) {
  n_t <- sum(D_mask); n_c <- length(delta) - n_t
  if (n_t == 0L || n_c == 0L) return(.cs_na_cell(n_t, n_c, n_total))
  n  <- length(delta)
  D  <- as.numeric(D_mask)
  iw <- rep(1, n)
  int.cov <- .cs_int_cov(X, n)

  ps.fit <- .cs_pscore(int.cov, D, iw)
  if (is.null(ps.fit)) return(.cs_na_cell(n_t, n_c, n_total))
  W <- ps.fit * (1 - ps.fit) * iw

  trim.ps <- (ps.fit < 1.01)
  trim.ps[D == 0] <- (ps.fit[D == 0] < trim.level)

  w.treat <- trim.ps * iw * D
  w.cont  <- trim.ps * iw * ps.fit * (1 - D) / (1 - ps.fit)
  if (mean(w.treat) == 0 || mean(w.cont) == 0) {
    return(.cs_na_cell(n_t, n_c, n_total))
  }
  att.treat <- w.treat * delta
  att.cont  <- w.cont * delta
  eta.treat <- mean(att.treat) / mean(w.treat)
  eta.cont  <- mean(att.cont)  / mean(w.cont)
  att <- eta.treat - eta.cont

  score.ps <- iw * (D - ps.fit) * int.cov
  Hinv <- .cs_safe_inv(crossprod(int.cov, W * int.cov))
  if (is.null(Hinv)) return(.cs_na_cell(n_t, n_c, n_total))
  asy.lin.rep.ps <- score.ps %*% (Hinv * n)

  inf.treat  <- (att.treat - w.treat * eta.treat) / mean(w.treat)
  inf.cont.1 <- (att.cont - w.cont * eta.cont)
  M2 <- colMeans(w.cont * (delta - eta.cont) * int.cov)
  inf.cont.2 <- asy.lin.rep.ps %*% M2
  inf.control <- (inf.cont.1 + inf.cont.2) / mean(w.cont)

  list(att = att, IF = as.numeric(inf.treat - inf.control),
       n_treated = as.integer(n_t), n_control = as.integer(n_c))
}


# DR: mirrors DRDID::drdid_panel.
#' @keywords internal
#' @noRd
.cs_inner_dr <- function(delta, D_mask, X, n_total, trim.level = 0.995) {
  n_t <- sum(D_mask); n_c <- length(delta) - n_t
  if (n_t == 0L || n_c == 0L) return(.cs_na_cell(n_t, n_c, n_total))
  n  <- length(delta)
  D  <- as.numeric(D_mask)
  iw <- rep(1, n)
  int.cov <- .cs_int_cov(X, n)

  ps.fit <- .cs_pscore(int.cov, D, iw)
  if (is.null(ps.fit)) {
    return(.cs_inner_or(delta, D_mask, NULL, NULL, n_total, X = X))
  }
  trim.ps <- (ps.fit < 1.01)
  trim.ps[D == 0] <- (ps.fit[D == 0] < trim.level)
  W <- ps.fit * (1 - ps.fit) * iw

  reg.coeff <- .cs_wls(int.cov, delta, iw, keep = (D == 0))
  if (is.null(reg.coeff)) return(.cs_na_cell(n_t, n_c, n_total))
  out.delta <- as.numeric(int.cov %*% reg.coeff)

  w.treat <- trim.ps * iw * D
  w.cont  <- trim.ps * iw * ps.fit * (1 - D) / (1 - ps.fit)
  if (mean(w.treat) == 0 || mean(w.cont) == 0) {
    return(.cs_na_cell(n_t, n_c, n_total))
  }
  dr.att.treat <- w.treat * (delta - out.delta)
  dr.att.cont  <- w.cont * (delta - out.delta)
  eta.treat <- mean(dr.att.treat) / mean(w.treat)
  eta.cont  <- mean(dr.att.cont)  / mean(w.cont)
  att <- eta.treat - eta.cont

  weights.ols <- iw * (1 - D)
  wols.x  <- weights.ols * int.cov
  wols.eX <- weights.ols * (delta - out.delta) * int.cov
  XpXinv <- .cs_safe_inv(crossprod(wols.x, int.cov) / n)
  if (is.null(XpXinv)) return(.cs_na_cell(n_t, n_c, n_total))
  asy.lin.rep.wols <- wols.eX %*% XpXinv

  score.ps <- iw * (D - ps.fit) * int.cov
  Hinv <- .cs_safe_inv(crossprod(int.cov, W * int.cov))
  if (is.null(Hinv)) return(.cs_na_cell(n_t, n_c, n_total))
  asy.lin.rep.ps <- score.ps %*% (Hinv * n)

  inf.treat.1 <- (dr.att.treat - w.treat * eta.treat)
  M1 <- colMeans(w.treat * int.cov)
  inf.treat.2 <- asy.lin.rep.wols %*% M1
  inf.treat <- (inf.treat.1 - inf.treat.2) / mean(w.treat)

  inf.cont.1 <- (dr.att.cont - w.cont * eta.cont)
  M2 <- colMeans(w.cont * (delta - out.delta - eta.cont) * int.cov)
  inf.cont.2 <- asy.lin.rep.ps %*% M2
  M3 <- colMeans(w.cont * int.cov)
  inf.cont.3 <- asy.lin.rep.wols %*% M3
  inf.control <- (inf.cont.1 + inf.cont.2 - inf.cont.3) / mean(w.cont)

  list(att = att, IF = as.numeric(inf.treat - inf.control),
       n_treated = as.integer(n_t), n_control = as.integer(n_c))
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
