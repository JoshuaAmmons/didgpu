# Repeated-cross-section Callaway-Sant'Anna cell estimators.
#
# did::att_gt uses these when a panel is unbalanced and
# allow_unbalanced_panel = TRUE (compute.att_gt, the `!panel` branch): a
# cell is then every row observed at the two periods, not one row per
# unit, and the three estimators are DRDID's repeated-cross-section ones.
# Each function below mirrors its DRDID counterpart statement for
# statement, including the influence function, so ATT(g,t) and its
# per-row influence contribution match to machine precision.
#
#   .cs_rc_or   <- DRDID::reg_did_rc
#   .cs_rc_ipw  <- DRDID::std_ipw_did_rc
#   .cs_rc_dr   <- DRDID::drdid_rc      (the locally efficient DR estimator,
#                                        which is the one did calls)
#
# Arguments, per cell: y (outcome, one per ROW), post (1 at the current
# period, 0 at the base period), D (1 if the row's unit is in cohort g),
# X (covariate matrix without intercept, or NULL), w (row weights).
# Each returns list(att, IF) with IF one value per row; the caller sums it
# by unit, as did does (stats::aggregate(att.inf.func, list(rightids), sum)).

#' @keywords internal
#' @noRd
.cs_rc_or <- function(y, post, D, X = NULL, w = NULL) {
  n <- length(D)
  int.cov <- .cs_int_cov(X, n)
  i.weights <- if (is.null(w)) rep(1, n) else w
  i.weights <- i.weights / mean(i.weights)

  pre_filter  <- (D == 0) & (post == 0)
  post_filter <- (D == 0) & (post == 1)
  reg.coeff.pre  <- .cs_wls(int.cov, y, i.weights, keep = pre_filter)
  reg.coeff.post <- .cs_wls(int.cov, y, i.weights, keep = post_filter)
  if (is.null(reg.coeff.pre) || is.null(reg.coeff.post)) return(NULL)
  out.y.pre  <- as.numeric(int.cov %*% reg.coeff.pre)
  out.y.post <- as.numeric(int.cov %*% reg.coeff.post)

  w.treat.pre  <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont       <- i.weights * D
  reg.att.treat.pre  <- w.treat.pre * y
  reg.att.treat.post <- w.treat.post * y
  reg.att.cont       <- w.cont * (out.y.post - out.y.pre)
  eta.treat.pre  <- mean(reg.att.treat.pre)  / mean(w.treat.pre)
  eta.treat.post <- mean(reg.att.treat.post) / mean(w.treat.post)
  eta.cont       <- mean(reg.att.cont)       / mean(w.cont)
  att <- (eta.treat.post - eta.treat.pre) - eta.cont

  weights.ols.pre <- i.weights * (1 - D) * (1 - post)
  wols.x.pre  <- weights.ols.pre * int.cov
  wols.eX.pre <- weights.ols.pre * (y - out.y.pre) * int.cov
  XpX.inv.pre <- .cs_safe_inv(crossprod(wols.x.pre, int.cov) / n)
  weights.ols.post <- i.weights * (1 - D) * post
  wols.x.post  <- weights.ols.post * int.cov
  wols.eX.post <- weights.ols.post * (y - out.y.post) * int.cov
  XpX.inv.post <- .cs_safe_inv(crossprod(wols.x.post, int.cov) / n)
  if (is.null(XpX.inv.pre) || is.null(XpX.inv.post)) return(NULL)
  asy.lin.rep.ols.pre  <- wols.eX.pre  %*% XpX.inv.pre
  asy.lin.rep.ols.post <- wols.eX.post %*% XpX.inv.post

  inf.treat.pre  <- (reg.att.treat.pre  - w.treat.pre  * eta.treat.pre)  / mean(w.treat.pre)
  inf.treat.post <- (reg.att.treat.post - w.treat.post * eta.treat.post) / mean(w.treat.post)
  inf.treat <- inf.treat.post - inf.treat.pre
  inf.cont.1 <- (reg.att.cont - w.cont * eta.cont)
  M1 <- colMeans(w.cont * int.cov)
  inf.cont.2.post <- asy.lin.rep.ols.post %*% M1
  inf.cont.2.pre  <- asy.lin.rep.ols.pre  %*% M1
  inf.control <- (inf.cont.1 + inf.cont.2.post - inf.cont.2.pre) / mean(w.cont)
  list(att = att, IF = as.numeric(inf.treat - inf.control))
}

# DRDID fits its propensity score with fastglm(..., intercept = FALSE) on
# int.cov, i.e. a logit of D on the intercept-augmented covariates.
#' @keywords internal
#' @noRd
.cs_rc_ps <- function(int.cov, D, i.weights) {
  ps.fit <- .cs_pscore(int.cov, D, i.weights)
  if (is.null(ps.fit)) return(NULL)
  pmin(ps.fit, 1 - 1e-6)
}

#' @keywords internal
#' @noRd
.cs_rc_ipw <- function(y, post, D, X = NULL, w = NULL, trim.level = 0.995) {
  n <- length(D)
  int.cov <- .cs_int_cov(X, n)
  i.weights <- if (is.null(w)) rep(1, n) else w
  i.weights <- i.weights / mean(i.weights)
  ps.fit <- .cs_rc_ps(int.cov, D, i.weights)
  if (is.null(ps.fit)) return(NULL)
  W <- ps.fit * (1 - ps.fit) * i.weights
  trim.ps <- (ps.fit < 1.01)
  trim.ps[D == 0] <- (ps.fit[D == 0] < trim.level)

  w.treat.pre  <- trim.ps * i.weights * D * (1 - post)
  w.treat.post <- trim.ps * i.weights * D * post
  w.cont.pre   <- trim.ps * i.weights * ps.fit * (1 - D) * (1 - post) / (1 - ps.fit)
  w.cont.post  <- trim.ps * i.weights * ps.fit * (1 - D) * post / (1 - ps.fit)
  eta.treat.pre  <- w.treat.pre  * y / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * y / mean(w.treat.post)
  eta.cont.pre   <- w.cont.pre   * y / mean(w.cont.pre)
  eta.cont.post  <- w.cont.post  * y / mean(w.cont.post)
  att.treat.pre  <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre   <- mean(eta.cont.pre)
  att.cont.post  <- mean(eta.cont.post)
  att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)

  score.ps <- i.weights * (D - ps.fit) * int.cov
  Hessian.ps <- tryCatch(chol2inv(chol(t(int.cov) %*% (W * int.cov))) * n,
                         error = function(e) NULL)
  if (is.null(Hessian.ps)) return(NULL)
  asy.lin.rep.ps <- score.ps %*% Hessian.ps
  inf.treat.pre  <- eta.treat.pre  - w.treat.pre  * att.treat.pre  / mean(w.treat.pre)
  inf.treat.post <- eta.treat.post - w.treat.post * att.treat.post / mean(w.treat.post)
  inf.treat <- inf.treat.post - inf.treat.pre
  inf.cont.pre  <- eta.cont.pre  - w.cont.pre  * att.cont.pre  / mean(w.cont.pre)
  inf.cont.post <- eta.cont.post - w.cont.post * att.cont.post / mean(w.cont.post)
  inf.cont <- inf.cont.post - inf.cont.pre
  M2.pre  <- colMeans(w.cont.pre  * (y - att.cont.pre)  * int.cov) / mean(w.cont.pre)
  M2.post <- colMeans(w.cont.post * (y - att.cont.post) * int.cov) / mean(w.cont.post)
  inf.cont <- inf.cont + asy.lin.rep.ps %*% (M2.post - M2.pre)
  list(att = att, IF = as.numeric(inf.treat - inf.cont))
}

#' @keywords internal
#' @noRd
.cs_rc_dr <- function(y, post, D, X = NULL, w = NULL, trim.level = 0.995) {
  n <- length(D)
  int.cov <- .cs_int_cov(X, n)
  i.weights <- if (is.null(w)) rep(1, n) else w
  i.weights <- i.weights / mean(i.weights)
  ps.fit <- .cs_rc_ps(int.cov, D, i.weights)
  if (is.null(ps.fit)) return(NULL)
  W <- ps.fit * (1 - ps.fit) * i.weights
  trim.ps <- (ps.fit < 1.01)
  trim.ps[D == 0] <- (ps.fit[D == 0] < trim.level)

  fit_on <- function(keep) .cs_wls(int.cov, y, i.weights, keep)
  b.cont.pre   <- fit_on((D == 0) & (post == 0))
  b.cont.post  <- fit_on((D == 0) & (post == 1))
  b.treat.pre  <- fit_on((D == 1) & (post == 0))
  b.treat.post <- fit_on((D == 1) & (post == 1))
  if (is.null(b.cont.pre) || is.null(b.cont.post) ||
      is.null(b.treat.pre) || is.null(b.treat.post)) return(NULL)
  out.y.cont.pre   <- as.numeric(int.cov %*% b.cont.pre)
  out.y.cont.post  <- as.numeric(int.cov %*% b.cont.post)
  out.y.cont       <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  out.y.treat.pre  <- as.numeric(int.cov %*% b.treat.pre)
  out.y.treat.post <- as.numeric(int.cov %*% b.treat.post)

  w.treat.pre  <- trim.ps * i.weights * D * (1 - post)
  w.treat.post <- trim.ps * i.weights * D * post
  w.cont.pre   <- trim.ps * i.weights * ps.fit * (1 - D) * (1 - post) / (1 - ps.fit)
  w.cont.post  <- trim.ps * i.weights * ps.fit * (1 - D) * post / (1 - ps.fit)
  w.d   <- trim.ps * i.weights * D
  w.dt1 <- trim.ps * i.weights * D * post
  w.dt0 <- trim.ps * i.weights * D * (1 - post)

  eta.treat.pre  <- w.treat.pre  * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont) / mean(w.treat.post)
  eta.cont.pre   <- w.cont.pre   * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post  <- w.cont.post  * (y - out.y.cont) / mean(w.cont.post)
  eta.d.post   <- w.d   * (out.y.treat.post - out.y.cont.post) / mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post) / mean(w.dt1)
  eta.d.pre    <- w.d   * (out.y.treat.pre  - out.y.cont.pre)  / mean(w.d)
  eta.dt0.pre  <- w.dt0 * (out.y.treat.pre  - out.y.cont.pre)  / mean(w.dt0)
  att.treat.pre  <- mean(eta.treat.pre);  att.treat.post <- mean(eta.treat.post)
  att.cont.pre   <- mean(eta.cont.pre);   att.cont.post  <- mean(eta.cont.post)
  att.d.post     <- mean(eta.d.post);     att.dt1.post   <- mean(eta.dt1.post)
  att.d.pre      <- mean(eta.d.pre);      att.dt0.pre    <- mean(eta.dt0.pre)
  att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
         (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)

  alr <- function(wts, resid) {
    inv <- .cs_safe_inv(crossprod(wts * int.cov, int.cov) / n)
    if (is.null(inv)) return(NULL)
    (wts * resid * int.cov) %*% inv
  }
  alr.pre   <- alr(i.weights * (1 - D) * (1 - post), y - out.y.cont.pre)
  alr.post  <- alr(i.weights * (1 - D) * post,       y - out.y.cont.post)
  alr.pre.t <- alr(i.weights * D * (1 - post),       y - out.y.treat.pre)
  alr.post.t<- alr(i.weights * D * post,             y - out.y.treat.post)
  if (is.null(alr.pre) || is.null(alr.post) ||
      is.null(alr.pre.t) || is.null(alr.post.t)) return(NULL)
  score.ps <- i.weights * (D - ps.fit) * int.cov
  Hessian.ps <- tryCatch(chol2inv(chol(t(int.cov) %*% (W * int.cov))) * n,
                         error = function(e) NULL)
  if (is.null(Hessian.ps)) return(NULL)
  asy.lin.rep.ps <- score.ps %*% Hessian.ps

  inf.treat.pre  <- eta.treat.pre  - w.treat.pre  * att.treat.pre  / mean(w.treat.pre)
  inf.treat.post <- eta.treat.post - w.treat.post * att.treat.post / mean(w.treat.post)
  M1.post <- -colMeans(w.treat.post * post * int.cov) / mean(w.treat.post)
  M1.pre  <- -colMeans(w.treat.pre * (1 - post) * int.cov) / mean(w.treat.pre)
  inf.treat.or <- alr.post %*% M1.post + alr.pre %*% M1.pre
  inf.cont.pre  <- eta.cont.pre  - w.cont.pre  * att.cont.pre  / mean(w.cont.pre)
  inf.cont.post <- eta.cont.post - w.cont.post * att.cont.post / mean(w.cont.post)
  M2.pre  <- colMeans(w.cont.pre  * (y - out.y.cont - att.cont.pre)  * int.cov) / mean(w.cont.pre)
  M2.post <- colMeans(w.cont.post * (y - out.y.cont - att.cont.post) * int.cov) / mean(w.cont.post)
  inf.cont.ps <- asy.lin.rep.ps %*% (M2.post - M2.pre)
  M3.post <- -colMeans(w.cont.post * post * int.cov) / mean(w.cont.post)
  M3.pre  <- -colMeans(w.cont.pre * (1 - post) * int.cov) / mean(w.cont.pre)
  inf.cont.or <- alr.post %*% M3.post + alr.pre %*% M3.pre
  inf.eff1 <- eta.d.post   - w.d   * att.d.post   / mean(w.d)
  inf.eff2 <- eta.dt1.post - w.dt1 * att.dt1.post / mean(w.dt1)
  inf.eff3 <- eta.d.pre    - w.d   * att.d.pre    / mean(w.d)
  inf.eff4 <- eta.dt0.pre  - w.dt0 * att.dt0.pre  / mean(w.dt0)
  inf.eff <- (inf.eff1 - inf.eff2) - (inf.eff3 - inf.eff4)
  mom.post <- colMeans((w.d / mean(w.d) - w.dt1 / mean(w.dt1)) * int.cov)
  mom.pre  <- colMeans((w.d / mean(w.d) - w.dt0 / mean(w.dt0)) * int.cov)
  inf.or <- (alr.post.t - alr.post) %*% mom.post - (alr.pre.t - alr.pre) %*% mom.pre
  inf.treat <- inf.treat.post - inf.treat.pre + inf.treat.or
  inf.cont  <- inf.cont.post - inf.cont.pre + inf.cont.ps + inf.cont.or
  list(att = att, IF = as.numeric(inf.treat - inf.cont + inf.eff + inf.or))
}


# One repeated-cross-section cell, did's `!panel` branch
# (compute.att_gt, lines 158-226). Rows are every observation of a
# treated or control unit at the base or current period. Returns the ATT
# and the unit-level influence function ON didgpu's STORAGE SCALE: did
# multiplies the per-row influence function by n / n1 (n = units in the
# data, n1 = rows in the cell) and sums it by unit; every downstream
# consumer here multiplies a cell's stored IF by n / (units in the cell),
# so the stored value is did's times (units in the cell) / n.
#' @keywords internal
#' @noRd
.cs_rc_cell <- function(method, y, post, D, C, X, row_unit, n_units) {
  G <- D
  n_t_post <- sum(G * post); n_t_pre <- sum(G * (1 - post))
  n_c_post <- sum(C * post); n_c_pre <- sum(C * (1 - post))
  cell_units <- unique(row_unit)
  empty <- list(att = NA_real_, IF = rep(0, length(cell_units)),
                units = cell_units)
  if (n_t_post == 0 || n_t_pre == 0 || n_c_post == 0 || n_c_pre == 0) {
    return(empty)
  }
  fun <- switch(toupper(method), OR = .cs_rc_or, IPW = .cs_rc_ipw,
                DR = .cs_rc_dr, stop("unknown est_method: ", method))
  est <- tryCatch(fun(y = y, post = post, D = D, X = X), error = function(e) NULL)
  if (is.null(est) || !is.finite(est$att)) return(empty)
  n1 <- length(y)
  per_unit <- rowsum((n_units / n1) * est$IF, row_unit, reorder = FALSE)
  u <- rownames(per_unit)
  list(att = est$att,
       IF = as.numeric(per_unit[, 1]) * length(u) / n_units,
       units = utils::type.convert(u, as.is = TRUE))
}
