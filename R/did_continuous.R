# ============================================================================
# didgpu_did_continuous(): de Chaisemartin & D'Haultfoeuille (2024)-style
# difference-in-differences for a CONTINUOUS treatment with NO STAYERS.
#
# When the continuous treatment changes for (almost) every unit, there are no
# pure "stayers" to serve as controls. This estimator works in first
# differences: with dY = Y_post - Y_pre and dD = D_post - D_pre per unit,
# under parallel trends E[dY | dD = 0] is the common trend (identified from
# "quasi-stayers" with dD ~ 0), and
#   effect(d) = E[dY | dD = d] - E[dY | dD = 0]      (effect of a dose change d)
#   ACR(d)    = d/dd E[dY | dD = d]                   (avg causal response)
# Two estimators of E[dY | dD = d]:
#   * estimator = "parametric"   : a degree-`degree` polynomial in dD (OLS).
#                                  sqrt(n) rate. The author's parametric
#                                  alternative.
#   * estimator = "nonparametric": local-linear (kernel) regression in dD,
#                                  using quasi-stayers near each evaluation
#                                  point. n^{2/5} rate (slower; it estimates a
#                                  derivative). EXPERIMENTAL.
# Multiplier-bootstrap SEs. Both estimators are validated by SIMULATION
# (recovering a known dose-response); unlike didgpu's other new estimators,
# there is no maintained R reference package to cross-check bit-for-bit, so
# the nonparametric path in particular is flagged experimental.
#
# Note: didgpu()'s `continuous = k` argument provides the related
# DIDmultiplegtDYN-equivalent parametric continuous estimator within the
# dynamic staggered framework; this function targets the no-stayers
# first-difference design.
#
# Reference: de Chaisemartin, D'Haultfoeuille, Pasquier & Vazquez-Bare
# (2024), "Difference-in-Differences Estimators for Treatments Continuously
# Distributed at Every Period".
# ============================================================================


# Per-unit first differences dY, dD (uses the earliest and latest observed
# period per unit).
#' @keywords internal
#' @noRd
.didc_first_diff <- function(d, idv, tv, yv, dv) {
  data.table::setorderv(d, c(idv, tv))
  d[, list(dY = get(yv)[.N] - get(yv)[1L],
           dD = get(dv)[.N] - get(dv)[1L]), by = c(idv)]
}

# Local-linear regression slope (ACR) and level at points `at`, bandwidth h,
# Gaussian kernel. Returns list(level, slope) at each `at`.
#' @keywords internal
#' @noRd
.didc_loclin <- function(x, y, at, h) {
  lev <- numeric(length(at)); slp <- numeric(length(at))
  for (j in seq_along(at)) {
    u <- (x - at[j]) / h
    w <- exp(-0.5 * u * u)
    sw <- sum(w)
    if (sw <= 0) { lev[j] <- NA; slp[j] <- NA; next }
    xc <- x - at[j]
    # weighted local linear: solve [1 xc] b = y with weights w
    s0 <- sum(w); s1 <- sum(w * xc); s2 <- sum(w * xc * xc)
    t0 <- sum(w * y); t1 <- sum(w * xc * y)
    det <- s0 * s2 - s1 * s1
    if (abs(det) < .Machine$double.eps) { lev[j] <- t0 / s0; slp[j] <- NA; next }
    a <- (s2 * t0 - s1 * t1) / det      # intercept = level at at[j]
    b <- (s0 * t1 - s1 * t0) / det      # slope = ACR at at[j]
    lev[j] <- a; slp[j] <- b
  }
  list(level = lev, slope = slp)
}


#' Continuous-treatment DiD with no stayers (de Chaisemartin-D'Haultfoeuille 2024)
#'
#' First-difference difference-in-differences for a continuous treatment that
#' changes for (almost) all units, so there are no pure stayers. Estimates the
#' level effect of a dose change, \eqn{effect(d) = E[dY|dD=d] - E[dY|dD=0]},
#' and the average causal response \eqn{ACR(d) = d/dd\,E[dY|dD=d]}, where
#' \eqn{dY, dD} are within-unit first differences and \eqn{E[dY|dD=0]} is the
#' common trend identified from quasi-stayers (units with \eqn{dD \approx 0}).
#'
#' @param df A data.frame / data.table panel.
#' @param outcome,treatment,id,time Character column names: outcome, the
#'   continuous treatment, unit id, and time.
#' @param estimator `"parametric"` (degree-`degree` polynomial in dD; sqrt(n))
#'   or `"nonparametric"` (local-linear in dD; n^{2/5}; EXPERIMENTAL).
#' @param degree Polynomial degree for the parametric estimator (default 2).
#' @param dvals Dose-change values at which to report effect(d)/ACR(d).
#'   Default: quantiles 0.1..0.9 of the nonzero dD.
#' @param bandwidth Local-linear bandwidth (nonparametric). Default: a
#'   Silverman-type rule, `1.06 * sd(dD) * n^(-1/5)`.
#' @param bootstrap_reps Multiplier-bootstrap replicates for SEs (default 200).
#' @param ci_level Confidence level in percent (default 95).
#' @param seed Bootstrap RNG seed.
#' @param verbose Logical.
#' @return A `didgpu_did_continuous_result`: a list with `dose` (dvals),
#'   `effect.d`, `acr.d` (+ SEs/CIs), `overall_acr`, `estimator`, and `args`.
#' @references de Chaisemartin, C., D'Haultfoeuille, X., Pasquier, F. &
#'   Vazquez-Bare, G. (2024). Difference-in-Differences Estimators for
#'   Treatments Continuously Distributed at Every Period.
#' @seealso [didgpu_cs_continuous()] (Callaway et al. 2024 dose-response);
#'   [didgpu()] with `continuous=` for the DIDmultiplegtDYN-equivalent
#'   parametric continuous estimator in the dynamic framework.
#' @examples
#' \donttest{
#' set.seed(1); nU <- 800L
#' dD <- stats::rnorm(nU, 0, 1)
#' dY <- 0.3 + 2 * dD - 0.5 * dD^2 + stats::rnorm(nU, 0, 0.5)
#' df <- data.frame(id = rep(seq_len(nU), each = 2L),
#'                  t = rep(1:2, nU),
#'                  D = as.numeric(rbind(0, dD)),
#'                  Y = as.numeric(rbind(0, dY)))
#' didgpu_did_continuous(df, "Y", "D", "id", "t", estimator = "parametric",
#'                       degree = 2, bootstrap_reps = 0)
#' }
#' @export
didgpu_did_continuous <- function(df, outcome, treatment, id, time,
                                  estimator = c("parametric", "nonparametric"),
                                  degree = 2L, dvals = NULL, bandwidth = NULL,
                                  bootstrap_reps = 200L, ci_level = 95,
                                  seed = 1L, verbose = TRUE) {
  estimator <- match.arg(estimator)
  for (nm in c("outcome", "treatment", "id", "time")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v))
      stop("`", nm, "` must be a single non-empty column name.")
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  degree <- as.integer(degree)
  d <- data.table::as.data.table(df)
  d <- d[!is.na(get(outcome)) & !is.na(get(treatment)) & !is.na(get(time))]
  fd <- .didc_first_diff(d, id, time, outcome, treatment)
  fd <- fd[is.finite(dD) & is.finite(dY)]
  if (nrow(fd) < degree + 2L) stop("too few units with first differences.")
  dD <- fd$dD; dY <- fd$dY
  nz <- dD[abs(dD) > 0]
  if (is.null(dvals)) {
    base <- if (length(nz) > 0) nz else dD
    dvals <- as.numeric(stats::quantile(base, probs = seq(0.1, 0.9, length.out = 5L)))
  }
  if (is.null(bandwidth))
    bandwidth <- 1.06 * stats::sd(dD) * length(dD)^(-1 / 5)

  # Point estimator: returns effect(d), acr(d), overall ACR, given (x=dD,y=dY).
  fit_fun <- function(x, y) {
    if (estimator == "parametric") {
      X <- outer(x, seq_len(degree), `^`)             # [x, x^2, ..., x^degree]
      fit <- stats::lm.fit(cbind(1, X), y)
      bet <- fit$coefficients                          # [intercept, b1..bdeg]
      slope_coef <- bet[-1L]
      polyval <- function(dd) as.numeric(outer(dd, seq_len(degree), `^`) %*% slope_coef)
      eff <- polyval(dvals)                            # E[dY|d]-E[dY|0] (intercept cancels)
      acr <- vapply(dvals, function(dd) {
        sum(slope_coef * seq_len(degree) * dd^(seq_len(degree) - 1L)) }, numeric(1))
      acr_overall <- mean(vapply(x, function(dd) {
        sum(slope_coef * seq_len(degree) * dd^(seq_len(degree) - 1L)) }, numeric(1)))
    } else {
      ll0 <- .didc_loclin(x, y, 0, bandwidth)$level    # trend (quasi-stayers)
      ll  <- .didc_loclin(x, y, dvals, bandwidth)
      eff <- ll$level - ll0
      acr <- ll$slope
      acr_overall <- mean(.didc_loclin(x, y, x, bandwidth)$slope, na.rm = TRUE)
    }
    list(eff = eff, acr = acr, acr_overall = acr_overall)
  }

  pt <- fit_fun(dD, dY)
  res <- list(dose = dvals, effect.d = pt$eff, acr.d = pt$acr,
              overall_acr = pt$acr_overall, estimator = estimator,
              effect.d_se = NA, acr.d_se = NA)

  bootstrap_reps <- as.integer(bootstrap_reps)
  if (bootstrap_reps > 0L) {
    n <- length(dD); set.seed(seed)
    beff <- matrix(NA_real_, bootstrap_reps, length(dvals))
    bacr <- matrix(NA_real_, bootstrap_reps, length(dvals))
    for (b in seq_len(bootstrap_reps)) {
      idx <- sample.int(n, n, replace = TRUE)
      fb <- tryCatch(fit_fun(dD[idx], dY[idx]), error = function(e) NULL)
      if (!is.null(fb)) { beff[b, ] <- fb$eff; bacr[b, ] <- fb$acr }
    }
    res$effect.d_se <- apply(beff, 2L, stats::sd, na.rm = TRUE)
    res$acr.d_se    <- apply(bacr, 2L, stats::sd, na.rm = TRUE)
    z <- stats::qnorm(1 - (1 - ci_level / 100) / 2)
    res$effect.d_lower <- pt$eff - z * res$effect.d_se
    res$effect.d_upper <- pt$eff + z * res$effect.d_se
    res$acr.d_lower <- pt$acr - z * res$acr.d_se
    res$acr.d_upper <- pt$acr + z * res$acr.d_se
  }

  res$args <- list(outcome = outcome, treatment = treatment, id = id,
                   time = time, degree = degree, bandwidth = bandwidth,
                   n_units = nrow(fd), ci_level = ci_level)
  class(res) <- c("didgpu_did_continuous_result", "list")
  if (verbose) print(res)
  invisible(res)
}


#' Print method for didgpu_did_continuous_result
#' @param x A `didgpu_did_continuous_result`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.didgpu_did_continuous_result <- function(x, ...) {
  cat(sprintf("de Chaisemartin-D'Haultfoeuille (2024) continuous DiD, no stayers [%s]\n",
              x$estimator))
  cat(sprintf("  %d units (first differences); ", x$args$n_units))
  if (x$estimator == "nonparametric")
    cat(sprintf("bandwidth %.4g  [EXPERIMENTAL: no reference cross-check]\n",
                x$args$bandwidth))
  else cat(sprintf("polynomial degree %d\n", x$args$degree))
  tab <- data.frame(dose_change = x$dose, effect = x$effect.d, ACR = x$acr.d)
  if (!all(is.na(x$effect.d_se))) { tab$effect_se <- x$effect.d_se; tab$ACR_se <- x$acr.d_se }
  cat("\n"); print(round(tab, 5), row.names = FALSE)
  cat(sprintf("\n  Overall ACR = %.5f\n", x$overall_acr))
  cat("  effect(d): effect of a dose change d vs quasi-stayers (dD~0).\n")
  cat("  ACR(d): average causal response (marginal effect) at change d.\n")
  invisible(x)
}
