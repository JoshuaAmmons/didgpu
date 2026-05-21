# ============================================================================
# didgpu_freyaldenhoven(): Freyaldenhoven, Hansen & Shapiro (2019) pre-event
# panel event-study estimator with an auxiliary "proxy" covariate.
#
# Standard event studies attribute confounding pre-trends to the policy. FHS
# add a covariate x_it that responds to the SAME confound but is NOT affected
# by the policy. Including x (instrumented, because x is itself endogenous)
# purges the confound, so the remaining event-study coefficients are
# consistent even when a confound generates pre-trends.
#
# We reimplement the eventstudyr parameterization exactly so coefficients
# match that reference:
#   * First-difference event study: z_fd = z - lag(z); interior leads/lags as
#     z_fd_lead{k} / z_fd_lag{k}; endpoint LEVELS z_lead{L} (entered as
#     1 - lead) and z_lag{F}; one term normalized out (omitted).
#       num_fd_leads = pre + overidpre,  num_fd_lags = post + overidpost - 1,
#       furthest_lag = num_fd_lags + 1.
#   * estimator = "OLS": TWFE (unit + time FE) on that design.
#   * estimator = "FHS": 2SLS with the proxy as an endogenous regressor,
#     instrumented by proxyIV (default: the z_fd_lead with the strongest
#     first-stage F-statistic for the proxy).
# Cluster-robust (CR1, on the unit) SEs. Coefficients match
# eventstudyr::EventStudy; SEs use CR1 (eventstudyr defaults to estimatr CR2).
#
# Reference: Freyaldenhoven, Hansen & Shapiro (2019), "Pre-event Trends in
# the Panel Event-Study Design", American Economic Review 109(9): 3307-3338.
# ============================================================================


# Build the FD-event-study design columns on a data.table (cols: id, t, z).
# Returns list(dt, terms) where terms is the ordered regressor set (minus the
# normalized column). Mirrors eventstudyr's EventStudy() construction.
#' @keywords internal
#' @noRd
.fhs_build_design <- function(dt, idvar, timevar, policyvar,
                              pre, post, overidpre, overidpost, normalize) {
  data.table::setorderv(dt, c(idvar, timevar))
  num_fd_lags  <- post + overidpost - 1L
  num_fd_leads <- pre + overidpre
  furthest_lag <- num_fd_lags + 1L
  pvar <- policyvar          # local name chosen NOT to collide with a column
  shift <- data.table::shift
  fd <- paste0(pvar, "_fd")

  # First difference and its leads/lags. (pvar/fd are read via get(); they do
  # not shadow a data column because no column is named "pvar"/"fd".)
  dt[, (fd) := get(pvar) - shift(get(pvar), 1L, type = "lag"), by = idvar]
  if (num_fd_leads >= 1L) for (k in seq_len(num_fd_leads))
    dt[, (paste0(pvar, "_fd_lead", k)) := shift(get(fd), k, type = "lead"), by = idvar]
  if (num_fd_lags >= 1L) for (k in seq_len(num_fd_lags))
    dt[, (paste0(pvar, "_fd_lag", k)) := shift(get(fd), k, type = "lag"), by = idvar]

  # Endpoint LEVELS: lead endpoint entered as (1 - lead); lag endpoint as level.
  lead_ep <- paste0(pvar, "_lead", num_fd_leads)
  lag_ep  <- paste0(pvar, "_lag", furthest_lag)
  dt[, (lead_ep) := 1 - shift(get(pvar), num_fd_leads, type = "lead"), by = idvar]
  dt[, (lag_ep)  := shift(get(pvar), furthest_lag, type = "lag"), by = idvar]

  # Normalized (omitted) column, per eventstudyr's rule.
  if (normalize < 0) {
    if (normalize == -(pre + overidpre + 1L))
      norm_col <- paste0(pvar, "_lead", -(normalize + 1L))
    else
      norm_col <- paste0(pvar, "_fd_lead", -normalize)
  } else if (normalize == 0) {
    norm_col <- if (normalize == post + overidpost) paste0(pvar, "_lag", normalize)
                else fd
  } else {
    norm_col <- if (normalize == post + overidpost) paste0(pvar, "_lag", normalize)
                else paste0(pvar, "_fd_lag", normalize)
  }

  lead_fd <- if (num_fd_leads >= 1L) paste0(pvar, "_fd_lead", num_fd_leads:1L) else character(0)
  lag_fd  <- if (num_fd_lags  >= 1L) paste0(pvar, "_fd_lag",  1:num_fd_lags)   else character(0)
  terms <- c(lead_ep, lead_fd, fd, lag_fd, lag_ep)
  terms <- terms[terms != norm_col]
  list(dt = dt, terms = terms, norm_col = norm_col)
}


# Cluster-robust (CR1) sandwich for an estimator with "hat" regressors Xhat
# (= X for OLS; = projected X for 2SLS), bread (X'Xhat)^-1, residuals u.
#' @keywords internal
#' @noRd
.fhs_cluster_vcov <- function(bread, Xhat, resid, cl, n_obs, k) {
  cli <- as.integer(factor(cl))
  G <- length(unique(cli))
  Xu <- Xhat * resid
  sg <- rowsum(Xu, cli, reorder = FALSE)
  meat <- crossprod(sg)                      # sum_g (Xhat_g' u_g)(.)'
  cc <- (G / (G - 1)) * ((n_obs - 1) / (n_obs - k))
  V <- bread %*% meat %*% t(bread) * cc
  (V + t(V)) / 2
}


#' Freyaldenhoven-Hansen-Shapiro (2019) pre-event proxy event study
#'
#' Panel event-study estimator that uses an auxiliary covariate (a "proxy"
#' affected by the confound but not the policy) to correct for confounding
#' pre-trends. `estimator = "OLS"` is the plain two-way FE event study;
#' `estimator = "FHS"` adds the proxy as an endogenous regressor and
#' instruments it (2SLS) with a far policy lead, purging the confound.
#' Reimplements `eventstudyr::EventStudy`'s first-difference parameterization,
#' so the event-study coefficients match that reference.
#'
#' @param df A data.frame / data.table panel.
#' @param outcome,policy,id,time Character column names: outcome, the (binary
#'   or continuous) policy variable, unit id, and integer time.
#' @param estimator `"OLS"` (default) or `"FHS"`.
#' @param proxy For `"FHS"`, the character name of the proxy covariate.
#' @param proxyIV For `"FHS"`, the instrument column. Default: the
#'   first-differenced policy lead with the strongest first-stage F.
#' @param pre,post Non-negative integers: anticipation leads (`pre`) and
#'   dynamic lags (`post`).
#' @param overidpre,overidpost Extra leads/lags (over-identification /
#'   endpoints). Defaults mirror eventstudyr (`overidpost = 1`,
#'   `overidpre = pre + post`).
#' @param normalize Event-time coefficient to omit (normalized to 0). Default
#'   `-(pre + 1)`.
#' @param cluster Optional unit-level cluster column for SEs. Defaults to `id`.
#' @param tol,max_iter Two-way demeaning controls.
#' @param verbose Logical.
#' @return A `didgpu_freyaldenhoven_result`: a list with `coefficients` (a
#'   matrix: Estimate / SE / LB.CI / UB.CI per event-study term, plus the
#'   proxy for FHS), `estimator`, `proxyIV`, and `args`.
#' @references Freyaldenhoven, S., Hansen, C. & Shapiro, J.M. (2019).
#'   Pre-event Trends in the Panel Event-Study Design. \emph{American Economic
#'   Review} 109(9): 3307-3338.
#' @examples
#' \donttest{
#' # eventstudyr's example data: policy z, outcome y_base, proxy x_r.
#' if (requireNamespace("eventstudyr", quietly = TRUE)) {
#'   d <- eventstudyr::example_data
#'   didgpu_freyaldenhoven(d, "y_base", "z", "id", "t",
#'                         estimator = "FHS", proxy = "x_r", pre = 0, post = 3)
#' }
#' }
#' @export
didgpu_freyaldenhoven <- function(df, outcome, policy, id, time,
                                  estimator = c("OLS", "FHS"),
                                  proxy = NULL, proxyIV = NULL,
                                  pre = 0L, post = 1L,
                                  overidpre = pre + post, overidpost = 1L,
                                  normalize = -(pre + 1L),
                                  cluster = NULL, tol = 1e-10,
                                  max_iter = 1000L, verbose = TRUE) {
  estimator <- match.arg(estimator)
  for (nm in c("outcome", "policy", "id", "time")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v))
      stop("`", nm, "` must be a single non-empty column name.")
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  if (estimator == "FHS" && is.null(proxy))
    stop("`proxy` must be supplied when estimator = 'FHS'.")
  pre <- as.integer(pre); post <- as.integer(post)
  overidpre <- as.integer(overidpre); overidpost <- as.integer(overidpost)
  normalize <- as.integer(normalize)
  cluster_col <- cluster %||% id

  dt <- data.table::as.data.table(df)
  dt <- dt[!is.na(get(outcome)) & !is.na(get(policy)) & !is.na(get(time))]
  des <- .fhs_build_design(dt, id, time, policy,
                           pre, post, overidpre, overidpost, normalize)
  dt <- des$dt; terms <- des$terms

  # FHS: pick the proxyIV (strongest first-stage lead) if unspecified.
  if (estimator == "FHS" && is.null(proxyIV)) {
    fd_leads <- grep(paste0("^", policy, "_fd_lead"), names(dt), value = TRUE)
    fd_leads <- fd_leads[fd_leads %in% terms]
    bestF <- -Inf
    for (v in fd_leads) {
      sub <- dt[!is.na(get(proxy)) & !is.na(get(v))]
      fm <- stats::lm(stats::reformulate(v, proxy), data = sub)
      Fv <- tryCatch(summary(fm)$fstatistic[["value"]], error = function(e) NA_real_)
      if (is.finite(Fv) && Fv > bestF) { bestF <- Fv; proxyIV <- v }
    }
    if (is.null(proxyIV)) stop("could not auto-select a proxyIV; specify one.")
  }

  # Assemble the modelling columns and drop rows with any NA (matches the
  # reference's listwise deletion at panel edges).
  exog <- terms
  if (estimator == "FHS") exog <- terms[terms != proxyIV]
  used <- c(outcome, exog, if (estimator == "FHS") c(proxy, proxyIV))
  used <- unique(used)
  mdt <- dt[stats::complete.cases(dt[, used, with = FALSE])]
  if (nrow(mdt) == 0L) stop("no complete rows after building leads/lags.")

  g_idx <- as.integer(factor(mdt[[id]]))
  t_idx <- as.integer(factor(mdt[[time]]))
  cl    <- mdt[[cluster_col]]

  # Two-way FE absorption via the shared iterative demeaner.
  Mcols <- c(outcome, exog, if (estimator == "FHS") c(proxy, proxyIV))
  M <- as.matrix(mdt[, Mcols, with = FALSE])
  Md <- .twfe_demean(M, g_idx, t_idx, tol = tol, max_iter = max_iter)
  yd <- Md[, 1L]
  n_obs <- nrow(Md)
  n_units <- length(unique(g_idx)); n_times <- length(unique(t_idx))

  rankfail <- function(e) stop(
    "FHS/OLS design is rank-deficient (collinear event-study terms after ",
    "FE absorption). The panel is likely too short for the requested ",
    "pre/post/overid window, or the policy lacks variation. Reduce pre/post ",
    "or widen the panel.")
  if (estimator == "OLS") {
    X <- Md[, exog, drop = FALSE]
    XtXi <- tryCatch(solve(crossprod(X)), error = rankfail)
    beta <- as.numeric(XtXi %*% crossprod(X, yd)); names(beta) <- exog
    resid <- as.numeric(yd - X %*% beta)
    kpar <- ncol(X) + n_units + (n_times - 1L)
    V <- .fhs_cluster_vcov(XtXi, X, resid, cl, n_obs, kpar)
    coef_names <- exog
  } else {
    Xe <- Md[, exog, drop = FALSE]
    p  <- Md[, proxy]; w <- Md[, proxyIV]
    X  <- cbind(Xe, proxy = p)                 # regressors (proxy endogenous)
    Z  <- cbind(Xe, proxyIV = w)               # instruments (proxyIV excluded)
    ZtZi <- tryCatch(solve(crossprod(Z)), error = rankfail)
    PZX  <- Z %*% (ZtZi %*% crossprod(Z, X))   # projection of X onto Z
    bread <- tryCatch(solve(crossprod(X, PZX)), error = rankfail)  # (X'P_Z X)^-1
    beta <- as.numeric(bread %*% crossprod(PZX, yd))
    names(beta) <- c(exog, proxy)
    resid <- as.numeric(yd - X %*% beta)
    kpar <- ncol(X) + n_units + (n_times - 1L)
    V <- .fhs_cluster_vcov(bread, PZX, resid, cl, n_obs, kpar)
    coef_names <- c(exog, proxy)
  }

  se <- sqrt(pmax(diag(V), 0)); names(se) <- coef_names
  z975 <- stats::qnorm(0.975)
  cmat <- cbind(Estimate = beta, SE = se,
                LB.CI = beta - z975 * se, UB.CI = beta + z975 * se)
  rownames(cmat) <- coef_names

  out <- structure(list(
    coefficients = cmat, estimator = estimator, proxyIV = proxyIV,
    norm_col = des$norm_col,
    args = list(outcome = outcome, policy = policy, id = id, time = time,
                proxy = proxy, pre = pre, post = post, overidpre = overidpre,
                overidpost = overidpost, normalize = normalize,
                cluster = cluster_col, n_obs = n_obs, n_units = n_units)
  ), class = c("didgpu_freyaldenhoven_result", "list"))
  if (verbose) print(out)
  invisible(out)
}


#' Print method for didgpu_freyaldenhoven_result
#' @param x A `didgpu_freyaldenhoven_result`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.didgpu_freyaldenhoven_result <- function(x, ...) {
  cat(sprintf("Freyaldenhoven-Hansen-Shapiro (2019) event study [%s]\n",
              x$estimator))
  cat(sprintf("  %d obs, %d units; normalized term: %s",
              x$args$n_obs, x$args$n_units, x$norm_col))
  if (x$estimator == "FHS")
    cat(sprintf("; proxy = %s, proxyIV = %s", x$args$proxy, x$proxyIV))
  cat("\n\n")
  print(round(x$coefficients, 5))
  if (x$estimator == "OLS")
    cat("\n  OLS = plain TWFE event study (no confound correction).\n")
  else
    cat("\n  FHS: proxy instrumented to purge a pre-trend confound.\n")
  invisible(x)
}
