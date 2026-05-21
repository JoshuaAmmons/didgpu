# ============================================================================
# didgpu_twfe(): the naive two-way fixed-effects (TWFE) dynamic event study.
#
# The universal "naive" baseline reported alongside the modern robust
# estimators (didgpu / didgpu_cs / didgpu_fect) in virtually every
# applied DiD paper. For a NON-ABSORBING binary treatment it takes the
# distributed-lag form
#
#   Y_it = alpha_i + lambda_t
#          + sum_{k=1}^{effects} beta_k  D_{i, t-(k-1)}   (Effect_k)
#          + sum_{j=1}^{placebo} gamma_j D_{i, t+j}       (Placebo_j)
#          + eps_it
#
# where Effect_1 is the contemporaneous (same-period) treatment
# coefficient and Effect_k the coefficient on the (k-1)-period lag;
# Placebo_j is the j-period LEAD (a pre-trend / anticipation check).
#
# IMPORTANT: this is a deliberately naive estimator. Under
# heterogeneous treatment effects with staggered timing, TWFE
# event-study coefficients are contaminated by "forbidden
# comparisons" (Sun & Abraham 2021; de Chaisemartin & D'Haultfoeuille
# 2020). Report it as a baseline and compare it to didgpu() /
# didgpu_cs(), which are robust to that bias.
#
# Output mirrors didgpu()'s Effects / Placebos matrices so the two
# slot into the same tables and plots.
# ============================================================================


#' Naive two-way fixed-effects (TWFE) dynamic event study
#'
#' Fits the distributed-lag TWFE specification with unit and time
#' fixed effects and cluster-robust standard errors. Intended as the
#' "naive" baseline to report next to the bias-robust estimators
#' [didgpu()] and [didgpu_cs()].
#'
#' @param df A data.frame / data.table panel.
#' @param outcome,group,time,treatment Character column names: the
#'   outcome, unit id, time id, and (binary, possibly non-absorbing)
#'   treatment indicator.
#' @param effects Integer >= 1. Number of post/contemporaneous lag
#'   terms (Effect_1 = contemporaneous, Effect_k = lag k-1).
#' @param placebo Integer >= 0. Number of lead terms (Placebo_j =
#'   lead j), the pre-trend / anticipation checks.
#' @param cluster Optional character column name to cluster SEs on.
#'   Defaults to `group`.
#' @param tol Convergence tolerance for the iterative two-way
#'   fixed-effect demeaning. Default `1e-10`.
#' @param max_iter Max demeaning iterations. Default `1000`.
#' @param verbose Logical.
#' @return An object of class `didgpu_twfe_result`: a list with
#'   `coef` (named vector), `results` (with `Effects` and `Placebos`
#'   matrices: Estimate / SE / LB.CI / UB.CI / N), and `args`.
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
#'                             tau_profile = c(0.5, 1.0), seed = 17L)
#' fit <- didgpu_twfe(p, "Y", "unit", "period", "D",
#'                     effects = 3L, placebo = 2L, verbose = FALSE)
#' print(fit)
#' }
#' @export
didgpu_twfe <- function(df, outcome, group, time, treatment,
                        effects = 1L, placebo = 0L,
                        cluster = NULL, tol = 1e-10,
                        max_iter = 1000L, verbose = TRUE) {
  stopifnot(is.data.frame(df) || data.table::is.data.table(df))
  for (nm in c("outcome", "group", "time", "treatment")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v))
      stop("`", nm, "` must be a single non-empty column name.")
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  effects <- as.integer(effects)
  placebo <- as.integer(placebo)
  if (effects < 1L) stop("`effects` must be >= 1.")
  if (placebo < 0L) stop("`placebo` must be >= 0.")
  cluster_col <- cluster %||% group
  if (!cluster_col %in% names(df)) stop("cluster column not in df: ", cluster_col)

  # Build the working table by direct extraction so a cluster column
  # that coincides with group/time/treatment doesn't create duplicate
  # names (the common cluster = group default).
  src_df <- data.table::as.data.table(df)
  d <- data.table::data.table(
    Y_XX  = src_df[[outcome]],
    G_XX  = src_df[[group]],
    T_XX  = src_df[[time]],
    D_XX  = src_df[[treatment]],
    CL_XX = src_df[[cluster_col]])
  d <- d[!is.na(Y_XX) & !is.na(D_XX)]
  data.table::setkey(d, G_XX, T_XX)

  # --- Build the lead/lag treatment regressors by TIME VALUE (gap-safe). ---
  # Effect_k regressor = D_{i, t-(k-1)};  Placebo_j regressor = D_{i, t+j}.
  reg_names <- character(0)
  for (k in seq_len(effects)) {
    lagamt <- k - 1L
    nm <- paste0("Effect_", k)
    # src maps (group, original_t + lagamt) -> D at original_t, so a
    # join on (group, t) attaches D_{i, t - lagamt} to each row.
    src <- d[, list(G_XX, T_XX = T_XX + lagamt, vv = D_XX)]
    d[src, (nm) := i.vv, on = c("G_XX", "T_XX")]
    d[is.na(get(nm)), (nm) := 0]
    reg_names <- c(reg_names, nm)
  }
  for (j in seq_len(placebo)) {
    nm <- paste0("Placebo_", j)
    src <- d[, list(G_XX, T_XX = T_XX - j, vv = D_XX)]
    d[src, (nm) := i.vv, on = c("G_XX", "T_XX")]
    d[is.na(get(nm)), (nm) := 0]
    reg_names <- c(reg_names, nm)
  }

  # --- Two-way FE absorption via iterative demeaning (Gaure 2013). ---
  cols <- c("Y_XX", reg_names)
  M <- as.matrix(d[, cols, with = FALSE])
  g_idx <- as.integer(factor(d$G_XX))
  t_idx <- as.integer(factor(d$T_XX))
  Md <- .twfe_demean(M, g_idx, t_idx, tol = tol, max_iter = max_iter)

  Yd <- Md[, 1L]
  Xd <- Md[, -1L, drop = FALSE]
  colnames(Xd) <- reg_names

  # --- OLS on the demeaned data. ---
  XtX <- crossprod(Xd)
  XtX_inv <- tryCatch(solve(XtX), error = function(e)
    stop("TWFE design is rank-deficient (collinear lead/lag terms). ",
         "Reduce `effects`/`placebo` or check the panel."))
  beta <- as.numeric(XtX_inv %*% crossprod(Xd, Yd))
  names(beta) <- reg_names
  resid <- as.numeric(Yd - Xd %*% beta)

  # --- Cluster-robust (CR1) variance. ---
  n_obs   <- nrow(Xd)
  n_units <- length(unique(g_idx))
  n_times <- length(unique(t_idx))
  k_reg   <- ncol(Xd)
  # Params: regressors + unit FE + time FE - 1 (shared intercept).
  k_params <- k_reg + n_units + (n_times - 1L)
  vcov <- .twfe_cluster_vcov(Xd, resid, d$CL_XX, XtX_inv,
                             n_obs, k_params)
  se <- sqrt(pmax(diag(vcov), 0))
  names(se) <- reg_names

  z <- stats::qnorm(0.975)
  mk_mat <- function(nms) {
    if (length(nms) == 0L)
      return(matrix(numeric(0), nrow = 0, ncol = 5,
                    dimnames = list(NULL,
                      c("Estimate", "SE", "LB.CI", "UB.CI", "N"))))
    est <- beta[nms]; s <- se[nms]
    m <- cbind(Estimate = est, SE = s,
               LB.CI = est - z * s, UB.CI = est + z * s,
               N = n_obs)
    rownames(m) <- nms
    m
  }
  eff_nms <- paste0("Effect_", seq_len(effects))
  pl_nms  <- if (placebo > 0L) paste0("Placebo_", seq_len(placebo)) else character(0)

  out <- structure(list(
    coef    = beta,
    results = list(Effects = mk_mat(eff_nms),
                   Placebos = mk_mat(pl_nms),
                   vcov = vcov),
    args = list(outcome = outcome, group = group, time = time,
                treatment = treatment, effects = effects,
                placebo = placebo, cluster = cluster_col,
                n_obs = n_obs, n_units = n_units, n_times = n_times)
  ), class = c("didgpu_twfe_result", "list"))
  if (verbose) print(out)
  invisible(out)
}


# Iterative two-way demeaning (alternating-projections / Gaure 2013).
# Removes unit means then time means from every column of M until the
# max column-wise change falls below tol. Matches the fixed-effect
# absorption used by lm(y ~ . + factor(unit) + factor(time)) to
# numerical precision.
#' @keywords internal
#' @noRd
.twfe_demean <- function(M, g_idx, t_idx, tol = 1e-10, max_iter = 1000L) {
  ng <- max(g_idx); nt <- max(t_idx)
  gcount <- tabulate(g_idx, ng); tcount <- tabulate(t_idx, nt)
  for (iter in seq_len(max_iter)) {
    prev <- M
    # subtract unit means (reorder=TRUE -> rows 1..ng, so direct index)
    gm <- rowsum(M, g_idx, reorder = TRUE) / gcount
    M <- M - gm[g_idx, , drop = FALSE]
    # subtract time means
    tm <- rowsum(M, t_idx, reorder = TRUE) / tcount
    M <- M - tm[t_idx, , drop = FALSE]
    if (max(abs(M - prev)) < tol) break
  }
  M
}


# CR1 cluster-robust sandwich variance:
#   V = (X'X)^{-1} [ sum_g (X_g' u_g)(X_g' u_g)' ] (X'X)^{-1} * c
# with the standard small-sample correction
#   c = (G/(G-1)) * ((N-1)/(N-K)).
#' @keywords internal
#' @noRd
.twfe_cluster_vcov <- function(Xd, resid, cluster, XtX_inv, n_obs, k_params) {
  cl <- as.integer(factor(cluster))
  G <- length(unique(cl))
  k <- ncol(Xd)
  meat <- matrix(0, k, k)
  Xu <- Xd * resid                       # row-wise X_i * u_i
  sg <- rowsum(Xu, cl, reorder = FALSE)  # per-cluster sum of X_i u_i (G x k)
  for (gi in seq_len(nrow(sg))) {
    s <- sg[gi, ]
    meat <- meat + tcrossprod(s)
  }
  cc <- (G / (G - 1)) * ((n_obs - 1) / (n_obs - k_params))
  V <- XtX_inv %*% meat %*% XtX_inv * cc
  (V + t(V)) / 2                          # symmetrize
}


#' Print method for didgpu_twfe_result
#' @param x A `didgpu_twfe_result`.
#' @param ... Unused.
#' @return The input invisibly.
#' @export
print.didgpu_twfe_result <- function(x, ...) {
  cat("Two-way FE dynamic event study (naive TWFE baseline)\n")
  cat(sprintf("  %d obs, %d units, %d periods; clustered on '%s'\n",
              x$args$n_obs, x$args$n_units, x$args$n_times, x$args$cluster))
  cat("\nEffects (Effect_1 = contemporaneous):\n")
  print(round(x$results$Effects[, c("Estimate", "SE", "LB.CI", "UB.CI"),
                                drop = FALSE], 4))
  if (nrow(x$results$Placebos) > 0L) {
    cat("\nPlacebos (leads; pre-trend check):\n")
    print(round(x$results$Placebos[, c("Estimate", "SE", "LB.CI", "UB.CI"),
                                   drop = FALSE], 4))
  }
  cat("\nNOTE: TWFE is a naive baseline, biased under heterogeneous\n")
  cat("effects. Compare to didgpu() / didgpu_cs() for robust estimates.\n")
  invisible(x)
}
