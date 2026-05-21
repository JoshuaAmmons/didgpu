# ============================================================================
# Callaway & Sant'Anna (2021): "Difference-in-Differences with Multiple Time
# Periods" — the dominant DiD framework when treatment is staggered across
# cohorts.
#
# Reference: Callaway, B. and Sant'Anna, P. (2021), "Difference-in-Differences
# with multiple time periods", Journal of Econometrics 225(2): 200-230.
# R package: `did` by Brantly Callaway (https://github.com/bcallaway11/did).
#
# Core object: ATT(g, t) = the average treatment effect on units first
# treated at time g, evaluated at calendar time t (t >= g). The CS framework
# decomposes the panel into these (g, t) cells, estimates each one robustly
# (using only never-treated or not-yet-treated units as controls — never
# cohorts that switched at different times), and then aggregates ATT(g, t)
# into more interpretable summaries:
#
#   - Event-study:    ATT^es(e) = weighted avg of ATT(g, g+e) over g
#   - Group-specific: ATT^g(g) = weighted avg of ATT(g, t) for t >= g
#   - Calendar time:  ATT^c(t) = weighted avg of ATT(g, t) for g <= t
#   - Overall:        ATT^O = weighted avg over all (g, t) cells
#
# Three inner estimators for each ATT(g, t):
#   - "OR"  (outcome regression): regress Y change on covariates, controls only
#   - "IPW" (inverse probability weighting): weight by Pr(treated | X)
#   - "DR"  (doubly robust): combines OR + IPW; preferred in practice
#
# Control groups:
#   - "never": never-treated only (simplest, requires never-treated to exist)
#   - "notyet": not-yet-treated as of time t (more flexible)
#
# GPU acceleration targets:
#   1. The per-(g, t) inner regressions — embarrassingly parallel across cells
#   2. The bootstrap loop (B=1000 typical for inference)
#   3. The aggregation step (linear combinations across (g, t) cells)
#
# STATUS: scaffolded; OR + never-treated implementation is real (see
# R/cs_or.R), other methods stub out with a clear roadmap.
# ============================================================================


#' Estimate Callaway-Sant'Anna (2021) staggered-treatment DiD
#'
#' Estimates group-time average treatment effects ATT(g, t) for each
#' cohort g (defined by first-treatment period) and each post-treatment
#' time t >= g, then aggregates into the requested summary.
#'
#' Three inner estimators (`est_method`):
#' \itemize{
#'   \item `"OR"` (outcome regression, default): fit a linear model of
#'         the outcome change on covariates among controls only;
#'         predict counterfactual for treated; ATT = mean of
#'         (observed - predicted) over treated cohort cells.
#'   \item `"IPW"` (inverse-probability weighting): weight treated and
#'         control cells by inverse propensity score.
#'   \item `"DR"` (doubly robust): combines OR and IPW; consistent if
#'         either model is correct. Preferred in practice.
#' }
#'
#' Four aggregation schemes (`aggregation`):
#' \itemize{
#'   \item `"event"` (default): event-study, indexed by event-time e = t - g.
#'   \item `"group"`: per-cohort average effect over post-treatment periods.
#'   \item `"calendar"`: per-calendar-time average over cohorts already treated.
#'   \item `"overall"`: single scalar summary, the weighted average of
#'         all post-treatment (g, t) cells.
#' }
#'
#' @param df A panel data.frame.
#' @param outcome,group,time,treatment Column names. `treatment` must be
#'   binary (0/1).
#' @param control_group Either `"never"` (never-treated; default) or
#'   `"notyet"` (not-yet-treated).
#' @param est_method One of `"OR"`, `"IPW"`, `"DR"`. Default `"OR"`.
#' @param aggregation One of `"event"`, `"group"`, `"calendar"`,
#'   `"overall"`. Default `"event"`. Use [didgpu_cs_aggregate()] to
#'   compute additional aggregations from a single fit.
#' @param covariates Optional character vector of time-invariant
#'   covariate column names for OR / IPW / DR adjustment.
#' @param bootstrap_reps Integer. Number of bootstrap reps for SE
#'   estimation. Default `0L` (point estimate only).
#' @param bootstrap_kind One of `"cluster"` (resample units, refit per
#'   rep) or `"multiplier"` (multiplier wild bootstrap on the
#'   influence functions; much faster for large B). Default `"cluster"`.
#' @param ci_level Numeric in (0, 100). Default `95`.
#' @param seed Integer.
#' @param backend One of `"auto"`, `"r"`, `"cuda"`.
#' @param verbose Logical. Print progress per (g, t).
#' @return An object of class `didgpu_cs_result`:
#'   \itemize{
#'     \item `att_gt`: long-form data.frame of ATT(g, t) estimates.
#'     \item `aggregation`: the chosen aggregation summary
#'           (event-study / group / calendar / overall).
#'     \item `args`: the canonical args bundle (for re-aggregation).
#'   }
#'
#' @references
#' Callaway, B. and Sant'Anna, P. (2021). "Difference-in-Differences
#' with multiple time periods." *Journal of Econometrics* 225(2): 200-230.
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 100L, n_periods = 12L,
#'                             tau_profile = c(0.5, 1.0), seed = 17L)
#' fit <- didgpu_cs(p, "Y", "unit", "period", "D",
#'                   est_method = "OR", aggregation = "event",
#'                   bootstrap_reps = 0L, backend = "r", verbose = FALSE)
#' print(fit)
#' }
#' @export
didgpu_cs <- function(
    df, outcome, group, time, treatment,
    control_group = c("never", "notyet"),
    est_method    = c("OR", "IPW", "DR"),
    aggregation   = c("event", "group", "calendar", "overall"),
    covariates    = NULL,
    bootstrap_reps = 0L,
    bootstrap_kind = c("cluster", "multiplier"),
    ci_level      = 95,
    seed          = 1L,
    backend       = "auto",
    verbose       = TRUE) {
  bootstrap_kind <- match.arg(bootstrap_kind)

  control_group <- match.arg(control_group)
  est_method    <- match.arg(est_method)
  aggregation   <- match.arg(aggregation)

  # Arg validation.
  stopifnot(is.data.frame(df) || data.table::is.data.table(df))
  for (nm in c("outcome", "group", "time", "treatment")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v)) {
      stop("`", nm, "` must be a single non-empty character column name.")
    }
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  if (!is.null(covariates)) {
    miss <- setdiff(covariates, names(df))
    if (length(miss)) stop("covariates not in df: ", paste(miss, collapse = ", "))
  }
  bootstrap_reps <- as.integer(bootstrap_reps)
  seed <- as.integer(seed)

  # Resolve backend.
  cuda_ok <- isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE))
  resolved_backend <- if (backend %in% c("auto", "cuda") && cuda_ok) "cuda"
                       else "r"

  args <- list(
    outcome = outcome, group = group, time = time, treatment = treatment,
    control_group = control_group, est_method = est_method,
    aggregation = aggregation, covariates = covariates,
    bootstrap_reps = bootstrap_reps,
    bootstrap_kind = bootstrap_kind,
    ci_level = ci_level,
    seed = seed, backend = resolved_backend
  )

  # All three methods (OR / IPW / DR) and both control groups
  # (never / notyet) now go through the unified .cs_compute_att_gt
  # dispatcher. Covariates supported across all three.
  att_gt <- .cs_compute_att_gt(df, args, verbose = verbose)
  if (bootstrap_reps > 0L) {
    boot_kind <- args$bootstrap_kind %||% "cluster"
    if (boot_kind == "multiplier") {
      att_gt <- .cs_multiplier_bootstrap_se(att_gt, args, verbose = verbose)
    } else if (identical(args$backend, "cuda")) {
      # GPU IF-shortcut cluster bootstrap. Falls back to the per-rep
      # R recomputation if CUDA returns NULL.
      att_gt_cuda <- .cs_cluster_bootstrap_cuda(att_gt, args)
      att_gt <- if (!is.null(att_gt_cuda)) att_gt_cuda
                else .cs_bootstrap_se(df, args, att_gt, verbose = verbose)
    } else {
      att_gt <- .cs_bootstrap_se(df, args, att_gt, verbose = verbose)
    }
  }
  agg <- .cs_aggregate(att_gt, aggregation, args)
  placebo <- .cs_placebo_test(att_gt, args)
  out <- structure(list(
    att_gt      = att_gt,
    aggregation = agg,
    placebo     = placebo,
    args        = args
  ), class = c("didgpu_cs_result", "list"))
  return(out)
}


#' Re-aggregate a fitted didgpu_cs result with a different scheme
#'
#' @param fit A `didgpu_cs_result` object.
#' @param aggregation One of `"event"`, `"group"`, `"calendar"`, `"overall"`.
#' @return The same `didgpu_cs_result` with a new `aggregation` slot.
#' @export
didgpu_cs_aggregate <- function(fit, aggregation = c("event", "group",
                                                       "calendar", "overall")) {
  stopifnot(inherits(fit, "didgpu_cs_result"))
  aggregation <- match.arg(aggregation)
  agg <- .cs_aggregate(fit$att_gt, aggregation, fit$args)
  fit$aggregation <- agg
  fit$args$aggregation <- aggregation
  fit
}


#' Print method for didgpu_cs_result
#' @param x A `didgpu_cs_result`.
#' @param ... Unused.
#' @return The input invisibly.
#' @export
print.didgpu_cs_result <- function(x, ...) {
  cat(sprintf("Callaway-Sant'Anna estimate (method = '%s', controls = '%s')\n",
              x$args$est_method, x$args$control_group))
  cat(sprintf("  %d ATT(g, t) cells across %d cohort(s)\n",
              nrow(x$att_gt), length(unique(x$att_gt$g))))
  cat(sprintf("  bootstrap reps  : %d\n", x$args$bootstrap_reps %||% 0L))
  cat(sprintf("\nAggregation: %s\n", x$args$aggregation))
  print(x$aggregation)
  cat("\nFor a different aggregation: didgpu_cs_aggregate(fit, 'event'|'group'|'calendar'|'overall')\n")
  invisible(x)
}


.cs_not_implemented <- function(est_method, control_group) {
  stop(sprintf(
    paste("didgpu_cs(est_method = '%s', control_group = '%s') is",
          "scaffolded but not yet implemented.\nv1 supports OR + never-",
          "treated controls + no covariates. Roadmap:\n",
          "  1. OR with covariates (linear regression on Y change ~ X)\n",
          "  2. IPW (inverse propensity weight)\n",
          "  3. DR (doubly-robust, the preferred CS estimator)\n",
          "  4. control_group = 'notyet' (use not-yet-treated as controls)\n",
          "Reference: did package by Brantly Callaway."),
    est_method, control_group), call. = FALSE)
}
