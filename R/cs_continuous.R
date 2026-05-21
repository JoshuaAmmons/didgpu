# ============================================================================
# didgpu_cs_continuous(): Callaway, Goodman-Bacon & Sant'Anna (2024)
# difference-in-differences with a CONTINUOUS treatment (dose).
#
# With a continuous dose D and a comparison group (never- / not-yet-treated,
# D = 0), the method estimates the dose-response curve nonparametrically:
#   ATT(d)  = E[dY | D = d] - E[dY | D = 0]        (level effect at dose d)
#   ACRT(d) = d/dd ATT(d)                          (causal response / slope)
# where dY is the within-unit outcome change from the base to the treatment
# period. E[dY | D = d] is fit by a B-spline regression of dY on a spline
# basis of the (positive) dose; ATT(d) subtracts the comparison group's mean
# change; ACRT(d) is the analytic spline derivative times the slope
# coefficients. This mirrors contdid::cont_did_acrt exactly.
#
# Native reimplementation: the spline BASIS uses splines2 (the same library
# contdid uses, so the basis matches bit-for-bit), but the DiD logic,
# dose-response, ATT/ACRT extraction and inference are didgpu's own. SEs are
# from a unit-level multiplier (Rademacher) bootstrap.
#
# Scope: a single treated cohort vs a never-treated (D = 0) comparison, base
# period = (adoption - 1), treatment period = adoption. This is the canonical
# continuous-DiD design (and matches contdid on 2-period data). Cross-checked
# against contdid::cont_did.
#
# Reference: Callaway, Goodman-Bacon & Sant'Anna (2024), "Difference-in-
# Differences with a Continuous Treatment", NBER WP 32117.
# ============================================================================

#' Knot placement at interior quantiles of the (positive) dose.
#' @keywords internal
#' @noRd
.cd_knots_quantile <- function(x, num_knots) {
  if (num_knots <= 0L) return(numeric(0))
  stats::quantile(x, probs = seq(0, 1, length.out = num_knots + 2L))[
    -c(1L, num_knots + 2L)]
}


#' Callaway-Goodman-Bacon-Sant'Anna (2024) continuous-treatment DiD
#'
#' Estimates the continuous-treatment dose-response: the level effect
#' \eqn{ATT(d)} and the causal response (slope) \eqn{ACRT(d) = ATT'(d)} of a
#' continuous dose, comparing dose-\eqn{d} units to a never-treated (dose 0)
#' comparison group via a within-unit before/after change. The dose-response
#' is fit by a B-spline regression; ATT(d) subtracts the comparison mean
#' change and ACRT(d) is the analytic spline derivative. Standard errors come
#' from a unit-level multiplier bootstrap.
#'
#' @param df A data.frame / data.table panel.
#' @param yname,dname,gname,tname,idname Character column names: outcome,
#'   continuous dose, cohort (0 = never-treated; positive = adoption period),
#'   time, and unit id.
#' @param dvals Numeric dose values at which to report the curve. Default:
#'   quantiles 0.1..0.99 of the positive doses.
#' @param degree B-spline degree (default 3, cubic).
#' @param num_knots Number of interior knots (default 0). Knots are placed at
#'   interior quantiles of the positive dose.
#' @param control_group `"nevertreated"` (default).
#' @param bootstrap_reps Multiplier-bootstrap replicates for SEs (default
#'   200; 0 to skip).
#' @param ci_level Confidence level in percent (default 95).
#' @param seed Bootstrap RNG seed.
#' @param verbose Logical.
#' @return A `didgpu_cs_continuous_result`: a list with `dose` (the dvals),
#'   `att.d`, `acrt.d` (and their SEs / CIs when bootstrapped), `overall_att`,
#'   `overall_acrt`, and `args`.
#' @references Callaway, B., Goodman-Bacon, A. & Sant'Anna, P.H.C. (2024).
#'   Difference-in-Differences with a Continuous Treatment. NBER WP 32117.
#' @seealso [didgpu_cs()] for the binary staggered estimator.
#' @examples
#' \donttest{
#' if (requireNamespace("contdid", quietly = TRUE) &&
#'     requireNamespace("splines2", quietly = TRUE)) {
#'   d <- contdid::simulate_contdid_data(n = 400, num_time_periods = 2)
#'   didgpu_cs_continuous(d, "Y", "D", "G", "time_period", "id",
#'                        degree = 3, num_knots = 2, bootstrap_reps = 0)
#' }
#' }
#' @export
didgpu_cs_continuous <- function(df, yname, dname, gname, tname, idname,
                                 dvals = NULL, degree = 3L, num_knots = 0L,
                                 control_group = "nevertreated",
                                 bootstrap_reps = 200L, ci_level = 95,
                                 seed = 1L, verbose = TRUE) {
  if (!requireNamespace("splines2", quietly = TRUE))
    stop("didgpu_cs_continuous requires the 'splines2' package for the ",
         "B-spline basis. Install it with install.packages('splines2').")
  for (nm in c("yname", "dname", "gname", "tname", "idname")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v))
      stop("`", nm, "` must be a single non-empty column name.")
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  degree <- as.integer(degree); num_knots <- as.integer(num_knots)
  d <- data.table::as.data.table(df)
  d <- d[!is.na(get(yname)) & !is.na(get(gname)) & !is.na(get(tname))]

  # Treated cohort(s): positive gname. Canonical case: a single cohort.
  cohorts <- sort(unique(d[[gname]][d[[gname]] > 0]))
  if (length(cohorts) == 0L) stop("no treated cohort (all gname <= 0).")
  if (length(cohorts) > 1L)
    stop("didgpu_cs_continuous currently supports a single treated cohort ",
         "vs a never-treated comparison. Got cohorts: ",
         paste(cohorts, collapse = ", "), ".")
  g <- cohorts[1L]
  base_period <- g - 1L

  # Per-unit before/after change dY = Y(post=g) - Y(base=g-1); dose zeroed
  # for the comparison group (gname == 0), matching contdid's 2x2 subset.
  wide <- d[get(tname) %in% c(base_period, g),
            list(Y = get(yname)[1L], G = get(gname)[1L], D = get(dname)[1L],
                 dy = get(yname)[match(g, get(tname))] -
                      get(yname)[match(base_period, get(tname))]),
            by = c(idname)]
  wide <- wide[!is.na(dy)]
  wide[, dose := D * (G == g)]
  treated <- wide[dose > 0]
  comp    <- wide[G == 0]
  if (nrow(treated) < degree + num_knots + 1L)
    stop("too few treated units for the requested spline (degree + knots).")
  if (nrow(comp) == 0L) stop("no never-treated (dose 0) comparison units.")

  knots <- .cd_knots_quantile(treated$dose, num_knots)
  if (is.null(dvals))
    dvals <- as.numeric(stats::quantile(treated$dose,
                                        probs = seq(0.1, 0.99, length.out = 5L)))
  comp_mean <- mean(comp$dy)

  # Core dose-response estimator (returns att.d, acrt.d, overall, the spline
  # fit). Reused by the bootstrap with reweighted units.
  fit_doseresp <- function(tdose, tdy, cmean, w_t = NULL) {
    B  <- splines2::bSpline(tdose, degree = degree, knots = knots, intercept = FALSE)
    Bg <- splines2::bSpline(dvals, degree = degree, knots = knots, intercept = FALSE)
    Bd <- splines2::dbs(dvals,   degree = degree, knots = knots, intercept = FALSE)
    X <- cbind(1, B)
    if (is.null(w_t)) {
      bet <- qr.solve(X, tdy)
    } else {
      sw <- sqrt(w_t)
      bet <- qr.solve(X * sw, tdy * sw)
    }
    slope <- bet[-1L]
    att <- as.numeric(cbind(1, Bg) %*% bet) - cmean
    acrt <- as.numeric(Bd %*% slope)
    Bt_all <- splines2::bSpline(tdose, degree = degree, knots = knots, intercept = FALSE)
    Bd_all <- splines2::dbs(tdose, degree = degree, knots = knots, intercept = FALSE)
    att_overall  <- mean(as.numeric(cbind(1, Bt_all) %*% bet)) - cmean
    acrt_overall <- mean(as.numeric(Bd_all %*% slope))
    list(att = att, acrt = acrt, att_overall = att_overall,
         acrt_overall = acrt_overall)
  }

  pt <- fit_doseresp(treated$dose, treated$dy, comp_mean)

  res <- list(dose = dvals, att.d = pt$att, acrt.d = pt$acrt,
              overall_att = pt$att_overall, overall_acrt = pt$acrt_overall,
              att.d_se = NA, acrt.d_se = NA)

  bootstrap_reps <- as.integer(bootstrap_reps)
  if (bootstrap_reps > 0L) {
    nt <- nrow(treated); nc <- nrow(comp)
    set.seed(seed)
    batt <- matrix(NA_real_, bootstrap_reps, length(dvals))
    bacr <- matrix(NA_real_, bootstrap_reps, length(dvals))
    for (b in seq_len(bootstrap_reps)) {
      # Rademacher multipliers reweight treated units; resample the
      # comparison mean via its own multipliers.
      wt <- sample(c(0.5, 1.5), nt, replace = TRUE)   # mean 1, positive weights
      cm <- stats::weighted.mean(comp$dy, sample(c(0.5, 1.5), nc, replace = TRUE))
      fb <- tryCatch(fit_doseresp(treated$dose, treated$dy, cm, w_t = wt),
                     error = function(e) NULL)
      if (!is.null(fb)) { batt[b, ] <- fb$att; bacr[b, ] <- fb$acrt }
    }
    res$att.d_se  <- apply(batt, 2L, stats::sd, na.rm = TRUE)
    res$acrt.d_se <- apply(bacr, 2L, stats::sd, na.rm = TRUE)
    z <- stats::qnorm(1 - (1 - ci_level / 100) / 2)
    res$att.d_lower  <- pt$att  - z * res$att.d_se
    res$att.d_upper  <- pt$att  + z * res$att.d_se
    res$acrt.d_lower <- pt$acrt - z * res$acrt.d_se
    res$acrt.d_upper <- pt$acrt + z * res$acrt.d_se
  }

  res$args <- list(yname = yname, dname = dname, gname = gname, tname = tname,
                   idname = idname, degree = degree, num_knots = num_knots,
                   knots = knots, cohort = g, n_treated = nrow(treated),
                   n_comparison = nrow(comp), ci_level = ci_level)
  class(res) <- c("didgpu_cs_continuous_result", "list")
  if (verbose) print(res)
  invisible(res)
}


#' Print method for didgpu_cs_continuous_result
#' @param x A `didgpu_cs_continuous_result`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.didgpu_cs_continuous_result <- function(x, ...) {
  cat("Callaway-Goodman-Bacon-Sant'Anna (2024) continuous-treatment DiD\n")
  cat(sprintf("  cohort %s vs never-treated; %d treated, %d comparison units\n",
              x$args$cohort, x$args$n_treated, x$args$n_comparison))
  cat(sprintf("  B-spline degree %d, %d interior knots\n",
              x$args$degree, x$args$num_knots))
  tab <- data.frame(dose = x$dose, ATT = x$att.d, ACRT = x$acrt.d)
  if (!all(is.na(x$att.d_se))) { tab$ATT_se <- x$att.d_se; tab$ACRT_se <- x$acrt.d_se }
  cat("\nDose-response:\n"); print(round(tab, 5), row.names = FALSE)
  cat(sprintf("\n  Overall ATT = %.5f   Overall ACRT = %.5f\n",
              x$overall_att, x$overall_acrt))
  cat("  ATT(d): level effect at dose d.  ACRT(d): marginal causal response.\n")
  invisible(x)
}
