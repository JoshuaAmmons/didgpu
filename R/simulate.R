# ============================================================================
# Simulated DGP for validation
#
# Generates a balanced panel with a known event-time treatment profile.
# Used by tests to verify that:
#   (a) didgpu recovers the truth within Monte Carlo error
#   (b) didgpu and DIDmultiplegtDYN agree on the same panel
#       (max-abs-diff < 1e-6 on event-time coefficients)
#
# Design choices, all aimed at the binary non-absorbing case:
#   - Units are randomly assigned to cohorts (defined by first-treatment
#     date `F_g`), including a never-treated cohort `F_g = Inf`.
#   - Treatment is staggered: each treated unit turns on at its `F_g`
#     and stays on (absorbing within this simulator; the user-facing
#     estimator does not require absorption, but the simulator is
#     simpler this way and the estimator handles both).
#   - The outcome is unit FE + time FE + a known event-time effect
#     profile `tau(k)` for k = 0, 1, ..., E, plus Gaussian noise.
#   - Pre-treatment effects (placebos) are zero by construction.
#
# To stress-test edge cases (never-treated only, always-treated, single
# cohort, etc.) the test files compose this generator with overrides.
# ============================================================================


#' Generate a simulated DiD panel with known event-time profile
#'
#' @param n_units integer. Number of units.
#' @param n_periods integer. Number of time periods (1..n_periods).
#' @param frac_treated numeric between 0 and 1. Share of units that ever get
#'   treated. The rest are never-treated controls.
#' @param min_treat_period integer. Earliest treatment-on period.
#' @param max_treat_period integer. Latest treatment-on period. Must be
#'   `<= n_periods`. Treated units' `F_g` is drawn uniformly between
#'   `min_treat_period` and `max_treat_period`.
#' @param tau_profile numeric vector. Event-time effects:
#'   `tau_profile[k+1]` is the effect at event time `k = 0, 1, ...`.
#'   Length `length(tau_profile) - 1` upper bound on k.
#' @param sigma numeric. Idiosyncratic noise SD.
#' @param unit_fe_sd numeric. Unit fixed-effect SD.
#' @param time_fe_sd numeric. Time fixed-effect SD.
#' @param seed integer. RNG seed.
#'
#' @return A data.frame with columns: `unit` (int), `period` (int),
#'   `D` (binary treatment indicator), `Y` (outcome). Sorted by
#'   (unit, period). Also has an attribute `"truth"`: a list with
#'   `F_g` (named numeric per unit; Inf = never-treated),
#'   `tau_profile`, `unit_fe`, `time_fe`.
#'
#' @examples
#' p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L,
#'                             tau_profile = c(0.5, 1.0, 1.2),
#'                             seed = 17L)
#' head(p)
#' # Inspect the underlying DGP:
#' truth <- attr(p, "truth")
#' table(truth$F_g)            # cohort sizes (Inf = never-treated)
#' truth$tau_profile
#' @export
didgpu_simulate_panel <- function(
    n_units          = 100L,
    n_periods        = 20L,
    frac_treated     = 0.5,
    min_treat_period = NULL,
    max_treat_period = NULL,
    tau_profile      = c(0.2, 0.4, 0.6, 0.5, 0.4),
    sigma            = 0.5,
    unit_fe_sd       = 1.0,
    time_fe_sd       = 0.3,
    seed             = 1L) {

  # Default treatment window: middle ~half of the panel. Scales with n_periods
  # so callers can pass tiny panels without having to override these.
  if (is.null(min_treat_period)) min_treat_period <- max(2L, as.integer(n_periods * 0.25))
  if (is.null(max_treat_period)) max_treat_period <- max(min_treat_period,
                                                          as.integer(n_periods * 0.75))

  stopifnot(
    n_units >= 2, n_periods >= 2,
    frac_treated >= 0, frac_treated <= 1,
    min_treat_period >= 1, max_treat_period <= n_periods,
    min_treat_period <= max_treat_period,
    length(tau_profile) >= 1, all(is.finite(tau_profile))
  )

  set.seed(seed)
  units   <- seq_len(n_units)
  periods <- seq_len(n_periods)

  # Assign first-treatment period F_g per unit. Inf == never-treated.
  n_treated <- round(n_units * frac_treated)
  treated_units <- if (n_treated > 0L) {
    sort(sample(units, n_treated))
  } else integer(0)
  F_g <- rep(Inf, n_units)
  names(F_g) <- as.character(units)
  if (n_treated > 0L) {
    F_g[as.character(treated_units)] <- sample(
      min_treat_period:max_treat_period, n_treated, replace = TRUE
    )
  }

  # Fixed effects.
  unit_fe <- stats::rnorm(n_units, 0, unit_fe_sd); names(unit_fe) <- as.character(units)
  time_fe <- stats::rnorm(n_periods, 0, time_fe_sd); names(time_fe) <- as.character(periods)

  # Materialise panel.
  panel <- data.table::CJ(unit = units, period = periods)
  panel[, F_g := F_g[as.character(unit)]]
  panel[, D := as.integer(period >= F_g)]

  # Event time k. Defined for treated units only; -Inf for never-treated.
  panel[, k_evt := ifelse(is.finite(F_g), period - F_g, -Inf)]

  # Treatment effect: zero before k=0; tau_profile[k+1] for k>=0; for
  # k beyond profile, hold the last value (so figure k = E persists).
  max_k <- length(tau_profile) - 1L
  tau_val <- function(k) {
    if (!is.finite(k) || k < 0) return(0)
    idx <- min(as.integer(k) + 1L, length(tau_profile))
    tau_profile[idx]
  }
  panel[, tau_k := vapply(k_evt, tau_val, numeric(1))]

  # Build Y.
  panel[, Y := unit_fe[as.character(unit)] +
              time_fe[as.character(period)] +
              tau_k +
              stats::rnorm(.N, 0, sigma)]

  data.table::setkeyv(panel, c("unit", "period"))
  out <- as.data.frame(panel[, list(unit, period, D, Y)])

  attr(out, "truth") <- list(
    F_g = F_g,
    tau_profile = tau_profile,
    unit_fe = unit_fe,
    time_fe = time_fe,
    sigma = sigma,
    seed = seed
  )
  out
}


#' Generate a panel with BOTH switcher-in and switcher-out units
#'
#' Useful for testing the cross-direction Neyman pooling. The first
#' `frac_in` of treated units start at d=0 and turn on (switcher-in);
#' the remaining `1 - frac_in` start at d=1 and turn off (switcher-out).
#' Never-treated units stay at d=0.
#'
#' @param n_units integer. Number of units.
#' @param n_periods integer. Number of time periods.
#' @param frac_treated numeric between 0 and 1. Share that ever switch.
#' @param frac_in numeric between 0 and 1. Of switchers, share going
#'   in (rest go out). 0 = all out-switchers, 1 = all in-switchers.
#' @param min_treat_period,max_treat_period Treatment window.
#' @param tau_in,tau_out numeric vectors. Event-time profiles for
#'   in-switchers and out-switchers respectively.
#' @param sigma,unit_fe_sd,time_fe_sd,seed As in `didgpu_simulate_panel`.
#'
#' @return Data.frame with (unit, period, D, Y) plus a `truth` attribute.
#' @export
didgpu_simulate_panel_bidir <- function(
    n_units          = 100L,
    n_periods        = 20L,
    frac_treated     = 0.6,
    frac_in          = 0.5,
    min_treat_period = NULL,
    max_treat_period = NULL,
    tau_in           = c(0.5, 1.0, 1.2),
    tau_out          = c(-0.5, -0.8, -1.0),
    sigma            = 0.4,
    unit_fe_sd       = 1.0,
    time_fe_sd       = 0.3,
    seed             = 1L) {

  if (is.null(min_treat_period)) min_treat_period <- max(2L, as.integer(n_periods * 0.3))
  if (is.null(max_treat_period)) max_treat_period <- max(min_treat_period,
                                                          as.integer(n_periods * 0.7))

  set.seed(seed)
  units <- seq_len(n_units)
  periods <- seq_len(n_periods)

  n_treated <- round(n_units * frac_treated)
  n_in <- round(n_treated * frac_in)
  n_out <- n_treated - n_in
  treated_units <- if (n_treated > 0L) sort(sample(units, n_treated)) else integer(0)
  in_units  <- if (n_in  > 0L) treated_units[seq_len(n_in)] else integer(0)
  out_units <- if (n_out > 0L) treated_units[(n_in + 1L):n_treated] else integer(0)

  F_g <- rep(Inf, n_units); names(F_g) <- as.character(units)
  baseline_d <- rep(0L, n_units); names(baseline_d) <- as.character(units)
  direction <- rep(NA_integer_, n_units); names(direction) <- as.character(units)
  if (n_in > 0L) {
    F_g[as.character(in_units)] <- sample(min_treat_period:max_treat_period, n_in, replace = TRUE)
    direction[as.character(in_units)] <- 1L
  }
  if (n_out > 0L) {
    F_g[as.character(out_units)] <- sample(min_treat_period:max_treat_period, n_out, replace = TRUE)
    baseline_d[as.character(out_units)] <- 1L
    direction[as.character(out_units)] <- 0L
  }

  unit_fe <- stats::rnorm(n_units, 0, unit_fe_sd); names(unit_fe) <- as.character(units)
  time_fe <- stats::rnorm(n_periods, 0, time_fe_sd); names(time_fe) <- as.character(periods)

  panel <- data.table::CJ(unit = units, period = periods)
  panel[, baseline_d := baseline_d[as.character(unit)]]
  panel[, F_g := F_g[as.character(unit)]]
  panel[, direction := direction[as.character(unit)]]
  # D: starts at baseline_d, switches when t >= F_g. For in-switchers
  # (baseline 0): D = 1 when t >= F_g. For out-switchers (baseline 1):
  # D = 0 when t >= F_g.
  panel[, D := ifelse(period < F_g, baseline_d,
                       ifelse(direction == 1L, 1L, 0L))]
  panel[is.na(D), D := baseline_d]

  # Event-time effect: tau_in for in-switchers, tau_out for out.
  tau_val <- function(k, dir) {
    if (!is.finite(k) || k < 0) return(0)
    prof <- if (isTRUE(dir == 1L)) tau_in else if (isTRUE(dir == 0L)) tau_out else numeric(0)
    if (length(prof) == 0L) return(0)
    prof[min(as.integer(k) + 1L, length(prof))]
  }
  panel[, k_evt := ifelse(is.finite(F_g), period - F_g, -Inf)]
  panel[, tau_k := mapply(tau_val, k_evt, direction)]

  panel[, Y := unit_fe[as.character(unit)] +
              time_fe[as.character(period)] +
              tau_k +
              stats::rnorm(.N, 0, sigma)]

  data.table::setkeyv(panel, c("unit", "period"))
  out <- as.data.frame(panel[, list(unit, period, D, Y)])
  attr(out, "truth") <- list(
    F_g = F_g, baseline_d = baseline_d, direction = direction,
    tau_in = tau_in, tau_out = tau_out, sigma = sigma, seed = seed
  )
  out
}
