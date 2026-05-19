# Helper: a deterministic small panel used by multiple test files.
small_panel <- function(seed = 17L) {
  didgpu_simulate_panel(
    n_units = 60L, n_periods = 12L, frac_treated = 0.6,
    min_treat_period = 4L, max_treat_period = 9L,
    tau_profile = c(0.5, 1.0, 1.2, 1.0),
    sigma = 0.4, seed = seed
  )
}

# Skip a test when DIDmultiplegtDYN is not installed.
skip_if_no_reference <- function() {
  testthat::skip_if_not_installed("DIDmultiplegtDYN")
}
