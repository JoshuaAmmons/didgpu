# Helper: a deterministic small panel used by multiple test files.
small_panel <- function(seed = 17L) {
  didgpu_simulate_panel(
    n_units = 60L, n_periods = 12L, frac_treated = 0.6,
    min_treat_period = 4L, max_treat_period = 9L,
    tau_profile = c(0.5, 1.0, 1.2, 1.0),
    sigma = 0.4, seed = seed
  )
}

# Is the DIDmultiplegtDYN reference oracle actually RUNNABLE here?
#
# DIDmultiplegtDYN (>= 2.x) computes via bare `pl$...` calls and relies on
# the 'polars' package being ATTACHED — but polars is only in its Suggests,
# so the reference package never attaches it itself and errors with
# "object 'pl' not found" when a caller hasn't run library(polars). polars
# is also not on CRAN (it ships from r-universe), so it may be absent or an
# incompatible version entirely.
#
# So "installed" is not enough. We (a) attach polars when it's available
# so the reference can find `pl`, then (b) smoke-test the reference once.
# If it still can't run, the parity tests SKIP instead of erroring the
# whole suite. Result cached so we pay the smoke test at most once.
.reference_runnable <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    cached <<- isTRUE(tryCatch({
      if (!requireNamespace("DIDmultiplegtDYN", quietly = TRUE)) return(FALSE)
      # DIDmultiplegtDYN needs polars attached (bare `pl$...`); attach it.
      if (requireNamespace("polars", quietly = TRUE) &&
          !"package:polars" %in% search()) {
        suppressMessages(suppressWarnings(
          library("polars", character.only = TRUE,
                  quietly = TRUE, warn.conflicts = FALSE)))
      }
      df <- didgpu_simulate_panel(n_units = 20L, n_periods = 6L, seed = 1L)
      df$D <- as.integer(df$D >= 0.5)
      suppressMessages(suppressWarnings(
        DIDmultiplegtDYN::did_multiplegt_dyn(
          df = as.data.frame(df), outcome = "Y", group = "unit",
          time = "period", treatment = "D", effects = 1, graph_off = TRUE)))
      TRUE
    }, error = function(e) FALSE))
    cached
  }
})

# Skip a reference-parity test when the DIDmultiplegtDYN oracle is missing
# OR installed-but-not-runnable (e.g. its 'polars' backend is absent or an
# incompatible version). Keeps the suite green on machines/CI where the
# reference cannot run, while still executing the comparison where it can.
skip_if_no_reference <- function() {
  testthat::skip_if_not_installed("DIDmultiplegtDYN")
  if (!.reference_runnable()) {
    testthat::skip(
      "DIDmultiplegtDYN installed but not runnable here (its optional 'polars' backend is missing or incompatible); skipping reference-parity check")
  }
}
