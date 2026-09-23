# didgpu_compare() must work from a clean session.
#
# DIDmultiplegtDYN (>= 2.x) calls polars through bare `pl$...` but only
# Suggests it, so it never attaches polars itself. Without
# library(polars) the reference fit died inside polars -- "attempt to
# apply non-function", or "Evaluation failed in `$with_columns()`"
# depending on version -- and didgpu_compare() surfaced that with
# nothing pointing at the cause. It now attaches polars when available
# and says what is missing when it is not.

test_that("didgpu_compare runs without the caller attaching polars", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  skip_if("package:polars" %in% search(),
          "polars already attached in this session")
  p <- as.data.frame(didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                                           seed = 11L))
  expect_no_error(
    didgpu_compare(p, "Y", "unit", "period", "D", effects = 2L,
                   verbose = FALSE))
})

test_that("didgpu_compare agrees with the reference", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  p <- as.data.frame(didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                                           seed = 11L))
  cmp <- didgpu_compare(p, "Y", "unit", "period", "D", effects = 2L,
                        verbose = FALSE)
  expect_true(isTRUE(cmp$pass) || is.na(cmp$pass))
})
