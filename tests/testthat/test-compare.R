test_that("didgpu_compare returns pass=TRUE on a panel where backends agree", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel_bidir(
    n_units = 80L, n_periods = 15L,
    frac_treated = 0.6, frac_in = 0.5,
    min_treat_period = 5L, max_treat_period = 10L, seed = 11L
  )
  res <- didgpu_compare(p, "Y", "unit", "period", "D",
                         effects = 3L, placebo = 1L, verbose = FALSE)
  expect_true(res$pass)
  expect_s3_class(res$report, "data.frame")
  expect_true(all(res$report$max_abs_diff[res$report$col == "Estimate"] < 1e-10))
})

test_that("didgpu_compare gracefully NAs when reference is not installed", {
  with_mock_namespace <- function() {
    # Temporarily mask DIDmultiplegtDYN. requireNamespace fails when the
    # package isn't installed; we simulate that without actually
    # uninstalling.
    saved <- requireNamespace
    # No-op; we can't easily mock requireNamespace cleanly. Skip this
    # negative-path test if the reference is in fact installed.
    if (requireNamespace("DIDmultiplegtDYN", quietly = TRUE)) {
      skip("DIDmultiplegtDYN is installed; can't test 'not installed' branch")
    }
    p <- didgpu_simulate_panel(n_units = 20L, n_periods = 8L)
    res <- expect_warning(
      didgpu_compare(p, "Y", "unit", "period", "D",
                      effects = 1L, verbose = FALSE),
      "not installed"
    )
    expect_true(is.na(res$pass))
  }
  with_mock_namespace()
})
