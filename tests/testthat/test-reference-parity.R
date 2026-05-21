# Numerical equivalence: didgpu (reference backend) must produce point
# estimates identical to DIDmultiplegtDYN. SEs differ by construction
# (bootstrap vs. analytic) and are not compared here -- see
# test-bootstrap-converges.R for the SE convergence test.

test_that("reference backend point estimate matches DIDmultiplegtDYN exactly", {
  skip_if_no_reference()
  p <- small_panel()

  fit_ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 4, placebo = 2, graph_off = TRUE
    )
  ))

  fit_us <- didgpu(
    df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
    effects = 4L, placebo = 2L,
    bootstrap_reps = 0L, backend = "reference",
    verbose = FALSE
  )

  diff_e <- max(abs(fit_us$results$Effects[, "Estimate"] -
                    fit_ref$results$Effects[, 1]))
  diff_p <- max(abs(fit_us$results$Placebos[, "Estimate"] -
                    fit_ref$results$Placebos[, 1]))

  expect_lt(diff_e, 1e-12)
  expect_lt(diff_p, 1e-12)
})

test_that("resume yields identical aggregate to a fresh run with same config", {
  skip_if_no_reference()
  p <- small_panel()
  cdir <- tempfile("didgpu_test_")
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)

  # First run: 5 bootstrap reps to completion.
  fit1 <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 1L,
                  bootstrap_reps = 5L, seed = 42L,
                  checkpoint_dir = cdir, backend = "reference",
                  verbose = FALSE)

  # Second invocation on the same dir: should add no cells, return same.
  fit2 <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 1L,
                  bootstrap_reps = 5L, seed = 42L,
                  checkpoint_dir = cdir, backend = "reference",
                  verbose = FALSE)

  expect_equal(fit1$results$Effects[, "Estimate"],
               fit2$results$Effects[, "Estimate"])
  expect_equal(fit1$results$Effects[, "SE"],
               fit2$results$Effects[, "SE"])
  expect_equal(fit1$results$p_jointeffects,
               fit2$results$p_jointeffects)
})

test_that("config mismatch on resume is rejected", {
  skip_if_no_reference()
  p <- small_panel()
  cdir <- tempfile("didgpu_test_")
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)

  didgpu(p, "Y", "unit", "period", "D",
         effects = 3L, placebo = 1L, bootstrap_reps = 2L, seed = 1L,
         checkpoint_dir = cdir, backend = "reference", verbose = FALSE)

  # Changing effects mid-resume should error, not silently overwrite.
  expect_error(
    didgpu(p, "Y", "unit", "period", "D",
           effects = 4L, placebo = 1L, bootstrap_reps = 2L, seed = 1L,
           checkpoint_dir = cdir, backend = "reference", verbose = FALSE),
    "config mismatch"
  )
})
