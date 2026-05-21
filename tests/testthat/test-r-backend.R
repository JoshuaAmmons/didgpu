# The r-backend currently supports only effects = 1, placebo = 0,
# no controls. Within that subset, its output must match the reference
# to floating-point precision (the entire point of "r" is being a
# native, independent implementation that agrees with DIDmultiplegtDYN).

test_that("r-backend point estimate matches reference at effects=1, placebo=0", {
  skip_if_no_reference()

  for (seed in c(17L, 42L, 101L)) {
    p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                                frac_treated = 0.6,
                                min_treat_period = 4L, max_treat_period = 9L,
                                tau_profile = c(0.5, 1.0),
                                sigma = 0.4, seed = seed)

    ref <- suppressMessages(suppressWarnings(
      DIDmultiplegtDYN::did_multiplegt_dyn(
        df = as.data.frame(p), outcome = "Y", group = "unit",
        time = "period", treatment = "D",
        effects = 1, placebo = 0, graph_off = TRUE
      )
    ))

    us <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 1L, placebo = 0L,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)

    diff <- abs(us$results$Effects[1, "Estimate"] - ref$results$Effects[1, 1])
    expect_lt(diff, 1e-10,
              label = sprintf("seed=%d diff=%.2e", seed, diff))
  }
})

test_that("r-backend matches reference at effects in {2, 3, 5}", {
  skip_if_no_reference()
  for (eff in c(2L, 3L, 5L)) {
    for (seed in c(17L, 42L)) {
      p <- didgpu_simulate_panel(n_units = 80L, n_periods = 18L,
                                  frac_treated = 0.6,
                                  min_treat_period = 4L, max_treat_period = 9L,
                                  tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8, 0.6),
                                  sigma = 0.4, seed = seed)
      ref <- suppressMessages(suppressWarnings(
        DIDmultiplegtDYN::did_multiplegt_dyn(
          df = as.data.frame(p), outcome = "Y", group = "unit",
          time = "period", treatment = "D",
          effects = as.double(eff), placebo = 0, graph_off = TRUE
        )
      ))
      us <- didgpu(p, "Y", "unit", "period", "D",
                    effects = eff, placebo = 0L,
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
      ref_e <- as.numeric(ref$results$Effects[, 1])
      us_e  <- as.numeric(us$results$Effects[, "Estimate"])
      diff  <- max(abs(ref_e - us_e))
      expect_lt(diff, 1e-10,
                label = sprintf("effects=%d seed=%d maxdiff=%.2e",
                                eff, seed, diff))
    }
  }
})

test_that("Effects/Placebos sample-size columns match reference", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel_bidir(
    n_units = 80L, n_periods = 15L,
    frac_treated = 0.6, frac_in = 0.5,
    min_treat_period = 5L, max_treat_period = 10L, seed = 11L
  )
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 1, graph_off = TRUE
    )
  ))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 1L,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)

  expect_equal(ncol(us$results$Effects), 8L)
  expect_equal(colnames(us$results$Effects),
               c("Estimate", "SE", "LB.CI", "UB.CI",
                 "N", "Switchers", "N.w", "Switchers.w"))
  # Numerical match against the reference's N and Switchers columns.
  expect_equal(as.numeric(us$results$Effects[, "N"]),
               as.numeric(ref$results$Effects[, 5]))
  expect_equal(as.numeric(us$results$Effects[, "Switchers"]),
               as.numeric(ref$results$Effects[, 6]))
  expect_equal(as.numeric(us$results$Placebos[, "N"]),
               as.numeric(ref$results$Placebos[, 5]))
  expect_equal(as.numeric(us$results$Placebos[, "Switchers"]),
               as.numeric(ref$results$Placebos[, 6]))
})

test_that("r-backend ATE matches reference for effects >= 1", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel(
    n_units = 80L, n_periods = 18L, frac_treated = 0.6,
    min_treat_period = 4L, max_treat_period = 9L,
    tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8),
    sigma = 0.4, seed = 17L
  )
  for (eff in c(1L, 2L, 3L, 5L)) {
    ref <- suppressMessages(suppressWarnings(
      DIDmultiplegtDYN::did_multiplegt_dyn(
        df = as.data.frame(p), outcome = "Y", group = "unit",
        time = "period", treatment = "D",
        effects = as.double(eff), placebo = 0, graph_off = TRUE
      )
    ))
    us <- didgpu(p, "Y", "unit", "period", "D",
                  effects = eff, placebo = 0L,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
    diff <- abs(as.numeric(us$results$ATE[1, "Estimate"]) -
                as.numeric(ref$results$ATE[1, 1]))
    # Allow up to 4 ulps to absorb the order-of-summation difference
    # between the per-group then per-k aggregation paths.
    expect_lt(diff, 4 * .Machine$double.eps,
              label = sprintf("effects=%d diff=%.2e", eff, diff))
  }
})

test_that("r-backend through full bootstrap orchestrator matches reference", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel(
    n_units = 60L, n_periods = 14L, frac_treated = 0.6,
    min_treat_period = 4L, max_treat_period = 9L,
    tau_profile = c(0.5, 1.0, 1.2),
    sigma = 0.4, seed = 19L
  )
  cdir_r   <- tempfile("didgpu_r_e2e_")
  cdir_ref <- tempfile("didgpu_ref_e2e_")
  on.exit(unlink(c(cdir_r, cdir_ref), recursive = TRUE), add = TRUE)

  fit_r <- didgpu(p, "Y", "unit", "period", "D",
                   effects = 3L, placebo = 1L,
                   bootstrap_reps = 10L, seed = 1L,
                   checkpoint_dir = cdir_r,
                   backend = "r", verbose = FALSE)
  fit_ref <- didgpu(p, "Y", "unit", "period", "D",
                     effects = 3L, placebo = 1L,
                     bootstrap_reps = 10L, seed = 1L,
                     checkpoint_dir = cdir_ref,
                     backend = "reference", verbose = FALSE)

  # Point estimates: bit-identical (both backends call the same kernel
  # on the same cluster-resampled panels with the same seed).
  expect_lt(max(abs(fit_r$results$Effects[, "Estimate"] -
                    fit_ref$results$Effects[, "Estimate"])), 1e-10)
  expect_lt(max(abs(fit_r$results$Placebos[, "Estimate"] -
                    fit_ref$results$Placebos[, "Estimate"])), 1e-10)
  # SEs: should match within FP tolerance.
  expect_lt(max(abs(fit_r$results$Effects[, "SE"] -
                    fit_ref$results$Effects[, "SE"])), 1e-10)

  # Resume on r-backend dir produces identical output.
  fit_r2 <- didgpu(p, "Y", "unit", "period", "D",
                    effects = 3L, placebo = 1L,
                    bootstrap_reps = 10L, seed = 1L,
                    checkpoint_dir = cdir_r,
                    backend = "r", verbose = FALSE)
  expect_equal(fit_r$results$Effects[, "Estimate"],
               fit_r2$results$Effects[, "Estimate"])
  expect_equal(fit_r$results$Effects[, "SE"],
               fit_r2$results$Effects[, "SE"])
})

test_that("r-backend matches reference with placebos", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel(
    n_units = 80L, n_periods = 18L, frac_treated = 0.6,
    min_treat_period = 7L, max_treat_period = 12L,
    tau_profile = c(0.5, 1.0, 1.2),
    sigma = 0.4, seed = 17L
  )
  cfgs <- list(c(2L, 1L), c(3L, 2L), c(1L, 3L))
  for (cfg in cfgs) {
    eff <- cfg[1]; pl <- cfg[2]
    ref <- suppressMessages(suppressWarnings(
      DIDmultiplegtDYN::did_multiplegt_dyn(
        df = as.data.frame(p), outcome = "Y", group = "unit",
        time = "period", treatment = "D",
        effects = as.double(eff), placebo = as.double(pl), graph_off = TRUE
      )
    ))
    us <- didgpu(p, "Y", "unit", "period", "D",
                  effects = eff, placebo = pl,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
    expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                      as.numeric(ref$results$Effects[, 1]))), 1e-10)
    expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                      as.numeric(ref$results$Placebos[, 1]))), 1e-10)
    # Auto-clamp: didgpu should report the same n_placebos as the reference.
    expect_equal(us$results$N_Placebos, nrow(ref$results$Placebos),
                 label = sprintf("effects=%d placebo=%d", eff, pl))
  }
})
