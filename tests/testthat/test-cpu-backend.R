# CPU backend (Rcpp): didgpu_cpp_core_one_event_time in src/cpu_core.cpp.
# Currently supports the binary-no-controls case. Falls back to r for
# unsupported feature combinations.

test_that("cpu backend matches r backend on basic binary panel", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              seed = 17L)
  us_r <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 1L,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  us_c <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 1L,
                  bootstrap_reps = 0L, backend = "cpu", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us_r$results$Effects[, "Estimate"]) -
                    as.numeric(us_c$results$Effects[, "Estimate"]))),
            1e-10)
  expect_lt(max(abs(as.numeric(us_r$results$Placebos[, "Estimate"]) -
                    as.numeric(us_c$results$Placebos[, "Estimate"]))),
            1e-10)
  # Sample-size columns must match exactly.
  for (col in c("N", "Switchers", "N.w", "Switchers.w")) {
    expect_equal(as.numeric(us_c$results$Effects[, col]),
                 as.numeric(us_r$results$Effects[, col]))
  }
})

test_that("cpu backend matches r backend with switchers = 'in'", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0),
                              seed = 17L)
  us_r <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 2L, switchers = "in",
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  us_c <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 2L, switchers = "in",
                  bootstrap_reps = 0L, backend = "cpu", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us_r$results$Effects[, "Estimate"]) -
                    as.numeric(us_c$results$Effects[, "Estimate"]))),
            1e-10)
})

test_that("cpu backend falls back to r for unsupported features", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0),
                              seed = 17L)
  set.seed(1L)
  p$w  <- runif(nrow(p), 0.5, 2)
  p$X1 <- rnorm(nrow(p))
  # Each of these should fall back to r-backend; the cpu backend's
  # public report should still say backend = "r" via the fallback.
  for (extra in list(
    list(weight = "w"),
    list(controls = "X1"),
    list(normalized = TRUE),
    list(trends_lin = TRUE),
    list(same_switchers = TRUE)
  )) {
    args <- c(list(p, "Y", "unit", "period", "D",
                    effects = 2L, bootstrap_reps = 0L,
                    backend = "cpu", verbose = FALSE),
              extra)
    fit <- do.call(didgpu, args)
    # The fall-through path goes through .backend_r_impl which tags the
    # cell with backend = "r" — so the public report shows the fall-back.
    # We check estimates equal r-backend directly.
    args_r <- c(list(p, "Y", "unit", "period", "D",
                      effects = 2L, bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE),
                extra)
    fit_r <- do.call(didgpu, args_r)
    expect_equal(as.numeric(fit$results$Effects[, "Estimate"]),
                 as.numeric(fit_r$results$Effects[, "Estimate"]),
                 tolerance = 1e-12,
                 info = sprintf("cpu fallback for %s",
                                 paste(names(extra), collapse = ", ")))
  }
})

test_that("cpu backend matches reference end-to-end on a clean panel", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel(n_units = 50L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0),
                              seed = 17L)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = as.double(2), placebo = 0, graph_off = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 2L, bootstrap_reps = 0L,
                backend = "cpu", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            1e-10)
})

test_that("didgpu_backend_info shows cpu as available", {
  info <- didgpu_backend_info()
  expect_true("cpu" %in% info$backend)
  expect_true(info$available[info$backend == "cpu"])
})

test_that("cpu backend works in bootstrap orchestrator", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0),
                              seed = 17L)
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 2L, bootstrap_reps = 5L, seed = 1L,
                 backend = "cpu", verbose = FALSE)
  expect_equal(nrow(fit$results$Effects), 2L)
  # Bootstrap SE should be finite and positive.
  se <- as.numeric(fit$results$Effects[, "SE"])
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
})
