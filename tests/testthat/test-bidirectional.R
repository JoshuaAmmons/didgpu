# Bidirectional panels: both switcher-in and switcher-out units present.
# This is where the cross-direction sign convention and the
# (time, d_sq) cohort matching matter.

test_that("r-backend matches reference on bidirectional panel (effects)", {
  skip_if_no_reference()

  for (seed in c(11L, 17L, 42L)) {
    p <- didgpu_simulate_panel_bidir(
      n_units = 80L, n_periods = 15L,
      frac_treated = 0.6, frac_in = 0.5,
      min_treat_period = 5L, max_treat_period = 10L,
      seed = seed
    )

    for (eff in c(1L, 2L, 3L)) {
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
      diff <- max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                      as.numeric(ref$results$Effects[, 1])))
      expect_lt(diff, 1e-10,
                label = sprintf("bidir seed=%d effects=%d diff=%.2e",
                                seed, eff, diff))
    }
  }
})

test_that("r-backend matches reference on bidirectional panel with placebos", {
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
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))), 1e-10)
  expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                    as.numeric(ref$results$Placebos[, 1]))), 1e-10)
  expect_lt(abs(as.numeric(us$results$ATE[1, "Estimate"]) -
                as.numeric(ref$results$ATE[1, 1])), 4 * .Machine$double.eps)
})

test_that("r-backend matches reference for an out-only panel with always-treated controls", {
  skip_if_no_reference()
  set.seed(11L)
  n_units <- 80L; n_periods <- 15L; n_treated <- 48L
  treated_units <- sort(sample(seq_len(n_units), n_treated))
  F_g <- rep(Inf, n_units)
  F_g[treated_units] <- sample(5L:10L, n_treated, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  tau_out <- c(-0.5, -0.8, -1.0, -1.1, -1.2)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- ifelse(panel$period >= panel$F_g, 0L, 1L)
  panel$k_evt <- panel$period - panel$F_g
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  panel$tau_k <- 0
  panel$tau_k[post] <- tau_out[pmin(panel$k_evt[post] + 1L, length(tau_out))]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             panel$tau_k + rnorm(nrow(panel), 0, 0.4)
  panel <- panel[order(panel$unit, panel$period), c("unit", "period", "D", "Y")]

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = panel, outcome = "Y", group = "unit", time = "period",
      treatment = "D", effects = 3, placebo = 0, graph_off = TRUE
    )
  ))
  us <- didgpu(panel, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))), 1e-10)
})

test_that("r-backend matches reference for out-only with no always-treated (only switchers as controls)", {
  skip_if_no_reference()
  # 24 out-switchers, no never-treated. Each switcher uses other unswitched
  # out-switchers as its controls. This is the edge case where some
  # switcher cells have no concurrent controls (and must be excluded).
  set.seed(11L)
  n_units <- 24L; n_periods <- 15L
  F_g <- sample(5L:10L, n_units, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  tau_out <- c(-0.5, -0.8, -1.0, -1.1, -1.2)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- ifelse(panel$period >= panel$F_g, 0L, 1L)
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  post <- panel$k_evt >= 0
  panel$tau_k[post] <- tau_out[pmin(panel$k_evt[post] + 1L, length(tau_out))]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             panel$tau_k + rnorm(nrow(panel), 0, 0.4)
  panel <- panel[order(panel$unit, panel$period), c("unit", "period", "D", "Y")]

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = panel, outcome = "Y", group = "unit", time = "period",
      treatment = "D", effects = 3, placebo = 0, graph_off = TRUE
    )
  ))
  us <- didgpu(panel, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))), 1e-10)
})

test_that("cluster argument: point estimate invariant, bootstrap SEs change", {
  set.seed(11L)
  n_units <- 80L; n_periods <- 15L
  clust <- sample.int(20L, n_units, replace = TRUE)
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), 48L))
  F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$cluster <- clust[panel$unit]
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g & is.finite(panel$F_g))
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  panel$tau_k[post] <- c(0.5, 1.0, 1.2, 1.0)[pmin(panel$k_evt[post] + 1L, 4L)]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] + panel$tau_k +
             rnorm(nrow(panel), 0, 0.4)
  panel <- panel[order(panel$unit, panel$period), ]

  # Point estimates must match regardless of cluster arg.
  fit_no <- didgpu(panel, "Y", "unit", "period", "D",
                    effects = 3L, placebo = 0L, bootstrap_reps = 0L,
                    backend = "r", verbose = FALSE)
  fit_c <- didgpu(panel, "Y", "unit", "period", "D",
                   effects = 3L, placebo = 0L, cluster = "cluster",
                   bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_equal(fit_no$results$Effects[, "Estimate"],
               fit_c$results$Effects[, "Estimate"])

  # Bootstrap SEs must differ when clustering at different levels.
  fit_bg <- didgpu(panel, "Y", "unit", "period", "D",
                    effects = 3L, placebo = 0L,
                    bootstrap_reps = 20L, seed = 1L,
                    backend = "r", verbose = FALSE)
  fit_bc <- didgpu(panel, "Y", "unit", "period", "D",
                    effects = 3L, placebo = 0L, cluster = "cluster",
                    bootstrap_reps = 20L, seed = 1L,
                    backend = "r", verbose = FALSE)
  expect_true(any(abs(fit_bg$results$Effects[, "SE"] -
                      fit_bc$results$Effects[, "SE"]) > 1e-4))
})

test_that("switchers argument matches reference at 'in' and 'out'", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel_bidir(
    n_units = 80L, n_periods = 15L,
    frac_treated = 0.6, frac_in = 0.5,
    min_treat_period = 5L, max_treat_period = 10L, seed = 11L
  )
  for (sw in c("in", "out")) {
    ref <- suppressMessages(suppressWarnings(
      DIDmultiplegtDYN::did_multiplegt_dyn(
        df = as.data.frame(p), outcome = "Y", group = "unit",
        time = "period", treatment = "D",
        effects = 3, placebo = 0, graph_off = TRUE,
        switchers = sw
      )
    ))
    us <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 0L, switchers = sw,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
    diff <- max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1])))
    expect_lt(diff, 1e-10, label = sprintf("switchers=%s diff=%.2e", sw, diff))
  }
})

test_that("r-backend through full bootstrap matches reference on bidirectional", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel_bidir(
    n_units = 60L, n_periods = 13L,
    frac_treated = 0.6, frac_in = 0.5,
    min_treat_period = 5L, max_treat_period = 9L, seed = 23L
  )
  cdir_r <- tempfile("didgpu_bidir_e2e_r_")
  cdir_ref <- tempfile("didgpu_bidir_e2e_ref_")
  on.exit(unlink(c(cdir_r, cdir_ref), recursive = TRUE), add = TRUE)

  fit_r <- didgpu(p, "Y", "unit", "period", "D",
                   effects = 3L, placebo = 1L,
                   bootstrap_reps = 8L, seed = 1L,
                   checkpoint_dir = cdir_r,
                   backend = "r", verbose = FALSE)
  fit_ref <- didgpu(p, "Y", "unit", "period", "D",
                     effects = 3L, placebo = 1L,
                     bootstrap_reps = 8L, seed = 1L,
                     checkpoint_dir = cdir_ref,
                     backend = "reference", verbose = FALSE)
  expect_lt(max(abs(fit_r$results$Effects[, "Estimate"] -
                    fit_ref$results$Effects[, "Estimate"])), 1e-10)
  expect_lt(max(abs(fit_r$results$Effects[, "SE"] -
                    fit_ref$results$Effects[, "SE"])), 1e-10)
})
