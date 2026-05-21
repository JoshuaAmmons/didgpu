# normalized = TRUE: divide each per-event-time DID by its pooled
# cumulative treatment-change magnitude delta_D_k. Reference:
# did_multiplegt_main.R:1086-1090 (delta_D pooling) + 1124-1126 (DID /=).
# For binary on/off, delta_D_k = k, so this divides the k-th cumulative
# effect by k, producing the average per-period effect.

test_that("binary, normalized = TRUE matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                              frac_treated = 0.6,
                              min_treat_period = 5L, max_treat_period = 10L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 17L)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE, normalized = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("binary normalized = TRUE divides cumulative DID by k", {
  # For binary treatment, delta_D_k = k (because each in-switcher gains
  # exactly 1 unit of D for k periods, normalized by N_inc = N_inc).
  # So the k-th normalized effect == (unnormalized k-th effect) / k.
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                              frac_treated = 0.6,
                              min_treat_period = 5L, max_treat_period = 10L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 17L)
  us_n <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 0L, normalized = FALSE,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  us_y <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 0L, normalized = TRUE,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  cumulative <- as.numeric(us_n$results$Effects[, "Estimate"])
  normalized <- as.numeric(us_y$results$Effects[, "Estimate"])
  expect_lt(max(abs(normalized - cumulative / seq_along(cumulative))),
            1e-10)
})

test_that("binary, switchers = in + normalized = TRUE matches reference", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                              frac_treated = 0.6,
                              min_treat_period = 5L, max_treat_period = 10L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 17L)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE,
      switchers = "in", normalized = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L,
                switchers = "in", normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("binary placebos with normalized = TRUE match reference", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                              frac_treated = 0.6,
                              min_treat_period = 5L, max_treat_period = 10L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 17L)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 2, graph_off = TRUE, normalized = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 2L, normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                    as.numeric(ref$results$Placebos[, 1]))),
            1e-10)
})

test_that("continuous = 1 + normalized = TRUE matches reference", {
  skip_if_no_reference()
  # Use the same DGP as test-continuous.R for cross-comparison.
  set.seed(11L)
  n_units <- 80L; n_periods <- 15L
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), as.integer(n_units * 0.6)))
  F_g[treated] <- sample(5L:10L, length(treated), replace = TRUE)
  baseline_d <- runif(n_units, 0, 1)
  post_d <- baseline_d + rnorm(n_units, 0, 0.5)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- ifelse(panel$period >= panel$F_g,
                    post_d[panel$unit], baseline_d[panel$unit])
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  prof <- c(0.3, 0.5, 0.6, 0.5)
  panel$tau_k[post] <- (post_d[panel$unit[post]] - baseline_d[panel$unit[post]]) *
                       prof[pmin(panel$k_evt[post] + 1L, length(prof))]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             panel$tau_k + rnorm(nrow(panel), 0, 0.4)
  p <- panel[order(panel$unit, panel$period), c("unit", "period", "D", "Y")]

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE,
      continuous = 1, normalized = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L,
                continuous = 1L, normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("multivalued (D in {0,1,2,3}) + normalized = TRUE matches reference", {
  skip_if_no_reference()
  set.seed(23L)
  n_units <- 80L; n_periods <- 15L
  fg <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), 48L))
  fg[treated] <- sample(5L:10L, 48L, replace = TRUE)
  doses <- sample(c(1L, 2L, 3L), n_units, replace = TRUE)
  mv <- data.table::data.table(
    unit = rep(seq_len(n_units), each = n_periods),
    period = rep(seq_len(n_periods), n_units))
  mv[, F_g := fg[unit]]
  mv[, dose := doses[unit]]
  mv[, D := as.integer(ifelse(period >= F_g, dose, 0L))]
  mv[, k_evt := period - F_g]
  mv[, tau_k := 0]
  mv[is.finite(F_g) & k_evt >= 0,
     tau_k := dose * c(0.3, 0.5, 0.6, 0.5)[pmin(k_evt + 1L, 4L)]]
  mv_unit_fe <- rnorm(n_units, 0, 1)
  mv_time_fe <- rnorm(n_periods, 0, 0.3)
  mv[, Y := mv_unit_fe[unit] + mv_time_fe[period] + tau_k + rnorm(.N, 0, 0.4)]
  p <- as.data.frame(mv[order(unit, period), .(unit, period, D, Y)])

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE, normalized = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("ATE under normalized = TRUE matches reference (unnormalized)", {
  # The reference's Av_tot_eff is built from per-row U_Gg contributions
  # summed across event-times and is NEVER divided by delta_D, regardless
  # of the normalized= flag. So our ATE under normalized=TRUE must match
  # the ATE under normalized=FALSE (and both must match the reference).
  skip_if_no_reference()
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              frac_treated = 0.6,
                              min_treat_period = 4L, max_treat_period = 9L,
                              tau_profile = c(0.5, 1.0, 1.2, 1.0),
                              sigma = 0.4, seed = 17L)
  us_y <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 0L, normalized = TRUE,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  us_n <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 0L, normalized = FALSE,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  ref_y <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE, normalized = TRUE)))
  # ATE under normalized=TRUE matches the unnormalized ATE.
  expect_lt(abs(as.numeric(us_y$results$ATE[1, "Estimate"]) -
                as.numeric(us_n$results$ATE[1, "Estimate"])),
            8 * .Machine$double.eps)
  # And both match the reference's ATE.
  expect_lt(abs(as.numeric(us_y$results$ATE[1, "Estimate"]) -
                as.numeric(ref_y$results$ATE[1, "Estimate"])),
            8 * .Machine$double.eps)
})

test_that("normalized = TRUE composes correctly with weight", {
  skip_if_no_reference()
  set.seed(31L)
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
                              frac_treated = 0.6,
                              min_treat_period = 4L, max_treat_period = 8L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 31L)
  p$w <- runif(nrow(p), 0.5, 2.5)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE,
      weight = "w", normalized = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, weight = "w",
                normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            1e-10)
})

test_that("normalized = TRUE composes correctly with controls", {
  skip_if_no_reference()
  set.seed(33L)
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
                              frac_treated = 0.6,
                              min_treat_period = 4L, max_treat_period = 8L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 33L)
  p$X <- rnorm(nrow(p))
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE,
      controls = "X", normalized = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, controls = "X",
                normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            1e-8)
})

test_that("normalized = TRUE composes with only_never_switchers + same_switchers", {
  skip_if_no_reference()
  set.seed(37L)
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
                              frac_treated = 0.6,
                              min_treat_period = 4L, max_treat_period = 8L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 37L)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 2, placebo = 0, graph_off = TRUE,
      only_never_switchers = TRUE, same_switchers = TRUE,
      normalized = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 2L, placebo = 0L,
                only_never_switchers = TRUE, same_switchers = TRUE,
                normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("normalized = FALSE (default) gives unnormalized cumulative effects", {
  # Regression guard: passing normalized=FALSE (or omitting it) must
  # NOT change existing behaviour. Mirrors test-r-backend.R.
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                              frac_treated = 0.6,
                              min_treat_period = 5L, max_treat_period = 10L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 17L)
  us_default <- didgpu(p, "Y", "unit", "period", "D",
                       effects = 3L, placebo = 0L,
                       bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  us_explicit <- didgpu(p, "Y", "unit", "period", "D",
                        effects = 3L, placebo = 0L, normalized = FALSE,
                        bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_identical(as.numeric(us_default$results$Effects[, "Estimate"]),
                   as.numeric(us_explicit$results$Effects[, "Estimate"]))
})
