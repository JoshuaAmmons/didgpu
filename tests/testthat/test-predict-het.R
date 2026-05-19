# predict_het: regress per-group ATT contribution on time-invariant
# covariates with HC1 robust SEs and a joint F-test. Reference:
# did_multiplegt_main.R:1641-1745.

build_het_panel <- function(seed = 17L, n_units = 100L) {
  p <- didgpu_simulate_panel(n_units = n_units, n_periods = 15L,
                              frac_treated = 0.6,
                              min_treat_period = 5L, max_treat_period = 10L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = seed)
  set.seed(seed + 100L)
  unit_size <- runif(length(unique(p$unit)), 0, 5)
  p$size <- unit_size[p$unit]
  unit_age <- rnorm(length(unique(p$unit)), 50, 10)
  p$age <- unit_age[p$unit]
  p
}

test_that("predict_het single covariate matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- build_het_panel()
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = as.double(3), placebo = 0, graph_off = TRUE,
      predict_het = list(c("size"), c(-1)))))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L,
                predict_het = list("size", -1L),
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_true(!is.null(us$results$predict_het))
  expect_equal(nrow(us$results$predict_het), 3L)
  expect_equal(us$results$predict_het$effect, c(1L, 2L, 3L))
  expect_equal(us$results$predict_het$covariate, c("size", "size", "size"))
  for (col in c("Estimate", "SE", "t", "LB", "UB", "pF")) {
    expect_lt(max(abs(us$results$predict_het[[col]] -
                      ref$results$predict_het[[col]])),
              1e-10,
              label = sprintf("predict_het %s column", col))
  }
  expect_equal(us$results$predict_het$N, ref$results$predict_het$N)
})

test_that("predict_het multiple covariates matches reference", {
  skip_if_no_reference()
  p <- build_het_panel()
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = as.double(3), placebo = 0, graph_off = TRUE,
      predict_het = list(c("size", "age"), c(-1)))))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L,
                predict_het = list(c("size", "age"), -1L),
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_equal(nrow(us$results$predict_het), 6L)
  # Reference sorts by (covariate, effect); we sort by (effect, covariate).
  # Align both by sorting on the same composite key before comparing.
  key <- function(d) order(d$effect, d$covariate)
  ref_s <- ref$results$predict_het[key(ref$results$predict_het), ]
  us_s  <- us$results$predict_het[key(us$results$predict_het), ]
  for (col in c("Estimate", "SE", "t", "LB", "UB")) {
    expect_lt(max(abs(us_s[[col]] - ref_s[[col]])),
              1e-10,
              label = sprintf("multi-cov predict_het %s column", col))
  }
})

test_that("predict_het subset event-times matches reference", {
  skip_if_no_reference()
  p <- build_het_panel()
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = as.double(3), placebo = 0, graph_off = TRUE,
      predict_het = list(c("size"), c(1, 3)))))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L,
                predict_het = list("size", c(1L, 3L)),
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_equal(us$results$predict_het$effect, c(1L, 3L))
  for (col in c("Estimate", "SE")) {
    expect_lt(max(abs(us$results$predict_het[[col]] -
                      ref$results$predict_het[[col]])),
              1e-10,
              label = sprintf("subset predict_het %s column", col))
  }
})

test_that("predict_het with bad covariates errors clearly", {
  p <- build_het_panel()
  expect_error(
    didgpu(p, "Y", "unit", "period", "D",
            effects = 2L, predict_het = list("nope", -1L),
            bootstrap_reps = 0L, backend = "r", verbose = FALSE),
    "predict_het covariates not in df"
  )
})

test_that("predict_het with normalized = TRUE is ignored (returns NULL)", {
  p <- build_het_panel()
  # Capture the validation message but don't require its exact text.
  fit <- suppressMessages(
    didgpu(p, "Y", "unit", "period", "D",
            effects = 2L, normalized = TRUE,
            predict_het = list("size", -1L),
            bootstrap_reps = 0L, backend = "r", verbose = FALSE))
  expect_null(fit$results$predict_het)
})

test_that("predict_het wrong shape errors", {
  p <- build_het_panel()
  expect_error(
    didgpu(p, "Y", "unit", "period", "D",
            effects = 2L, predict_het = "size",
            bootstrap_reps = 0L, backend = "r", verbose = FALSE),
    "must be a list of length 2"
  )
})
