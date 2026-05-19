# same_switchers_pl: additionally restrict placebo dist masks on a
# still_switcher_pl_XX indicator that requires the switcher to have
# valid pre-period diff_y at every placebo horizon q in 1..placebo.
# Reference: did_multiplegt_dyn_core.R:177-215.
#
# The reference requires that same_switchers = TRUE is also set when
# same_switchers_pl = TRUE; didgpu mirrors that constraint.

build_pl_panel <- function(seed = 17L) {
  didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                         frac_treated = 0.6,
                         min_treat_period = 6L, max_treat_period = 10L,
                         tau_profile = c(0.5, 1.0, 1.2),
                         sigma = 0.4, seed = seed)
}

test_that("same_switchers + same_switchers_pl matches reference (placebo = 2)", {
  skip_if_no_reference()
  p <- build_pl_panel()
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = as.double(2), placebo = as.double(2), graph_off = TRUE,
      same_switchers = TRUE, same_switchers_pl = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 2L, placebo = 2L,
                same_switchers = TRUE, same_switchers_pl = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
  expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                    as.numeric(ref$results$Placebos[, 1]))),
            1e-10)
})

test_that("same_switchers + same_switchers_pl matches at placebo = 1, 3", {
  skip_if_no_reference()
  p <- build_pl_panel()
  for (npl in c(1L, 3L)) {
    ref <- suppressMessages(suppressWarnings(
      DIDmultiplegtDYN::did_multiplegt_dyn(
        df = as.data.frame(p), outcome = "Y", group = "unit",
        time = "period", treatment = "D",
        effects = as.double(2), placebo = as.double(npl), graph_off = TRUE,
        same_switchers = TRUE, same_switchers_pl = TRUE)))
    us <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 2L, placebo = npl,
                  same_switchers = TRUE, same_switchers_pl = TRUE,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
    expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                      as.numeric(ref$results$Placebos[, 1]))),
              1e-10,
              label = sprintf("placebo = %d", npl))
  }
})

test_that("same_switchers_pl without same_switchers errors with clear message", {
  p <- build_pl_panel()
  expect_error(
    didgpu(p, "Y", "unit", "period", "D",
            effects = 2L, placebo = 1L,
            same_switchers_pl = TRUE,
            bootstrap_reps = 0L, backend = "r", verbose = FALSE),
    "same_switchers_pl.*requires.*same_switchers"
  )
})

test_that("same_switchers_pl with placebo = 0 errors", {
  p <- build_pl_panel()
  expect_error(
    didgpu(p, "Y", "unit", "period", "D",
            effects = 2L, placebo = 0L,
            same_switchers = TRUE, same_switchers_pl = TRUE,
            bootstrap_reps = 0L, backend = "r", verbose = FALSE),
    "requires .placebo > 0"
  )
})

test_that("same_switchers_pl = FALSE (default) is unchanged", {
  p <- build_pl_panel()
  us_default <- didgpu(p, "Y", "unit", "period", "D",
                       effects = 2L, placebo = 2L,
                       same_switchers = TRUE,
                       bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  us_explicit <- didgpu(p, "Y", "unit", "period", "D",
                        effects = 2L, placebo = 2L,
                        same_switchers = TRUE, same_switchers_pl = FALSE,
                        bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_identical(as.numeric(us_default$results$Placebos[, "Estimate"]),
                   as.numeric(us_explicit$results$Placebos[, "Estimate"]))
})
