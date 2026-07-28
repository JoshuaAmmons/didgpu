# Regression test: bootstrap aggregation must survive degenerate resamples.
#
# Under sparse-switching treatments, a cluster-bootstrap resample can contain
# no valid switcher cell at some horizon, yielding a zero-length effects
# vector for that iteration. Before the 0.1.2 fix, .aggregate_to_result()'s
# vapply() calls hard-required full-length vectors, so a single such
# iteration aborted the whole estimation with
#   "values must be length K ... but FUN(X[[i]]) result is length 0"
# and the failure probability GREW with bootstrap_reps: exactly the
# large-rep final runs users need for publication inference were the ones
# crashing (observed in the wild: a 78-group / 4-switcher panel had
# 103/2000 degenerate resamples).
#
# The fix drops degenerate iterations with a warning; SEs and the bootstrap
# covariance use the surviving iterations, reported via results$n_boot and
# results$n_boot_dropped.

# Panel built so degenerate resamples are near-certain at 400 reps:
# only TWO switcher units among 30, both adopting at t = 10. A unit
# resample (30 draws with replacement) excludes any given unit with
# probability (29/30)^30 ~ 0.36, hence excludes BOTH switchers with
# probability ~ 0.13 => ~50 degenerate iterations expected at 400 reps.
make_sparse_panel <- function(seed = 42L) {
  set.seed(seed)
  nU <- 30L; Tn <- 20L
  unit <- rep(seq_len(nU), each = Tn)
  period <- rep(seq_len(Tn), times = nU)
  D <- as.integer(unit <= 2L & period >= 10L)   # units 1-2 switch at t=10
  alpha <- rnorm(nU)[unit]
  Y <- alpha + 0.1 * period + 0.4 * D + rnorm(nU * Tn, 0, 0.5)
  data.frame(unit, period, Y, D)
}

test_that("sparse-switcher panels estimate at high reps instead of crashing", {
  p <- make_sparse_panel()
  expect_warning(
    fit <- didgpu(df = p, outcome = "Y", group = "unit", time = "period",
                  treatment = "D", effects = 5L, placebo = 3L,
                  cluster = "unit", bootstrap_reps = 400L, seed = 1L,
                  backend = "cpu", verbose = FALSE),
    regexp = "degenerate cells"
  )
  E <- fit$results$Effects
  expect_true(is.matrix(E) && nrow(E) >= 1L)
  expect_true(is.finite(E[1, "Estimate"]))
  expect_true(is.finite(E[1, "SE"]) && E[1, "SE"] > 0)

  # accounting: dropped + surviving == requested reps
  expect_true(fit$results$n_boot_dropped > 0L)
  expect_identical(fit$results$n_boot + fit$results$n_boot_dropped, 400L)
})

test_that("panels with no degenerate resamples report n_boot_dropped == 0", {
  # dense switching: half the units switch -> excluding all of them in a
  # resample is essentially impossible, so nothing should be dropped.
  set.seed(7)
  nU <- 30L; Tn <- 20L
  unit <- rep(seq_len(nU), each = Tn)
  period <- rep(seq_len(Tn), times = nU)
  D <- as.integer(unit <= 15L & period >= 10L)
  Y <- rnorm(nU)[unit] + 0.1 * period + 0.4 * D + rnorm(nU * Tn, 0, 0.5)
  p <- data.frame(unit, period, Y, D)
  expect_no_warning(
    fit <- didgpu(df = p, outcome = "Y", group = "unit", time = "period",
                  treatment = "D", effects = 5L, placebo = 3L,
                  cluster = "unit", bootstrap_reps = 200L, seed = 1L,
                  backend = "cpu", verbose = FALSE)
  )
  expect_identical(fit$results$n_boot_dropped, 0L)
  expect_identical(fit$results$n_boot, 200L)
})

test_that("point estimates are unaffected by the degenerate-drop (iter 0 untouched)", {
  # the drop only filters bootstrap iterations; the point estimate comes
  # from the full-sample cell b=0 and must be identical at any rep count.
  p <- make_sparse_panel()
  f_lo <- suppressWarnings(
    didgpu(df = p, outcome = "Y", group = "unit", time = "period",
           treatment = "D", effects = 5L, placebo = 3L, cluster = "unit",
           bootstrap_reps = 25L, seed = 1L, backend = "cpu", verbose = FALSE))
  f_hi <- suppressWarnings(
    didgpu(df = p, outcome = "Y", group = "unit", time = "period",
           treatment = "D", effects = 5L, placebo = 3L, cluster = "unit",
           bootstrap_reps = 400L, seed = 99L, backend = "cpu", verbose = FALSE))
  expect_equal(f_lo$results$Effects[, "Estimate"],
               f_hi$results$Effects[, "Estimate"], tolerance = 1e-12)
})
