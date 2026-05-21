# Tests for didgpu_did_static() — dCDH (2020) DID_M instantaneous estimator.

test_that("DID_M matches a hand-computed value", {
  # 4 units x 3 periods. At t3: switch-in dY=3 vs stay0 dY=1 -> DID+=2;
  # switch-out dY=0 vs stay1 dY=1 -> DID-=1. DID_M = (2+1)/2 = 1.5.
  d <- data.frame(
    Y = c(1, 2, 5,  1, 2, 3,  4, 5, 5,  4, 5, 6),
    unit = rep(1:4, each = 3), period = rep(1:3, times = 4),
    D = c(0, 0, 1,  0, 0, 0,  1, 1, 0,  1, 1, 1))
  r <- didgpu_did_static(d, "Y", "unit", "period", "D",
                         bootstrap_reps = 0L, verbose = FALSE)
  expect_s3_class(r, "didgpu_did_static_result")
  expect_equal(r$did, 1.5, tolerance = 1e-10)
  expect_equal(r$n_switchers, 2)
})

test_that("a single absorbing switch reduces to the 2x2 DiD", {
  d <- data.frame(Y = c(1, 4, 2, 3), unit = c(1, 1, 2, 2),
                  period = c(1, 2, 1, 2), D = c(0, 1, 0, 0))
  r <- didgpu_did_static(d, "Y", "unit", "period", "D",
                         bootstrap_reps = 0L, verbose = FALSE)
  expect_equal(r$did, (4 - 1) - (3 - 2), tolerance = 1e-10)   # 2
})

test_that("switch-OUT direction is handled (reversible treatment)", {
  # unit 1 stays treated (control), unit 2 switches out 1->0.
  # DID- = stay1 dY - switchout dY = (6-5) - (3-5) = 3.
  d <- data.frame(Y = c(5, 6, 5, 3), unit = c(1, 1, 2, 2),
                  period = c(1, 2, 1, 2), D = c(1, 1, 1, 0))
  r <- didgpu_did_static(d, "Y", "unit", "period", "D",
                         bootstrap_reps = 0L, verbose = FALSE)
  expect_equal(r$did, (6 - 5) - (3 - 5), tolerance = 1e-10)   # 3
})

test_that("result structure + cluster bootstrap SE", {
  set.seed(1); nU <- 60L; Tn <- 6L
  d <- data.frame(unit = rep(1:nU, each = Tn), period = rep(1:Tn, times = nU))
  d$D <- as.integer(stats::runif(nrow(d)) < 0.4)            # non-absorbing on/off
  d$Y <- rnorm(nU)[d$unit] + 0.1 * d$period + 0.7 * d$D + rnorm(nrow(d), 0, 0.3)
  r <- didgpu_did_static(d, "Y", "unit", "period", "D",
                         bootstrap_reps = 40L, seed = 2L, verbose = FALSE)
  expect_true(is.finite(r$did))
  expect_true(is.finite(r$se) && r$se > 0)
  expect_length(r$ci, 2L)
  expect_true(all(c("time", "n_switch_in", "n_switch_out", "did_in", "did_out")
                  %in% names(r$per_period)))
})

test_that("weighted DID_M respects observation weights", {
  # Two switch-in units at t2 with different dY; weight one heavily.
  d <- data.frame(
    Y = c(0, 2,  0, 6,  0, 0,  0, 0),         # u1 dY=2, u2 dY=6, stay0 u3/u4 dY=0
    unit = rep(1:4, each = 2), period = rep(1:2, times = 4),
    D = c(0, 1,  0, 1,  0, 0,  0, 0),
    w = c(3, 3,  1, 1,  1, 1,  1, 1))
  r <- didgpu_did_static(d, "Y", "unit", "period", "D", weight = "w",
                         bootstrap_reps = 0L, verbose = FALSE)
  # weighted mean dY of switchers = (3*2 + 1*6)/4 = 3; stay0 dY = 0 -> DID_M=3.
  expect_equal(r$did, 3, tolerance = 1e-10)
})

test_that("errors on non-binary treatment", {
  d <- data.frame(Y = 1:4, unit = c(1, 1, 2, 2), period = c(1, 2, 1, 2),
                  D = c(0, 2, 0, 1))
  expect_error(didgpu_did_static(d, "Y", "unit", "period", "D", verbose = FALSE),
               "binary")
})

test_that("print runs and returns invisibly", {
  d <- data.frame(Y = c(1, 4, 2, 3), unit = c(1, 1, 2, 2),
                  period = c(1, 2, 1, 2), D = c(0, 1, 0, 0))
  r <- didgpu_did_static(d, "Y", "unit", "period", "D",
                         bootstrap_reps = 0L, verbose = FALSE)
  expect_output(print(r), "DID_M")
  expect_identical(withVisible(print(r))$visible, FALSE)
})

# Cross-check against the reference DID_M when DIDmultiplegt is installed
# (skipped in CI / on machines without it).
test_that("matches DIDmultiplegt::did_multiplegt DID_M when available", {
  skip_if_not_installed("DIDmultiplegt")
  set.seed(11); nU <- 100L; Tn <- 6L
  d <- data.frame(unit = rep(1:nU, each = Tn), period = rep(1:Tn, times = nU))
  d$D <- as.integer(stats::runif(nrow(d)) < 0.45)
  d$Y <- rnorm(nU)[d$unit] + 0.1 * d$period + 0.8 * d$D + rnorm(nrow(d), 0, 0.4)
  ours <- didgpu_did_static(d, "Y", "unit", "period", "D",
                            bootstrap_reps = 0L, verbose = FALSE)$did
  # Modern DIDmultiplegt: did_multiplegt(mode = "old", df, Y, G, T, D);
  # the DID_M point estimate is $effect.
  ref <- tryCatch(
    as.numeric(DIDmultiplegt::did_multiplegt(
      mode = "old", df = d, Y = "Y", G = "unit", T = "period", D = "D")$effect),
    error = function(e) NA_real_)
  skip_if(length(ref) != 1L || is.na(ref), "did_multiplegt return shape not recognized")
  expect_equal(ours, ref, tolerance = 1e-6)
})
