# Tests for didgpu_bacon() — Goodman-Bacon decomposition (#101).

# Balanced, binary, absorbing staggered panel with never-treated controls
# and NO always-treated units (so the decomposition sample == full sample,
# letting us cross-check beta_twfe against didgpu_twfe()).
make_stag <- function(seed = 101L) {
  set.seed(seed)
  nU <- 90L; Tn <- 12L
  adopt <- sample(c(rep(Inf, 30), rep(4, 20), rep(7, 20), rep(10, 20)))
  unit <- rep(1:nU, each = Tn); period <- rep(1:Tn, times = nU)
  adopt_i <- adopt[unit]
  D <- as.integer(period >= adopt_i)
  alpha <- rnorm(nU)[unit]; lambda <- period / 5
  tsa <- ifelse(is.finite(adopt_i) & period >= adopt_i, period - adopt_i + 1L, 0L)
  te <- ifelse(tsa > 0, 0.5 * tsa + 0.3 * (adopt_i == 7), 0)   # heterogeneous
  Y <- alpha + lambda + te + rnorm(nU * Tn, 0, 0.5)
  data.frame(unit, period, Y, D)
}

test_that("the Goodman-Bacon identity holds: sum(weight * 2x2) == TWFE beta", {
  bd <- didgpu_bacon(make_stag(), "Y", "unit", "period", "D")
  expect_s3_class(bd, "didgpu_bacon")
  expect_equal(bd$beta_check, bd$beta_twfe, tolerance = 1e-9)
  expect_equal(sum(bd$comparisons$weight), 1, tolerance = 1e-10)
  expect_true(all(bd$comparisons$weight >= -1e-12))
})

test_that("beta_twfe matches the independent didgpu_twfe() coefficient", {
  p  <- make_stag()
  bd <- didgpu_bacon(p, "Y", "unit", "period", "D")
  tw <- didgpu_twfe(p, "Y", "unit", "period", "D", effects = 1L,
                    placebo = 0L, verbose = FALSE)
  expect_equal(bd$beta_twfe, as.numeric(tw$coef["Effect_1"]), tolerance = 1e-7)
})

test_that("decomposition has the three comparison types and a sane forbidden weight", {
  bd <- didgpu_bacon(make_stag(), "Y", "unit", "period", "D")
  expect_setequal(bd$summary$type,
                  c("Treated vs Untreated", "Earlier vs Later Treated",
                    "Later vs Earlier Treated"))
  expect_true(bd$forbidden_weight > 0 && bd$forbidden_weight < 1)
  # per-type weights sum to the total (== 1)
  expect_equal(sum(bd$summary$weight), 1, tolerance = 1e-10)
})

test_that("always-treated units are dropped and the identity still holds", {
  p <- make_stag()
  # Force units 1..5 to be always-treated (D == 1 in every period).
  p$D[p$unit %in% 1:5] <- 1L
  bd <- didgpu_bacon(p, "Y", "unit", "period", "D")
  expect_equal(bd$n_always_treated, 5L)
  expect_equal(bd$beta_check, bd$beta_twfe, tolerance = 1e-9)
})

test_that("no never-treated group: only timing-pair comparisons", {
  # All units adopt eventually (no Inf), none in period 1.
  set.seed(7); nU <- 60L; Tn <- 10L
  adopt <- sample(rep(c(3, 6, 9), each = 20))
  unit <- rep(1:nU, each = Tn); period <- rep(1:Tn, times = nU); a <- adopt[unit]
  D <- as.integer(period >= a)
  Y <- rnorm(nU)[unit] + period / 4 +
       ifelse(period >= a, 0.5 * (period - a + 1L), 0) + rnorm(nU * Tn, 0, 0.4)
  bd <- didgpu_bacon(data.frame(unit, period, Y, D), "Y", "unit", "period", "D")
  expect_false(bd$has_never_treated)
  expect_false("Treated vs Untreated" %in% bd$summary$type)
  expect_equal(bd$beta_check, bd$beta_twfe, tolerance = 1e-9)
})

test_that("informative errors on out-of-scope inputs", {
  base <- expand.grid(period = 1:4, unit = 1:6)
  base <- base[order(base$unit, base$period), ]
  set.seed(1); base$Y <- rnorm(nrow(base))

  # non-binary treatment
  b_bin <- base; b_bin$D <- 0L; b_bin$D[1] <- 2L
  expect_error(didgpu_bacon(b_bin, "Y", "unit", "period", "D"), "binary")

  # non-absorbing: unit 1 turns treatment off (0,1,0,0)
  b_abs <- base; b_abs$D <- 0L; b_abs$D[b_abs$unit == 1] <- c(0L, 1L, 0L, 0L)
  expect_error(didgpu_bacon(b_abs, "Y", "unit", "period", "D"), "staggered")

  # unbalanced (drop a row)
  b_unb <- base[-1, ]; b_unb$D <- 0L
  expect_error(didgpu_bacon(b_unb, "Y", "unit", "period", "D"), "unbalanced")
})

test_that("print runs and returns invisibly", {
  bd <- didgpu_bacon(make_stag(), "Y", "unit", "period", "D")
  expect_output(print(bd), "Goodman-Bacon")
  expect_identical(withVisible(print(bd))$visible, FALSE)
})
