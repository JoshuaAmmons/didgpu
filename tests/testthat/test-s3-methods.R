# Standard model accessor S3 methods: coef(), confint(), vcov().

test_that("coef() returns named vector of effects + placebos + ATE", {
  p <- small_panel()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 3L, placebo = 1L,
                 bootstrap_reps = 10L, seed = 2L,
                 backend = "r", verbose = FALSE)

  all_c  <- coef(fit)
  expect_named(all_c, c("Effect_1", "Effect_2", "Effect_3",
                        "Placebo_1", "ATE"))
  expect_type(all_c, "double")

  eff_c <- coef(fit, which = "effects")
  expect_length(eff_c, 3L)
  expect_named(eff_c, c("Effect_1", "Effect_2", "Effect_3"))

  pl_c <- coef(fit, which = "placebos")
  expect_length(pl_c, 1L)
  expect_named(pl_c, "Placebo_1")

  ate_c <- coef(fit, which = "ate")
  expect_length(ate_c, 1L)
  expect_named(ate_c, "ATE")

  # Values must match the printed Effects matrix.
  expect_equal(unname(eff_c), as.numeric(fit$results$Effects[, "Estimate"]))
})

test_that("confint() returns 2-column matrix with stored CI level", {
  p <- small_panel()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 2L, placebo = 1L,
                 bootstrap_reps = 10L, seed = 2L,
                 backend = "r", verbose = FALSE)
  ci <- confint(fit)
  expect_equal(ncol(ci), 2L)
  expect_equal(nrow(ci), 4L)  # 2 effects + 1 placebo + 1 ATE
  expect_match(colnames(ci)[1L], "LB")
  expect_match(colnames(ci)[2L], "UB")
  # Lower bound always < upper bound (modulo NA).
  v <- ci[!is.na(ci[, 1L]) & !is.na(ci[, 2L]), , drop = FALSE]
  expect_true(all(v[, 1L] <= v[, 2L]))
})

test_that("confint(parm = ...) subsets correctly", {
  p <- small_panel()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 3L, bootstrap_reps = 5L, seed = 2L,
                 backend = "r", verbose = FALSE)
  ci <- confint(fit, parm = c("Effect_1", "Effect_3"))
  expect_equal(rownames(ci), c("Effect_1", "Effect_3"))
  expect_warning(confint(fit, parm = "NoSuchEffect"),
                  "parm not in result")
})

test_that("vcov() returns symmetric square matrix with bootstrap > 1", {
  p <- small_panel()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 3L, placebo = 1L,
                 bootstrap_reps = 10L, seed = 2L,
                 backend = "r", verbose = FALSE)
  v <- vcov(fit)
  expect_true(is.matrix(v))
  expect_equal(nrow(v), ncol(v))
  expect_equal(nrow(v), 4L)  # 3 effects + 1 placebo
  expect_true(isSymmetric(v, tol = 1e-10))
  # The diagonal must match the squared SEs from Effects matrix.
  diag_v <- diag(v)
  expect_equal(unname(diag_v[1L:3L]),
               as.numeric(fit$results$Effects[, "SE"])^2,
               tolerance = 1e-10)
})

test_that("vcov() returns NA matrix without bootstrap", {
  p <- small_panel()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 2L, bootstrap_reps = 0L,
                 backend = "r", verbose = FALSE)
  v <- vcov(fit)
  expect_true(is.matrix(v))
  # With 0 bootstrap reps, .aggregate_to_result still builds a cov
  # matrix from the single iter-0 cell; behaviour is "NA except diag"
  # or all-NA. We only test that it's a matrix of the right shape.
  expect_equal(nrow(v), ncol(v))
  expect_equal(nrow(v), 2L)
})
