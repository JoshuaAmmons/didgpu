# Tests for didgpu_did_continuous() — de Chaisemartin & D'Haultfoeuille (2024)
# continuous treatment, no stayers (first-difference design).

# A 2-period panel with a continuous dose CHANGE dD ~ N(0,1) for every unit
# (no stayers) and a known concave dose-response:
#   dY = 0.3 + 2*dD - 0.5*dD^2 + noise  =>
#     effect(d) = E[dY|d] - E[dY|0] = 2d - 0.5 d^2
#     ACR(d)    = d/dd E[dY|d]      = 2 - d
make_nostayer_panel <- function(seed = 1L, nU = 3000L, sd_noise = 0.5) {
  set.seed(seed)
  dD <- stats::rnorm(nU, 0, 1)
  dY <- 0.3 + 2 * dD - 0.5 * dD^2 + stats::rnorm(nU, 0, sd_noise)
  data.frame(id = rep(seq_len(nU), each = 2L),
             time_period = rep(1:2, nU),
             D = as.numeric(rbind(0, dD)),
             Y = as.numeric(rbind(0, dY)))
}

test_that("result has the documented structure (parametric)", {
  d <- make_nostayer_panel()
  dv <- c(-1, -0.5, 0.5, 1)
  r <- didgpu_did_continuous(d, "Y", "D", "id", "time_period",
                             estimator = "parametric", degree = 2, dvals = dv,
                             bootstrap_reps = 0, verbose = FALSE)
  expect_s3_class(r, "didgpu_did_continuous_result")
  expect_length(r$effect.d, 4L); expect_length(r$acr.d, 4L)
  expect_true(all(is.finite(r$effect.d)) && all(is.finite(r$acr.d)))
  expect_true(is.finite(r$overall_acr))
  expect_identical(r$estimator, "parametric")
})

test_that("parametric recovers the known dose-response", {
  d <- make_nostayer_panel(seed = 3L, nU = 6000L)
  dv <- c(-1, -0.5, 0.5, 1, 1.5)
  r <- didgpu_did_continuous(d, "Y", "D", "id", "time_period",
                             estimator = "parametric", degree = 2, dvals = dv,
                             bootstrap_reps = 0, verbose = FALSE)
  true_eff <- 2 * dv - 0.5 * dv^2
  true_acr <- 2 - dv
  expect_lt(max(abs(r$effect.d - true_eff)), 0.05)
  expect_lt(max(abs(r$acr.d   - true_acr)), 0.08)
  # ACR is decreasing in d (concave response): 2 - d
  expect_gt(r$acr.d[1], r$acr.d[5])
})

test_that("nonparametric local-linear recovers the same shape", {
  skip_on_cran()
  d <- make_nostayer_panel(seed = 7L, nU = 8000L)
  dv <- c(-0.5, 0.5, 1)
  r <- didgpu_did_continuous(d, "Y", "D", "id", "time_period",
                             estimator = "nonparametric", dvals = dv,
                             bootstrap_reps = 0, verbose = FALSE)
  expect_identical(r$estimator, "nonparametric")
  true_eff <- 2 * dv - 0.5 * dv^2
  true_acr <- 2 - dv
  expect_lt(max(abs(r$effect.d - true_eff)), 0.10)
  expect_lt(max(abs(r$acr.d   - true_acr)), 0.15)
})

test_that("bootstrap SEs are finite and non-negative", {
  d <- make_nostayer_panel()
  dv <- c(-0.5, 0.5, 1)
  r <- didgpu_did_continuous(d, "Y", "D", "id", "time_period",
                             estimator = "parametric", degree = 2, dvals = dv,
                             bootstrap_reps = 50, seed = 2L, verbose = FALSE)
  expect_true(all(is.finite(r$effect.d_se)) && all(r$effect.d_se >= 0))
  expect_true(all(is.finite(r$acr.d_se)) && all(r$acr.d_se >= 0))
  expect_length(r$effect.d_lower, length(dv))
  expect_true(all(r$effect.d_lower <= r$effect.d_upper))
})

test_that("errors on misuse", {
  d <- make_nostayer_panel(nU = 50L)
  expect_error(didgpu_did_continuous(d, "nope", "D", "id", "time_period",
                                     verbose = FALSE), "column not in df")
  # too few units for the requested polynomial degree
  d1 <- make_nostayer_panel(nU = 3L)
  expect_error(didgpu_did_continuous(d1, "Y", "D", "id", "time_period",
                                     degree = 5, verbose = FALSE),
               "too few units")
})

test_that("print runs and returns invisibly", {
  d <- make_nostayer_panel(nU = 500L)
  r <- didgpu_did_continuous(d, "Y", "D", "id", "time_period",
                             estimator = "parametric", degree = 2,
                             bootstrap_reps = 0, verbose = FALSE)
  expect_output(print(r), "continuous DiD")
  expect_identical(withVisible(print(r))$visible, FALSE)
  # nonparametric print carries the EXPERIMENTAL marker
  rn <- didgpu_did_continuous(d, "Y", "D", "id", "time_period",
                              estimator = "nonparametric",
                              bootstrap_reps = 0, verbose = FALSE)
  expect_output(print(rn), "EXPERIMENTAL")
})
