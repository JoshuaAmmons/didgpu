# Tests for didgpu_cs_continuous() — Callaway-Goodman-Bacon-Sant'Anna (2024).

# A 2-period continuous-dose panel: half never-treated (G=0, D=0), half
# treated at period 2 with a continuous dose and a concave dose-response.
make_cont_panel <- function(seed = 1L, nU = 400L) {
  set.seed(seed)
  half <- nU %/% 2L
  G <- rep(c(0, 2), times = c(nU - half, half))
  D <- ifelse(G == 2, stats::runif(nU, 0.1, 1), 0)
  a <- stats::rnorm(nU)
  d <- data.frame(id = rep(seq_len(nU), each = 2L),
                  time_period = rep(1:2, times = nU),
                  G = rep(G, each = 2L), D = rep(D, each = 2L))
  d$Y <- a[d$id] + 0.3 * d$time_period +
         ifelse(d$time_period == 2L, 2 * d$D - 1.2 * d$D^2, 0) +
         stats::rnorm(nrow(d), 0, 0.3)
  d
}

test_that("dose-response has the documented structure", {
  skip_if_not_installed("splines2")
  d <- make_cont_panel()
  r <- didgpu_cs_continuous(d, "Y", "D", "G", "time_period", "id",
                            dvals = seq(0.2, 0.9, length.out = 5), degree = 3,
                            num_knots = 2, bootstrap_reps = 0, verbose = FALSE)
  expect_s3_class(r, "didgpu_cs_continuous_result")
  expect_length(r$att.d, 5L); expect_length(r$acrt.d, 5L)
  expect_true(all(is.finite(r$att.d)) && all(is.finite(r$acrt.d)))
  expect_true(is.finite(r$overall_att) && is.finite(r$overall_acrt))
})

test_that("recovers a roughly concave dose-response (ACRT decreasing)", {
  skip_if_not_installed("splines2")
  d <- make_cont_panel(seed = 5L, nU = 1500L)
  # degree >= 2 so the spline derivative (ACRT) can vary; degree 1 would give
  # a constant ACRT. The DGP marginal effect 2 - 2.4*d is decreasing.
  r <- didgpu_cs_continuous(d, "Y", "D", "G", "time_period", "id",
                            dvals = c(0.2, 0.8), degree = 2, num_knots = 0,
                            bootstrap_reps = 0, verbose = FALSE)
  expect_gt(r$acrt.d[1], r$acrt.d[2])
})

test_that("bootstrap SEs are finite and non-negative", {
  skip_if_not_installed("splines2")
  d <- make_cont_panel()
  r <- didgpu_cs_continuous(d, "Y", "D", "G", "time_period", "id",
                            degree = 3, num_knots = 1, bootstrap_reps = 50,
                            seed = 2L, verbose = FALSE)
  expect_true(all(is.finite(r$att.d_se)) && all(r$att.d_se >= 0))
  expect_true(all(is.finite(r$acrt.d_se)))
  expect_length(r$att.d_lower, length(r$dose))
})

test_that("errors on misuse", {
  skip_if_not_installed("splines2")
  d <- make_cont_panel()
  expect_error(didgpu_cs_continuous(d, "nope", "D", "G", "time_period", "id",
                                    verbose = FALSE), "column not in df")
  d2 <- d; d2$G[d2$id %in% 1:40] <- 1L     # inject a second treated cohort
  expect_error(didgpu_cs_continuous(d2, "Y", "D", "G", "time_period", "id",
                                    verbose = FALSE), "single treated cohort")
})

test_that("print runs and returns invisibly", {
  skip_if_not_installed("splines2")
  d <- make_cont_panel()
  r <- didgpu_cs_continuous(d, "Y", "D", "G", "time_period", "id",
                            degree = 3, num_knots = 1, bootstrap_reps = 0, verbose = FALSE)
  expect_output(print(r), "continuous-treatment")
  expect_identical(withVisible(print(r))$visible, FALSE)
})

# Primary correctness check: ATT(d)/ACRT(d) match contdid::cont_did exactly.
test_that("ATT(d) and ACRT(d) match contdid::cont_did", {
  skip_if_not_installed("splines2")
  skip_if_not_installed("contdid")
  set.seed(1); d <- contdid::simulate_contdid_data(n = 600, num_time_periods = 2)
  dv <- seq(0.2, 0.9, length.out = 5)
  ref <- suppressWarnings(contdid::cont_did(
    yname = "Y", dname = "D", gname = "G", tname = "time_period", idname = "id",
    data = d, target_parameter = "level", aggregation = "dose",
    treatment_type = "continuous", control_group = "nevertreated",
    degree = 3, num_knots = 2, dvals = dv, bstrap = FALSE, cband = FALSE))
  ours <- didgpu_cs_continuous(d, "Y", "D", "G", "time_period", "id",
                               dvals = dv, degree = 3, num_knots = 2,
                               bootstrap_reps = 0, verbose = FALSE)
  expect_equal(ours$att.d,  as.numeric(ref$att.d),  tolerance = 1e-9)
  expect_equal(ours$acrt.d, as.numeric(ref$acrt.d), tolerance = 1e-9)
  expect_equal(ours$overall_att,  ref$overall_att,  tolerance = 1e-9)
  expect_equal(ours$overall_acrt, ref$overall_acrt, tolerance = 1e-9)
})
