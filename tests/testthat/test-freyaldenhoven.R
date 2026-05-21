# Tests for didgpu_freyaldenhoven() — FHS (2019) pre-event proxy event study.

make_fhs_panel <- function(seed = 1L, nU = 80L, Tn = 24L) {
  set.seed(seed)
  d <- data.frame(id = rep(1:nU, each = Tn), t = rep(1:Tn, times = nU))
  # Adoption spread widely across interior periods (+ ~25% never-treated) so
  # the lead/lag endpoint design has ample variation after FE absorption.
  adopt <- sample(c(rep(Inf, nU %/% 4L),
                    sample(6:20, nU - nU %/% 4L, replace = TRUE)))
  d$z <- as.integer(d$t >= adopt[d$id])
  eta <- rnorm(nU)[d$id] + 0.2 * d$t                  # confound (trends)
  d$x <- eta + rnorm(nrow(d), 0, 0.3)                 # proxy responds to confound
  d$y <- rnorm(nU)[d$id] + 0.1 * d$t + 0.5 * d$z + 0.4 * eta + rnorm(nrow(d), 0, 0.3)
  d
}

test_that("OLS returns the expected first-difference event-study term set", {
  d <- make_fhs_panel()
  r <- didgpu_freyaldenhoven(d, "y", "z", "id", "t", estimator = "OLS",
                             pre = 0, post = 3, verbose = FALSE)
  expect_s3_class(r, "didgpu_freyaldenhoven_result")
  expect_true(all(c("Estimate", "SE", "LB.CI", "UB.CI") %in% colnames(r$coefficients)))
  nm <- rownames(r$coefficients)
  expect_true(all(c("z_lead3", "z_fd_lead2", "z_fd", "z_fd_lag1", "z_lag4") %in% nm))
  expect_false("z_fd_lead1" %in% nm)            # normalized out
  expect_true(all(is.finite(r$coefficients[, "Estimate"])))
})

test_that("FHS adds the proxy, auto-selects an instrument, drops it from regressors", {
  d <- make_fhs_panel()
  r <- didgpu_freyaldenhoven(d, "y", "z", "id", "t", estimator = "FHS",
                             proxy = "x", pre = 0, post = 3, verbose = FALSE)
  expect_true("x" %in% rownames(r$coefficients))
  expect_match(r$proxyIV, "^z_fd_lead")
  expect_false(r$proxyIV %in% rownames(r$coefficients))
})

test_that("FHS requires a proxy", {
  d <- make_fhs_panel()
  expect_error(
    didgpu_freyaldenhoven(d, "y", "z", "id", "t", estimator = "FHS", verbose = FALSE),
    "proxy")
})

test_that("print runs and returns invisibly", {
  d <- make_fhs_panel()
  r <- didgpu_freyaldenhoven(d, "y", "z", "id", "t", estimator = "OLS",
                             pre = 0, post = 3, verbose = FALSE)
  expect_output(print(r), "Freyaldenhoven")
  expect_identical(withVisible(print(r))$visible, FALSE)
})

# Primary correctness check: event-study coefficients match eventstudyr's
# OLS and FHS to numerical precision (skipped where eventstudyr is absent).
test_that("matches eventstudyr OLS + FHS coefficients when available", {
  skip_if_not_installed("eventstudyr")
  d <- as.data.frame(eventstudyr::example_data)
  for (est in c("OLS", "FHS")) {
    es_args <- list(estimator = est, data = d, outcomevar = "y_base",
                    policyvar = "z", idvar = "id", timevar = "t",
                    pre = 0, post = 3, kernel = "estimatr")
    if (est == "FHS") es_args$proxy <- "x_r"
    es <- suppressWarnings(do.call(eventstudyr::EventStudy, es_args))
    ref <- tryCatch(stats::coef(es$output), error = function(e) NULL)
    skip_if(is.null(ref), "could not extract eventstudyr coefficients")
    our_args <- list(df = d, outcome = "y_base", policy = "z", id = "id",
                     time = "t", estimator = est, pre = 0, post = 3, verbose = FALSE)
    if (est == "FHS") our_args$proxy <- "x_r"
    ours <- do.call(didgpu_freyaldenhoven, our_args)$coefficients[, "Estimate"]
    common <- intersect(names(ref), names(ours))
    expect_gt(length(common), 5L)
    expect_equal(unname(ours[common]), unname(ref[common]),
                 tolerance = 1e-6, info = paste("estimator", est))
  }
})
