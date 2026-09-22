# Always-treated units must be dropped from the fect estimation sample.
#
# A unit treated in every observed period has no control cell, so its unit
# fixed effect (fe) / factor loading (ife, mc) is unidentified.
# .fect_fe_fit() sets an unidentified alpha to 0, which makes the imputed
# counterfactual Y_hat = xi alone -- the unit's entire LEVEL then lands in
# the residual and is reported as treatment effect.
#
# Always-treated units are selected on level (they are precisely the units
# already treated before the window opened), so the bias does not average
# out. On the known-zero DGP below, retaining 10 such units reported
# ATE = +1.63 against a true effect of exactly 0, and it survived every
# factor count r = 1..4 as well as method = "fe".

# Known-zero DGP: true ATT is exactly 0, and the 10 always-treated units
# are given a high level (+3), mimicking units already treated before the
# sample window.
.always_treated_panel <- function(seed = 7L) {
  set.seed(seed)
  NU <- 60L; NT <- 12L
  ufe <- stats::rnorm(NU, 0, 1); tfe <- stats::rnorm(NT, 0, 0.3)
  ufe[1:10] <- ufe[1:10] + 3
  g <- expand.grid(period = 1:NT, unit = 1:NU)
  g$D <- 0L
  g$D[g$unit %in% 1:10] <- 1L                      # always-treated
  g$D[g$unit %in% 11:30 & g$period >= 7L] <- 1L    # switchers at t = 7
  g$Y <- ufe[g$unit] + tfe[g$period] +             # TRUE tau = 0
         stats::rnorm(nrow(g), 0, 0.4)
  g[order(g$unit, g$period), c("unit", "period", "D", "Y")]
}

.fit <- function(d, method) {
  didgpu_fect(df = d, outcome = "Y", group = "unit", time = "period",
              treatment = "D", method = method, bootstrap_reps = 0L,
              verbose = FALSE)
}

test_that("always-treated units are dropped, warned about, and counted", {
  p <- .always_treated_panel()
  # Always-treated units have zero untreated periods, so they fail any
  # min_T0 >= 1. The warning now names the min_T0 rule.
  expect_warning(f <- .fit(p, "fe"), "always-treated")
  expect_equal(f$n_always_treated_dropped, 10L)
})

test_that("retaining always-treated units no longer biases the ATT", {
  p <- .always_treated_panel()
  for (m in c("fe", "ife", "mc")) {
    suppressWarnings(f <- .fit(p, m))
    # Pre-fix this was ~+1.6 to +1.8 against a true effect of 0.
    expect_lt(abs(f$results$ATE[1L, 1L]), 0.5, label = paste("ATE for", m))
  }
})

test_that("dropping always-treated by hand changes nothing", {
  # The estimate must be invariant to whether the caller pre-filters, which
  # is the strongest statement that they contribute nothing.
  p <- .always_treated_panel()
  q <- p[!p$unit %in% 1:10, ]
  for (m in c("fe", "ife", "mc")) {
    suppressWarnings(a <- .fit(p, m))
    b <- .fit(q, m)
    expect_equal(a$results$ATE[1L, 1L], b$results$ATE[1L, 1L], tolerance = 1e-12)
    expect_equal(b$n_always_treated_dropped, 0L)
  }
})

test_that("a panel of only always-treated units is an error, not a number", {
  p <- .always_treated_panel()
  p <- p[p$unit %in% 1:10, ]
  expect_error(.fit(p, "fe"), "fewer than min_T0")
})

test_that("fe matches the fect package once always-treated are dropped", {
  skip_if_not_installed("fect")
  p <- .always_treated_panel()

  # NOTE on tolerance. didgpu's alternating-projections fit stops at
  # `tol` (default 1e-5), which leaves a real gap against fect's answer:
  #   tol = 1e-5  -> 5.97e-05
  #   tol = 1e-8  -> 2.39e-06
  #   tol = 1e-12 -> 1.91e-08
  # The two estimators therefore agree DEFINITIONALLY, and didgpu is not
  # bit-for-bit with fect at its default tolerance. Pin tol here so this
  # test measures the definition rather than the stopping rule.
  ours <- suppressWarnings(
    didgpu_fect(df = p, outcome = "Y", group = "unit", time = "period",
                treatment = "D", method = "fe", bootstrap_reps = 0L,
                verbose = FALSE, tol = 1e-12, max_iter = 100000L)
  )
  theirs <- suppressMessages(suppressWarnings(
    fect::fect(Y ~ D, data = p, index = c("unit", "period"),
               method = "fe", se = FALSE, force = "two-way")
  ))
  expect_equal(unname(ours$results$ATE[1L, 1L]), unname(theirs$att.avg),
               tolerance = 1e-6)
})

test_that("default tolerance leaves a documented gap against fect", {
  skip_if_not_installed("fect")
  # Guards the claim above: if the default fit ever becomes exact, or
  # drifts materially worse, this test should be revisited rather than
  # the docs quietly becoming wrong.
  p <- .always_treated_panel()
  suppressWarnings(ours <- .fit(p, "fe"))
  theirs <- suppressMessages(suppressWarnings(
    fect::fect(Y ~ D, data = p, index = c("unit", "period"),
               method = "fe", se = FALSE, force = "two-way")
  ))
  expect_lt(abs(ours$results$ATE[1L, 1L] - theirs$att.avg), 1e-3)
})
