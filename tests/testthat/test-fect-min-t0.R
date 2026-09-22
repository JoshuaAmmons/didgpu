# Units need enough UNTREATED periods for their counterfactual to be
# identified, and fect drops the ones that do not have them.
#
#     method "fe"                       -> min.T0 = 1
#     method "ife" / "mc" / "both" / ... -> min.T0 = 5
#
# didgpu kept every unit and extrapolated. With a strong factor
# structure that is catastrophic, not marginal: on a 100-unit panel with
# 2 latent factors and a TRUE ATT of +1.0, method "ife" returned
# -0.17445 -- the wrong sign -- while fect returned +1.00064. Dropping
# the 16 units with fewer than 5 untreated periods (the same 84 fect
# keeps) moved didgpu to +0.94980.
#
# The rule also subsumes the always-treated case: such units have zero
# untreated periods and fail any min_T0 >= 1.

.mt0_panel <- function(seed = 1L, lsd = 3, fsd = 2.5, nu = 100L,
                       np = 15L, r = 2L, att = 1.0) {
  set.seed(seed)
  ufe <- stats::rnorm(nu, 0, 0.5); tfe <- stats::rnorm(np, 0, 0.3)
  L <- matrix(stats::rnorm(nu * r, 0, lsd), nu, r)
  F <- matrix(stats::rnorm(r * np, 0, fsd), r, np)
  Fg <- rep(Inf, nu)
  tr <- sort(sample(seq_len(nu), nu / 2))
  Fg[tr] <- sample(3:(np - 1), length(tr), replace = TRUE)
  g <- expand.grid(period = 1:np, unit = 1:nu)
  g <- g[order(g$unit, g$period), ]
  g$D <- as.integer(g$period >= Fg[g$unit])
  g$Y <- ufe[g$unit] + tfe[g$period] +
         (L %*% F)[cbind(g$unit, g$period)] + att * g$D +
         stats::rnorm(nrow(g), 0, 0.2)
  g[, c("unit", "period", "D", "Y")]
}

.mt0_fit <- function(p, method, ...) {
  suppressWarnings(didgpu_fect(df = p, outcome = "Y", group = "unit",
    time = "period", treatment = "D", method = method,
    bootstrap_reps = 0L, verbose = FALSE, ...))
}

test_that("min_T0 defaults follow fect: 1 for fe, 5 for ife and mc", {
  p <- .mt0_panel()
  expect_equal(.mt0_fit(p, "fe")$min_T0, 1L)
  expect_equal(.mt0_fit(p, "ife", r = 2L)$min_T0, 5L)
  expect_equal(.mt0_fit(p, "mc")$min_T0, 5L)
})

test_that("short-history units are dropped and counted", {
  p <- .mt0_panel()
  npre <- tapply(p$D, p$unit, function(x) sum(x == 0))
  f <- .mt0_fit(p, "ife", r = 2L)
  expect_equal(f$n_units_dropped, sum(npre < 5L))
  expect_gt(f$n_units_dropped, 0L)
})

test_that("ife recovers the ATT on a strong factor structure", {
  # Pre-fix this returned -0.17445 against a true effect of +1.0.
  est <- .mt0_fit(.mt0_panel(), "ife", r = 2L)$results$ATE[1L, 1L]
  expect_gt(est, 0.8)
  expect_lt(est, 1.2)
})

test_that("ife tracks fect across factor strengths", {
  skip_if_not_installed("fect")
  for (lsd in c(0.5, 1.5, 3)) {
    p <- .mt0_panel(lsd = lsd, fsd = lsd * 0.83)
    ours <- .mt0_fit(p, "ife", r = 2L)$results$ATE[1L, 1L]
    theirs <- suppressMessages(suppressWarnings(
      fect::fect(Y ~ D, data = p, index = c("unit", "period"),
                 method = "ife", r = 2L, se = FALSE,
                 force = "two-way")))$att.avg
    # Pre-fix the gap reached 1.18 (sign inversion); this guards the
    # substantive agreement, not bit-for-bit parity.
    expect_lt(abs(ours - theirs), 0.1, label = sprintf("lsd = %.1f", lsd))
  }
})

test_that("min_T0 can be overridden and changes the sample", {
  p <- .mt0_panel()
  loose <- .mt0_fit(p, "ife", r = 2L, min_T0 = 1L)
  tight <- .mt0_fit(p, "ife", r = 2L, min_T0 = 8L)
  expect_lt(loose$n_units_dropped, tight$n_units_dropped)
  expect_error(.mt0_fit(p, "ife", r = 2L, min_T0 = 0L), "positive integer")
})
