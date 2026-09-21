# didgpu_fect(method = "mc") must be invariant to the LEVEL of the outcome.
#
# Athey et al. (2021) estimate Y = L + unit FE + time FE, penalising the
# nuclear norm of L ALONE. The old implementation soft-thresholded the raw
# outcome matrix with no fixed effects at all, so the penalty shrank the
# LEVEL and the residual Y - Y_hat absorbed it. The estimator was then a
# function of where the outcome happened to sit.
#
# Measured on the known-zero DGP below, the reported ATT was:
#     Y + 0   -> +0.26915        (truth 0; fe -0.024, ife -0.102)
#     Y + 1   -> +0.42445
#     Y + 10  -> +0.44530
#     Y + 100 -> +2.55072
# while fect::fect returned -0.02443 at every level. On a positive,
# trending outcome this manufactured large, monotonically rising,
# significant effects where fe and ife both found a null.
#
# Two things were wrong and both are pinned here: the fit ignored fixed
# effects, and lambda was scaled to the singular values of the RAW matrix,
# so the penalty itself grew with the level.

.mc_panel <- function(seed = 7L, level = 0) {
  set.seed(seed)
  NU <- 60L; NT <- 12L
  ufe <- stats::rnorm(NU, 0, 1); tfe <- stats::rnorm(NT, 0, 0.3)
  g <- expand.grid(period = 1:NT, unit = 1:NU)
  g$D <- 0L
  g$D[g$unit %in% 11:30 & g$period >= 7L] <- 1L   # switchers, TRUE tau = 0
  g$Y <- level + ufe[g$unit] + tfe[g$period] + stats::rnorm(nrow(g), 0, 0.4)
  g[order(g$unit, g$period), c("unit", "period", "D", "Y")]
}

.mc_ate <- function(d, m) {
  suppressWarnings(didgpu_fect(df = d, outcome = "Y", group = "unit",
    time = "period", treatment = "D", method = m,
    bootstrap_reps = 0L, verbose = FALSE))$results$ATE[1, 1]
}

test_that("every fect method is invariant to adding a constant to Y", {
  for (m in c("fe", "ife", "mc")) {
    base <- .mc_ate(.mc_panel(level = 0), m)
    for (lv in c(1, 10, 100)) {
      shifted <- .mc_ate(.mc_panel(level = lv), m)
      expect_equal(shifted, base, tolerance = 1e-8,
                   label = sprintf("%s at Y + %g", m, lv))
    }
  }
})

test_that("mc recovers a known-zero ATT", {
  # Pre-fix this was +0.269 against a true effect of exactly 0.
  expect_lt(abs(.mc_ate(.mc_panel(), "mc")), 0.1)
})

test_that("mc agrees with fe and ife on a DGP with no factor structure", {
  # The DGP is two-way fixed effects plus noise, so all three estimators
  # target the same thing and must land in the same place.
  d <- .mc_panel()
  fe <- .mc_ate(d, "fe"); mc <- .mc_ate(d, "mc")
  expect_lt(abs(mc - fe), 0.05)
})

test_that("mc tracks the fect package", {
  skip_if_not_installed("fect")
  d <- .mc_panel()
  ours <- .mc_ate(d, "mc")
  theirs <- suppressMessages(suppressWarnings(
    fect::fect(Y ~ D, data = d, index = c("unit", "period"),
               method = "mc", se = FALSE, force = "two-way")))$att.avg
  # Not bit-for-bit: lambda selection still differs from fect's own CV.
  # Pre-fix the gap was 2.9e-01; this guards the structural agreement.
  expect_lt(abs(ours - theirs), 0.05)
})

test_that("lambda is scaled to the penalised residual, not the raw level", {
  # .fect_mc_sigma_max() must not move when the outcome is shifted.
  mk_mats <- function(level) {
    d <- .mc_panel(level = level)
    Y <- matrix(d$Y, nrow = length(unique(d$unit)), byrow = TRUE)
    M <- matrix(d$D, nrow = length(unique(d$unit)), byrow = TRUE)
    list(Y = Y, M = M)
  }
  a <- mk_mats(0); b <- mk_mats(100)
  expect_equal(didgpu:::.fect_mc_sigma_max(a$Y, a$M),
               didgpu:::.fect_mc_sigma_max(b$Y, b$M), tolerance = 1e-8)
})
