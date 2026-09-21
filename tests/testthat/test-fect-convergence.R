# fect fits must converge at default settings, and must say so.
#
# Two regressions are pinned here.
#
# 1. .fect_fe_fit() had a second stopping rule:
#        if (iter > 1L && abs(loss - prev_loss) < tol) break
#    `loss` is a SUM of squared residuals, so its absolute change falls
#    below a tolerance meant for parameter units long before alpha and xi
#    have settled. It fired first: on a 60x10 panel the fit exited after
#    6 iterations with delta = 2.0e-04, twenty times the requested 1e-05,
#    and reported itself done. That early exit is what kept method "fe"
#    3.3e-05 away from fect::fect at default settings; removing it took
#    the gap to 7.8e-07.
#
# 2. The solver reported iter / delta / lambda per cell, but nothing
#    surfaced them, so a caller could not distinguish a converged fit
#    from one that had exhausted max_iter. didgpu_fect() now carries
#    $diagnostics and warns when the fit did not converge.

.conv_panel <- function(seed = 11L) {
  as.data.frame(didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                                      seed = seed))
}
.conv_fit <- function(p, method, ...) {
  suppressWarnings(didgpu_fect(df = p, outcome = "Y", group = "unit",
    time = "period", treatment = "D", method = method,
    bootstrap_reps = 0L, verbose = FALSE, ...))
}

test_that("didgpu_fect reports convergence diagnostics", {
  p <- .conv_panel()
  for (m in c("fe", "ife", "mc")) {
    d <- .conv_fit(p, m)$diagnostics
    expect_false(is.null(d), label = m)
    expect_true(all(c("iter", "delta", "converged", "tol", "max_iter") %in%
                      names(d)), label = m)
    expect_true(is.finite(d$iter), label = m)
  }
})

test_that("fe converges at default settings", {
  # Pre-fix: iter = 6, delta = 2.0e-04 against tol = 1e-05, reported done.
  d <- .conv_fit(.conv_panel(), "fe")$diagnostics
  expect_true(isTRUE(d$converged))
  expect_lt(d$delta, d$tol)
})

test_that("fe agrees closely with fect::fect once it converges", {
  skip_if_not_installed("fect")
  set.seed(7)
  NU <- 60L; NT <- 12L
  ufe <- stats::rnorm(NU, 0, 1); tfe <- stats::rnorm(NT, 0, 0.3)
  g <- expand.grid(period = 1:NT, unit = 1:NU)
  g$D <- 0L
  g$D[g$unit %in% 11:30 & g$period >= 7L] <- 1L
  g$Y <- ufe[g$unit] + tfe[g$period] + stats::rnorm(nrow(g), 0, 0.4)
  g <- g[order(g$unit, g$period), c("unit", "period", "D", "Y")]
  ours <- .conv_fit(g, "fe")$results$ATE[1, 1]
  theirs <- suppressMessages(suppressWarnings(
    fect::fect(Y ~ D, data = g, index = c("unit", "period"),
               method = "fe", se = FALSE, force = "two-way")))$att.avg
  # Pre-fix this gap was 2.81e-05.
  expect_lt(abs(ours - theirs), 1e-5)
})

test_that("a fit that cannot converge warns instead of reporting a number", {
  p <- .conv_panel()
  expect_warning(
    didgpu_fect(df = p, outcome = "Y", group = "unit", time = "period",
                treatment = "D", method = "mc", max_iter = 2L,
                bootstrap_reps = 0L, verbose = FALSE),
    "did not converge")
})

test_that("the non-convergence warning does not fire on a good fit", {
  p <- .conv_panel()
  expect_no_warning(
    didgpu_fect(df = p, outcome = "Y", group = "unit", time = "period",
                treatment = "D", method = "fe", bootstrap_reps = 0L,
                verbose = FALSE))
})
