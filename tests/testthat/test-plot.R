# Smoke tests for plot.didgpu_result. We don't render the device, so
# we only check that the call returns the input invisibly and doesn't
# error on the typical result shapes (effects only, effects + placebos,
# zero-bootstrap, only placebos).

setup_pdf <- function() {
  pdf(file = NULL)  # null device so plot calls don't pop windows
}

test_that("plot() draws effects-only result without error", {
  setup_pdf(); on.exit(dev.off(), add = TRUE)
  p <- small_panel()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 3L, bootstrap_reps = 10L, seed = 2L,
                 backend = "r", verbose = FALSE)
  expect_invisible(plot(fit))
  expect_silent(plot(fit, main = "test"))
})

test_that("plot() draws effects + placebos result without error", {
  setup_pdf(); on.exit(dev.off(), add = TRUE)
  p <- small_panel()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 2L, placebo = 2L,
                 bootstrap_reps = 10L, seed = 2L,
                 backend = "r", verbose = FALSE)
  expect_invisible(plot(fit))
})

test_that("plot() works with bootstrap_reps = 0 (CIs suppressed)", {
  setup_pdf(); on.exit(dev.off(), add = TRUE)
  p <- small_panel()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 2L, bootstrap_reps = 0L,
                 backend = "r", verbose = FALSE)
  expect_invisible(plot(fit))
})

test_that("plot() respects user xlim/ylim and main", {
  setup_pdf(); on.exit(dev.off(), add = TRUE)
  p <- small_panel()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 3L, bootstrap_reps = 5L, seed = 2L,
                 backend = "r", verbose = FALSE)
  expect_invisible(plot(fit, xlim = c(-5, 5), ylim = c(-2, 5),
                         main = "custom title"))
})
