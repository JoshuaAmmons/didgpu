# The ATE must not depend on which backend computed it.
#
# .backend_cuda() had
#     ate <- if (h$l_eff == 1L) ce$effects[1] else NA_real_
# so backend "cuda" returned NA for the ATE whenever effects > 1, while
# the per-horizon effects themselves were correct. backend = "auto"
# resolves to "cuda" on any machine with a GPU, so an ordinary
# didgpu(effects = 3) call returned an NA ATE there and a number on CPU.
#
# The weights it needed (ce$n_inc) were already being computed and
# returned as n_inc_effects; the CUDA branch simply never used them.

.ate_panel <- function(seed = 11L) {
  as.data.frame(didgpu_simulate_panel(n_units = 120L, n_periods = 12L,
                                      seed = seed))
}
.ate_fit <- function(p, bk, eff) {
  didgpu(df = p, outcome = "Y", group = "unit", time = "period",
         treatment = "D", effects = eff, placebo = 0L, bootstrap_reps = 0L,
         backend = bk, verbose = FALSE)
}
.backends <- function() {
  bk <- c("r", "cpu")
  if (isTRUE(tryCatch(didgpu_has_cuda_support(), error = function(e) FALSE)))
    bk <- c(bk, "cuda", "auto")
  bk
}

test_that("the ATE is finite and backend-invariant for effects > 1", {
  p <- .ate_panel()
  for (eff in c(1L, 3L, 5L)) {
    ref <- .ate_fit(p, "r", eff)$results$ATE[1L, 1L]
    expect_true(is.finite(ref), label = paste("r backend, effects =", eff))
    for (bk in .backends()) {
      got <- .ate_fit(p, bk, eff)$results$ATE[1L, 1L]
      expect_true(is.finite(got),
                  label = sprintf("%s backend, effects = %d", bk, eff))
      expect_equal(got, ref, tolerance = 1e-12,
                   label = sprintf("%s backend, effects = %d", bk, eff))
    }
  }
})

test_that("on a binary absorbing panel the ATE is the switcher-weighted mean", {
  # This identity is a PROPERTY OF THIS PANEL, not the definition of the
  # ATE: didgpu_simulate_panel() is binary and absorbing, so every
  # switcher's dose change is exactly 1 and Av_tot_eff's denominator
  # sum_k N_k * delta_k collapses to sum_k N_k. See
  # test-ate-estimand.R for the general case, which this used to get
  # wrong by assuming delta_k == 1 always.
  p <- .ate_panel()
  expect_true(all(p$D %in% c(0, 1)))
  for (bk in .backends()) {
    f <- .ate_fit(p, bk, 4L)
    e <- as.numeric(f$results$Effects[, "Estimate"])
    n <- as.numeric(f$results$Effects[, "Switchers"])
    ok <- is.finite(e) & n > 0
    expect_equal(f$results$ATE[1L, 1L],
                 sum(e[ok] * n[ok]) / sum(n[ok]),
                 tolerance = 1e-12, label = paste("backend", bk))
  }
})
