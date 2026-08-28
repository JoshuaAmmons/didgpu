# didgpu_fect(method = "ife", r = 0) must run.
#
# An interactive-fixed-effects model with zero factors IS the plain
# two-way FE model, and sweeping r = 0..k is the standard way to ask
# whether a finding depends on the factor structure at all. Previously
# r = 0 errored: svd(M, nu = 0, nv = 0) omits the `u` component entirely
# (it is present only when nu > 0), so .fect_svd_r()'s
# `s$u %*% diag(...)` failed with "requires numeric/complex matrix/vector
# arguments".

.p10 <- function() {
  as.data.frame(didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                                      seed = 11L))
}
.ate <- function(method, r = 2L) {
  didgpu_fect(df = .p10(), outcome = "Y", group = "unit", time = "period",
              treatment = "D", method = method, r = r,
              bootstrap_reps = 0L, verbose = FALSE)$results$ATE[1L, 1L]
}

test_that("ife with r = 0 runs", {
  expect_no_error(.ate("ife", 0L))
  expect_true(is.finite(.ate("ife", 0L)))
})

test_that("ife with r = 0 reduces exactly to two-way FE", {
  # This is the semantic content of r = 0: no factor term at all.
  expect_equal(.ate("ife", 0L), .ate("fe"), tolerance = 1e-12)
})

test_that("the r sweep used for factor-structure diagnostics runs end to end", {
  for (r in 0:3) expect_true(is.finite(.ate("ife", as.integer(r))),
                             label = paste("r =", r))
})
