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
.ate <- function(method, r = 2L, ...) {
  suppressWarnings(
    didgpu_fect(df = .p10(), outcome = "Y", group = "unit", time = "period",
                treatment = "D", method = method, r = r,
                bootstrap_reps = 0L, verbose = FALSE, ...))$results$ATE[1L, 1L]
}

test_that("ife with r = 0 runs", {
  expect_no_error(.ate("ife", 0L))
  expect_true(is.finite(.ate("ife", 0L)))
})

test_that("ife with r = 0 reduces exactly to two-way FE at matched min_T0", {
  # r = 0 means no factor term, so the ESTIMATOR collapses to two-way FE.
  # The two calls do not agree by default, though, and should not: fect
  # sets min.T0 = 1 for method "fe" but 5 for "ife", so they fit
  # different samples. Verified against the reference on a 60-unit panel
  # where 22 units have fewer than 5 untreated periods:
  #     fect    fe +0.302820   ife(r = 0) +0.236833
  #     didgpu  fe +0.302820   ife(r = 0) +0.236833
  # Both packages separate them, and didgpu reproduces both numbers. Pin
  # the collapse at a matched min_T0, and pin the divergence at default.
  expect_equal(.ate("ife", 0L, min_T0 = 1L), .ate("fe"), tolerance = 1e-12)
})

test_that("ife and fe differ at their default min_T0, as in fect", {
  p <- .p10()
  npre <- tapply(p$D, p$unit, function(x) sum(x == 0))
  skip_if(all(npre >= 5), "no short-history units on this panel")
  expect_false(isTRUE(all.equal(.ate("ife", 0L), .ate("fe"),
                                tolerance = 1e-10)))
})

test_that("the r sweep used for factor-structure diagnostics runs end to end", {
  for (r in 0:3) expect_true(is.finite(.ate("ife", as.integer(r))),
                             label = paste("r =", r))
})
