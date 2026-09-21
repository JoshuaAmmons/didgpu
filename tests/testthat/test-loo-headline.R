# didgpu_loo() must report the OVERALL ATT, never the first row of the
# requested aggregation.
#
# Regression. .loo_extract_headline() took fit$aggregation$estimate[1] for
# a didgpu_cs_result, with a comment claiming that "works for all four
# aggregations". It does not. For the default aggregation = "event",
# row 1 is the MOST NEGATIVE event time -- the longest pre-treatment
# horizon -- so LOO reported a pre-treatment placebo as its headline.
#
# It also hid itself well: dropping a single cohort seldom changes which
# cells populate the earliest lead, so nearly every entity came back with
# an IDENTICAL estimate. That reads as "no entity is influential" rather
# than as a bug. On a 56-year panel with a 1985 cohort it returned the
# e = -42 cell (-0.006317) instead of the ATT (-0.032831). Affected both
# by = "cohort" and by = "unit".

.loo_panel <- function(seed = 11L) {
  as.data.frame(didgpu_simulate_panel(n_units = 120L, n_periods = 12L,
                                      frac_treated = 0.6, seed = seed))
}

.loo_fit <- function(p, agg = "event") {
  didgpu_cs(df = p, outcome = "Y", group = "unit", time = "period",
            treatment = "D", control_group = "never", est_method = "DR",
            aggregation = agg, bootstrap_reps = 0L, backend = "r",
            verbose = FALSE)
}

test_that("the LOO headline is the overall ATT, not the first aggregation row", {
  p <- .loo_panel()
  f <- .loo_fit(p, "event")
  a <- as.data.frame(f$aggregation)
  first_row <- a$estimate[1]
  overall <- didgpu:::.cs_aggregate(f$att_gt, "overall", f$args)$estimate[1]
  # The trap only exists when these differ, which they do for "event".
  expect_gt(abs(first_row - overall), 1e-6)
  expect_lt(a$event_time[1], 0)   # row 1 really is pre-treatment

  lo <- didgpu_loo(f, by = "cohort", df = p, verbose = FALSE)
  expect_gt(nrow(lo), 1L)
  expect_equal(sum(abs(lo$estimate - first_row) < 1e-10), 0L)
})

test_that("LOO estimates vary across entities for cohort and unit", {
  # The old behaviour returned a near-constant value; genuine LOO must move.
  p <- .loo_panel()
  f <- .loo_fit(p, "event")
  for (byv in c("cohort", "unit")) {
    lo <- didgpu_loo(f, by = byv, df = p, verbose = FALSE)
    est <- lo$estimate[is.finite(lo$estimate)]
    expect_gt(length(est), 1L)
    expect_gt(stats::sd(est), 0)
    expect_gt(length(unique(round(est, 10))), 1L)
  }
})

test_that("the headline is invariant to the fit's aggregation scheme", {
  # Whether the user asked for "event" or "overall", LOO answers the same
  # question: how does the overall ATT move when this entity is dropped.
  p <- .loo_panel()
  le <- didgpu_loo(.loo_fit(p, "event"),   by = "cohort", df = p, verbose = FALSE)
  lo <- didgpu_loo(.loo_fit(p, "overall"), by = "cohort", df = p, verbose = FALSE)
  le <- le[order(le$leave_out), ]; lo <- lo[order(lo$leave_out), ]
  expect_equal(le$estimate, lo$estimate, tolerance = 1e-10)
})

test_that("delta is measured against the overall ATT", {
  p <- .loo_panel()
  f <- .loo_fit(p, "event")
  overall <- didgpu:::.cs_aggregate(f$att_gt, "overall", f$args)$estimate[1]
  lo <- didgpu_loo(f, by = "cohort", df = p, verbose = FALSE)
  ok <- is.finite(lo$estimate) & is.finite(lo$delta)
  expect_equal(lo$estimate[ok] - lo$delta[ok], rep(overall, sum(ok)),
               tolerance = 1e-10)
})
