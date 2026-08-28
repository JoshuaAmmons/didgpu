# didgpu_tidy() must accept a didgpu_cs_result.
#
# didgpu_cs_result is a first-class result type elsewhere in the package
# (didgpu_loo() and didgpu_honest_did() both document and accept it), but
# didgpu_tidy() guarded with stopifnot(inherits(x, "didgpu_result")) and
# the CS class does not carry that parent. Callers got the opaque
# `inherits(x, "didgpu_result") is not TRUE` rather than either a tidy
# frame or a usable message.

.cs_fit <- function(reps = 0L) {
  p <- as.data.frame(didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                                           seed = 11L))
  didgpu_cs(df = p, outcome = "Y", group = "unit", time = "period",
            treatment = "D", bootstrap_reps = reps, backend = "r",
            verbose = FALSE)
}

test_that("didgpu_tidy accepts a didgpu_cs_result", {
  t <- didgpu_tidy(.cs_fit())
  expect_s3_class(t, "data.frame")
  expect_gt(nrow(t), 0L)
  expect_true(all(c("term", "estimate", "std.error", "statistic",
                    "p.value", "conf.low", "conf.high", "kind") %in% names(t)))
  expect_setequal(unique(t$kind), c("aggregate", "att_gt"))
})

test_that("tidy works with bootstrap_reps = 0, where CIs are absent", {
  # att_gt only gains ci_low / ci_high on the bootstrap path; indexing the
  # missing columns previously produced a rows-mismatch in data.frame().
  t <- didgpu_tidy(.cs_fit(0L))
  cells <- t[t$kind == "att_gt", ]
  expect_gt(nrow(cells), 0L)
  expect_true(all(is.na(cells$conf.low)))
})

test_that("bootstrap SEs and CIs reach the tidy frame", {
  t <- didgpu_tidy(.cs_fit(40L))
  cells <- t[t$kind == "att_gt", ]
  expect_true(any(is.finite(cells$std.error)))
  expect_true(any(is.finite(cells$conf.low)))
})

test_that("conf.int = FALSE drops the CI columns", {
  t <- didgpu_tidy(.cs_fit(), conf.int = FALSE)
  expect_false(any(c("conf.low", "conf.high") %in% names(t)))
})

test_that("an unsupported object gets a message naming both classes", {
  expect_error(didgpu_tidy(list(1)), "didgpu_cs_result")
})
