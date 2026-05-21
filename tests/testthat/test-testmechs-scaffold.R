# TestMechs sharp mediation tests — scaffold phase. Verifies the
# public API exists, validates args, and the partial-density helper
# works on toy data.

test_that("didgpu_test_sharp_null exists and runs for every method", {
  testthat::skip_if_not_installed("quadprog")
  set.seed(1L)
  df <- data.frame(
    D = sample(c(0L, 1L), 200L, replace = TRUE),
    M = sample(seq_len(3L), 200L, replace = TRUE),
    Y = rnorm(200L)
  )
  for (method in c("CS", "ARP", "FSST")) {
    res <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                    method = method, B = 30L, seed = 1L)
    expect_equal(res$method, method)
    expect_true(is.finite(res$test_stat))
  }
})

test_that("didgpu_test_sharp_null validates required args", {
  df <- data.frame(D = 0:1, M = 1:2, Y = c(1.0, 2.0))
  expect_error(
    didgpu_test_sharp_null(df, "Dnope", "M", "Y"),
    "column not in df: Dnope"
  )
  expect_error(
    didgpu_test_sharp_null(df, "D", "M", "Y", method = "BAD"),
    "should be one of"
  )
})

test_that("didgpu_test_sharp_null rejects non-binary D", {
  df <- data.frame(D = c(0L, 1L, 2L, 0L),
                    M = c(1L, 1L, 2L, 2L),
                    Y = c(1.0, 2.0, 3.0, 4.0))
  expect_error(
    didgpu_test_sharp_null(df, "D", "M", "Y"),
    "in.*0.*1"
  )
})

test_that("didgpu_lb_frac_affected exists and stubs out cleanly", {
  set.seed(1L)
  df <- data.frame(
    D = sample(c(0L, 1L), 100L, replace = TRUE),
    M = sample(seq_len(3L), 100L, replace = TRUE),
    Y = rnorm(100L))
  expect_error(
    didgpu_lb_frac_affected(df, "D", "M", "Y", B = 20L, seed = 1L),
    "not yet implemented"
  )
})

test_that(".testmechs_bin_y bins continuous Y into the requested quantile bins", {
  set.seed(1L)
  y <- rnorm(500L)
  bins <- didgpu:::.testmechs_bin_y(y, n_bins = 5L)
  expect_true(all(bins %in% 1L:5L))
  # Roughly equal-frequency bins.
  tab <- table(bins)
  expect_true(min(tab) > 50L)  # at least 50 per bin out of 500
})

test_that(".testmechs_bin_y handles already-discrete Y without changing semantics", {
  y <- c(1L, 1L, 2L, 2L, 3L, 3L)
  bins <- didgpu:::.testmechs_bin_y(y, n_bins = 5L)
  expect_equal(bins, c(1L, 1L, 2L, 2L, 3L, 3L))
})

test_that(".testmechs_partial_density computes a correctly-sized beta vector", {
  set.seed(1L)
  n <- 1000L
  d <- sample(c(0L, 1L), n, replace = TRUE)
  m <- sample(seq_len(3L), n, replace = TRUE)
  y <- sample(seq_len(4L), n, replace = TRUE)
  out <- didgpu:::.testmechs_partial_density(d, m, y)
  # Length = 2 * K * d_y = 2 * 3 * 4 = 24.
  expect_equal(length(out$beta), 24L)
  expect_equal(out$K, 3L)
  expect_equal(out$d_y, 4L)
  # Each per-d block sums to 1 (it's a joint distribution over Y, M
  # conditional on D).
  block_d0 <- out$beta[1:(3L * 4L)]
  block_d1 <- out$beta[(3L * 4L + 1L):24L]
  expect_lt(abs(sum(block_d0) - 1), 1e-10)
  expect_lt(abs(sum(block_d1) - 1), 1e-10)
})
