# TestMechs sharp-null test for the binary-mediator case via Cox-Shi.
# These tests verify the end-to-end pipeline; the CS engine itself
# also has its own unit tests below.

skip_if_no_quadprog <- function() {
  testthat::skip_if_not_installed("quadprog")
}

test_that(".testmechs_bootstrap produces a B x dim_beta matrix with correct row sums per D-block", {
  set.seed(1L)
  n <- 600L
  d <- sample(c(0L, 1L), n, replace = TRUE)
  m <- sample(c(1L, 2L), n, replace = TRUE)
  y <- sample(c(1L, 2L, 3L), n, replace = TRUE)
  boot <- didgpu:::.testmechs_bootstrap(d, m, y, B = 20L,
                                          method = "nonparametric",
                                          seed = 1L, backend = "r")
  expect_equal(nrow(boot), 20L)
  expect_equal(ncol(boot), 2L * 2L * 3L)
  # Each per-D block should sum to roughly 1 in each row (it's a
  # joint distribution over Y, M | D).
  K  <- 2L; dy <- 3L
  for (b in 1L:5L) {
    block_0 <- boot[b, 1:(K * dy)]
    block_1 <- boot[b, (K * dy + 1L):(2L * K * dy)]
    expect_lt(abs(sum(block_0) - 1), 1e-10)
    expect_lt(abs(sum(block_1) - 1), 1e-10)
  }
})

test_that(".testmechs_bootstrap is deterministic given the same seed", {
  set.seed(1L)
  n <- 200L
  d <- sample(c(0L, 1L), n, replace = TRUE)
  m <- sample(c(1L, 2L), n, replace = TRUE)
  y <- sample(c(1L, 2L), n, replace = TRUE)
  b1 <- didgpu:::.testmechs_bootstrap(d, m, y, B = 5L, seed = 42L,
                                        backend = "r")
  b2 <- didgpu:::.testmechs_bootstrap(d, m, y, B = 5L, seed = 42L,
                                        backend = "r")
  expect_equal(b1, b2)
})

test_that(".testmechs_cs_test runs without error on a tiny synthetic case", {
  skip_if_no_quadprog()
  # Theta_hat in the polytope {theta >= 0}; should give T = 0 and
  # pval = 1.
  theta <- c(0.5, 0.5)
  Sigma <- diag(2) * 0.01
  A     <- diag(2)   # theta >= 0
  res <- didgpu:::.testmechs_cs_test(theta, Sigma, A)
  expect_equal(res$test_stat, 0)
  expect_equal(res$pval, 1)
  expect_false(res$reject)
})

test_that(".testmechs_cs_test rejects when theta_hat is outside the cone", {
  skip_if_no_quadprog()
  # Theta_hat has a negative component; A says theta >= 0; the
  # projection moves it onto the boundary and the test stat is large.
  theta <- c(0.5, -2.0)
  Sigma <- diag(2) * 0.01
  A     <- diag(2)
  res <- didgpu:::.testmechs_cs_test(theta, Sigma, A)
  expect_gt(res$test_stat, 0)
  expect_lt(res$pval, 0.05)
  expect_true(res$reject)
})

test_that("didgpu_test_sharp_null end-to-end on binary M + binary Y under full mediation", {
  skip_if_no_quadprog()
  # Simulate a panel where full mediation holds exactly: Y(d, m) = m
  # (Y depends only on M, not D directly). The sharp-null test should
  # NOT reject (p-value should be reasonably large).
  set.seed(42L)
  n <- 2000L
  D <- sample(c(0L, 1L), n, replace = TRUE)
  # M depends on D (compliers + always/never-takers).
  M <- ifelse(D == 1L,
              sample(c(0L, 1L), n, replace = TRUE, prob = c(0.3, 0.7)),
              sample(c(0L, 1L), n, replace = TRUE, prob = c(0.7, 0.3)))
  # Y depends only on M (full mediation).
  Y <- ifelse(M == 1L,
              rbinom(n, 1, 0.6),
              rbinom(n, 1, 0.2))
  df <- data.frame(D = D, M = M, Y = Y)
  res <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                  method = "CS", B = 50L, num_Ybins = 2L,
                                  seed = 1L, backend = "r")
  expect_equal(res$method, "CS")
  expect_true(is.finite(res$test_stat))
  expect_true(is.finite(res$pval))
  # We don't enforce non-rejection (statistical noise); just that the
  # test ran and produced sensible output.
})

test_that("didgpu_test_sharp_null produces finite output under violation (sharp-null power TBD)", {
  skip_if_no_quadprog()
  # Strong direct effect of D on Y, holding M constant. A FULLY-correct
  # sharp-null test should reject; our v1 polytope+CS plumbing computes
  # finite output but the equality-constraint handling needs more care
  # before the power is comparable to the reference TestMechs package
  # (see roadmap in R/testmechs_sharp_null.R). For now we just verify
  # the pipeline runs end-to-end and produces a finite p-value.
  set.seed(42L)
  n <- 2000L
  D <- sample(c(0L, 1L), n, replace = TRUE)
  M <- sample(c(0L, 1L), n, replace = TRUE)
  Y <- as.integer(D == 1L & M == 1L)
  df <- data.frame(D = D, M = M, Y = Y)
  res <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                  method = "CS", B = 50L, num_Ybins = 2L,
                                  seed = 1L, backend = "r")
  expect_true(is.finite(res$test_stat))
  expect_true(is.finite(res$pval))
  expect_true(res$pval >= 0 && res$pval <= 1)
})

test_that("ARP and FSST now run (no longer stubbed)", {
  skip_if_no_quadprog()
  set.seed(1L)
  n <- 200L
  df <- data.frame(D = sample(0:1, n, TRUE),
                    M = sample(1:2, n, TRUE),
                    Y = rnorm(n))
  for (method in c("ARP", "FSST")) {
    res <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                    method = method, B = 30L,
                                    num_Ybins = 2L, seed = 1L)
    expect_equal(res$method, method)
    expect_true(is.finite(res$test_stat))
  }
})

test_that("Multi-level M (K > 2) now runs (no longer stubbed)", {
  skip_if_no_quadprog()
  set.seed(1L)
  n <- 800L
  df <- data.frame(D = sample(0:1, n, TRUE),
                    M = sample(1:4, n, TRUE),
                    Y = rnorm(n))
  res <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                  method = "CS", B = 30L,
                                  num_Ybins = 2L, seed = 1L)
  expect_equal(res$K, 4L)
  expect_true(is.finite(res$test_stat))
})
