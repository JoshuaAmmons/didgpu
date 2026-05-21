# TestMechs ARP and FSST methods (alternatives to CS).

skip_if_no_qp <- function() {
  testthat::skip_if_not_installed("quadprog")
}

build_tm_panel <- function(n = 600L, K = 2L, seed = 1L) {
  set.seed(seed)
  D <- sample(c(0L, 1L), n, replace = TRUE)
  M <- sample(seq_len(K), n, replace = TRUE)
  Y <- rnorm(n) + 0.3 * (M == 2L) + 0.1 * D
  data.frame(D = D, M = M, Y = Y, stringsAsFactors = FALSE)
}

test_that("ARP method runs end-to-end on binary-M case", {
  skip_if_no_qp()
  df <- build_tm_panel(n = 800L, K = 2L, seed = 1L)
  res <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                  method = "ARP", B = 30L,
                                  num_Ybins = 3L, seed = 1L)
  expect_equal(res$method, "ARP")
  expect_true(is.finite(res$test_stat))
  expect_true(is.finite(res$cv) || is.na(res$cv))
})

test_that("FSST method runs end-to-end on binary-M case", {
  skip_if_no_qp()
  df <- build_tm_panel(n = 800L, K = 2L, seed = 1L)
  res <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                  method = "FSST", B = 30L,
                                  num_Ybins = 3L, seed = 1L)
  expect_equal(res$method, "FSST")
  expect_true(is.finite(res$test_stat))
  expect_true(is.finite(res$cv) || is.na(res$cv))
  expect_true(is.finite(res$pval) || is.na(res$pval))
})

test_that("All three TestMechs methods produce results on the same panel", {
  skip_if_no_qp()
  df <- build_tm_panel(n = 1000L, K = 2L, seed = 23L)
  results <- list()
  for (method in c("CS", "ARP", "FSST")) {
    results[[method]] <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                                  method = method,
                                                  B = 30L, num_Ybins = 3L,
                                                  seed = 1L)
  }
  for (method in c("CS", "ARP", "FSST")) {
    expect_true(is.finite(results[[method]]$test_stat))
  }
})

test_that("Multi-level M (K = 3) works with CS method", {
  skip_if_no_qp()
  df <- build_tm_panel(n = 1500L, K = 3L, seed = 1L)
  res <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                  method = "CS", B = 30L,
                                  num_Ybins = 3L, seed = 1L)
  expect_equal(res$K, 3L)
  expect_true(is.finite(res$test_stat))
})

test_that("Multi-level M (K = 3) works with ARP method", {
  skip_if_no_qp()
  df <- build_tm_panel(n = 1500L, K = 3L, seed = 1L)
  res <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                  method = "ARP", B = 30L,
                                  num_Ybins = 3L, seed = 1L)
  expect_equal(res$K, 3L)
  expect_true(is.finite(res$test_stat))
})

test_that("Multi-level M (K = 4) works with all three methods", {
  skip_if_no_qp()
  df <- build_tm_panel(n = 2000L, K = 4L, seed = 7L)
  for (method in c("CS", "ARP", "FSST")) {
    res <- didgpu_test_sharp_null(df, "D", "M", "Y",
                                    method = method, B = 20L,
                                    num_Ybins = 3L, seed = 1L)
    expect_equal(res$K, 4L)
  }
})
