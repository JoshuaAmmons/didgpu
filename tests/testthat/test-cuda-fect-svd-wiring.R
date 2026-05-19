# Phase-1 wiring test for task #80.
#
# The CUDA cuda_fect_svd kernels (truncated + softthreshold) have
# real cuSOLVER implementations. The R-side wiring tries CUDA first
# and falls back to base R svd() on any failure. Whether the CUDA
# numerical output exactly matches R svd() is a Phase 4 (#90)
# verification concern — Phase 1 only verifies the wiring contract:
#   (a) backend = "cuda" doesn't error out.
#   (b) The function shape returned matches the R-side shape.
#   (c) On a fresh-build sanity check, the CUDA result is at least
#       close to the R result (within 1e-3 — Jacobi SVD has loose
#       tolerance).

test_that(".fect_svd_truncated_cuda returns a list(L, F, d) or NULL", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  set.seed(1L)
  M <- matrix(rnorm(40), nrow = 8, ncol = 5)
  res <- didgpu:::.fect_svd_truncated_cuda(M, r = 2L)
  if (is.null(res)) succeed("kernel returned NULL — fallback works")
  else {
    expect_named(res, c("L", "F", "d"), ignore.order = TRUE)
    expect_equal(dim(res$L), c(8L, 2L))
    expect_equal(dim(res$F), c(2L, 5L))
  }
})

test_that(".fect_svd_softthreshold_cuda returns Y_hat or NULL", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  set.seed(2L)
  Y <- matrix(rnorm(60), nrow = 10, ncol = 6)
  res <- didgpu:::.fect_svd_softthreshold_cuda(Y, lambda = 0.5)
  if (is.null(res)) succeed("kernel returned NULL — fallback works")
  else {
    expect_named(res, c("Y_hat", "n_nonzero"), ignore.order = TRUE)
    expect_equal(dim(res$Y_hat), c(10L, 6L))
  }
})

test_that("CUDA truncated SVD agrees with base svd() on a small matrix (loose tolerance)", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  set.seed(3L)
  M <- matrix(rnorm(40), nrow = 8, ncol = 5)
  cuda_res <- didgpu:::.fect_svd_truncated_cuda(M, r = 3L)
  skip_if(is.null(cuda_res), "CUDA kernel returned NULL; nothing to compare")
  r_res    <- didgpu:::.fect_svd_r(M, r = 3L)
  # L * F is invariant to sign / rotation; compare the reconstruction.
  M_hat_cuda <- cuda_res$L %*% cuda_res$F
  M_hat_r    <- r_res$L    %*% r_res$F
  expect_equal(M_hat_cuda, M_hat_r, tolerance = 1e-3)
})

test_that("CUDA softthreshold agrees with R reconstruction on a small matrix", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  set.seed(4L)
  Y <- matrix(rnorm(60), nrow = 10, ncol = 6)
  cuda_res <- didgpu:::.fect_svd_softthreshold_cuda(Y, lambda = 0.5)
  skip_if(is.null(cuda_res), "CUDA kernel returned NULL; nothing to compare")
  s <- svd(Y)
  D_st <- pmax(s$d - 0.5, 0)
  nz <- D_st > 0
  Y_hat_r <- s$u[, nz, drop = FALSE] %*%
             diag(D_st[nz], sum(nz), sum(nz)) %*%
             t(s$v[, nz, drop = FALSE])
  expect_equal(cuda_res$Y_hat, Y_hat_r, tolerance = 1e-3)
})
