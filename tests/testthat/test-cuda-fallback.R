# A CUDA kernel that fails on the device -- out of memory when other jobs
# share the GPU ("CUDA DID kernel failed with code 2") -- used to fail the
# whole fit. The backend now computes that fit on the CPU and says so.

test_that("a failing CUDA DID kernel falls back to the CPU kernel", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(), error = function(e) FALSE)),
              "didgpu built without CUDA")
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                             tau_profile = c(0.5, 1.0, 1.2), seed = 3L)
  p$D <- as.integer(p$D >= 0.5)
  cpu <- didgpu(p, "Y", "unit", "period", "D", effects = 3L, placebo = 1L,
                backend = "cpu", graph_off = TRUE)
  local_mocked_bindings(
    .cuda_one_event_time = function(...) stop("CUDA DID kernel failed with code 2"),
    .package = "didgpu")
  expect_message(
    gpu <- didgpu(p, "Y", "unit", "period", "D", effects = 3L, placebo = 1L,
                  backend = "cuda", graph_off = TRUE),
    "out of memory.*CPU backend instead")
  expect_equal(gpu$results$Effects, cpu$results$Effects, tolerance = 1e-12)
  expect_equal(gpu$results$ATE, cpu$results$ATE, tolerance = 1e-12)
})

test_that("other CUDA backend errors still stop the fit", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(), error = function(e) FALSE)),
              "didgpu built without CUDA")
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 8L,
                             tau_profile = c(0.5, 1.0), seed = 3L)
  p$D <- as.integer(p$D >= 0.5)
  local_mocked_bindings(
    .cuda_one_event_time = function(...) stop("something else went wrong"),
    .package = "didgpu")
  expect_error(didgpu(p, "Y", "unit", "period", "D", effects = 2L,
                      backend = "cuda", graph_off = TRUE),
               "something else went wrong")
})
