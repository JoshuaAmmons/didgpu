# Overnight stress-test driver. Runs the full test suite + the overnight
# stress harness at N = 100 per scenario, captures every result, writes
# a morning-readable report to inst/overnight_reports/.

Sys.setenv(DIDGPU_OVERNIGHT = "1")
Sys.setenv(DIDGPU_OVERNIGHT_N = "100")
Sys.setenv(DIDGPU_FUZZ_N      = "100")

library(didgpu)
library(testthat)

setwd("C:/Users/ammonsj/DID GPU/didgpu")

started_at <- Sys.time()
cat("=========================================================\n")
cat(" didgpu overnight stress run\n")
cat(sprintf(" started: %s\n", format(started_at)))
cat(sprintf(" DIDGPU_OVERNIGHT_N: %s\n", Sys.getenv("DIDGPU_OVERNIGHT_N")))
cat(sprintf(" DIDGPU_FUZZ_N:      %s\n", Sys.getenv("DIDGPU_FUZZ_N")))
cat("=========================================================\n\n")

# ---- 1. Full standard test suite ----
cat("[1/2] Full standard test suite (testthat::test_dir)\n\n")
suite_t0 <- Sys.time()
suite_res <- tryCatch(
  testthat::test_dir("tests/testthat",
                     reporter = testthat::SummaryReporter$new(),
                     stop_on_failure = FALSE),
  error = function(e) {
    cat("ERROR in suite:", conditionMessage(e), "\n")
    NULL
  })
suite_elapsed <- as.numeric(difftime(Sys.time(), suite_t0, units = "mins"))
cat(sprintf("\n[1/2] suite elapsed: %.1f min\n\n", suite_elapsed))

# ---- 2. Overnight stress tests at N=100 per scenario ----
cat("[2/2] Overnight stress tests (DIDGPU_OVERNIGHT=1)\n\n")
overnight_t0 <- Sys.time()
overnight_res <- tryCatch(
  testthat::test_file("tests/testthat/test-overnight-stress.R",
                      reporter = testthat::SummaryReporter$new(),
                      stop_on_failure = FALSE),
  error = function(e) {
    cat("ERROR in overnight:", conditionMessage(e), "\n")
    NULL
  })
overnight_elapsed <- as.numeric(difftime(Sys.time(), overnight_t0, units = "mins"))
cat(sprintf("\n[2/2] overnight elapsed: %.1f min\n\n", overnight_elapsed))

total_elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "mins"))
cat("=========================================================\n")
cat(sprintf(" TOTAL elapsed: %.1f min\n", total_elapsed))
cat(sprintf(" finished: %s\n", format(Sys.time())))
cat("=========================================================\n")
