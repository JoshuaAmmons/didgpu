# didgpu_by(): subgroup wrapper.

test_that("didgpu_by splits panel by levels and returns per-level results", {
  p <- small_panel()
  # Add a grouping column with two levels.
  p$region <- ifelse(p$unit %% 2L == 0L, "north", "south")
  fit_by <- didgpu_by(p, "region",
                      outcome = "Y", group = "unit",
                      time = "period", treatment = "D",
                      effects = 2L, bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE)
  expect_s3_class(fit_by, "didgpu_by_result")
  expect_named(fit_by, c("north", "south"))
  expect_s3_class(fit_by[["north"]], "didgpu_result")
  expect_s3_class(fit_by[["south"]], "didgpu_result")
  # Each subgroup's estimates differ from the other (with very high
  # probability on a noisy panel).
  e_n <- as.numeric(fit_by[["north"]]$results$Effects[, "Estimate"])
  e_s <- as.numeric(fit_by[["south"]]$results$Effects[, "Estimate"])
  expect_false(all(is.na(e_n) | is.na(e_s) | e_n == e_s))
})

test_that("didgpu_by reconciles with full-panel didgpu for a single-level by", {
  # Trivial case: by_var with only one level should give a one-element
  # list whose single fit equals the full-panel fit.
  p <- small_panel()
  p$all <- "x"
  expect_warning(fit_by <- didgpu_by(p, "all",
                                      outcome = "Y", group = "unit",
                                      time = "period", treatment = "D",
                                      effects = 2L, bootstrap_reps = 0L,
                                      backend = "r", verbose = FALSE),
                  "only 1 distinct level")
  fit_full <- didgpu(p, "Y", "unit", "period", "D",
                      effects = 2L, bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE)
  expect_equal(length(fit_by), 1L)
  expect_equal(as.numeric(fit_by[["x"]]$results$Effects[, "Estimate"]),
               as.numeric(fit_full$results$Effects[, "Estimate"]))
})

test_that("didgpu_by writes per-subgroup checkpoint subdirectories", {
  p <- small_panel()
  p$region <- ifelse(p$unit %% 2L == 0L, "north", "south")
  cdir <- tempfile("didgpu_by_")
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)
  fit_by <- didgpu_by(p, "region",
                      outcome = "Y", group = "unit",
                      time = "period", treatment = "D",
                      effects = 2L, bootstrap_reps = 2L, seed = 1L,
                      checkpoint_dir = cdir,
                      backend = "r", verbose = FALSE)
  expect_true(dir.exists(file.path(cdir, "north")))
  expect_true(dir.exists(file.path(cdir, "south")))
  expect_true(file.exists(file.path(cdir, "north", "meta.json")))
  expect_true(file.exists(file.path(cdir, "south", "manifest.csv")))
})

test_that("print.didgpu_by_result runs without error", {
  p <- small_panel()
  p$region <- ifelse(p$unit %% 2L == 0L, "n", "s")
  fit_by <- didgpu_by(p, "region",
                      outcome = "Y", group = "unit",
                      time = "period", treatment = "D",
                      effects = 2L, bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE)
  expect_output(print(fit_by), "didgpu_by result")
  expect_output(print(fit_by), "region = n")
  expect_output(print(fit_by), "region = s")
})
