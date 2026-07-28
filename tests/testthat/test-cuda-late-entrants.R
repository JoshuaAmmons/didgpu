# Regression test: CUDA backend must survive groups with NA baseline
# treatment (late entrants / gaps at the global first period).
#
# core_r.R sets d_sq_XX (baseline treatment) to the group's treatment at
# the GLOBAL first period; a group unobserved there gets d_sq_XX == NA.
# On the CPU path such groups fall out of every cohort mask and
# contribute zero. On the CUDA path, before the fix in cuda_glue.R, the
# NA flowed through match() into cohort_key as NA_integer_, which
# reaches the kernels as INT_MIN and triggered an illegal memory access:
#   "CUDA DID kernel failed with code 700"
# on EVERY subsequent call in the process (error 700 poisons the CUDA
# context). Any real-world unbalanced panel -- e.g. tokens/firms entering
# the sample over time -- hit this immediately.
#
# The fix parks NA-key rows in a padding cohort whose kernel contribution
# is identically zero, matching CPU semantics.

make_gappy_panel <- function(seed = 12L, nU = 30L, Tn = 60L) {
  set.seed(seed)
  unit <- rep(seq_len(nU), each = Tn)
  period <- rep(seq_len(Tn), times = nU)
  D <- integer(nU * Tn)
  for (u in 1:(nU %/% 2)) {
    on <- sort(sample(3:(Tn - 2), sample(2:5, 1)))
    for (o in on) D[unit == u & period %in% o:min(o + 1, Tn)] <- 1L
  }
  Y <- rnorm(nU)[unit] + 0.05 * period + 0.4 * D + rnorm(nU * Tn, 0, 0.5)
  d <- data.frame(unit, period, Y, D)
  # knock out the first periods for a few units => NA baseline treatment
  drop <- (d$unit %in% c(4L, 9L, 21L) & d$period <= 3L)
  d[!drop, ]
}

test_that("CUDA survives late-entrant groups (NA baseline) and matches CPU", {
  skip_if_no_cuda()
  p <- make_gappy_panel()
  f_cpu <- didgpu(df = p, outcome = "Y", group = "unit", time = "period",
                  treatment = "D", effects = 5L, placebo = 3L,
                  cluster = "unit", bootstrap_reps = 0L, seed = 1L,
                  backend = "cpu", verbose = FALSE)
  f_gpu <- didgpu(df = p, outcome = "Y", group = "unit", time = "period",
                  treatment = "D", effects = 5L, placebo = 3L,
                  cluster = "unit", bootstrap_reps = 0L, seed = 1L,
                  backend = "cuda", verbose = FALSE)
  expect_equal(f_gpu$results$Effects[, "Estimate"],
               f_cpu$results$Effects[, "Estimate"], tolerance = 1e-9)
  expect_equal(f_gpu$results$Placebos[, "Estimate"],
               f_cpu$results$Placebos[, "Estimate"], tolerance = 1e-9)
})

test_that("CUDA context is not poisoned after a gappy-panel estimation", {
  skip_if_no_cuda()
  # before the fix, the illegal access left the CUDA context unusable, so
  # even a clean balanced panel failed afterwards in the same process.
  p_gap <- make_gappy_panel()
  invisible(didgpu(df = p_gap, outcome = "Y", group = "unit",
                   time = "period", treatment = "D", effects = 3L,
                   placebo = 2L, cluster = "unit", bootstrap_reps = 0L,
                   seed = 1L, backend = "cuda", verbose = FALSE))
  set.seed(3)
  nU <- 20L; Tn <- 20L
  unit <- rep(seq_len(nU), each = Tn); period <- rep(seq_len(Tn), nU)
  D <- as.integer(unit <= 8L & period >= 10L)
  Y <- rnorm(nU)[unit] + 0.4 * D + rnorm(nU * Tn, 0, 0.5)
  f <- didgpu(df = data.frame(unit, period, Y, D), outcome = "Y",
              group = "unit", time = "period", treatment = "D",
              effects = 3L, placebo = 2L, cluster = "unit",
              bootstrap_reps = 0L, seed = 1L, backend = "cuda",
              verbose = FALSE)
  expect_true(is.finite(f$results$Effects[1, "Estimate"]))
})
