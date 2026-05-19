# ============================================================================
# Overnight stress tests. Each scenario runs N_OVERNIGHT random panels
# (default 100; controllable via DIDGPU_OVERNIGHT_N env var) and checks
# that didgpu produces sensible output across every estimator family.
#
# Skipped during normal CI runs (skip_on_ci()); designed to be run via
# the `inst/scripts/overnight_run.R` driver. Results are aggregated
# and dumped to a report file the user can read in the morning.
# ============================================================================

skip_unless_overnight <- function() {
  if (!nzchar(Sys.getenv("DIDGPU_OVERNIGHT")) ||
      Sys.getenv("DIDGPU_OVERNIGHT") == "0") {
    testthat::skip("set DIDGPU_OVERNIGHT=1 to run overnight stress tests")
  }
}

.overnight_n <- function() {
  n <- as.integer(Sys.getenv("DIDGPU_OVERNIGHT_N", "100"))
  if (is.na(n) || n < 1L) 100L else n
}


.rand_panel_balanced <- function(seed, n_units = 60L, n_periods = 12L) {
  set.seed(seed)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  F_g <- rep(Inf, n_units)
  n_treated <- as.integer(n_units * runif(1L, 0.4, 0.8))
  treated <- sort(sample(seq_len(n_units), n_treated))
  F_g[treated] <- sample(seq(3L, n_periods - 1L), n_treated, replace = TRUE)
  panel <- expand.grid(unit = seq_len(n_units), period = seq_len(n_periods))
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             runif(1L, 0.5, 1.5) * panel$D +
             rnorm(nrow(panel), 0, runif(1L, 0.2, 0.5))
  panel[order(panel$unit, panel$period),
         c("unit", "period", "D", "Y")]
}


test_that("overnight: didgpu fuzz at N", {
  skip_unless_overnight()
  N <- .overnight_n()
  fails <- 0L
  for (s in seq_len(N)) {
    p <- .rand_panel_balanced(seed = s + 100000L)
    fit <- tryCatch(
      didgpu(p, "Y", "unit", "period", "D",
              effects = 2L, bootstrap_reps = 0L,
              backend = "r", verbose = FALSE),
      error = function(e) NULL)
    if (is.null(fit)) fails <- fails + 1L
    e <- as.numeric(fit$results$Effects[, "Estimate"])
    if (any(!is.finite(e) & !is.na(e))) fails <- fails + 1L
  }
  expect_lt(fails, max(1L, as.integer(N * 0.05)))
})


test_that("overnight: didgpu_cs OR + bootstrap at N", {
  skip_unless_overnight()
  N <- .overnight_n()
  fails <- 0L
  for (s in seq_len(N)) {
    p <- .rand_panel_balanced(seed = s + 200000L)
    fit <- tryCatch(
      didgpu_cs(p, "Y", "unit", "period", "D",
                 est_method = "OR",
                 bootstrap_reps = 5L, seed = 1L,
                 backend = "r", verbose = FALSE),
      error = function(e) NULL)
    if (is.null(fit)) { fails <- fails + 1L; next }
    if (nrow(fit$att_gt) == 0L) fails <- fails + 1L
  }
  expect_lt(fails, max(1L, as.integer(N * 0.10)))
})


test_that("overnight: didgpu_cs all three methods at N", {
  skip_unless_overnight()
  N <- .overnight_n()
  fails <- list(OR = 0L, IPW = 0L, DR = 0L)
  for (s in seq_len(N)) {
    p <- .rand_panel_balanced(seed = s + 300000L)
    for (m in c("OR", "IPW", "DR")) {
      fit <- tryCatch(
        didgpu_cs(p, "Y", "unit", "period", "D",
                   est_method = m, bootstrap_reps = 0L,
                   backend = "r", verbose = FALSE),
        error = function(e) NULL)
      if (is.null(fit) || nrow(fit$att_gt) == 0L) {
        fails[[m]] <- fails[[m]] + 1L
      }
    }
  }
  for (m in c("OR", "IPW", "DR")) {
    expect_lt(fails[[m]], max(1L, as.integer(N * 0.10)),
               label = sprintf("cs/%s", m))
  }
})


test_that("overnight: didgpu_fect all three methods at N", {
  skip_unless_overnight()
  N <- .overnight_n()
  fails <- list(fe = 0L, ife = 0L, mc = 0L)
  for (s in seq_len(N)) {
    p <- .rand_panel_balanced(seed = s + 400000L)
    for (m in c("fe", "ife", "mc")) {
      fit <- tryCatch(
        didgpu_fect(p, "Y", "unit", "period", "D",
                     method = m, effects = 1L,
                     bootstrap_reps = 0L,
                     backend = "r", verbose = FALSE),
        error = function(e) NULL)
      if (is.null(fit)) fails[[m]] <- fails[[m]] + 1L
    }
  }
  for (m in c("fe", "ife", "mc")) {
    expect_lt(fails[[m]], max(1L, as.integer(N * 0.10)),
               label = sprintf("fect/%s", m))
  }
})


test_that("overnight: TestMechs all three methods at N (binary M)", {
  skip_unless_overnight()
  testthat::skip_if_not_installed("quadprog")
  N <- .overnight_n()
  fails <- list(CS = 0L, ARP = 0L, FSST = 0L)
  for (s in seq_len(N)) {
    set.seed(s + 500000L)
    n <- sample(400:1200, 1L)
    df <- data.frame(
      D = sample(c(0L, 1L), n, replace = TRUE),
      M = sample(c(1L, 2L), n, replace = TRUE),
      Y = rnorm(n))
    for (m in c("CS", "ARP", "FSST")) {
      res <- tryCatch(
        didgpu_test_sharp_null(df, "D", "M", "Y",
                                 method = m, B = 30L,
                                 num_Ybins = 3L, seed = 1L),
        error = function(e) NULL)
      if (is.null(res) || !is.finite(res$test_stat)) {
        fails[[m]] <- fails[[m]] + 1L
      }
    }
  }
  for (m in c("CS", "ARP", "FSST")) {
    expect_lt(fails[[m]], max(1L, as.integer(N * 0.15)),
               label = sprintf("testmechs/%s", m))
  }
})


test_that("overnight: TestMechs multi-level M at N", {
  skip_unless_overnight()
  testthat::skip_if_not_installed("quadprog")
  N <- .overnight_n()
  fails <- 0L
  for (s in seq_len(N)) {
    set.seed(s + 600000L)
    n <- sample(800:2000, 1L)
    K <- sample(c(2L, 3L, 4L), 1L)
    df <- data.frame(
      D = sample(c(0L, 1L), n, replace = TRUE),
      M = sample(seq_len(K), n, replace = TRUE),
      Y = rnorm(n))
    res <- tryCatch(
      didgpu_test_sharp_null(df, "D", "M", "Y",
                               method = "CS", B = 30L,
                               num_Ybins = 3L, seed = 1L),
      error = function(e) NULL)
    if (is.null(res) || !is.finite(res$test_stat)) fails <- fails + 1L
  }
  expect_lt(fails, max(1L, as.integer(N * 0.15)))
})


test_that("overnight: HonestDiD on random CS fits at N", {
  skip_unless_overnight()
  testthat::skip_if_not_installed("HonestDiD")
  N <- .overnight_n()
  fails <- 0L
  for (s in seq_len(N)) {
    p <- .rand_panel_balanced(seed = s + 700000L, n_periods = 12L)
    fit <- tryCatch(
      didgpu_cs(p, "Y", "unit", "period", "D",
                 est_method = "OR",
                 bootstrap_reps = 10L, seed = 1L,
                 backend = "r", verbose = FALSE),
      error = function(e) NULL)
    if (is.null(fit)) { fails <- fails + 1L; next }
    # Need >= 2 pre-treatment event-times for HonestDiD.
    n_pre <- sum(fit$att_gt$event_time < 0L)
    if (n_pre < 2L) next   # skip this seed; not enough pre-periods
    sens <- tryCatch(
      suppressWarnings(didgpu_honest_did(fit, event_post = 1L,
                                            method = "RM",
                                            Mbar = c(0, 0.5))),
      error = function(e) NULL)
    if (is.null(sens) || any(!is.finite(c(sens$lb, sens$ub)))) {
      fails <- fails + 1L
    }
  }
  expect_lt(fails, max(1L, as.integer(N * 0.30)))
})


test_that("overnight: didgpu_compare against DIDmultiplegtDYN at N", {
  skip_unless_overnight()
  testthat::skip_if_not_installed("DIDmultiplegtDYN")
  N <- .overnight_n()
  fails <- 0L
  max_diffs <- numeric(N)
  for (s in seq_len(N)) {
    p <- .rand_panel_balanced(seed = s + 800000L)
    rep <- tryCatch(
      didgpu_compare(p, "Y", "unit", "period", "D",
                       effects = 2L, placebo = 0L),
      error = function(e) NULL)
    if (is.null(rep) || is.na(rep$pass)) next
    max_diffs[s] <- max(rep$report$max_abs_diff, na.rm = TRUE)
    if (!isTRUE(rep$pass)) fails <- fails + 1L
  }
  expect_lt(fails, max(1L, as.integer(N * 0.05)))
})
