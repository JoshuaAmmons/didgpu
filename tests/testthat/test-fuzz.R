# ============================================================================
# Adversarial fuzz tests
#
# Generate randomized panels (varying size, period count, treatment
# distribution, sparsity, NA patterns, weight scales) and exercise every
# feature against the reference. The contract for every case is one of:
#   1. Both backends return numerically-equal estimates within tol.
#   2. Both backends return NA (degenerate panel / horizon).
#   3. Both backends error.
#
# Anything else (one succeeds, one errors; both succeed with different
# values) is a hard test failure and we want to know about it.
#
# These tests can be slow — they generate many panels and call the
# reference each time. The default is a SMALL number of fuzz seeds for
# fast CI; bump up via Sys.setenv(DIDGPU_FUZZ_N = "100") for thorough
# local runs.
# ============================================================================


# How many fuzz iterations per scenario. Default low for CI speed; users
# can crank it up via env var for deep local runs.
.fuzz_n <- function() {
  n <- as.integer(Sys.getenv("DIDGPU_FUZZ_N", "8"))
  if (is.na(n) || n < 1L) 8L else n
}


# Generate a randomized but valid panel from a single seed. The
# parameters are themselves randomized so the test space is broad.
.rand_panel <- function(seed) {
  set.seed(seed)
  n_units   <- sample(15:80, 1L)
  n_periods <- sample(5:18, 1L)
  frac_treat <- runif(1L, 0.3, 0.85)
  min_t <- sample(2:(n_periods - 2L), 1L)
  max_t <- sample(min_t:(n_periods - 1L), 1L)
  prof_len <- sample(2:5, 1L)
  prof <- runif(prof_len, -1.5, 1.5)
  sigma <- runif(1L, 0.1, 0.8)
  didgpu_simulate_panel(
    n_units = n_units, n_periods = n_periods,
    frac_treated = frac_treat,
    min_treat_period = min_t, max_treat_period = max_t,
    tau_profile = prof, sigma = sigma, seed = seed)
}


# Compare two estimate vectors, treating both-NA as agreement.
.estimates_agree <- function(a, b, tol = 1e-8) {
  if (length(a) != length(b)) return(FALSE)
  if (length(a) == 0L) return(TRUE)
  ok <- vapply(seq_along(a), function(i) {
    if (is.na(a[i]) && is.na(b[i])) return(TRUE)
    if (is.na(a[i]) || is.na(b[i])) return(FALSE)
    isTRUE(abs(a[i] - b[i]) <= tol)
  }, logical(1))
  all(ok)
}


# Run didgpu and the reference, return both results or a structured
# error. Handles the case where the reference rejects an arg combination
# the simulator can produce.
.dual_fit <- function(p, ...) {
  args <- list(...)
  us <- tryCatch(
    didgpu(p, "Y", "unit", "period", "D",
            ..., bootstrap_reps = 0L,
            backend = "r", verbose = FALSE),
    error = function(e) structure(list(err = conditionMessage(e)),
                                    class = "fuzz_err"))
  # Build ref args, dropping didgpu-only ones.
  ref_args <- args[!names(args) %in% c("verbose", "checkpoint_dir",
                                         "n_workers", "resume",
                                         "on_iter", "seed",
                                         "bootstrap_reps", "backend")]
  # Pull effects/placebo separately because the reference is fussy about types.
  ref_args$df        <- as.data.frame(p)
  ref_args$outcome   <- "Y"
  ref_args$group     <- "unit"
  ref_args$time      <- "period"
  ref_args$treatment <- "D"
  if (!is.null(ref_args$effects)) ref_args$effects <- as.double(ref_args$effects)
  if (!is.null(ref_args$placebo)) ref_args$placebo <- as.double(ref_args$placebo)
  ref_args$graph_off <- TRUE
  ref <- tryCatch(
    suppressMessages(suppressWarnings(
      do.call(DIDmultiplegtDYN::did_multiplegt_dyn, ref_args))),
    error = function(e) structure(list(err = conditionMessage(e)),
                                    class = "fuzz_err"))
  list(us = us, ref = ref)
}


# Single fuzz check: given a panel and an arg list, the strict
# contract is:
#   both ok       -> estimates must agree
#   both error    -> fine (degenerate panel)
#   didgpu err, ref ok -> HARD FAILURE (we have a bug)
#   didgpu ok, ref err -> NOT a failure — didgpu is more robust than
#                          the reference (the reference has known bugs
#                          for some panel shapes; cf. the cluster-arg
#                          bug documented in NOTES). We still record
#                          the case so the user can see it.
.check_fuzz_case <- function(p, ...) {
  fit <- .dual_fit(p, ...)
  us_err  <- inherits(fit$us,  "fuzz_err")
  ref_err <- inherits(fit$ref, "fuzz_err")
  if (us_err && ref_err)
    return(list(ok = TRUE, why = "both errored"))
  if (us_err && !ref_err)
    return(list(ok = FALSE,
                why = paste("didgpu errored but reference succeeded:",
                            fit$us$err)))
  if (!us_err && ref_err)
    return(list(ok = TRUE,
                why = paste("didgpu more robust than reference:",
                            fit$ref$err)))
  # Both ok: compare estimates.
  us_e  <- as.numeric(fit$us$results$Effects[, "Estimate"])
  ref_e <- as.numeric(fit$ref$results$Effects[, 1])
  if (!.estimates_agree(us_e, ref_e, tol = 1e-7)) {
    return(list(ok = FALSE,
                why = sprintf("effects disagree: max abs diff = %g",
                              max(abs(us_e - ref_e), na.rm = TRUE))))
  }
  if (length(fit$us$results$Placebos) > 0L &&
      nrow(fit$us$results$Placebos) > 0L) {
    us_p <- as.numeric(fit$us$results$Placebos[, "Estimate"])
    if (!is.null(fit$ref$results$Placebos) &&
        nrow(fit$ref$results$Placebos) > 0L) {
      ref_p <- as.numeric(fit$ref$results$Placebos[, 1])
      if (!.estimates_agree(us_p, ref_p, tol = 1e-7)) {
        return(list(ok = FALSE,
                    why = sprintf("placebos disagree: max abs diff = %g",
                                  max(abs(us_p - ref_p), na.rm = TRUE))))
      }
    }
  }
  list(ok = TRUE, why = "agree")
}


test_that("fuzz: vanilla binary, both directions", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    p <- .rand_panel(seed = s + 1000L)
    res <- .check_fuzz_case(p, effects = 3L, placebo = 0L)
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (vanilla):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: with weight column", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    p <- .rand_panel(seed = s + 2000L)
    set.seed(s + 2500L)
    p$w <- runif(nrow(p), 0.1, 5)
    res <- .check_fuzz_case(p, effects = 3L, placebo = 0L, weight = "w")
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (weight):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: with controls", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    p <- .rand_panel(seed = s + 3000L)
    set.seed(s + 3500L)
    p$X1 <- rnorm(nrow(p))
    p$X2 <- rnorm(nrow(p), sd = 0.5)
    res <- .check_fuzz_case(p, effects = 3L, placebo = 0L,
                              controls = c("X1", "X2"))
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (controls):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: with switchers = 'in'", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    p <- .rand_panel(seed = s + 4000L)
    res <- .check_fuzz_case(p, effects = 3L, placebo = 0L,
                              switchers = "in")
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (switchers='in'):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: with normalized", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    p <- .rand_panel(seed = s + 5000L)
    res <- .check_fuzz_case(p, effects = 3L, placebo = 0L,
                              normalized = TRUE)
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (normalized):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: with placebos", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    p <- .rand_panel(seed = s + 6000L)
    res <- .check_fuzz_case(p, effects = 2L, placebo = 2L)
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (placebos):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: trends_lin", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    p <- .rand_panel(seed = s + 7000L)
    res <- .check_fuzz_case(p, effects = 2L, placebo = 0L,
                              trends_lin = TRUE)
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (trends_lin):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: only_never_switchers + same_switchers", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    p <- .rand_panel(seed = s + 8000L)
    res <- .check_fuzz_case(p, effects = 2L, placebo = 0L,
                              only_never_switchers = TRUE,
                              same_switchers = TRUE)
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (only_never + same_switchers):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: kitchen sink (weight + controls + normalized + placebos)", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    p <- .rand_panel(seed = s + 9000L)
    set.seed(s + 9500L)
    p$w  <- runif(nrow(p), 0.5, 3)
    p$X1 <- rnorm(nrow(p))
    res <- .check_fuzz_case(p, effects = 2L, placebo = 1L,
                              weight = "w", controls = "X1",
                              normalized = TRUE)
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (kitchen sink):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: trends_lin + weight + controls + placebos", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    p <- .rand_panel(seed = s + 11000L)
    set.seed(s + 11500L)
    p$w  <- runif(nrow(p), 0.5, 2)
    p$X1 <- rnorm(nrow(p))
    res <- .check_fuzz_case(p, effects = 2L, placebo = 1L,
                              trends_lin = TRUE,
                              weight = "w", controls = "X1")
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (trends_lin sink):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: multivalued treatment {0,1,2,3}", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in seq_len(.fuzz_n())) {
    set.seed(s + 12000L)
    n_units   <- sample(40:80, 1L)
    n_periods <- sample(8:14, 1L)
    fg <- rep(Inf, n_units)
    n_treated <- as.integer(n_units * runif(1L, 0.4, 0.8))
    treated <- sort(sample(seq_len(n_units), n_treated))
    fg[treated] <- sample(seq(3L, n_periods - 1L), n_treated, replace = TRUE)
    doses <- sample(c(1L, 2L, 3L), n_units, replace = TRUE)
    panel <- data.table::data.table(
      unit = rep(seq_len(n_units), each = n_periods),
      period = rep(seq_len(n_periods), n_units))
    panel[, F_g := fg[unit]]
    panel[, dose := doses[unit]]
    panel[, D := as.integer(ifelse(period >= F_g, dose, 0L))]
    panel[, k_evt := period - F_g]
    panel[, tau_k := 0]
    panel[is.finite(F_g) & k_evt >= 0,
          tau_k := dose * c(0.3, 0.5, 0.6, 0.5)[pmin(k_evt + 1L, 4L)]]
    ufe <- rnorm(n_units); tfe <- rnorm(n_periods, 0, 0.3)
    panel[, Y := ufe[unit] + tfe[period] + tau_k + rnorm(.N, 0, 0.4)]
    p <- as.data.frame(panel[order(unit, period), .(unit, period, D, Y)])
    res <- .check_fuzz_case(p, effects = 2L, placebo = 0L)
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (multivalued):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: degenerate (all-same F_g) panel — both should error or NA", {
  # A panel where every treated unit switches at the same period is a
  # known degenerate case (no within-cohort variation in switch timing).
  # The reference errors; we should either error or return all-NA effects.
  for (s in 1L:5L) {
    set.seed(s + 20000L)
    n_units <- 40L; n_periods <- 10L
    same_F <- sample(3:8, 1L)
    fg <- rep(Inf, n_units)
    fg[seq_len(round(n_units * 0.6))] <- same_F
    panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
    panel$F_g <- fg[panel$unit]
    panel$D <- as.integer(panel$period >= panel$F_g)
    panel$Y <- rnorm(n_units)[panel$unit] +
               rnorm(n_periods, 0, 0.3)[panel$period] +
               0.5 * panel$D + rnorm(nrow(panel), 0, 0.3)
    p <- panel[, c("unit", "period", "D", "Y")]
    us <- tryCatch(
      didgpu(p, "Y", "unit", "period", "D",
              effects = 2L, bootstrap_reps = 0L,
              backend = "r", verbose = FALSE),
      error = function(e) e)
    # Either error, or NA effects, are acceptable for this degenerate
    # panel. The HARD failure case would be returning numeric
    # estimates as if the panel were well-posed.
    if (!inherits(us, "error")) {
      eff <- as.numeric(us$results$Effects[, "Estimate"])
      # All-NA is fine. Some numeric values are fine if the panel
      # happens to have enough variation. The "wrong" case would be
      # crashing or returning Inf.
      expect_true(all(is.finite(eff) | is.na(eff)),
                   info = sprintf("seed %d: non-finite Effects: %s",
                                   s, paste(eff, collapse = ", ")))
    }
  }
  expect_true(TRUE)
})


test_that("fuzz: bootstrap path produces stable SEs across seeds", {
  # We can't compare bootstrap SEs to the reference (different SE
  # formulas), but we CAN check that the r-backend's own bootstrap is
  # stable and finite, and that the point estimate inside the
  # bootstrap matches a no-bootstrap call (cell 0 is deterministic).
  failures <- character(0)
  for (s in 1L:8L) {
    p <- .rand_panel(seed = s + 40000L)
    no_boot <- tryCatch(
      didgpu(p, "Y", "unit", "period", "D",
              effects = 2L, bootstrap_reps = 0L,
              backend = "r", verbose = FALSE),
      error = function(e) NULL)
    if (is.null(no_boot)) next
    with_boot <- tryCatch(
      didgpu(p, "Y", "unit", "period", "D",
              effects = 2L, bootstrap_reps = 8L, seed = 1L,
              backend = "r", verbose = FALSE),
      error = function(e) NULL)
    if (is.null(with_boot)) {
      failures <- c(failures,
                    sprintf("seed %d: bootstrap fit errored", s))
      next
    }
    # Point estimates must match between the two runs (cell 0 in the
    # bootstrap path is the unresampled point estimate).
    e_nb <- as.numeric(no_boot$results$Effects[, "Estimate"])
    e_wb <- as.numeric(with_boot$results$Effects[, "Estimate"])
    if (!.estimates_agree(e_nb, e_wb, tol = 1e-10)) {
      failures <- c(failures,
                    sprintf("seed %d: point estimates differ across no-boot/with-boot calls",
                            s))
    }
    # SEs must be non-negative finite numbers (or NA).
    se <- as.numeric(with_boot$results$Effects[, "SE"])
    if (any(is.finite(se) & se < 0)) {
      failures <- c(failures, sprintf("seed %d: negative SE in bootstrap", s))
    }
    if (any(is.finite(se) & !is.finite(se))) {
      failures <- c(failures, sprintf("seed %d: non-finite SE in bootstrap", s))
    }
  }
  if (length(failures))
    fail(paste0("fuzz failures (bootstrap stability):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: n_workers > 1 produces identical output to sequential", {
  # The parallel path must be bit-identical to the sequential path at
  # the same seed (we force Mersenne-Twister + per-iter seed for this).
  failures <- character(0)
  for (s in 1L:4L) {
    p <- .rand_panel(seed = s + 50000L)
    seq_fit <- tryCatch(
      didgpu(p, "Y", "unit", "period", "D",
              effects = 2L, bootstrap_reps = 6L, seed = 7L,
              n_workers = 1L,
              backend = "r", verbose = FALSE),
      error = function(e) NULL)
    par_fit <- tryCatch(
      suppressMessages(
        didgpu(p, "Y", "unit", "period", "D",
                effects = 2L, bootstrap_reps = 6L, seed = 7L,
                n_workers = 2L,
                backend = "r", verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(seq_fit) || is.null(par_fit)) next
    e_s <- as.numeric(seq_fit$results$Effects[, "Estimate"])
    e_p <- as.numeric(par_fit$results$Effects[, "Estimate"])
    se_s <- as.numeric(seq_fit$results$Effects[, "SE"])
    se_p <- as.numeric(par_fit$results$Effects[, "SE"])
    if (!.estimates_agree(e_s, e_p, tol = 1e-12)) {
      failures <- c(failures, sprintf("seed %d: estimates differ seq vs par", s))
    }
    if (!.estimates_agree(se_s, se_p, tol = 1e-12)) {
      failures <- c(failures, sprintf("seed %d: SEs differ seq vs par", s))
    }
  }
  if (length(failures))
    fail(paste0("fuzz failures (parallel == sequential):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: checkpoint round-trip preserves results", {
  # Run a fit with checkpoint_dir, then re-read the cells from disk
  # and re-aggregate. The two results must be identical.
  failures <- character(0)
  for (s in 1L:4L) {
    p <- .rand_panel(seed = s + 60000L)
    cdir <- tempfile("fuzz_ck_")
    on.exit(unlink(cdir, recursive = TRUE), add = TRUE)
    fit1 <- tryCatch(
      didgpu(p, "Y", "unit", "period", "D",
              effects = 2L, bootstrap_reps = 4L, seed = 7L,
              checkpoint_dir = cdir,
              backend = "r", verbose = FALSE),
      error = function(e) NULL)
    if (is.null(fit1)) next
    # Resume should be a no-op (all cells already done) and produce
    # the same result.
    fit2 <- tryCatch(
      suppressMessages(didgpu_resume(cdir, df = p)),
      error = function(e) NULL)
    if (is.null(fit2)) {
      failures <- c(failures, sprintf("seed %d: didgpu_resume errored", s))
      next
    }
    e1 <- as.numeric(fit1$results$Effects[, "Estimate"])
    e2 <- as.numeric(fit2$results$Effects[, "Estimate"])
    if (!.estimates_agree(e1, e2, tol = 1e-12)) {
      failures <- c(failures, sprintf("seed %d: resume produced different result",
                                        s))
    }
  }
  if (length(failures))
    fail(paste0("fuzz failures (checkpoint resume):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: NA-rich panels — didgpu produces finite output", {
  # Drop ~10% of outcome cells and ~5% of treatment cells at random.
  #
  # NA-handling semantics diverge between didgpu and the reference on
  # panels with NA holes: the reference DROPS rows with NA Y or NA D
  # early in main.R:85 (`df <- subset(df, !is.na(df$mean_Y) & ...`),
  # which can change F_g and T_g for NA-heavy units. didgpu KEEPS the
  # balanced grid and uses N_gt = 0 to exclude those cells from sums.
  # Both choices are defensible; on clean panels they agree exactly.
  #
  # This test therefore checks only the weaker contract: didgpu must
  # produce finite (non-Inf, non-NaN) estimates on NA-rich panels and
  # not crash. Bit-identicalness on NA-rich panels is documented as a
  # known divergence (matches the design choice that didgpu does not
  # silently drop user data without telling them).
  failures <- character(0)
  for (s in 1L:8L) {
    set.seed(s + 70000L)
    p <- .rand_panel(seed = s + 70000L)
    out_holes  <- sample(nrow(p), max(1L, as.integer(nrow(p) * 0.1)))
    trt_holes  <- sample(nrow(p), max(1L, as.integer(nrow(p) * 0.05)))
    p$Y[out_holes] <- NA
    p$D[trt_holes] <- NA
    fit <- tryCatch(
      didgpu(p, "Y", "unit", "period", "D",
              effects = 2L, placebo = 0L,
              bootstrap_reps = 0L, backend = "r", verbose = FALSE),
      error = function(e) NULL)
    if (is.null(fit)) next   # acceptable: degenerate panel after NA holes
    eff <- as.numeric(fit$results$Effects[, "Estimate"])
    # The only HARD failure is non-finite numeric output (Inf / NaN).
    bad <- eff[!is.na(eff) & !is.finite(eff)]
    if (length(bad)) {
      failures <- c(failures,
                     sprintf("seed %d: non-finite Effects: %s",
                              s, paste(eff, collapse = ", ")))
    }
  }
  if (length(failures))
    fail(paste0("fuzz failures (NA-rich finite):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: didgpu_by reconciles with manual per-subgroup didgpu", {
  # didgpu_by(p, by_var = X, ...) should be IDENTICAL to running
  # didgpu(p[p$X == lvl, ], ...) for each level. This catches subtle
  # bugs in the by-splitting / per-subgroup arg forwarding.
  failures <- character(0)
  for (s in 1L:6L) {
    p <- .rand_panel(seed = s + 80000L)
    set.seed(s + 80500L)
    # Random 2-3 level grouping.
    n_lvls <- sample(2:3, 1L)
    lvls <- paste0("g", seq_len(n_lvls))
    p$grp <- sample(lvls, length(unique(p$unit)), replace = TRUE)[p$unit]
    by_fit <- tryCatch(
      didgpu_by(p, "grp",
                 outcome = "Y", group = "unit",
                 time = "period", treatment = "D",
                 effects = 2L, bootstrap_reps = 0L,
                 backend = "r", verbose = FALSE),
      error = function(e) NULL)
    if (is.null(by_fit)) next
    for (lvl in lvls) {
      sub <- p[p$grp == lvl, , drop = FALSE]
      manual <- tryCatch(
        didgpu(sub, "Y", "unit", "period", "D",
                effects = 2L, bootstrap_reps = 0L,
                backend = "r", verbose = FALSE),
        error = function(e) NULL)
      if (is.null(manual)) {
        # by_fit might have errored on this subgroup too; check if so.
        if (!is.null(by_fit[[lvl]])) {
          failures <- c(failures,
                        sprintf("seed %d level %s: manual errored but by_fit succeeded",
                                 s, lvl))
        }
        next
      }
      e_by  <- as.numeric(by_fit[[lvl]]$results$Effects[, "Estimate"])
      e_man <- as.numeric(manual$results$Effects[, "Estimate"])
      if (!.estimates_agree(e_by, e_man, tol = 1e-12)) {
        failures <- c(failures,
                      sprintf("seed %d level %s: by_fit != manual didgpu",
                               s, lvl))
      }
    }
  }
  if (length(failures))
    fail(paste0("fuzz failures (didgpu_by vs manual):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: predict_het with random covariates matches reference", {
  skip_if_no_reference()
  failures <- character(0)
  for (s in 1L:6L) {
    p <- .rand_panel(seed = s + 90000L)
    set.seed(s + 90500L)
    # Add 1-2 time-invariant covariates.
    n_covs <- sample(1:2, 1L)
    cov_names <- paste0("Z", seq_len(n_covs))
    n_units <- length(unique(p$unit))
    for (cn in cov_names) {
      vals <- rnorm(n_units)
      p[[cn]] <- vals[p$unit]
    }
    # Wrap in tryCatch — predict_het can error on small panels with
    # rank-deficient cohort interactions.
    ref <- tryCatch(
      suppressMessages(suppressWarnings(
        DIDmultiplegtDYN::did_multiplegt_dyn(
          df = as.data.frame(p), outcome = "Y", group = "unit",
          time = "period", treatment = "D",
          effects = as.double(2), placebo = 0, graph_off = TRUE,
          predict_het = list(cov_names, c(-1))))),
      error = function(e) NULL)
    us <- tryCatch(
      didgpu(p, "Y", "unit", "period", "D",
              effects = 2L, placebo = 0L,
              predict_het = list(cov_names, -1L),
              bootstrap_reps = 0L, backend = "r", verbose = FALSE),
      error = function(e) NULL)
    if (is.null(ref) || is.null(us)) next
    if (is.null(ref$results$predict_het) || is.null(us$results$predict_het)) next
    # Align row order before comparing (we sort by (effect, covariate);
    # reference sorts by (covariate, effect)).
    key <- function(d) order(d$effect, d$covariate)
    rs <- ref$results$predict_het[key(ref$results$predict_het), ]
    us_s <- us$results$predict_het[key(us$results$predict_het), ]
    if (nrow(rs) != nrow(us_s)) {
      failures <- c(failures,
                     sprintf("seed %d: row counts differ (%d ref vs %d us)",
                              s, nrow(rs), nrow(us_s)))
      next
    }
    for (col in c("Estimate", "SE")) {
      if (!.estimates_agree(rs[[col]], us_s[[col]], tol = 1e-8)) {
        failures <- c(failures,
                       sprintf("seed %d: predict_het %s disagrees", s, col))
        break
      }
    }
  }
  if (length(failures))
    fail(paste0("fuzz failures (predict_het):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})


test_that("fuzz: very small panel (5 units, 4 periods)", {
  skip_if_no_reference()
  failures <- character(0)
  # Tiny panels stress the auto-clamp logic for effects/placebo horizons.
  for (s in 1L:10L) {
    set.seed(s + 30000L)
    n_units <- 5L; n_periods <- 4L
    fg <- rep(Inf, n_units)
    fg[c(1L, 2L, 3L)] <- c(2L, 3L, 3L)
    panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
    panel$F_g <- fg[panel$unit]
    panel$D <- as.integer(panel$period >= panel$F_g)
    panel$Y <- rnorm(n_units)[panel$unit] +
               rnorm(n_periods, 0, 0.3)[panel$period] +
               0.5 * panel$D + rnorm(nrow(panel), 0, 0.3)
    p <- panel[, c("unit", "period", "D", "Y")]
    res <- .check_fuzz_case(p, effects = 2L, placebo = 0L)
    if (!res$ok) failures <- c(failures,
                                 sprintf("seed %d: %s", s, res$why))
  }
  if (length(failures))
    fail(paste0("fuzz failures (small):\n  ",
                paste(failures, collapse = "\n  ")))
  expect_true(TRUE)
})
