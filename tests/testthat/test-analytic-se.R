# SEs must be DIDmultiplegtDYN's analytic SEs, not a bootstrap
# approximation of them.
#
# Regression. didgpu had no analytic SE path at all:
#   * at bootstrap_reps = 0 the SE column came back all NA;
#   * at bootstrap_reps > 0 it reported the bootstrap SD, which is a
#     DIFFERENT estimator from the reference's, e.g.
#       didgpu (50 reps)  0.056174 0.051302 0.054590 0.059531 0.047412
#       reference         0.055504 0.054278 0.057839 0.052243 0.055922
#     test-reference-parity.R said so outright ("SEs come from different
#     estimators (bootstrap vs. analytic) and are not compared here"),
#     so the gap was known and simply never closed.
#
# It was also the whole performance story. The reference computes its
# SEs in one pass; didgpu's default of 100 bootstrap reps meant 101 full
# fits, so didgpu came out ~5x slower overall despite being ~10-25x
# faster per fit. The default is now bootstrap_reps = 0.
#
# DESCRIPTION promises bit-for-bit equivalence, so these are equality
# tests. See R/analytic_se.R for the derivation.

.ase_panel <- function(seed = 7L, nu = 200L, np = 20L, kind = "binary") {
  set.seed(seed)
  ufe <- stats::rnorm(nu, 0, 1); tfe <- stats::rnorm(np, 0, 0.3)
  Fg <- rep(Inf, nu)
  tr <- sort(sample(seq_len(nu), round(nu * 0.6)))
  Fg[tr] <- sample(4:(np - 4), length(tr), replace = TRUE)
  g <- expand.grid(period = seq_len(np), unit = seq_len(nu))
  g <- g[order(g$unit, g$period), ]
  on <- g$period >= Fg[g$unit]
  g$D <- switch(kind,
    binary = as.numeric(on),
    dose   = { set.seed(seed + 1L); dd <- sample(c(1, 3, 5), nu, TRUE)
               ifelse(on, dd[g$unit], 0) },
    nonabs = as.numeric(on & g$period < Fg[g$unit] + 3),
    bothdir = { set.seed(seed + 2L); dd <- sample(c(2, -2), nu, TRUE)
                ifelse(on, 1 + dd[g$unit], 1) })
  g$Y <- ufe[g$unit] + tfe[g$period] + 0.3 * g$D +
         stats::rnorm(nrow(g), 0, 0.4)
  g$state <- ((g$unit - 1L) %/% 5L) + 1L
  g[, c("unit", "period", "D", "Y", "state")]
}
.ase_fit <- function(d, bk = "r", ...) {
  suppressMessages(suppressWarnings(
    didgpu(df = d, outcome = "Y", group = "unit", time = "period",
           treatment = "D", backend = bk, verbose = FALSE, ...)))
}
.ase_backends <- function() {
  bk <- c("r", "cpu")
  if (isTRUE(tryCatch(didgpu_has_cuda_support(), error = function(e) FALSE)))
    bk <- c(bk, "cuda", "auto")
  bk
}
.ase_ref <- function(d, ...) {
  if (!"package:polars" %in% search()) {
    suppressMessages(attachNamespace("polars"))
  }
  suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = d, outcome = "Y", group = "unit", time = "period",
      treatment = "D", graph_off = TRUE, ...)))
}

test_that("SEs are reported without any bootstrap", {
  # Pre-fix this column was entirely NA at bootstrap_reps = 0.
  f <- .ase_fit(.ase_panel(), effects = 4L, placebo = 2L, bootstrap_reps = 0L)
  expect_true(all(is.finite(f$results$Effects[, "SE"])))
  expect_true(all(is.finite(f$results$Placebos[, "SE"])))
  expect_true(is.finite(f$results$ATE[1L, "SE"]))
})

test_that("bootstrap_reps defaults to 0", {
  # The reference needs one pass; didgpu used to take 101 by default.
  expect_identical(eval(formals(didgpu)$bootstrap_reps), 0L)
})

test_that("effect, placebo and ATE SEs match DIDmultiplegtDYN", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  cases <- list(
    list(lbl = "binary",        d = .ase_panel(),                     a = list()),
    list(lbl = "dose",          d = .ase_panel(11L, kind = "dose"),   a = list()),
    list(lbl = "non-absorbing", d = .ase_panel(5L, kind = "nonabs"),  a = list()),
    list(lbl = "in and out",    d = .ase_panel(9L, kind = "bothdir"), a = list()),
    list(lbl = "normalized",    d = .ase_panel(11L, kind = "dose"),
         a = list(normalized = TRUE)),
    list(lbl = "switchers=in",  d = .ase_panel(9L, kind = "bothdir"),
         a = list(switchers = "in")),
    list(lbl = "switchers=out", d = .ase_panel(9L, kind = "bothdir"),
         a = list(switchers = "out"))
  )
  for (cs in cases) {
    ours <- do.call(.ase_fit, c(list(cs$d, "r", effects = 4L, placebo = 2L,
                                     bootstrap_reps = 0L), cs$a))
    theirs <- do.call(.ase_ref, c(list(cs$d, effects = 4, placebo = 2), cs$a))
    expect_equal(as.numeric(ours$results$Effects[, "SE"]),
                 as.numeric(theirs$results$Effects[, "SE"]),
                 tolerance = 1e-12, label = paste("effects SE,", cs$lbl))
    expect_equal(as.numeric(ours$results$Placebos[, "SE"]),
                 as.numeric(theirs$results$Placebos[, "SE"]),
                 tolerance = 1e-12, label = paste("placebo SE,", cs$lbl))
    expect_equal(ours$results$ATE[1L, "SE"],
                 as.numeric(theirs$results$ATE[1L, "SE"]),
                 tolerance = 1e-12, label = paste("ATE SE,", cs$lbl))
  }
})

test_that("clustered SEs match DIDmultiplegtDYN", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  d <- .ase_panel()
  ours <- .ase_fit(d, "r", effects = 4L, placebo = 2L, bootstrap_reps = 0L,
                   cluster = "state")
  theirs <- .ase_ref(d, effects = 4, placebo = 2, cluster = "state")
  expect_equal(as.numeric(ours$results$Effects[, "SE"]),
               as.numeric(theirs$results$Effects[, "SE"]), tolerance = 1e-12)
  expect_equal(as.numeric(ours$results$Placebos[, "SE"]),
               as.numeric(theirs$results$Placebos[, "SE"]), tolerance = 1e-12)
  expect_equal(ours$results$ATE[1L, "SE"],
               as.numeric(theirs$results$ATE[1L, "SE"]), tolerance = 1e-12)
  # Guard: clustering must actually move the SEs, or this proves nothing.
  unc <- .ase_fit(d, "r", effects = 4L, placebo = 2L, bootstrap_reps = 0L)
  expect_gt(max(abs(as.numeric(ours$results$Effects[, "SE"]) -
                    as.numeric(unc$results$Effects[, "SE"]))), 1e-6)
})

test_that("joint nullity tests match DIDmultiplegtDYN", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  # These come from the same influence functions, via the covariance
  # identity cov(i,j) = (var(U_i + U_j) - var(U_i) - var(U_j)) / 2.
  d <- .ase_panel()
  ours <- .ase_fit(d, "r", effects = 5L, placebo = 3L, bootstrap_reps = 0L)
  theirs <- .ase_ref(d, effects = 5, placebo = 3)
  expect_equal(ours$results$p_jointeffects,
               as.numeric(theirs$results$p_jointeffects), tolerance = 1e-10)
  expect_equal(ours$results$p_jointplacebo,
               as.numeric(theirs$results$p_jointplacebo), tolerance = 1e-10)
})

test_that("SEs do not depend on which backend computed them", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  d <- .ase_panel()
  theirs <- .ase_ref(d, effects = 5, placebo = 3)
  ref_e <- as.numeric(theirs$results$Effects[, "SE"])
  ref_p <- as.numeric(theirs$results$Placebos[, "SE"])
  ref_a <- as.numeric(theirs$results$ATE[1L, "SE"])
  for (bk in .ase_backends()) {
    f <- .ase_fit(d, bk, effects = 5L, placebo = 3L, bootstrap_reps = 0L)
    expect_equal(as.numeric(f$results$Effects[, "SE"]), ref_e,
                 tolerance = 1e-12, label = paste("effects SE on", bk))
    expect_equal(as.numeric(f$results$Placebos[, "SE"]), ref_p,
                 tolerance = 1e-12, label = paste("placebo SE on", bk))
    expect_equal(f$results$ATE[1L, "SE"], ref_a,
                 tolerance = 1e-12, label = paste("ATE SE on", bk))
  }
})

test_that("asking for a bootstrap does not change the reported SE", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  # The bootstrap still runs (the joint tests fall back to its covariance
  # when the influence matrix is unavailable), but the SE column stays
  # analytic.
  d <- .ase_panel()
  a <- .ase_fit(d, "r", effects = 4L, placebo = 2L, bootstrap_reps = 0L)
  b <- .ase_fit(d, "r", effects = 4L, placebo = 2L, bootstrap_reps = 20L,
                seed = 1L)
  expect_equal(as.numeric(a$results$Effects[, "SE"]),
               as.numeric(b$results$Effects[, "SE"]), tolerance = 1e-12)
})

test_that("a checkpoint from another version is refused, not resumed", {
  # Cells written by an older didgpu may come from a different estimator
  # (the ATE changed when it became Av_tot_eff; the SEs changed when they
  # became analytic). Resuming would mix them silently.
  d <- .ase_panel(nu = 40L, np = 10L)
  cdir <- file.path(tempdir(), paste0("ckpt_", as.integer(runif(1, 1, 1e8))))
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)
  suppressMessages(suppressWarnings(
    didgpu(df = d, outcome = "Y", group = "unit", time = "period",
           treatment = "D", effects = 2L, placebo = 0L, bootstrap_reps = 2L,
           backend = "r", verbose = FALSE, checkpoint_dir = cdir)))
  meta_path <- file.path(cdir, "meta.json")
  expect_true(file.exists(meta_path))
  m0 <- jsonlite::fromJSON(meta_path)
  expect_identical(as.integer(m0$cell_rev), didgpu:::.didgpu_cell_rev)
  rerun <- function() suppressMessages(suppressWarnings(
    didgpu(df = d, outcome = "Y", group = "unit", time = "period",
           treatment = "D", effects = 2L, placebo = 0L,
           bootstrap_reps = 2L, backend = "r", verbose = FALSE,
           checkpoint_dir = cdir)))
  put <- function(m) writeLines(jsonlite::toJSON(m, auto_unbox = TRUE,
                                pretty = TRUE, null = "null"),
                                meta_path, useBytes = TRUE)

  # An unchanged checkpoint resumes.
  expect_no_error(rerun())

  # An older estimator revision is refused -- even under the SAME
  # package version string, which is exactly the case a version-only
  # check missed (every build this release reports 0.1.2).
  m <- m0; m$cell_rev <- didgpu:::.didgpu_cell_rev - 1L; put(m)
  expect_error(rerun(), "different didgpu build")

  # A checkpoint from before the stamp existed is refused too.
  m <- m0; m$cell_rev <- NULL; put(m)
  expect_error(rerun(), "different didgpu build")
})

test_that("controls fall back to the bootstrap rather than a near-miss SE", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  # With `controls` the reference subtracts a control-estimation
  # correction from the influence function (part2_switch,
  # did_multiplegt_dyn_core.R:502-535) that didgpu does not compute.
  # Reporting the uncorrected number would be wrong by ~9e-05, so the
  # analytic path is switched off there: SEs are NA unless a bootstrap
  # is requested, and the user is told.
  set.seed(13); nu <- 150L; np <- 16L
  ufe <- stats::rnorm(nu, 0, 1); tfe <- stats::rnorm(np, 0, 0.3)
  Fg <- rep(Inf, nu); tr <- sort(sample(seq_len(nu), 90L))
  Fg[tr] <- sample(4:12, 90L, replace = TRUE)
  g <- expand.grid(period = seq_len(np), unit = seq_len(nu))
  g <- g[order(g$unit, g$period), ]
  g$D <- as.numeric(g$period >= Fg[g$unit])
  g$X <- stats::rnorm(nrow(g))
  g$Y <- ufe[g$unit] + tfe[g$period] + 0.3 * g$D + 0.5 * g$X +
         stats::rnorm(nrow(g), 0, 0.4)
  g <- g[, c("unit", "period", "D", "Y", "X")]

  expect_message(
    didgpu(df = g, outcome = "Y", group = "unit", time = "period",
           treatment = "D", effects = 3L, placebo = 0L, bootstrap_reps = 0L,
           backend = "r", verbose = FALSE, controls = "X"),
    "analytic standard errors are not available")

  f <- .ase_fit(g, "r", effects = 3L, placebo = 0L, bootstrap_reps = 0L,
                controls = "X")
  expect_true(all(is.na(f$results$Effects[, "SE"])))
  # The point estimates are still exact against the reference.
  theirs <- .ase_ref(g, effects = 3, placebo = 0, controls = "X")
  expect_equal(as.numeric(f$results$Effects[, "Estimate"]),
               as.numeric(theirs$results$Effects[, "Estimate"]),
               tolerance = 1e-12)
  # And a bootstrap still gives a usable SE.
  fb <- .ase_fit(g, "r", effects = 3L, placebo = 0L, bootstrap_reps = 10L,
                 seed = 1L, controls = "X")
  expect_true(all(is.finite(fb$results$Effects[, "SE"])))
})

test_that("vcov(), the joint tests and didgpu_joint_placebo() agree", {
  # Once the SEs became analytic, vcov() was still the bootstrap
  # covariance: its diagonal stopped matching SE^2 (0.0039 vs 0.0132 in
  # the R CMD check failure), and didgpu_joint_placebo() -- which reads
  # vcov -- stopped reproducing p_jointplacebo (0.641 vs 0.678). All
  # three now come from one analytic matrix.
  f <- .ase_fit(.ase_panel(), "r", effects = 3L, placebo = 3L,
                bootstrap_reps = 0L)
  v <- vcov(f)
  expect_true(all(is.finite(v)))
  expect_true(isSymmetric(unname(v), tol = 1e-12))
  expect_equal(unname(diag(v)),
               c(as.numeric(f$results$Effects[, "SE"]),
                 as.numeric(f$results$Placebos[, "SE"]))^2,
               tolerance = 1e-12)
  expect_equal(didgpu_joint_placebo(f)$p_value, f$results$p_jointplacebo,
               tolerance = 1e-8)
})

test_that("joint tests under normalized = TRUE match DIDmultiplegtDYN", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  # The reference divides each influence vector by delta_k BEFORE
  # polarising (did_multiplegt_main.R:1170-1175). Polarising the raw
  # vectors against the normalised SEs gives a wrong covariance, so this
  # is pinned on a design where delta_k != 1.
  d <- .ase_panel(11L, kind = "dose")
  ours <- .ase_fit(d, "r", effects = 4L, placebo = 3L, bootstrap_reps = 0L,
                   normalized = TRUE)
  theirs <- .ase_ref(d, effects = 4, placebo = 3, normalized = TRUE)
  expect_equal(ours$results$p_jointeffects,
               as.numeric(theirs$results$p_jointeffects), tolerance = 1e-10)
  expect_equal(ours$results$p_jointplacebo,
               as.numeric(theirs$results$p_jointplacebo), tolerance = 1e-10)
})
