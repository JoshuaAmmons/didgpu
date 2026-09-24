# `ATE` must be DIDmultiplegtDYN's Av_tot_eff: the average total effect
# PER UNIT OF TREATMENT.
#
# Regression. didgpu formed the ATE as the switcher-weighted mean of the
# per-horizon DIDs:
#     ATE = sum_k N_k * DID_k / sum_k N_k
# The reference instead carries its own denominator (U_Gg_den_XX, built
# from delta_D_i_XX in did_multiplegt_dyn_core.R:912-923):
#     ATE = sum_k N_k * DID_k / sum_k N_k * delta_k
# where delta_k is the average current-period treatment change among
# event-time-k switchers.
#
# The two coincide exactly when treatment is binary and absorbing, because
# then delta_k == 1 for every k -- which is why the whole test suite, and
# every worked example in the docs, missed it. On any other design the
# reported ATE was off by the average dose: on a dose-5 panel it was 5x
# the reference's, and on a non-monotone 0 -> 2 -> 1 panel 1.5x.
#
# DESCRIPTION promises "bit-for-bit numerical equivalence" with the
# reference, so these are equality tests, not tolerance tests.

.est_base <- function(seed = 11L, nu = 90L, np = 12L) {
  set.seed(seed)
  ufe <- stats::rnorm(nu, 0, 1); tfe <- stats::rnorm(np, 0, 0.3)
  Fg <- rep(Inf, nu)
  tr <- sort(sample(seq_len(nu), round(nu * 0.6)))
  Fg[tr] <- sample(4:9, length(tr), replace = TRUE)
  g <- expand.grid(period = seq_len(np), unit = seq_len(nu))
  g <- g[order(g$unit, g$period), ]
  g$Fg <- Fg[g$unit]
  g$Y  <- ufe[g$unit] + tfe[g$period] +
          0.3 * as.numeric(g$period >= g$Fg) +
          stats::rnorm(nrow(g), 0, 0.4)
  g
}

.est_designs <- function() {
  g <- .est_base()
  bin <- g; bin$D <- as.numeric(bin$period >= bin$Fg)
  nab <- g; nab$D <- as.numeric(nab$period >= nab$Fg &
                                 nab$period < nab$Fg + 3)
  set.seed(3); dose <- sample(c(1, 3, 5), 90L, replace = TRUE)
  mv  <- g; mv$D <- ifelse(mv$period >= mv$Fg, dose[mv$unit], 0)
  nm  <- g; nm$D <- ifelse(nm$period >= nm$Fg,
                           ifelse(nm$period < nm$Fg + 2, 2, 1), 0)
  keep <- c("unit", "period", "D", "Y")
  list("binary absorbing"   = bin[, keep],
       "non-absorbing"      = nab[, keep],
       "multivalued dose"   = mv[, keep],
       "non-monotone dose"  = nm[, keep])
}

.est_fit <- function(d, bk = "r", ...) {
  suppressMessages(suppressWarnings(
    didgpu(df = d, outcome = "Y", group = "unit", time = "period",
           treatment = "D", effects = 4L, placebo = 0L,
           bootstrap_reps = 0L, backend = bk, verbose = FALSE, ...)))
}
.est_backends <- function() {
  bk <- c("r", "cpu")
  if (isTRUE(tryCatch(didgpu_has_cuda_support(), error = function(e) FALSE)))
    bk <- c(bk, "cuda", "auto")
  bk
}

test_that("scaling the dose by c divides the ATE by c", {
  # Closed form, no reference package needed. Multiplying the treatment
  # by a constant leaves every per-horizon effect untouched (the dose
  # enters only through F_g and the baseline, both unchanged) but scales
  # the per-unit-of-treatment denominator by exactly c. Under the old
  # definition the ATE did not move at all.
  g <- .est_base()
  on <- g$period >= g$Fg
  base <- g[, c("unit", "period", "Y")]; base$D <- as.numeric(on)
  ref <- .est_fit(base)
  ref_ate <- ref$results$ATE[1L, 1L]
  expect_true(is.finite(ref_ate))
  for (cc in c(2, 5, 0.5)) {
    d <- base; d$D <- cc * as.numeric(on)
    f <- .est_fit(d)
    # Per-horizon effects are untouched ...
    expect_equal(as.numeric(f$results$Effects[, "Estimate"]),
                 as.numeric(ref$results$Effects[, "Estimate"]),
                 tolerance = 1e-12, label = paste("effects at c =", cc))
    # ... while the ATE is exactly 1/c of the binary one.
    expect_equal(cc * f$results$ATE[1L, 1L], ref_ate, tolerance = 1e-12,
                 label = paste("ATE at c =", cc))
  }
})

test_that("`normalized` does not move the ATE", {
  # The reference normalises only the per-horizon effects; Av_tot_eff
  # carries its own denominator and is unaffected.
  d <- .est_designs()[["multivalued dose"]]
  expect_equal(.est_fit(d, normalized = TRUE)$results$ATE[1L, 1L],
               .est_fit(d)$results$ATE[1L, 1L], tolerance = 1e-12)
})

test_that("the ATE is backend-invariant on non-absorbing designs", {
  # The old formula lived separately in .backend_cpu and .backend_cuda,
  # and backend = "auto" resolves to cuda on a GPU machine, so this was
  # the default path.
  for (nm in names(.est_designs())) {
    d <- .est_designs()[[nm]]
    ref <- .est_fit(d, "r")$results$ATE[1L, 1L]
    expect_true(is.finite(ref), label = nm)
    for (bk in .est_backends()) {
      expect_equal(.est_fit(d, bk)$results$ATE[1L, 1L], ref,
                   tolerance = 1e-12, label = paste(nm, "on backend", bk))
    }
  }
})

test_that("the ATE matches DIDmultiplegtDYN's Av_tot_eff", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  # DIDmultiplegtDYN calls bare `pl$...`, so polars must be attached.
  if (!"package:polars" %in% search()) {
    suppressMessages(attachNamespace("polars"))
  }
  for (nm in names(.est_designs())) {
    d <- .est_designs()[[nm]]
    theirs <- suppressMessages(suppressWarnings(
      DIDmultiplegtDYN::did_multiplegt_dyn(
        df = d, outcome = "Y", group = "unit", time = "period",
        treatment = "D", effects = 4, placebo = 0, graph_off = TRUE)))
    expect_equal(.est_fit(d)$results$ATE[1L, 1L],
                 as.numeric(theirs$results$ATE[1L, 1L]),
                 tolerance = 1e-12, label = nm)
  }
})

test_that("the non-absorbing designs would have caught the old formula", {
  # Guard the tests above: if these designs happened to have delta_k == 1
  # they would prove nothing. Under the old definition the ATE equalled
  # the plain switcher-weighted mean; check that it no longer does.
  for (nm in c("multivalued dose", "non-monotone dose")) {
    f <- .est_fit(.est_designs()[[nm]])
    e <- as.numeric(f$results$Effects[, "Estimate"])
    n <- as.numeric(f$results$Effects[, "Switchers"])
    ok <- is.finite(e) & n > 0
    old <- sum(e[ok] * n[ok]) / sum(n[ok])
    expect_gt(abs(f$results$ATE[1L, 1L] - old), 1e-3)
  }
})

test_that("the ATE matches the reference with a non-zero baseline dose", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  if (!"package:polars" %in% search()) suppressMessages(attachNamespace("polars"))
  # Baseline treatment is 1, not 0, and units switch both up and down.
  # delta_k is then neither 1 nor equal across directions, so the pooled
  # denominator has to be right in both arms.
  set.seed(31); nu <- 90L; np <- 12L
  ufe <- stats::rnorm(nu, 0, 1); tfe <- stats::rnorm(np, 0, 0.3)
  Fg <- rep(Inf, nu); tr <- sort(sample(seq_len(nu), 60L))
  Fg[tr] <- sample(4:9, 60L, replace = TRUE)
  g <- expand.grid(period = seq_len(np), unit = seq_len(nu))
  g <- g[order(g$unit, g$period), ]
  set.seed(32); dir <- sample(c(2, -2, 3), nu, replace = TRUE)
  g$D <- ifelse(g$period >= Fg[g$unit], 1 + dir[g$unit], 1)
  g$Y <- ufe[g$unit] + tfe[g$period] + 0.3 * g$D +
         stats::rnorm(nrow(g), 0, 0.4)
  g <- g[, c("unit", "period", "D", "Y")]
  for (sw in c("", "in", "out")) {
    ours <- suppressMessages(suppressWarnings(
      didgpu(df = g, outcome = "Y", group = "unit", time = "period",
             treatment = "D", effects = 3L, placebo = 0L,
             bootstrap_reps = 0L, backend = "r", verbose = FALSE,
             switchers = sw)))
    theirs <- suppressMessages(suppressWarnings(
      DIDmultiplegtDYN::did_multiplegt_dyn(
        df = g, outcome = "Y", group = "unit", time = "period",
        treatment = "D", effects = 3, placebo = 0, graph_off = TRUE,
        switchers = sw)))
    expect_equal(ours$results$ATE[1L, 1L],
                 as.numeric(theirs$results$ATE[1L, 1L]),
                 tolerance = 1e-12, label = paste("switchers =", sw))
  }
})

test_that("the ATE matches the reference under `continuous`", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  if (!"package:polars" %in% search()) suppressMessages(attachNamespace("polars"))
  # Under `continuous` the treatment is binarised internally and the
  # original values are kept as treatment_XX_orig / d_sq_XX_orig. The
  # denominator has to read the ORIGINAL dose, as the reference does.
  set.seed(21); nu <- 80L; np <- 10L
  ufe <- stats::rnorm(nu, 0, 1); tfe <- stats::rnorm(np, 0, 0.3)
  Fg <- rep(Inf, nu); tr <- sort(sample(seq_len(nu), 48L))
  Fg[tr] <- sample(3:8, 48L, replace = TRUE)
  g <- expand.grid(period = seq_len(np), unit = seq_len(nu))
  g <- g[order(g$unit, g$period), ]
  set.seed(22)
  base <- stats::runif(nu, 0, 2); jump <- stats::runif(nu, -1.5, 1.5)
  g$D <- ifelse(g$period >= Fg[g$unit], base[g$unit] + jump[g$unit],
                base[g$unit])
  g$Y <- ufe[g$unit] + tfe[g$period] + 0.3 * g$D +
         stats::rnorm(nrow(g), 0, 0.4)
  g <- g[, c("unit", "period", "D", "Y")]
  for (cont in c(1, 2)) {
    ours <- suppressMessages(suppressWarnings(
      didgpu(df = g, outcome = "Y", group = "unit", time = "period",
             treatment = "D", effects = 3L, placebo = 0L,
             bootstrap_reps = 0L, backend = "r", verbose = FALSE,
             continuous = cont)))
    theirs <- suppressMessages(suppressWarnings(
      DIDmultiplegtDYN::did_multiplegt_dyn(
        df = g, outcome = "Y", group = "unit", time = "period",
        treatment = "D", effects = 3, placebo = 0, graph_off = TRUE,
        continuous = cont)))
    expect_equal(ours$results$ATE[1L, 1L],
                 as.numeric(theirs$results$ATE[1L, 1L]),
                 tolerance = 1e-12, label = paste("continuous =", cont))
  }
})
