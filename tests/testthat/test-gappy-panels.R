# Panels with holes and repeats must match DIDmultiplegtDYN exactly.
#
# Regression, found on the annual tax panels of a real application (which
# drop loss-making years, leaving gaps in 40-50% of firms' histories):
#
#   tric_tax, only_never_switchers = TRUE   Effect_2 -0.0156 vs -0.0145
#   ntr_tax,  only_never_switchers = TRUE   effects off by up to 5e-4,
#                                           placebos by up to 7e-4
#
# while every balanced-panel parity test passed. Three things were
# missing from .prep_panel, all of which the reference does before it
# balances the panel (did_multiplegt_main.R:119-300):
#
#   1. Collapsing un-aggregated data. A repeated (group, time) becomes one
#      cell with weighted-mean outcome/treatment and N_gt = the summed
#      weight. didgpu kept both rows, double-counting them and -- since
#      lags are taken by row -- misaligning that group's differences.
#   2. The missing-treatment rules. A group whose last observation before
#      its switch is not the period right before it has an unknown switch
#      date; the reference demotes it to a control truncated at its last
#      clean period. Under only_never_switchers = TRUE that turns it INTO
#      a never-switcher, which is why the gap showed there most.
#   3. The sample restrictions and T_g, in order: drop cohorts with no
#      variation in switch dates, fix G, drop (period, cohort) cells with
#      no control, and take T_g per cohort -- not per group.
#
# The simulators here are built to hit each of those.

.gp_panel <- function(seed, nu = 160L, np = 18L, gap_frac = 0.25,
                      dup_n = 0L, multi = FALSE, onoff = TRUE) {
  set.seed(seed)
  ufe <- stats::rnorm(nu, 0, 1); tfe <- stats::rnorm(np, 0, 0.3)
  Fg <- rep(Inf, nu)
  tr <- sort(sample(seq_len(nu), round(nu * 0.55)))
  Fg[tr] <- sample(3:(np - 3), length(tr), replace = TRUE)
  off <- rep(Inf, nu)
  if (onoff) off[tr] <- Fg[tr] + sample(c(2:6, Inf), length(tr), TRUE)
  dose <- if (multi) sample(1:3, nu, TRUE) else rep(1L, nu)
  # A few groups start treated, so a second baseline cohort exists.
  base1 <- sample(setdiff(seq_len(nu), tr), 8L)
  g <- expand.grid(t = seq_len(np), unit = seq_len(nu))
  g <- g[order(g$unit, g$t), ]
  on <- g$t >= Fg[g$unit] & g$t < off[g$unit]
  g$D <- ifelse(on, dose[g$unit], 0L)
  g$D[g$unit %in% base1] <- 1L
  g$Y <- ufe[g$unit] + tfe[g$t] + 0.3 * g$D + stats::rnorm(nrow(g), 0, 0.4)
  # Holes: drop a share of rows at random, never the first period, so a
  # switch is often preceded by a gap.
  drop <- which(g$t > 1L & stats::runif(nrow(g)) < gap_frac)
  if (length(drop)) g <- g[-drop, ]   # g[-integer(0), ] would drop EVERY row
  if (dup_n > 0L) {
    dup <- g[sample(nrow(g), dup_n), ]
    dup$Y <- dup$Y + stats::rnorm(dup_n, 0, 0.2)
    g <- rbind(g, dup)
  }
  g$cl <- ((g$unit - 1L) %/% 4L) + 1L
  g[order(g$unit, g$t), c("unit", "t", "D", "Y", "cl")]
}

.gp_ours <- function(d, bk = "r", ...) {
  suppressMessages(suppressWarnings(
    didgpu(df = d, outcome = "Y", group = "unit", time = "t", treatment = "D",
           backend = bk, verbose = FALSE, ...)))
}
.gp_ref <- function(d, ...) {
  if (!"package:polars" %in% search()) suppressMessages(attachNamespace("polars"))
  suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = d, outcome = "Y", group = "unit", time = "t", treatment = "D",
      graph_off = TRUE, ...)))
}
.gp_expect_same <- function(o, r, label) {
  ce <- function(m, j) as.numeric(m[, j])
  expect_equal(ce(o$results$Effects, 1), ce(r$results$Effects, 1),
               tolerance = 1e-10, label = paste(label, "effects"))
  expect_equal(ce(o$results$Effects, "SE"), ce(r$results$Effects, 2),
               tolerance = 1e-10, label = paste(label, "effect SEs"))
  if (!is.null(r$results$Placebos) && nrow(r$results$Placebos) > 0L) {
    expect_equal(ce(o$results$Placebos, 1), ce(r$results$Placebos, 1),
                 tolerance = 1e-10, label = paste(label, "placebos"))
  }
  expect_equal(o$results$ATE[1L, 1L], as.numeric(r$results$ATE[1L, 1L]),
               tolerance = 1e-10, label = paste(label, "ATE"))
  # The sample-size columns: N and Switchers. These are what exposed the
  # collapse bug -- didgpu counted duplicated firm-years twice.
  expect_equal(ce(o$results$Effects, "N"), ce(r$results$Effects, 5),
               label = paste(label, "N"))
  expect_equal(ce(o$results$Effects, "Switchers"), ce(r$results$Effects, 6),
               label = paste(label, "Switchers"))
}

test_that("gappy panels match, with and without only_never_switchers", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  for (seed in c(3L, 11L)) {
    d <- .gp_panel(seed)
    for (ons in c(FALSE, TRUE)) {
      .gp_expect_same(.gp_ours(d, effects = 4L, placebo = 2L, cluster = "cl",
                               only_never_switchers = ons),
                      .gp_ref(d, effects = 4, placebo = 2, cluster = "cl",
                              only_never_switchers = ons),
                      sprintf("seed %d, only_never_switchers = %s", seed, ons))
    }
  }
})

test_that("a switch right after a gap is demoted to a truncated control", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  # Guard: the simulator must actually produce ambiguous switch dates,
  # or this test proves nothing.
  d <- .gp_panel(7L, gap_frac = 0.35)
  dt <- data.table::as.data.table(d)[order(unit, t)]
  dt[, d_sq := D[1L], by = unit]
  dt[, Fg := { x <- t[D != d_sq]; if (length(x)) min(x) else NA_integer_ }, by = unit]
  dt[, last_obs := { x <- t[t < Fg[1L]]; if (length(x)) max(x) else NA_integer_ }, by = unit]
  n_amb <- dt[, .(amb = !is.na(Fg[1L]) && last_obs[1L] < Fg[1L] - 1L), by = unit][, sum(amb)]
  expect_gt(n_amb, 5L)
  .gp_expect_same(.gp_ours(d, effects = 3L, placebo = 1L,
                           only_never_switchers = TRUE),
                  .gp_ref(d, effects = 3, placebo = 1,
                          only_never_switchers = TRUE),
                  "ambiguous switch dates")
})

test_that("duplicated (group, time) rows are collapsed into one cell", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  d <- .gp_panel(5L, gap_frac = 0.15, dup_n = 25L)
  expect_gt(anyDuplicated(d[, c("unit", "t")]), 0L)
  for (ons in c(FALSE, TRUE)) {
    .gp_expect_same(.gp_ours(d, effects = 4L, placebo = 2L, cluster = "cl",
                             only_never_switchers = ons),
                    .gp_ref(d, effects = 4, placebo = 2, cluster = "cl",
                            only_never_switchers = ons),
                    sprintf("duplicates, only_never_switchers = %s", ons))
  }
})

test_that("multivalued, normalized, on a gappy panel with duplicates", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  d <- .gp_panel(9L, multi = TRUE, dup_n = 15L)
  .gp_expect_same(.gp_ours(d, effects = 4L, placebo = 2L, cluster = "cl",
                           only_never_switchers = TRUE, normalized = TRUE),
                  .gp_ref(d, effects = 4, placebo = 2, cluster = "cl",
                          only_never_switchers = TRUE, normalized = TRUE),
                  "multivalued normalized")
})

test_that("every backend agrees on a gappy panel", {
  # The fast kernels read the same prepped panel, so the gap rules must
  # reach them too.
  d <- .gp_panel(3L)
  ref <- .gp_ours(d, "r", effects = 4L, placebo = 2L)
  bks <- "cpu"
  if (isTRUE(tryCatch(didgpu_has_cuda_support(), error = function(e) FALSE)))
    bks <- c(bks, "cuda")
  for (bk in bks) {
    o <- .gp_ours(d, bk, effects = 4L, placebo = 2L)
    expect_equal(as.numeric(o$results$Effects[, 1]),
                 as.numeric(ref$results$Effects[, 1]), tolerance = 1e-12,
                 label = paste("effects on", bk))
    expect_equal(o$results$ATE[1L, 1L], ref$results$ATE[1L, 1L],
                 tolerance = 1e-12, label = paste("ATE on", bk))
  }
})
