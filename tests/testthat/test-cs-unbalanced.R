# Callaway-Sant'Anna on unbalanced panels, and the base period, must
# match did::att_gt + aggte.
#
# Regressions:
#   * On an unbalanced panel didgpu_cs() used, cell by cell, whatever
#     units happened to be observed in both periods -- matching NEITHER of
#     did's modes. did either balances the panel first (the default,
#     allow_unbalanced_panel = FALSE) or switches to DRDID's
#     repeated-cross-section estimators (allow_unbalanced_panel = TRUE).
#     On a real application's non-trade sales-growth panel the overall
#     ATT was -0.060 against did's -0.096.
#   * didgpu_cs() always used the universal base period (g - 1 for every
#     cell) while did defaults to "varying", so every pre-treatment cell
#     of the event study differed under defaults.
#   * With control_group = "notyet" and a universal base, did picks
#     not-yet-treated controls by the LATER of the two periods compared;
#     didgpu used the current one.
#   * Cohorts were always inferred from the treatment column. When a
#     unit's adoption-period row is missing it lands in a later cohort;
#     did reads the cohort from gname. `first_treat` now takes it.
#
# did 2.3.0 has two bugs this file steers around, both fixed in 2.5.0:
#   * an INTEGER gname silently removes the never-treated group
#     (bcallaway11/did#264) -- the simulators here store it as double;
#   * with allow_unbalanced_panel = TRUE its influence functions land on
#     the wrong units, so every SE that combines cells depends on how the
#     units happen to be numbered. Those SEs are only compared against
#     did >= 2.5.0. Point estimates and per-cell SEs are unaffected and
#     are compared against any version.

.cu_panel <- function(seed, nu = 300L, np = 10L, drop = 0) {
  set.seed(seed)
  g <- sample(c(0, 4:8), nu, TRUE, prob = c(0.35, rep(0.13, 5)))
  d <- data.table::CJ(id = seq_len(nu), t = seq_len(np))
  d[, gv := as.numeric(g[id])]            # double: see did#264 above
  d[, D := as.integer(gv > 0 & t >= gv)]
  ufe <- stats::rnorm(nu); tfe <- stats::rnorm(np, 0, 0.3)
  d[, y := ufe[id] + tfe[t] + 0.5 * D * (1 + 0.1 * (t - gv)) +
         stats::rnorm(.N, 0, 0.5)]
  if (drop > 0) d <- d[stats::runif(.N) > drop | t == 1L]
  as.data.frame(d)
}
.cu_ours <- function(d, m, bp, unb, cg = "never") {
  suppressMessages(suppressWarnings(didgpu_cs(
    d, "y", "id", "t", "D", control_group = cg, est_method = m,
    aggregation = "event", base_period = bp, allow_unbalanced_panel = unb,
    first_treat = "gv", backend = "r", verbose = FALSE)))
}
.cu_ref <- function(d, m, bp, unb, cg = "never") {
  suppressMessages(suppressWarnings(did::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "gv", data = d,
    control_group = c(never = "nevertreated", notyet = "notyettreated")[[cg]],
    est_method = c(OR = "reg", IPW = "ipw", DR = "dr")[[m]],
    base_period = bp, allow_unbalanced_panel = unb, bstrap = FALSE,
    cband = FALSE, print_details = FALSE)))
}
.cu_check <- function(d, m, bp, unb, cg = "never", combined_se = TRUE) {
  lbl <- sprintf("%s, %s, base %s, unbalanced %s", m, cg, bp, unb)
  ours <- .cu_ours(d, m, bp, unb, cg); ref <- .cu_ref(d, m, bp, unb, cg)
  a <- data.table::as.data.table(ours$att_gt)
  j <- merge(data.table::data.table(g = ref$group, t = ref$t, r = ref$att,
                                    r_se = ref$se),
             a[, list(g, t, att)], by = c("g", "t"), all = TRUE)
  expect_equal(is.na(j$att), is.na(j$r), label = paste(lbl, "NA cells"))
  expect_equal(j$att, j$r, tolerance = 1e-10, label = paste(lbl, "ATT(g,t)"))
  ev_r <- suppressWarnings(did::aggte(ref, type = "dynamic", na.rm = TRUE, bstrap = FALSE))
  s_r  <- suppressWarnings(did::aggte(ref, type = "simple",  na.rm = TRUE, bstrap = FALSE))
  ev_o <- data.table::as.data.table(ours$aggregation)
  s_o  <- data.table::as.data.table(didgpu_cs_aggregate(ours, "overall")$aggregation)
  ej <- merge(data.table::data.table(e = ev_r$egt, est = ev_r$att.egt, se = ev_r$se.egt),
              ev_o[, list(e = event_time, est_o = estimate, se_o = se)], by = "e")
  expect_equal(nrow(ej), length(ev_r$egt), label = paste(lbl, "event rows"))
  expect_equal(ej$est_o, ej$est, tolerance = 1e-10, label = paste(lbl, "event est"))
  expect_equal(s_o$estimate, s_r$overall.att, tolerance = 1e-10,
               label = paste(lbl, "overall ATT"))
  if (combined_se) {
    expect_equal(ej$se_o, ej$se, tolerance = 1e-10, label = paste(lbl, "event SE"))
    expect_equal(s_o$se, s_r$overall.se, tolerance = 1e-10,
                 label = paste(lbl, "overall SE"))
  }
}
.did_fixed_rc <- function() {
  requireNamespace("did", quietly = TRUE) &&
    utils::packageVersion("did") >= "2.5.0"
}

test_that("balanced panels match did under both base periods", {
  skip_if_not_installed("did")
  d <- .cu_panel(1)
  for (m in c("OR", "IPW", "DR")) for (bp in c("varying", "universal")) {
    .cu_check(d, m, bp, unb = FALSE)
  }
})

test_that("the not-yet-treated control group matches did, incl. universal base", {
  skip_if_not_installed("did")
  d <- .cu_panel(4)
  for (bp in c("varying", "universal")) .cu_check(d, "OR", bp, unb = FALSE, cg = "notyet")
})

test_that("unbalanced panel, default: balanced first, as did does", {
  skip_if_not_installed("did")
  d <- .cu_panel(2, drop = 0.2)
  for (m in c("OR", "IPW", "DR")) .cu_check(d, m, "varying", unb = FALSE)
})

test_that("allow_unbalanced_panel = TRUE: did's repeated-cross-section path", {
  skip_if_not_installed("did")
  d <- .cu_panel(2, drop = 0.2)
  fixed <- .did_fixed_rc()
  for (m in c("OR", "IPW", "DR")) for (bp in c("varying", "universal")) {
    .cu_check(d, m, bp, unb = TRUE, combined_se = fixed)
  }
})

test_that("repeated-cross-section SEs do not depend on unit numbering", {
  # The property did 2.3.0 lacked. Holds for didgpu regardless of the did
  # version installed.
  d <- .cu_panel(2, drop = 0.2)
  relab <- d
  set.seed(9); perm <- sample(unique(d$id))
  relab$id <- perm[match(d$id, sort(unique(d$id)))]
  se_of <- function(x) {
    f <- .cu_ours(x, "OR", "varying", TRUE)
    c(data.table::as.data.table(didgpu_cs_aggregate(f, "overall")$aggregation)$se,
      f$aggregation$se)
  }
  expect_equal(se_of(relab), se_of(d), tolerance = 1e-12)
})

test_that("without first_treat, a missing adoption row is warned about", {
  d <- .cu_panel(2, drop = 0.2)
  expect_warning(
    didgpu_cs(d, "y", "id", "t", "D", control_group = "never",
              est_method = "OR", allow_unbalanced_panel = TRUE,
              backend = "r", verbose = FALSE),
    "first_treat")
})

test_that("the default base period is did's 'varying'", {
  expect_identical(eval(formals(didgpu_cs)$base_period)[1], "varying")
  expect_false(eval(formals(didgpu_cs)$allow_unbalanced_panel))
})
