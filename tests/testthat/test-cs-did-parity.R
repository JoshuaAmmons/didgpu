# Bit-for-bit parity with the `did` package: ATT(g, t), its SE, and the
# dynamic (event-study) aggregation with ITS SE.
#
# didgpu_cs() delegates nothing -- it reimplements the estimator -- so the
# only way to keep the equivalence claim honest is to check it against the
# oracle. The per-cell estimators mirror DRDID function for function,
# because that is what did::att_gt() calls:
#     est_method = "reg" -> DRDID::reg_did_panel
#     est_method = "ipw" -> DRDID::std_ipw_did_panel
#     est_method = "dr"  -> DRDID::drdid_panel
#
# Why this file exists. The influence functions used to omit both the
# nuisance-estimation terms and the correct normalisers, and the
# aggregation had no SEs at all. Against did on a 200-unit panel the
# multiplier-bootstrap SEs came out at ~0.13x (OR) and ~0.50x (IPW/DR) of
# the truth -- confidence intervals two to eight times too narrow -- while
# every point estimate matched exactly. Point-estimate parity alone would
# not have caught it, so this file pins the SEs too.
#
# NOTE on pre-treatment cells: didgpu uses a universal base period and did
# defaults to a varying one, so ATT(g, t) for t < g legitimately differs.
# Parity is asserted on post-treatment cells and non-negative event times.

skip_if_no_did <- function() testthat::skip_if_not_installed("did")

.parity_panel <- function(seed = 17L, n_units = 200L, with_X = FALSE) {
  set.seed(seed)
  n_periods <- 12L
  ufe <- stats::rnorm(n_units, 0, 1)
  tfe <- stats::rnorm(n_periods, 0, 0.3)
  Fg <- rep(Inf, n_units)
  tr <- sort(sample(seq_len(n_units), as.integer(n_units * 0.6)))
  Fg[tr] <- sample(seq(3L, n_periods - 1L), length(tr), replace = TRUE)
  p <- expand.grid(unit = seq_len(n_units), period = seq_len(n_periods))
  p$D <- as.integer(p$period >= Fg[p$unit])
  p$Y <- ufe[p$unit] + tfe[p$period] + 1.0 * p$D +
         stats::rnorm(nrow(p), 0, 0.4)
  # did wants gname numeric with 0 for never-treated. An INTEGER column
  # makes did coerce Inf to NA internally, silently destroying the
  # never-treated group -- keep this double.
  p$G <- as.numeric(ifelse(is.finite(Fg[p$unit]), Fg[p$unit], 0))
  if (with_X) {
    xu <- stats::rnorm(n_units)
    p$x1 <- xu[p$unit]
    p$Y <- p$Y + 0.5 * p$x1
  }
  p[order(p$unit, p$period), ]
}

.did_map <- c(OR = "reg", IPW = "ipw", DR = "dr")

.fit_pair <- function(p, em, with_X) {
  xf <- if (with_X) ~x1 else ~1
  ref <- suppressMessages(suppressWarnings(did::att_gt(
    yname = "Y", tname = "period", idname = "unit", gname = "G",
    xformla = xf, data = p, control_group = "nevertreated",
    est_method = .did_map[[em]], bstrap = FALSE, cband = FALSE)))
  ours <- didgpu_cs(df = p, outcome = "Y", group = "unit", time = "period",
                    treatment = "D", control_group = "never",
                    est_method = em,
                    covariates = if (with_X) "x1" else NULL,
                    aggregation = "event", bootstrap_reps = 0L,
                    backend = "r", verbose = FALSE)
  list(ref = ref, ours = ours)
}

# Per-cell SE implied by the influence function, in DRDID's convention.
.cell_se <- function(fit) {
  IF <- attr(fit$att_gt, "IF_per_cell")
  vapply(seq_along(IF), function(k) {
    cl <- IF[[k]]
    if (is.null(cl) || !length(cl$IF)) return(NA_real_)
    nn <- length(cl$IF)
    stats::sd(cl$IF) * sqrt(nn - 1) / nn
  }, numeric(1))
}

test_that("ATT(g,t) and its SE match did::att_gt for every est_method", {
  skip_if_no_did()
  for (with_X in c(FALSE, TRUE)) {
    p <- .parity_panel(with_X = with_X)
    for (em in c("OR", "IPW", "DR")) {
      fp <- .fit_pair(p, em, with_X)
      o <- data.frame(g = fp$ours$att_gt$g, t = fp$ours$att_gt$t,
                      att = fp$ours$att_gt$att, se = .cell_se(fp$ours))
      d <- data.frame(g = fp$ref$group, t = fp$ref$t,
                      att_d = fp$ref$att, se_d = fp$ref$se)
      m <- merge(o, d, by = c("g", "t"))
      m <- m[m$t >= m$g & is.finite(m$se) & is.finite(m$se_d), ]
      expect_gt(nrow(m), 20L)
      lbl <- paste(em, "X =", with_X)
      expect_lt(max(abs(m$att - m$att_d)), 1e-10, label = paste("ATT", lbl))
      expect_lt(max(abs(m$se  - m$se_d)),  1e-10, label = paste("SE",  lbl))
    }
  }
})

test_that("event-study aggregation and its SE match did::aggte(dynamic)", {
  skip_if_no_did()
  for (with_X in c(FALSE, TRUE)) {
    p <- .parity_panel(with_X = with_X)
    for (em in c("OR", "IPW", "DR")) {
      fp <- .fit_pair(p, em, with_X)
      ag <- suppressMessages(suppressWarnings(
        did::aggte(fp$ref, type = "dynamic", bstrap = FALSE, cband = FALSE)))
      d1 <- data.frame(event_time = ag$egt, est_d = ag$att.egt,
                       se_d = ag$se.egt)
      m <- merge(fp$ours$aggregation, d1, by = "event_time")
      m <- m[m$event_time >= 0 & is.finite(m$se) & is.finite(m$se_d), ]
      expect_gt(nrow(m), 5L)
      lbl <- paste(em, "X =", with_X)
      expect_lt(max(abs(m$estimate - m$est_d)), 1e-10,
                label = paste("agg est", lbl))
      expect_lt(max(abs(m$se - m$se_d)), 1e-10,
                label = paste("agg SE", lbl))
    }
  }
})

test_that("the aggregation table carries SEs and CIs", {
  # Regression: .cs_aggregate() used to propagate point estimates only, so
  # every aggregate SE was absent and didgpu_tidy() reported NA for them.
  p <- .parity_panel()
  f <- didgpu_cs(df = p, outcome = "Y", group = "unit", time = "period",
                 treatment = "D", aggregation = "event",
                 bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_true(all(c("se", "ci_low", "ci_high") %in% names(f$aggregation)))
  expect_true(any(is.finite(f$aggregation$se)))
  fin <- f$aggregation[is.finite(f$aggregation$se), ]
  expect_true(all(fin$ci_low < fin$estimate & fin$estimate < fin$ci_high))
})

test_that("multiplier and cluster bootstrap SEs now agree with each other", {
  # The multiplier bootstrap consumes the influence functions; the cluster
  # bootstrap recomputes from scratch. They disagreed by a factor of two
  # to eight while the IFs were wrong.
  p <- .parity_panel(n_units = 120L)
  f1 <- didgpu_cs(df = p, outcome = "Y", group = "unit", time = "period",
                  treatment = "D", est_method = "DR", bootstrap_reps = 400L,
                  bootstrap_kind = "multiplier", seed = 1, backend = "r",
                  verbose = FALSE)
  f2 <- didgpu_cs(df = p, outcome = "Y", group = "unit", time = "period",
                  treatment = "D", est_method = "DR", bootstrap_reps = 400L,
                  bootstrap_kind = "cluster", seed = 1, backend = "r",
                  verbose = FALSE)
  m <- merge(data.frame(g = f1$att_gt$g, t = f1$att_gt$t, se1 = f1$att_gt$se),
             data.frame(g = f2$att_gt$g, t = f2$att_gt$t, se2 = f2$att_gt$se),
             by = c("g", "t"))
  m <- m[m$t >= m$g & is.finite(m$se1) & is.finite(m$se2), ]
  expect_gt(nrow(m), 10L)
  # Both are noisy at B = 400, so this is a sanity band, not equality.
  expect_lt(abs(stats::median(m$se1 / m$se2) - 1), 0.15)
})
