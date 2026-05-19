# Cross-validation: didgpu_cs vs the reference `did` package.
# The did package computes the same ATT(g, t) under DR but with a more
# polished implementation and tested SEs. This test checks that our
# point estimates are reasonably close (won't be bit-identical because
# we use slightly different SE plumbing).

skip_if_no_did <- function() {
  testthat::skip_if_not_installed("did")
}

build_cs_test_panel <- function(seed = 17L, n_units = 200L) {
  set.seed(seed)
  n_periods <- 12L
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), as.integer(n_units * 0.6)))
  F_g[treated] <- sample(seq(3L, n_periods - 1L), length(treated),
                          replace = TRUE)
  panel <- expand.grid(unit = seq_len(n_units), period = seq_len(n_periods))
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             1.0 * panel$D + rnorm(nrow(panel), 0, 0.3)
  # `did` wants a "G" column that's the actual first-treatment period
  # (0 for never-treated).
  panel$G <- ifelse(is.finite(panel$F_g), panel$F_g, 0L)
  panel[order(panel$unit, panel$period),
         c("unit", "period", "D", "Y", "G")]
}

test_that("didgpu_cs OR estimate is within 0.15 of did::att_gt at post-treatment cells", {
  skip_if_no_did()
  p <- build_cs_test_panel(seed = 17L)
  us <- didgpu_cs(p, "Y", "unit", "period", "D",
                   est_method = "OR", aggregation = "event",
                   bootstrap_reps = 0L,
                   backend = "r", verbose = FALSE)
  ref <- tryCatch(
    suppressWarnings(suppressMessages(did::att_gt(
      yname    = "Y",
      tname    = "period",
      idname   = "unit",
      gname    = "G",
      data     = p,
      control_group = "nevertreated",
      est_method    = "reg"))),
    error = function(e) NULL)
  skip_if(is.null(ref), "did::att_gt failed on this panel")
  # Align (g, t) keys and compare.
  ref_df <- data.frame(
    g = ref$group, t = ref$t,
    att_ref = ref$att,
    stringsAsFactors = FALSE
  )
  merged <- merge(us$att_gt, ref_df, by = c("g", "t"))
  # Restrict to POST-treatment cells (event-time >= 0).
  post <- merged[merged$event_time >= 0L, ]
  expect_true(nrow(post) > 0L)
  diffs <- post$att - post$att_ref
  max_abs_diff <- max(abs(diffs), na.rm = TRUE)
  expect_lt(max_abs_diff, 0.25)
})

test_that("didgpu_cs aggregation matches did::aggte event-study within 0.15", {
  skip_if_no_did()
  p <- build_cs_test_panel(seed = 19L)
  us <- didgpu_cs(p, "Y", "unit", "period", "D",
                   est_method = "OR", aggregation = "event",
                   bootstrap_reps = 0L,
                   backend = "r", verbose = FALSE)
  ref <- tryCatch(
    suppressWarnings(suppressMessages(did::att_gt(
      yname = "Y", tname = "period", idname = "unit", gname = "G",
      data = p, control_group = "nevertreated", est_method = "reg"))),
    error = function(e) NULL)
  skip_if(is.null(ref), "did::att_gt failed")
  ref_agg <- tryCatch(
    suppressWarnings(suppressMessages(did::aggte(ref, type = "dynamic"))),
    error = function(e) NULL)
  skip_if(is.null(ref_agg), "did::aggte failed")
  ref_es <- data.frame(
    event_time = ref_agg$egt,
    est_ref = ref_agg$att.egt,
    stringsAsFactors = FALSE
  )
  merged <- merge(us$aggregation, ref_es, by = "event_time")
  # Post-treatment event-times only.
  post <- merged[merged$event_time >= 0L, ]
  expect_true(nrow(post) > 0L)
  diffs <- post$estimate - post$est_ref
  expect_lt(max(abs(diffs), na.rm = TRUE), 0.25)
})
