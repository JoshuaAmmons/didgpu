# didgpu() must behave as DIDmultiplegtDYN::did_multiplegt_dyn() does: the
# same errors, the same messages, the same printed tables, the same result
# object, and the same numbers for every option the two share.

skip_if_not_installed("DIDmultiplegtDYN")
skip_if_not_installed("polars")
suppressMessages(library(polars))

.ux_panel <- function(G = 80, T = 9, seed = 11, multi = FALSE) {
  set.seed(seed)
  F <- sample(c(3:7, 99), G, replace = TRUE)
  base <- if (multi) sample(0:1, G, replace = TRUE) else rep(0, G)
  dose <- if (multi) sample(1:2, G, replace = TRUE) else rep(1, G)
  df <- data.frame(g = rep(seq_len(G), each = T), t = rep(seq_len(T), G))
  df$d <- ifelse(df$t >= F[df$g], base[df$g] + dose[df$g], base[df$g])
  df$x1 <- stats::rnorm(nrow(df))
  df$y <- stats::rnorm(nrow(df)) + 0.4 * df$d + 0.3 * df$x1
  df$het <- rep(stats::rnorm(G), each = T)
  df$cl <- (df$g - 1) %/% 4 + 1
  df$reg <- rep(sample(1:2, G, replace = TRUE), each = T)
  df
}
.capture <- function(f) {
  out <- character(0)
  r <- tryCatch(withCallingHandlers(f(),
         message = function(m) {
           out <<- c(out, paste0("M: ", sub("\n$", "", conditionMessage(m))))
           invokeRestart("muffleMessage")
         },
         warning = function(w) invokeRestart("muffleWarning")),
       error = function(e) { out <<- c(out, paste0("E: ", conditionMessage(e))); NULL })
  list(r = r, out = out)
}
.ref <- function(df, ...) .capture(function()
  DIDmultiplegtDYN::did_multiplegt_dyn(df, "y", "g", "t", "d", graph_off = TRUE, ...))
.gpu <- function(df, ..., backend = "r") .capture(function()
  didgpu(df, "y", "g", "t", "d", backend = backend, graph_off = TRUE, ...))
# The reference closes its print with two lines acknowledging its funding;
# didgpu leaves those out.
.print_ref <- function(x) utils::head(utils::capture.output(print(x)), -3)
.print_gpu <- function(x) utils::head(utils::capture.output(print(x)), -1)

test_that("nothing estimable: the reference's error, nothing else", {
  df <- .ux_panel(); df$d <- 0
  a <- .ref(df, effects = 2); g <- .gpu(df, effects = 2)
  expect_identical(g$out, a$out)
  expect_match(g$out[length(g$out)], "No treatment effect can be estimated", fixed = TRUE)
})

test_that("messages and printed tables match across options and backends", {
  df <- .ux_panel()
  specs <- list(list(effects = 4, placebo = 3), list(effects = 9, placebo = 7),
                list(effects = 2, placebo = 4), list(effects = 3, placebo = 2, controls = "x1"),
                list(effects = 3, placebo = 2, cluster = "cl"),
                list(effects = 3, placebo = 2, predict_het = list("het", -1)))
  for (sp in specs) {
    a <- do.call(.ref, c(list(df), sp))
    for (b in c("r", "cpu")) {
      g <- do.call(.gpu, c(list(df), sp, list(backend = b)))
      expect_identical(g$out, a$out, label = paste(b, deparse(sp)))
      expect_identical(.print_gpu(g$r), .print_ref(a$r), label = paste(b, deparse(sp)))
    }
  }
})

test_that("the result object carries the reference's fields", {
  df <- .ux_panel()
  a <- .ref(df, effects = 3, placebo = 2)$r
  g <- .gpu(df, effects = 3, placebo = 2)$r
  for (f in c("N_Effects", "N_Placebos", "delta_D_avg_total", "max_pl", "max_pl_gap",
              "p_jointeffects", "p_jointplacebo")) {
    expect_equal(g$results[[f]], a$results[[f]], tolerance = 1e-10, label = f)
  }
  for (m in c("Effects", "ATE", "Placebos")) {
    expect_identical(dimnames(g$results[[m]]), dimnames(a$results[[m]]))
    expect_equal(unname(g$results[[m]]), unname(a$results[[m]]), tolerance = 1e-10)
  }
  expect_identical(names(g$coef$b), names(a$coef$b))
})

test_that("each effect's N counts a control serving both directions once", {
  set.seed(5); G <- 90; T <- 8
  F <- sample(c(3:7, 99), G, replace = TRUE); to <- sample(c(0, 2), G, replace = TRUE)
  df <- data.frame(g = rep(seq_len(G), each = T), t = rep(seq_len(T), G))
  df$d <- ifelse(df$t >= F[df$g], to[df$g], 1)
  df$y <- stats::rnorm(nrow(df))
  a <- .ref(df, effects = 3)$r
  for (b in c("r", "cpu")) {
    g <- .gpu(df, effects = 3, backend = b)$r
    expect_equal(unname(g$results$Effects[, 5:8]), unname(a$results$Effects[, 5:8]))
    expect_equal(unname(g$results$ATE[, 5:8]), unname(a$results$ATE[, 5:8]))
  }
})

test_that("a missing control drops the row, as in the reference", {
  df <- .ux_panel(); df$x1[c(5, 40, 77, 300)] <- NA
  a <- .ref(df, effects = 3, controls = "x1")
  g <- .gpu(df, effects = 3, controls = "x1")
  expect_equal(unname(g$r$results$Effects[, 1:2]), unname(a$r$results$Effects[, 1:2]),
               tolerance = 1e-10)
  expect_identical(g$out, a$out)
})

test_that("the remaining reference options give the reference's numbers", {
  df <- .ux_panel(multi = TRUE)
  specs <- list(list(less_conservative_se = TRUE), list(more_granular_demeaning = TRUE),
                list(drop_if_d_miss_before_first_switch = TRUE),
                list(effects_equal = TRUE), list(effects_equal = "2, 3"),
                list(normalized = TRUE, normalized_weights = TRUE),
                list(predict_het = list("het", -1), predict_het_hc2bm = TRUE))
  for (sp in specs) {
    a <- do.call(.ref, c(list(df, effects = 3, placebo = 1), sp))
    g <- do.call(.gpu, c(list(df, effects = 3, placebo = 1), sp))
    lab <- deparse(sp)
    expect_equal(unname(g$r$results$Effects[, 1:2]), unname(a$r$results$Effects[, 1:2]),
                 tolerance = 1e-10, label = lab)
    expect_equal(g$r$results$p_equality_effects, a$r$results$p_equality_effects,
                 tolerance = 1e-10, label = lab)
    expect_identical(g$r$normalized_weights, a$r$normalized_weights, label = lab)
    if (!is.null(a$r$results$predict_het)) {
      expect_equal(as.matrix(g$r$results$predict_het[, 3:9]),
                   as.matrix(a$r$results$predict_het[, 3:9]), tolerance = 1e-10,
                   check.attributes = FALSE)
    }
    expect_identical(.print_gpu(g$r), .print_ref(a$r), label = lab)
  }
})

test_that("save_sample, design and date_first_switch are the reference's", {
  df <- .ux_panel(multi = TRUE)
  a <- .ref(df, effects = 2, save_sample = TRUE, design = c(0.8, "console"),
            date_first_switch = c("by_baseline_treat", "console"))$r
  g <- .gpu(df, effects = 2, save_sample = TRUE, design = c(0.8, "console"),
            date_first_switch = c("by_baseline_treat", "console"))$r
  expect_identical(g$save_sample, a$save_sample)
  expect_identical(g$design, a$design)
  expect_identical(g$date_first_switch, a$date_first_switch)
  expect_identical(.print_gpu(g), .print_ref(a))
})

test_that("by and by_path run per level as the reference does", {
  df <- .ux_panel(multi = TRUE)
  for (sp in list(list(by = "reg"), list(by_path = 2))) {
    a <- do.call(.ref, c(list(df, effects = 2), sp))
    g <- do.call(.gpu, c(list(df, effects = 2), sp))
    expect_identical(g$out, a$out)
    expect_identical(g$r$by_levels, a$r$by_levels)
    expect_identical(.print_gpu(g$r), .print_ref(a$r))
  }
})

test_that("trends_lin matches, except the reference's stale placebo Switchers", {
  # Under trends_lin DIDmultiplegtDYN fills every placebo row's
  # (unweighted) Switchers column with the LAST placebo's count -- a value
  # left over from its last per-placebo run; its own Switchers.w column,
  # which equals Switchers when there are no weights, has the right
  # counts. didgpu reports those, so that one column differs and nothing
  # else does.
  df <- .ux_panel()
  a <- .ref(df, effects = 3, placebo = 2, trends_lin = TRUE)
  g <- .gpu(df, effects = 3, placebo = 2, trends_lin = TRUE)
  expect_identical(g$out, a$out)
  expect_equal(unname(g$r$results$Effects), unname(a$r$results$Effects), tolerance = 1e-10)
  expect_equal(unname(g$r$results$Placebos[, -6]), unname(a$r$results$Placebos[, -6]),
               tolerance = 1e-10)
  expect_equal(unname(g$r$results$Placebos[, "Switchers"]),
               unname(a$r$results$Placebos[, "Switchers.w"]))
  expect_equal(g$r$results$p_jointplacebo, a$r$results$p_jointplacebo, tolerance = 1e-10)
})

test_that("clustered SEs count a group whose first period is missing", {
  # The cluster of a group was read from its first row; for a late
  # entrant that row is a balancing fill-in with no cluster, so the group
  # dropped out of the clustered variance.
  df <- .ux_panel(G = 100)
  late <- df$g %% 6 == 0 & df$t <= 2
  df <- df[!late, ]
  a <- .ref(df, effects = 3, placebo = 1, cluster = "cl")
  for (b in c("r", "cpu")) {
    g <- .gpu(df, effects = 3, placebo = 1, cluster = "cl", backend = b)
    expect_equal(unname(g$r$results$Effects[, 1:2]), unname(a$r$results$Effects[, 1:2]),
                 tolerance = 1e-10, label = b)
    expect_equal(unname(g$r$results$Placebos[, 1:2]), unname(a$r$results$Placebos[, 1:2]),
                 tolerance = 1e-10, label = b)
  }
})

test_that("reset and avg_time_periods match DIDmultiplegtDYN 2.4.0", {
  skip_if(utils::packageVersion("DIDmultiplegtDYN") < "2.4.0",
          "reset and avg_time_periods need DIDmultiplegtDYN >= 2.4.0")
  set.seed(21); G <- 100; T <- 12
  F <- sample(c(3:9, 99), G, replace = TRUE)
  df <- data.frame(g = rep(seq_len(G), each = T), t = rep(seq_len(T), G))
  off <- F + sample(2:4, G, TRUE)
  df$d <- as.numeric(df$t >= F[df$g] & df$t < off[df$g])
  df$y <- stats::rnorm(nrow(df)) + 0.4 * df$d
  for (sp in list(list(reset = 2), list(reset = 3, avg_time_periods = TRUE),
                  list(avg_time_periods = TRUE))) {
    a <- do.call(.ref, c(list(df, effects = 4, placebo = 2), sp))
    g <- do.call(.gpu, c(list(df, effects = 4, placebo = 2), sp))
    lab <- deparse(sp)
    expect_identical(g$out, a$out, label = lab)
    expect_equal(unname(g$r$results$Effects), unname(a$r$results$Effects),
                 tolerance = 1e-10, label = lab)
    expect_equal(g$r$avg_cumul, a$r$avg_cumul, tolerance = 1e-12, label = lab)
    expect_true(all(setdiff(names(a$r), "args") %in% names(g$r)), label = lab)
    expect_identical(.print_gpu(g$r), .print_ref(a$r), label = lab)
  }
})
