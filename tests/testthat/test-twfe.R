# didgpu_twfe(): naive TWFE dynamic event-study baseline.
#
# Correctness anchors:
#   - point estimates must match lm(Y ~ lead/lags + factor(unit) +
#     factor(time)) bit-for-bit (both are OLS with two-way FE),
#   - cluster-robust SEs must match a hand-computed CR1 sandwich.

# Build a panel with a known non-absorbing binary treatment.
make_twfe_panel <- function(n_units = 40L, n_periods = 12L, seed = 5L) {
  set.seed(seed)
  unit_fe <- stats::rnorm(n_units)
  time_fe <- stats::rnorm(n_periods, sd = 0.5)
  g <- rep(seq_len(n_units), each = n_periods)
  t <- rep(seq_len(n_periods), times = n_units)
  D <- as.integer(stats::runif(n_units * n_periods) < 0.25)  # non-absorbing
  Y <- unit_fe[g] + time_fe[t] + 1.0 * D + stats::rnorm(n_units * n_periods, sd = 0.3)
  data.frame(unit = g, period = t, D = D, Y = Y)
}

# Reference: lm with explicit lead/lag regressors + FE dummies.
ref_twfe <- function(p, effects, placebo) {
  d <- data.table::as.data.table(p)
  data.table::setkey(d, unit, period)
  regs <- character(0)
  for (k in seq_len(effects)) {
    nm <- paste0("Effect_", k); lag <- k - 1L
    src <- d[, list(unit, period = period + lag, vv = D)]
    d[src, (nm) := i.vv, on = c("unit", "period")]
    d[is.na(get(nm)), (nm) := 0]; regs <- c(regs, nm)
  }
  for (j in seq_len(placebo)) {
    nm <- paste0("Placebo_", j)
    src <- d[, list(unit, period = period - j, vv = D)]
    d[src, (nm) := i.vv, on = c("unit", "period")]
    d[is.na(get(nm)), (nm) := 0]; regs <- c(regs, nm)
  }
  form <- stats::as.formula(paste0(
    "Y ~ ", paste(regs, collapse = " + "),
    " + factor(unit) + factor(period)"))
  fit <- stats::lm(form, data = as.data.frame(d))
  list(coef = stats::coef(fit)[regs], dt = d, regs = regs, lm = fit)
}

test_that("didgpu_twfe point estimates match lm() with FE dummies (bit-for-bit)", {
  p <- make_twfe_panel()
  fit <- didgpu_twfe(p, "Y", "unit", "period", "D",
                      effects = 3L, placebo = 2L, verbose = FALSE)
  ref <- ref_twfe(p, 3L, 2L)
  got <- fit$coef[names(ref$coef)]
  expect_equal(unname(got), unname(ref$coef), tolerance = 1e-8)
})

test_that("didgpu_twfe cluster-robust SE matches a hand-computed CR1 sandwich", {
  p <- make_twfe_panel()
  fit <- didgpu_twfe(p, "Y", "unit", "period", "D",
                      effects = 3L, placebo = 1L, cluster = "unit",
                      verbose = FALSE)
  ref <- ref_twfe(p, 3L, 1L)

  # Hand CR1 on the lm: build the within-design and residuals from lm,
  # then sandwich with the same dof correction didgpu_twfe uses.
  X <- stats::model.matrix(ref$lm)[, ref$regs, drop = FALSE]
  # Residualize X and y on the FE columns so we work in the within space.
  fe_cols <- setdiff(colnames(stats::model.matrix(ref$lm)), c("(Intercept)", ref$regs))
  FE <- stats::model.matrix(ref$lm)[, fe_cols, drop = FALSE]
  resid_on_fe <- function(v) stats::lm.fit(cbind(1, FE), v)$residuals
  Xd <- apply(X, 2L, resid_on_fe)
  yd <- resid_on_fe(p$Y[as.integer(rownames(X))])
  b  <- solve(crossprod(Xd), crossprod(Xd, yd))
  u  <- as.numeric(yd - Xd %*% b)
  cl <- as.integer(factor(ref$dt$unit))
  sg <- rowsum(Xd * u, cl)
  meat <- Reduce(`+`, lapply(seq_len(nrow(sg)), function(i) tcrossprod(sg[i, ])))
  XtXi <- solve(crossprod(Xd))
  N <- nrow(Xd); G <- length(unique(cl))
  n_units <- length(unique(p$unit)); n_times <- length(unique(p$period))
  kpar <- ncol(Xd) + n_units + (n_times - 1L)
  cc <- (G / (G - 1)) * ((N - 1) / (N - kpar))
  V <- XtXi %*% meat %*% XtXi * cc
  ref_se <- sqrt(diag(V))

  got_se <- fit$results$Effects[, "SE"]
  expect_equal(unname(got_se), unname(ref_se[paste0("Effect_", 1:3)]),
               tolerance = 1e-6)
})

test_that("didgpu_twfe returns the expected structure and recovers ~true effect", {
  p <- make_twfe_panel(n_units = 120L, n_periods = 14L, seed = 9L)
  fit <- didgpu_twfe(p, "Y", "unit", "period", "D",
                      effects = 4L, placebo = 3L, verbose = FALSE)
  expect_s3_class(fit, "didgpu_twfe_result")
  expect_equal(nrow(fit$results$Effects), 4L)
  expect_equal(nrow(fit$results$Placebos), 3L)
  expect_true(all(c("Estimate", "SE", "LB.CI", "UB.CI", "N") %in%
                   colnames(fit$results$Effects)))
  # True contemporaneous effect is 1.0; Effect_1 should be close.
  expect_lt(abs(fit$results$Effects["Effect_1", "Estimate"] - 1.0), 0.15)
  # SEs positive and finite.
  expect_true(all(is.finite(fit$results$Effects[, "SE"])))
  expect_true(all(fit$results$Effects[, "SE"] > 0))
})

test_that("didgpu_twfe validates arguments", {
  p <- make_twfe_panel(n_units = 20L, n_periods = 6L)
  expect_error(didgpu_twfe(p, "nope", "unit", "period", "D", verbose = FALSE),
               "column not in df")
  expect_error(didgpu_twfe(p, "Y", "unit", "period", "D", effects = 0L,
                            verbose = FALSE), "effects")
})
