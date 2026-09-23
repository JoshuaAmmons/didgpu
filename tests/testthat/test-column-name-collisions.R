# A user column must never be mistaken for an argument naming a column.
#
# `[.data.table` evaluates its `i` (and `j`) expressions with the table's
# COLUMNS in scope. Several entry points filtered with
#
#     d <- d[!is.na(get(outcome)) & !is.na(get(group)) & ...]
#
# so a panel carrying a column literally called `outcome`, `group`,
# `time`, `treatment` or `weight` shadowed the argument holding that
# column's NAME, and get() was handed a vector instead of a string:
#
#     numeric column   -> "invalid first argument"
#     character column -> "first argument has length > 1"
#
# Those are ordinary names in a real analysis frame, so this fired on
# full panels while the same rows in a minimal four-column frame worked
# -- which reads as "extra columns break it" rather than as a name
# collision. Affected didgpu_bacon(), didgpu_fect(),
# didgpu_did_continuous() and didgpu_fhs(); didgpu() had the same issue
# on its `weight` argument.

.coll_panel <- function(seed = 11L, n_units = 120L, n_periods = 10L) {
  set.seed(seed)
  ufe <- stats::rnorm(n_units, 0, 1); tfe <- stats::rnorm(n_periods, 0, 0.3)
  Fg <- rep(Inf, n_units)
  tr <- sort(sample(seq_len(n_units), as.integer(n_units * 0.6)))
  Fg[tr] <- sample(3:(n_periods - 2L), length(tr), replace = TRUE)
  p <- expand.grid(period = 1:n_periods, unit = 1:n_units)
  p <- p[order(p$unit, p$period), ]
  p$D <- as.integer(p$period >= Fg[p$unit])
  p$Y <- ufe[p$unit] + tfe[p$period] - 0.02 * p$D +
         stats::rnorm(nrow(p), 0, 0.4)
  p[, c("unit", "period", "D", "Y")]
}

.hazard_names <- c("outcome", "group", "time", "treatment", "weight")

test_that("didgpu_bacon tolerates columns named after its arguments", {
  base <- .coll_panel()
  for (nm in .hazard_names) for (ty in c("chr", "num")) {
    p <- base
    p[[nm]] <- if (ty == "chr") paste0("v", p$unit) else stats::rnorm(nrow(p))
    expect_no_error(didgpu_bacon(p, "Y", "unit", "period", "D"))
  }
})

test_that("didgpu_fect tolerates them and returns unchanged estimates", {
  base <- .coll_panel()
  ref <- didgpu_fect(base, "Y", "unit", "period", "D",
                     bootstrap_reps = 0L, verbose = FALSE)$results$ATE[1L, 1L]
  for (nm in .hazard_names) for (ty in c("chr", "num")) {
    p <- base
    p[[nm]] <- if (ty == "chr") paste0("v", p$unit) else stats::rnorm(nrow(p))
    got <- didgpu_fect(p, "Y", "unit", "period", "D",
                       bootstrap_reps = 0L, verbose = FALSE)$results$ATE[1L, 1L]
    expect_equal(got, ref, tolerance = 1e-12,
                 label = paste(nm, ty))
  }
})

test_that("didgpu and didgpu_cs tolerate them too", {
  base <- .coll_panel()
  for (nm in .hazard_names) {
    p <- base
    p[[nm]] <- paste0("v", p$unit)
    expect_no_error(didgpu(p, "Y", "unit", "period", "D", effects = 2L,
                           placebo = 0L, bootstrap_reps = 0L,
                           backend = "r", verbose = FALSE))
    expect_no_error(didgpu_cs(p, "Y", "unit", "period", "D",
                              bootstrap_reps = 0L, backend = "r",
                              verbose = FALSE))
  }
})

test_that("a weight column named 'weight' is used, not shadowed", {
  # core_r.R resolved the weight column with get(weight) inside `j`.
  base <- .coll_panel()
  base$weight <- 1
  a <- didgpu(base, "Y", "unit", "period", "D", effects = 2L, placebo = 0L,
              bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  b <- didgpu(base, "Y", "unit", "period", "D", weight = "weight",
              effects = 2L, placebo = 0L, bootstrap_reps = 0L,
              backend = "r", verbose = FALSE)
  # Unit weights of 1 must reproduce the unweighted fit exactly.
  expect_equal(as.numeric(a$results$Effects[, "Estimate"]),
               as.numeric(b$results$Effects[, "Estimate"]),
               tolerance = 1e-12)
})
