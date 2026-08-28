# The CS cluster bootstrap must not care what the user's columns are named.
#
# Regression for a variable-shadowing crash. .cs_bootstrap_se() held the
# panel as a data.table and subset it with
#     block <- d[d[[args$group]] == u, , drop = FALSE]
# `[.data.table` evaluates its `i` expression with the table's COLUMNS in
# scope, so a panel carrying a column literally named `d` shadowed the
# local `d` with the treatment VECTOR; `d[[args$group]]` then became
# treatment[["g"]] and failed with "subscript out of bounds".
#
# `d` is an entirely ordinary name for a treatment indicator, so this hit
# real panels. It fired only with bootstrap_reps > 0 (point estimates go
# nowhere near this function), which made it look data-dependent rather
# than name-dependent: the same panel worked at reps = 0 and crashed at
# reps > 0, and didgpu's own simulated panels never tripped it because
# their treatment column is `D`.

.cs_panel <- function(treat_name = "D", outcome_name = "Y") {
  p <- as.data.frame(didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                                           seed = 11L))
  names(p)[names(p) == "D"] <- treat_name
  names(p)[names(p) == "Y"] <- outcome_name
  p
}

.cs_fit <- function(p, treat_name = "D", outcome_name = "Y", reps = 40L) {
  didgpu_cs(df = p, outcome = outcome_name, group = "unit", time = "period",
            treatment = treat_name, control_group = "never",
            est_method = "DR", aggregation = "event",
            bootstrap_reps = reps, seed = 1, backend = "r", verbose = FALSE)
}

test_that("bootstrap survives a treatment column named 'd'", {
  expect_no_error(.cs_fit(.cs_panel("d"), "d"))
})

test_that("bootstrap survives column names that collide with internals", {
  # Every local in .cs_bootstrap_se is a potential shadow if the panel is
  # ever held as a data.table again. Pin the whole set.
  for (nm in c("d", "panel", "units", "block", "picks", "key", "u", "b",
               "i", "grp_vec", "rows_by_unit", "pcount", "OFFSET")) {
    expect_no_error(.cs_fit(.cs_panel(nm), nm), message = nm)
  }
})

test_that("results are invariant to the treatment column's name", {
  a <- .cs_fit(.cs_panel("D"), "D")
  b <- .cs_fit(.cs_panel("d"), "d")
  expect_equal(a$att_gt$att, b$att_gt$att, tolerance = 0)
  expect_equal(a$att_gt$se,  b$att_gt$se,  tolerance = 0)
})

test_that("an outcome column named 'd' is also safe", {
  expect_no_error(.cs_fit(.cs_panel("trt", "d"), "trt", "d"))
})
