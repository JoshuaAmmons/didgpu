# ============================================================================
# broom-style methods for didgpu_result.
#
# These are registered as S3 methods on the broom::tidy and broom::glance
# generics, so the user pattern is:
#
#   library(broom); library(didgpu)
#   fit <- didgpu(...)
#   tidy(fit)             # one row per coefficient
#   glance(fit)           # one-row summary of the fit
#
# We provide direct didgpu_tidy() / didgpu_glance() entry points too,
# so the methods are usable without loading broom.
# ============================================================================


#' Tidy a didgpu_result into a one-row-per-coefficient data.frame
#'
#' Return shape matches the broom convention: columns `term`, `estimate`,
#' `std.error`, `statistic`, `p.value`, `conf.low`, `conf.high`. The
#' `term` column distinguishes effects (`Effect_1`, `Effect_2`, ...,
#' `ATE`) from placebos (`Placebo_1`, ...).
#'
#' @param x A `didgpu_result` or `didgpu_cs_result` object.
#'   For a `didgpu_cs_result` the return has one row per entry of
#'   the `aggregation` table (`kind = "aggregate"`) followed by one
#'   row per ATT(g, t) cell (`kind = "att_gt"`). `.cs_aggregate()`
#'   propagates point estimates only, so `std.error` and the CI
#'   columns are `NA` on the aggregate rows; the cell rows carry the
#'   real SEs and CIs.
#' @param conf.int Logical. Include confidence interval columns. Default TRUE.
#' @param conf.level Confidence level. Defaults to whatever the fit used.
#' @param ... Unused.
#' @return A data.frame.
#' @examples
#' p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 17L)
#' fit <- didgpu(p, "Y", "unit", "period", "D",
#'                effects = 2L, placebo = 1L,
#'                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
#' didgpu_tidy(fit)
#' @export
didgpu_tidy <- function(x, conf.int = TRUE, conf.level = NULL, ...) {
  # didgpu_cs_result is a first-class result type elsewhere in the
  # package (didgpu_loo(), didgpu_honest_did()), but it is shaped
  # differently from a didgpu_result -- $att_gt / $aggregation rather
  # than $results$Effects -- so it needs its own tidier. Previously it
  # fell through to the stopifnot() below and failed with the opaque
  # message `inherits(x, "didgpu_result") is not TRUE`.
  if (inherits(x, "didgpu_cs_result")) {
    return(.tidy_cs(x, conf.int = conf.int))
  }
  if (!inherits(x, "didgpu_result")) {
    stop("didgpu_tidy(): `x` must be a didgpu_result or a ",
         "didgpu_cs_result; got ",
         paste(class(x), collapse = "/"), ".", call. = FALSE)
  }

  rows <- list()
  if (!is.null(x$results$Effects) && nrow(x$results$Effects) > 0L) {
    eff <- x$results$Effects
    rows[[length(rows) + 1L]] <- .row_block(
      term       = rownames(eff) %||% paste0("Effect_", seq_len(nrow(eff))),
      estimate   = eff[, "Estimate"],
      std.error  = eff[, "SE"],
      conf.low   = eff[, "LB.CI"],
      conf.high  = eff[, "UB.CI"],
      kind       = "effect"
    )
  }
  if (!is.null(x$results$ATE)) {
    ate <- x$results$ATE
    rows[[length(rows) + 1L]] <- .row_block(
      term       = "ATE",
      estimate   = ate[1, "Estimate"],
      std.error  = ate[1, "SE"],
      conf.low   = ate[1, "LB.CI"],
      conf.high  = ate[1, "UB.CI"],
      kind       = "ate"
    )
  }
  if (!is.null(x$results$Placebos) && nrow(x$results$Placebos) > 0L) {
    pl <- x$results$Placebos
    rows[[length(rows) + 1L]] <- .row_block(
      term       = rownames(pl) %||% paste0("Placebo_", seq_len(nrow(pl))),
      estimate   = pl[, "Estimate"],
      std.error  = pl[, "SE"],
      conf.low   = pl[, "LB.CI"],
      conf.high  = pl[, "UB.CI"],
      kind       = "placebo"
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!isTRUE(conf.int)) {
    out$conf.low <- NULL; out$conf.high <- NULL
  }
  out
}

# Tidy a didgpu_cs_result.
#
# Two granularities are reported, because only one of them has SEs:
#   * $aggregation -- the requested scheme (event / group / calendar /
#     overall). .cs_aggregate() propagates POINT ESTIMATES ONLY, so
#     std.error and the CI columns are NA here. We do not invent one.
#   * $att_gt -- the per-(g, t) ATT cells, which do carry se / ci_low /
#     ci_high (analytic, or bootstrap when bootstrap_reps > 0).
#' @keywords internal
#' @noRd
.tidy_cs <- function(x, conf.int = TRUE) {
  rows <- list()

  ag <- x$aggregation
  if (!is.null(ag) && nrow(ag) > 0L) {
    lvl <- if ("event_time" %in% names(ag)) ag$event_time
           else if ("level" %in% names(ag)) ag$level
           else seq_len(nrow(ag))
    scheme <- x$args$aggregation %||% "agg"
    rows[[length(rows) + 1L]] <- .row_block(
      term      = paste0(scheme, "_", lvl),
      estimate  = ag$estimate,
      std.error = if ("se" %in% names(ag)) ag$se else rep(NA_real_, nrow(ag)),
      conf.low  = rep(NA_real_, nrow(ag)),
      conf.high = rep(NA_real_, nrow(ag)),
      kind      = "aggregate"
    )
  }

  gt <- x$att_gt
  if (!is.null(gt) && nrow(gt) > 0L) {
    # ci_low / ci_high are written by the bootstrap path only, so they are
    # absent entirely when bootstrap_reps = 0. Fill rather than index a
    # missing column (which yields NULL and a rows-mismatch in data.frame).
    .col <- function(nm) {
      if (nm %in% names(gt)) gt[[nm]] else rep(NA_real_, nrow(gt))
    }
    rows[[length(rows) + 1L]] <- .row_block(
      term      = sprintf("ATT_g%s_t%s", gt$g, gt$t),
      estimate  = gt$att,
      std.error = .col("se"),
      conf.low  = .col("ci_low"),
      conf.high = .col("ci_high"),
      kind      = "att_gt"
    )
  }

  if (!length(rows)) {
    return(data.frame(term = character(0), estimate = numeric(0),
                      std.error = numeric(0), statistic = numeric(0),
                      p.value = numeric(0), conf.low = numeric(0),
                      conf.high = numeric(0), kind = character(0),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!isTRUE(conf.int)) { out$conf.low <- NULL; out$conf.high <- NULL }
  out
}


.row_block <- function(term, estimate, std.error, conf.low, conf.high, kind) {
  statistic <- estimate / std.error
  p.value <- 2 * stats::pnorm(-abs(statistic))
  data.frame(
    term      = term,
    estimate  = as.numeric(estimate),
    std.error = as.numeric(std.error),
    statistic = as.numeric(statistic),
    p.value   = as.numeric(p.value),
    conf.low  = as.numeric(conf.low),
    conf.high = as.numeric(conf.high),
    kind      = kind,
    stringsAsFactors = FALSE
  )
}


#' One-row glance summary of a didgpu_result
#'
#' Columns: `n_effects`, `n_placebos`, `n_switchers`, `n_obs_effect_1`,
#' `p_jointeffects`, `p_jointplacebo`, `n_boot`, `backend`, `seed`.
#'
#' @param x A `didgpu_result`.
#' @param ... Unused.
#' @return A one-row data.frame.
#' @examples
#' p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 17L)
#' fit <- didgpu(p, "Y", "unit", "period", "D",
#'                effects = 2L, bootstrap_reps = 0L,
#'                backend = "r", verbose = FALSE)
#' didgpu_glance(fit)
#' @export
didgpu_glance <- function(x, ...) {
  stopifnot(inherits(x, "didgpu_result"))
  e <- x$results$Effects
  p <- x$results$Placebos
  n_switch_e1 <- if (!is.null(e) && nrow(e) > 0L) as.integer(e[1, "Switchers"]) else NA_integer_
  n_obs_e1    <- if (!is.null(e) && nrow(e) > 0L) as.integer(e[1, "N"]) else NA_integer_

  data.frame(
    n_effects       = as.integer(x$results$N_Effects),
    n_placebos      = as.integer(x$results$N_Placebos),
    n_switchers_e1  = n_switch_e1,
    n_obs_e1        = n_obs_e1,
    p_jointeffects  = as.numeric(x$results$p_jointeffects %||% NA),
    p_jointplacebo  = as.numeric(x$results$p_jointplacebo %||% NA),
    n_boot          = as.integer(x$results$n_boot %||% NA),
    backend         = as.character(x$args$backend %||% NA),
    seed            = as.integer(x$args$seed %||% NA),
    stringsAsFactors = FALSE
  )
}


# -------- Register against broom generics if broom is loaded --------

# broom integration: we don't depend on broom, but if the user loads
# it, broom::tidy(fit) and broom::glance(fit) will dispatch to the S3
# methods below (registered in NAMESPACE). Users without broom can
# call didgpu_tidy() / didgpu_glance() directly.

#' broom::tidy method for didgpu_result
#'
#' @param x A `didgpu_result`.
#' @param ... Passed to [didgpu_tidy()].
#' @return A data.frame.
#' @export
tidy.didgpu_result <- function(x, ...) didgpu_tidy(x, ...)

#' broom::glance method for didgpu_result
#'
#' @param x A `didgpu_result`.
#' @param ... Unused.
#' @return A one-row data.frame.
#' @export
glance.didgpu_result <- function(x, ...) didgpu_glance(x, ...)
