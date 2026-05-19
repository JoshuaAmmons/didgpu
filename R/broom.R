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
#' @param x A `didgpu_result` object.
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
  stopifnot(inherits(x, "didgpu_result"))

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
