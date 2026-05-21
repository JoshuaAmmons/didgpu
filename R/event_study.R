# ============================================================================
# Event-study helpers: extract the long-format data needed to plot the
# event-time profile (placebos to the left of zero, effects to the right).
# ============================================================================


#' Extract long-format event-study data from a didgpu_result
#'
#' Returns one row per event-time (placebos at negative k, effects at
#' k = 1..effects), with columns `event_time`, `estimate`, `std.error`,
#' `conf.low`, `conf.high`, `kind`. Suitable for direct plotting with
#' ggplot2 — e.g.:
#'
#' \preformatted{
#'   library(ggplot2)
#'   ggplot(didgpu_event_study_data(fit),
#'          aes(event_time, estimate)) +
#'     geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
#'     geom_hline(yintercept = 0, linetype = "dashed")
#' }
#'
#' Note: placebos are conventionally plotted at "negative" event time
#' (`-1, -2, ...`) even though they are computed at horizon `1, 2, ...`
#' from the switch. The `event_time` column follows the convention.
#'
#' @param x A `didgpu_result`.
#' @return A data.frame.
#' @examples
#' p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 17L)
#' fit <- didgpu(p, "Y", "unit", "period", "D",
#'                effects = 2L, placebo = 1L,
#'                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
#' didgpu_event_study_data(fit)
#' @export
didgpu_event_study_data <- function(x) {
  stopifnot(inherits(x, "didgpu_result"))

  rows <- list()

  if (!is.null(x$results$Placebos) && nrow(x$results$Placebos) > 0L) {
    pl <- x$results$Placebos
    rows[[length(rows) + 1L]] <- data.frame(
      event_time = -seq_len(nrow(pl)),
      estimate   = as.numeric(pl[, "Estimate"]),
      std.error  = as.numeric(pl[, "SE"]),
      conf.low   = as.numeric(pl[, "LB.CI"]),
      conf.high  = as.numeric(pl[, "UB.CI"]),
      kind       = "placebo",
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(x$results$Effects) && nrow(x$results$Effects) > 0L) {
    e <- x$results$Effects
    rows[[length(rows) + 1L]] <- data.frame(
      event_time = seq_len(nrow(e)),
      estimate   = as.numeric(e[, "Estimate"]),
      std.error  = as.numeric(e[, "SE"]),
      conf.low   = as.numeric(e[, "LB.CI"]),
      conf.high  = as.numeric(e[, "UB.CI"]),
      kind       = "effect",
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(data.frame(
      event_time = integer(0), estimate = numeric(0),
      std.error = numeric(0), conf.low = numeric(0),
      conf.high = numeric(0), kind = character(0)
    ))
  }
  out <- out[order(out$event_time), , drop = FALSE]
  rownames(out) <- NULL
  out
}
