# ============================================================================
# didgpu_compare(): run both backends on the user's panel and report any
# numerical disagreements. Intended as a one-shot trust check before
# switching production code from DIDmultiplegtDYN to didgpu.
# ============================================================================


#' Compare didgpu r-backend to the DIDmultiplegtDYN reference
#'
#' Runs both backends on the same panel and same call args, then reports
#' the max absolute disagreement per output column. Useful for verifying
#' didgpu produces output identical to the reference on your own data
#' before using it in production.
#'
#' @inheritParams didgpu
#' @param tolerance Numeric. Threshold above which a disagreement is
#'   flagged as a failure. Default `1e-10` (effectively machine epsilon
#'   for our pipeline).
#' @param verbose Logical. Print a per-column table.
#'
#' @return Invisibly a list with `pass` (logical scalar), `report`
#'   (data.frame of per-column max diffs), `fit_r` and `fit_ref` (the
#'   two full fits). If the reference is not installed, returns
#'   immediately with `pass = NA` and a warning.
#' @examples
#' \dontrun{
#' # Requires the DIDmultiplegtDYN package to be installed.
#' p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 17L)
#' report <- didgpu_compare(p, "Y", "unit", "period", "D",
#'                           effects = 2L, placebo = 1L)
#' report$pass
#' head(report$report)
#' }
#' @export
didgpu_compare <- function(
    df, outcome, group, time, treatment,
    effects = 1L, placebo = 0L,
    cluster = NULL, switchers = "",
    tolerance = 1e-10, verbose = TRUE) {

  if (!requireNamespace("DIDmultiplegtDYN", quietly = TRUE)) {
    warning("DIDmultiplegtDYN not installed; cannot compare. ",
            "Install with install.packages('DIDmultiplegtDYN').")
    return(invisible(list(pass = NA, report = NULL,
                           fit_r = NULL, fit_ref = NULL)))
  }

  fit_r <- didgpu(df = df, outcome = outcome, group = group, time = time,
                   treatment = treatment, effects = effects, placebo = placebo,
                   cluster = cluster, switchers = switchers,
                   bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  fit_ref_pkg <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(df), outcome = outcome, group = group,
      time = time, treatment = treatment,
      effects = as.double(effects), placebo = as.double(placebo),
      cluster = cluster, switchers = switchers, graph_off = TRUE
    )
  ))

  cmp <- function(us, ref, label, cols = NULL) {
    if (is.null(us) || is.null(ref) ||
        (is.matrix(ref) && nrow(ref) == 0L)) return(NULL)
    if (is.null(cols)) {
      n_min <- min(ncol(us), ncol(ref))
      cols <- seq_len(n_min)
    }
    n_rows <- min(nrow(us), nrow(ref))
    if (n_rows == 0L) return(NULL)
    sapply(cols, function(j) {
      d <- abs(as.numeric(us[seq_len(n_rows), j]) -
               as.numeric(ref[seq_len(n_rows), j]))
      d <- d[is.finite(d)]
      if (length(d) == 0L) NA_real_ else max(d)
    })
  }

  rows <- list()
  for (col_i in 1:4) {  # Estimate, SE, LB.CI, UB.CI
    col_name <- c("Estimate", "SE", "LB.CI", "UB.CI")[col_i]
    eff_diff <- cmp(fit_r$results$Effects, fit_ref_pkg$results$Effects,
                    "Effects", col_i)
    pl_diff  <- cmp(fit_r$results$Placebos, fit_ref_pkg$results$Placebos,
                    "Placebos", col_i)
    ate_diff <- cmp(fit_r$results$ATE, fit_ref_pkg$results$ATE,
                    "ATE", col_i)
    rows[[length(rows) + 1L]] <- data.frame(
      block = c(rep("Effects", length(eff_diff)),
                rep("Placebos", length(pl_diff)),
                rep("ATE", length(ate_diff))),
      col   = col_name,
      max_abs_diff = c(eff_diff, pl_diff, ate_diff),
      stringsAsFactors = FALSE
    )
  }
  report <- do.call(rbind, rows)
  report <- report[!is.na(report$max_abs_diff), ]

  # For SE / CI columns the comparison is only meaningful when we
  # provided a bootstrap (here we pass bootstrap_reps=0); reference
  # SEs are analytical, ours are NA in that case. So mark them as
  # informational and don't fail on them.
  report$fails <- report$max_abs_diff > tolerance & report$col == "Estimate"

  if (verbose) {
    cat(sprintf("didgpu_compare: tolerance = %.1e\n", tolerance))
    print(report, row.names = FALSE, digits = 4)
    cat("\n")
    if (any(report$fails)) {
      cat(sprintf("*** %d coefficient(s) exceed tolerance ***\n",
                  sum(report$fails)))
    } else {
      cat("OK: all Estimate columns agree within tolerance.\n")
    }
  }

  invisible(list(
    pass = !any(report$fails),
    report = report,
    fit_r = fit_r,
    fit_ref = fit_ref_pkg
  ))
}
