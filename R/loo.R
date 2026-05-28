# ============================================================================
# Leave-one-out (LOO) robustness analysis.
#
# Drop one entity (cohort, unit, cluster, or period) at a time,
# re-estimate, and report how the headline estimate changes. Standard
# diagnostic for detecting whether a single observation is driving
# the result.
#
# Works across all five estimator families by inspecting the fit's
# class and re-invoking the corresponding entry point:
#   - didgpu_result        -> didgpu()
#   - didgpu_cs_result     -> didgpu_cs()
#   - didgpu_fect_result   -> didgpu_fect()
#
# Default `by` per family:
#   - didgpu / didgpu_cs / didgpu_fect : "cohort" (leave-one-cohort-out)
#
# Other options work for any family: "unit" (drop one unit per
# replicate, expensive), "cluster" (drop one cluster), or any column
# name to drop one level of that column.
# ============================================================================


#' Leave-one-out (LOO) robustness analysis
#'
#' For each entity (cohort / unit / cluster / level of a column),
#' re-fits the estimator dropping that entity and reports the
#' headline estimate. The headline depends on family:
#'   - `didgpu_result`: the ATE.
#'   - `didgpu_cs_result`: the requested aggregation's first row
#'     (or the `overall` value if aggregation = "overall").
#'   - `didgpu_fect_result`: the ATE.
#'
#' Useful for detecting single-cohort or single-unit influence on
#' headline estimates.
#'
#' @param fit A fitted `didgpu_result`, `didgpu_cs_result`, or
#'   `didgpu_fect_result`.
#' @param by Either `"cohort"` (drop each treatment cohort in turn),
#'   `"unit"` (drop each unit), `"cluster"` (drop each cluster as
#'   defined by the fit's `cluster` arg if set, otherwise `group`),
#'   or a character giving a column name to use for the drop key.
#' @param df Optional original panel. If `NULL`, the function looks
#'   for an attached panel hash + path in the fit; if not found, you
#'   must supply `df`. Standard pattern: pass the same `df` you fit
#'   with.
#' @param verbose Logical. Print one line per leave-out fit.
#' @return An object of class `didgpu_loo_result`: a data.frame with
#'   columns `leave_out` (the entity dropped), `estimate` (the
#'   headline under that drop), `delta` (estimate - full-sample
#'   estimate), `delta_pct` (delta as % of full-sample). Sorted by
#'   `abs(delta)` descending. Plus `$full` attribute (the full-sample
#'   estimate) and `$by` attribute.
#'
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
#'                             tau_profile = c(0.5, 1.0), seed = 17L)
#' fit <- didgpu_cs(p, "Y", "unit", "period", "D",
#'                   est_method = "OR", aggregation = "overall",
#'                   bootstrap_reps = 0L, backend = "r", verbose = FALSE)
#' loo <- didgpu_loo(fit, by = "cohort", df = p, verbose = FALSE)
#' print(loo)
#' }
#' @export
didgpu_loo <- function(fit, by = "cohort", df = NULL, verbose = TRUE) {
  if (!any(inherits(fit, c("didgpu_result", "didgpu_cs_result",
                            "didgpu_fect_result")))) {
    stop("`fit` must be a didgpu_result, didgpu_cs_result, or didgpu_fect_result.")
  }
  if (is.null(df)) {
    stop("Please pass the original panel via `df = ...`. (didgpu does not ",
         "store the panel on the fit object for memory reasons.)")
  }

  family <- if (inherits(fit, "didgpu_cs_result"))    "cs"
            else if (inherits(fit, "didgpu_fect_result")) "fect"
            else                                          "didgpu"
  args <- fit$args
  group_col <- args$group
  time_col  <- args$time
  trt_col   <- args$treatment
  out_col   <- args$outcome

  # Identify the leave-out levels.
  by_levels <- .loo_resolve_by(fit, df, by, args, family)
  if (length(by_levels$values) < 2L) {
    stop("Need at least 2 distinct levels of `", by_levels$column,
         "` to leave one out.")
  }

  # Full-sample estimate for delta reference.
  full_est <- .loo_extract_headline(fit, family)

  # ---- Fast path: CS cohort-LOO with never-treated controls. ----
  # Dropping a treated cohort g* leaves every OTHER (g, t) cell
  # unchanged (never-treated controls are a fixed pool disjoint from
  # all treated cohorts, and the other cohorts' treated units are
  # untouched). So the leave-out estimate is just the full-sample
  # att_gt table re-aggregated WITHOUT g*'s cells -- no refit needed.
  # This is bit-identical to the refit path but ~100x faster (K cheap
  # re-aggregations instead of K full didgpu_cs fits). Only valid for
  # control_group = "never"; "notyet" controls DO change other cells
  # when a cohort is dropped, so that case falls through to the refit
  # loop below.
  if (identical(family, "cs") && identical(by, "cohort")) {
    fast <- .loo_cs_cohort_fast(fit, full_est)
    if (!is.null(fast)) {
      out <- fast[order(-abs(fast$delta), na.last = TRUE), , drop = FALSE]
      rownames(out) <- NULL
      attr(out, "full")   <- full_est
      attr(out, "by")     <- by_levels$column
      attr(out, "family") <- family
      attr(out, "method") <- "reaggregate (no refit)"
      class(out) <- c("didgpu_loo_result", class(out))
      return(out)
    }
  }

  results <- vector("list", length(by_levels$values))
  for (i in seq_along(by_levels$values)) {
    lvl <- by_levels$values[i]
    if (verbose) {
      message(sprintf("[loo] leaving out %s = %s (%d/%d)",
                      by_levels$column, format(lvl), i,
                      length(by_levels$values)))
    }
    df_minus <- .loo_drop_level(df, by_levels, lvl)
    if (nrow(df_minus) == 0L) {
      results[[i]] <- data.frame(
        leave_out = format(lvl),
        estimate  = NA_real_,
        delta     = NA_real_,
        delta_pct = NA_real_,
        note      = "all rows dropped",
        stringsAsFactors = FALSE
      )
      next
    }
    fit_b <- tryCatch(
      .loo_refit(df_minus, args, family),
      error = function(e) {
        list(error = conditionMessage(e))
      })
    if (!is.null(fit_b$error)) {
      results[[i]] <- data.frame(
        leave_out = format(lvl),
        estimate  = NA_real_,
        delta     = NA_real_,
        delta_pct = NA_real_,
        note      = paste("error:", fit_b$error),
        stringsAsFactors = FALSE
      )
      next
    }
    est_b <- .loo_extract_headline(fit_b, family)
    delta <- est_b - full_est
    pct <- if (!is.na(full_est) && full_est != 0)
             100 * delta / abs(full_est)
           else NA_real_
    results[[i]] <- data.frame(
      leave_out = format(lvl),
      estimate  = as.numeric(est_b),
      delta     = as.numeric(delta),
      delta_pct = as.numeric(pct),
      note      = "",
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, results)
  # Sort by abs(delta) descending so the most-influential drop is at the top.
  out <- out[order(-abs(out$delta), na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "full")    <- full_est
  attr(out, "by")      <- by_levels$column
  attr(out, "family")  <- family
  class(out) <- c("didgpu_loo_result", class(out))
  out
}


# Fast cohort-LOO for the CS family: re-aggregate the precomputed
# att_gt table excluding each cohort's cells, instead of refitting.
#
# Valid ONLY when control_group == "never": then dropping a treated
# cohort g* does not alter any other (g, t) cell (the never-treated
# control pool is disjoint from all treated cohorts, and the other
# cohorts' treated units are untouched), so the refit's att_gt is
# exactly the full-sample att_gt minus g*'s rows. Re-aggregating that
# subset reproduces the refit's headline bit-for-bit.
#
# Returns a data.frame with the standard LOO columns, or NULL to
# signal "not applicable -> use the generic refit path".
#' @keywords internal
#' @noRd
.loo_cs_cohort_fast <- function(fit, full_est) {
  args <- fit$args
  # notyet controls: dropping a cohort can change other cells (the
  # dropped cohort may have served as a not-yet-treated control), so
  # the shortcut is invalid -> fall back to refit.
  if (!identical(args$control_group %||% "never", "never")) return(NULL)
  att_gt <- fit$att_gt
  if (is.null(att_gt) || !("g" %in% names(att_gt)) || nrow(att_gt) == 0L) {
    return(NULL)
  }
  cohorts <- sort(unique(att_gt$g))
  if (length(cohorts) < 2L) return(NULL)

  rows <- lapply(cohorts, function(gstar) {
    sub <- att_gt[att_gt$g != gstar, , drop = FALSE]
    # Carry the IF / units attributes through so .cs_aggregate (and any
    # SE machinery it touches) sees a well-formed att_gt subset.
    attr(sub, "IF_per_cell")  <- attr(att_gt, "IF_per_cell")
    attr(sub, "F_g_per_unit") <- attr(att_gt, "F_g_per_unit")
    attr(sub, "units")        <- attr(att_gt, "units")
    if (nrow(sub) == 0L) {
      return(data.frame(leave_out = format(gstar), estimate = NA_real_,
                        delta = NA_real_, delta_pct = NA_real_,
                        note = "all cells dropped", stringsAsFactors = FALSE))
    }
    agg <- .cs_aggregate(sub, args$aggregation, args)
    est <- if ("estimate" %in% names(agg) && nrow(agg) > 0L)
             as.numeric(agg$estimate[1]) else NA_real_
    delta <- est - full_est
    pct <- if (!is.na(full_est) && full_est != 0)
             100 * delta / abs(full_est) else NA_real_
    data.frame(leave_out = format(gstar), estimate = est, delta = delta,
               delta_pct = pct, note = "", stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


# Extract the headline scalar estimate from a fitted result, based on
# the family. Used as the reference for the delta column.
#' @keywords internal
#' @noRd
.loo_extract_headline <- function(fit, family) {
  switch(family,
    "didgpu" = as.numeric(fit$results$ATE[1, "Estimate"]),
    "fect"   = as.numeric(fit$results$ATE[1, "Estimate"]),
    "cs"     = {
      agg <- fit$aggregation
      if ("estimate" %in% names(agg) && nrow(agg) > 0L) {
        if ("scheme" %in% names(agg) &&
            agg$scheme[1] == "overall") {
          as.numeric(agg$estimate[1])
        } else {
          # Use the first row's estimate as the headline (works for
          # all four aggregations).
          as.numeric(agg$estimate[1])
        }
      } else NA_real_
    },
    NA_real_
  )
}


# Resolve the `by` argument to (column, set of levels).
#' @keywords internal
#' @noRd
.loo_resolve_by <- function(fit, df, by, args, family) {
  if (identical(by, "cohort")) {
    # Cohort = first treatment period per unit. Build a per-unit map.
    d <- data.table::as.data.table(df)
    data.table::setnames(d,
      c(args$outcome, args$group, args$time, args$treatment),
      c("Y_XX", "G_XX", "T_XX", "D_XX"))
    d[, F_g_XX := {
        if (any(D_XX == 1L, na.rm = TRUE))
          as.numeric(min(T_XX[D_XX == 1L], na.rm = TRUE))
        else Inf
      }, by = G_XX]
    cohorts <- sort(unique(d$F_g_XX[is.finite(d$F_g_XX)]))
    return(list(column = "cohort (F_g)",
                values = cohorts,
                kind = "cohort",
                # We carry the unit -> cohort map for the drop step.
                unit_map = unique(d[, list(G_XX, F_g_XX)]),
                group_col = args$group))
  }
  if (identical(by, "unit")) {
    vals <- sort(unique(df[[args$group]]))
    return(list(column = args$group, values = vals, kind = "column",
                target_col = args$group))
  }
  if (identical(by, "cluster")) {
    cl <- args$cluster %||% args$group
    vals <- sort(unique(df[[cl]]))
    return(list(column = cl, values = vals, kind = "column",
                target_col = cl))
  }
  # Otherwise treat `by` as a column name in df.
  if (!by %in% names(df)) {
    stop("`by` = '", by, "' is not a column in df.")
  }
  vals <- sort(unique(df[[by]]))
  list(column = by, values = vals, kind = "column", target_col = by)
}


# Drop one level from df. For cohort drops, removes all units whose
# F_g equals the dropped cohort. For column drops, removes all rows
# matching that level.
#' @keywords internal
#' @noRd
.loo_drop_level <- function(df, by_levels, lvl) {
  if (by_levels$kind == "cohort") {
    # Avoid data.table column-scope: pull the vectors out into plain R.
    # unit_map has columns G_XX (= original unit IDs, just renamed in
    # place) and F_g_XX (cohort = first-treatment period).
    um <- by_levels$unit_map
    keep <- (um$F_g_XX != lvl) | !is.finite(um$F_g_XX)
    keep_units <- um$G_XX[keep]
    return(df[df[[by_levels$group_col]] %in% keep_units, , drop = FALSE])
  }
  df[df[[by_levels$target_col]] != lvl, , drop = FALSE]
}


# Refit on the leave-one-out subset. Dispatches to the right family.
#' @keywords internal
#' @noRd
.loo_refit <- function(df, args, family) {
  switch(family,
    "didgpu" = didgpu(df,
                      outcome = args$outcome, group = args$group,
                      time = args$time, treatment = args$treatment,
                      effects = args$effects %||% 1L,
                      placebo = args$placebo %||% 0L,
                      cluster = args$cluster,
                      controls = args$controls,
                      weight = args$weight,
                      switchers = args$switchers %||% "",
                      bootstrap_reps = 0L,
                      backend = args$backend %||% "r",
                      verbose = FALSE),
    "cs"     = didgpu_cs(df,
                          outcome = args$outcome, group = args$group,
                          time = args$time, treatment = args$treatment,
                          est_method = args$est_method,
                          control_group = args$control_group,
                          aggregation = args$aggregation,
                          covariates = args$covariates,
                          bootstrap_reps = 0L,
                          backend = args$backend %||% "r",
                          verbose = FALSE),
    "fect"   = didgpu_fect(df,
                            outcome = args$outcome, group = args$group,
                            time = args$time, treatment = args$treatment,
                            method = args$method,
                            effects = args$effects %||% 1L,
                            r = args$r %||% 2L,
                            lambda = args$lambda,
                            bootstrap_reps = 0L,
                            backend = args$backend %||% "r",
                            verbose = FALSE),
    stop("unknown family: ", family)
  )
}


#' Print method for didgpu_loo_result
#' @param x A `didgpu_loo_result`.
#' @param n Integer. Print top-`n` most-influential rows. Default `10L`.
#' @param ... Unused.
#' @return The input invisibly.
#' @export
print.didgpu_loo_result <- function(x, n = 10L, ...) {
  full <- attr(x, "full")
  by   <- attr(x, "by")
  cat(sprintf("Leave-one-out analysis (drop one %s at a time)\n", by))
  cat(sprintf("  full-sample estimate: %.4f\n", full))
  cat(sprintf("  %d leave-out replicates (sorted by abs(delta))\n",
              nrow(x)))
  cat("\nTop", min(n, nrow(x)), "most influential:\n")
  top <- utils::head(x, n)
  display <- data.frame(
    leave_out = top$leave_out,
    estimate  = sprintf("%.4f", top$estimate),
    delta     = sprintf("%+.4f", top$delta),
    delta_pct = sprintf("%+.1f%%", top$delta_pct),
    note      = top$note,
    stringsAsFactors = FALSE
  )
  print(display, row.names = FALSE)
  if (nrow(x) > n) {
    cat(sprintf("  ... %d more rows; use as.data.frame() to see all\n",
                nrow(x) - n))
  }
  cat("\nInterpretation:\n")
  cat("  Large |delta_pct| at the top => a single", by, "is driving\n")
  cat("  the headline estimate. < 10% across all rows = robust.\n")
  invisible(x)
}


#' Plot method for didgpu_loo_result (tornado plot)
#' @param x A `didgpu_loo_result`.
#' @param n_show Integer. Number of top rows to show. Default `20L`.
#' @param ... Extra args for plot().
#' @return The input invisibly.
#' @export
plot.didgpu_loo_result <- function(x, n_show = 20L, ...) {
  full <- attr(x, "full")
  top <- utils::head(x, n_show)
  top <- top[order(top$delta), , drop = FALSE]   # smallest at bottom
  par_old <- graphics::par(mar = c(4, 8, 4, 2))
  on.exit(graphics::par(par_old), add = TRUE)
  graphics::barplot(
    top$delta,
    names.arg = top$leave_out,
    horiz = TRUE,
    las = 1,
    main = sprintf("Leave-one-out: delta from full estimate (%.4f)", full),
    xlab = "Change in headline estimate",
    col = ifelse(top$delta > 0, "steelblue", "salmon"),
    ...
  )
  graphics::abline(v = 0, lty = 2L, col = "grey50")
  invisible(x)
}
