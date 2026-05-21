# ============================================================================
# didgpu_by(): run didgpu() separately for each level of a grouping
# variable. Mirrors the `by =` argument in DIDmultiplegtDYN but is
# implemented as a wrapper here so all the existing didgpu features
# (checkpoint, resume, parallel workers, normalized, ...) just work
# for each subset.
# ============================================================================


#' Estimate didgpu separately for each level of a grouping variable
#'
#' Wrapper that splits the panel by the unique values of `by_var` and
#' runs [didgpu()] on each subset, returning a named list of
#' `didgpu_result` objects. All other arguments are forwarded
#' identically to each per-level call.
#'
#' If `checkpoint_dir` is supplied, each subgroup writes into
#' `checkpoint_dir/<level>/`. Resume works the same way as a normal
#' [didgpu()] call inside each subgroup directory.
#'
#' @param df A panel data.frame.
#' @param by_var Character. Name of the grouping column. Each row in
#'   `df` must have exactly one value of this column. Levels must be
#'   coercible to character (used as subgroup labels and as
#'   subdirectory names under `checkpoint_dir`, if set).
#' @param outcome,group,time,treatment Forwarded to [didgpu()].
#' @param ... Additional arguments forwarded to [didgpu()] (e.g.
#'   `effects`, `placebo`, `bootstrap_reps`, `normalized`).
#' @param checkpoint_dir If non-NULL, each subgroup writes into a
#'   subdirectory of this path named after the level.
#' @param verbose Logical. Print one line per subgroup as it starts.
#' @return An object of class `didgpu_by_result` — a named list of
#'   `didgpu_result` objects, one per subgroup, plus an attribute
#'   `by_var` recording the grouping column name.
#' @examples
#' p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L,
#'                             tau_profile = c(0.5, 1.0),
#'                             seed = 17L)
#' p$region <- ifelse(p$unit %% 2L == 0L, "north", "south")
#' fit_by <- didgpu_by(p, "region",
#'                      outcome = "Y", group = "unit",
#'                      time = "period", treatment = "D",
#'                      effects = 2L, bootstrap_reps = 0L,
#'                      backend = "r", verbose = FALSE)
#' print(fit_by)
#' @export
didgpu_by <- function(df, by_var,
                       outcome, group, time, treatment,
                       ...,
                       checkpoint_dir = NULL,
                       verbose = TRUE) {
  stopifnot(is.character(by_var), length(by_var) == 1L,
            nzchar(by_var))
  if (!by_var %in% names(df)) stop("by_var column not in df: ", by_var)

  levels <- sort(unique(df[[by_var]]))
  if (length(levels) < 2L) {
    warning("by_var has only ", length(levels), " distinct level(s); ",
            "didgpu_by() returns a list of length ", length(levels), ".")
  }
  out <- vector("list", length(levels))
  names(out) <- as.character(levels)

  for (i in seq_along(levels)) {
    lvl <- levels[i]
    if (verbose) {
      message(sprintf("[didgpu_by] subgroup %d/%d: %s = %s",
                      i, length(levels), by_var, as.character(lvl)))
    }
    sub_df <- df[df[[by_var]] == lvl, , drop = FALSE]
    sub_cdir <- if (!is.null(checkpoint_dir)) {
      file.path(checkpoint_dir,
                gsub("[^A-Za-z0-9_.-]", "_", as.character(lvl)))
    } else NULL
    out[[i]] <- tryCatch(
      didgpu(
        sub_df, outcome = outcome, group = group,
        time = time, treatment = treatment,
        ...,
        checkpoint_dir = sub_cdir,
        verbose = FALSE),
      error = function(e) {
        warning(sprintf(
          "[didgpu_by] subgroup %s = %s failed: %s",
          by_var, as.character(lvl), conditionMessage(e)))
        # Return a sentinel so the user can see which subgroups failed.
        structure(list(error = conditionMessage(e),
                       subgroup = as.character(lvl)),
                   class = c("didgpu_by_failed", "list"))
      })
  }
  structure(out, class = c("didgpu_by_result", "list"),
            by_var = by_var)
}


#' Print method for didgpu_by_result
#'
#' Shows a compact table of per-subgroup point estimates.
#'
#' @param x A `didgpu_by_result` object.
#' @param ... Unused.
#' @return The input invisibly.
#' @export
print.didgpu_by_result <- function(x, ...) {
  by_var <- attr(x, "by_var")
  cat(sprintf("didgpu_by result (by = '%s'): %d subgroup(s)\n",
              by_var, length(x)))
  for (nm in names(x)) {
    r <- x[[nm]]
    if (inherits(r, "didgpu_by_failed")) {
      cat(sprintf("\n  [%s = %s]  FAILED: %s\n",
                  by_var, nm, r$error))
      next
    }
    eff <- as.numeric(r$results$Effects[, "Estimate"])
    se  <- as.numeric(r$results$Effects[, "SE"])
    if (length(eff) == 0L) {
      cat(sprintf("\n  [%s = %s]  no estimable effects (degenerate subgroup)\n",
                  by_var, nm))
      next
    }
    cat(sprintf("\n  [%s = %s]  effects (SE):\n", by_var, nm))
    out <- mapply(function(e, s, i) {
      sprintf("    Effect_%d: %s (%s)",
              i,
              if (is.na(e)) "      NA" else sprintf("%8.4f", e),
              if (is.na(s)) "    NA" else sprintf("%6.3f", s))
    }, eff, se, seq_along(eff))
    cat(paste(out, collapse = "\n"), "\n")
  }
  invisible(x)
}
