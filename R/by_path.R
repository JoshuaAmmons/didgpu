# ============================================================================
# didgpu_compute_paths / didgpu_by_path
#
# Treatment-trajectory subgroup analysis. For each group g, build a
# string identifier from the treatment values at F_g - 1 (the baseline
# period) and at F_g, F_g + 1, ..., F_g + effects - 1 (the first
# `effects` post-switch periods). Groups whose treatment trajectory is
# identical share the same path; the panel can then be analysed
# per-path with `didgpu_by()`.
#
# Mirrors DIDmultiplegtDYN's `by_path` option, which is itself a
# wrapper around `did_multiplegt_dyn(..., by = "<path-column>")`.
# Reference: did_multiplegt_by_path() in DIDmultiplegtDYN.
# ============================================================================


#' Add a per-group treatment-trajectory column to a panel
#'
#' For each group, builds a comma-separated string of the treatment
#' values observed at (F_g - 1, F_g, F_g + 1, ..., F_g + effects - 1):
#' the baseline period plus the first `effects` post-switch periods.
#' Groups with no switch (`F_g > T_max`) get a "no-switch" path string.
#'
#' Use the resulting `path` column as a `by_var` to estimate treatment
#' effects separately per trajectory (e.g. via `didgpu_by()`).
#'
#' @param df A panel data.frame.
#' @param outcome,group,time,treatment Column names (we only need them
#'   to compute F_g consistently with `didgpu()`).
#' @param effects Integer. Number of post-switch periods to encode in
#'   the path (i.e. the path string has length `effects + 1`).
#' @param top_n Integer or NULL. If non-NULL, keep only the `top_n`
#'   most-common paths; groups with rarer paths are tagged with
#'   `NA_character_` for the path column. NULL = keep every distinct
#'   path.
#' @return The input data.frame with a new column `path` of type
#'   character. Order of rows preserved. Groups with no observed
#'   treatment trajectory get `NA`.
#'
#' @examples
#' p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
#'                             tau_profile = c(0.5, 1.0),
#'                             seed = 17L)
#' aug <- didgpu_compute_paths(p, "Y", "unit", "period", "D",
#'                              effects = 2L)
#' table(aug$path, useNA = "ifany")
#' @export
didgpu_compute_paths <- function(df, outcome, group, time, treatment,
                                   effects = 1L, top_n = NULL) {
  stopifnot(is.data.frame(df) || data.table::is.data.table(df))
  for (nm in c("outcome", "group", "time", "treatment")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v)) {
      stop("`", nm, "` must be a single non-empty character column name.")
    }
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  effects <- as.integer(effects)
  stopifnot(effects >= 1L)
  if (!is.null(top_n)) {
    stopifnot(is.numeric(top_n), length(top_n) == 1L, top_n >= 0L)
    top_n <- as.integer(top_n)
  }
  original_class <- if (data.table::is.data.table(df)) "data.table" else "data.frame"
  d <- data.table::as.data.table(df)
  # Stash the original group/time/treatment columns under internal
  # names so we can rename them back exactly. Keep the ORIGINAL group
  # column on d so the eventual merge has matching column types.
  d[, group_path_XX := d[[group]]]
  d[, time_path_XX := d[[time]]]
  d[, treatment_path_XX := d[[treatment]]]
  # Stable integer encodings for sorting.
  d[, group_path_int_XX := as.integer(factor(group_path_XX,
                                               levels = sort(unique(group_path_XX))))]
  d[, time_path_int_XX := as.integer(factor(time_path_XX,
                                              levels = sort(unique(time_path_XX))))]

  # Per-group baseline treatment (value at t_min).
  t_min_int <- min(d$time_path_int_XX)
  d[, d_sq_path_XX := {
      idx <- which(time_path_int_XX == t_min_int)
      if (length(idx) > 0L) treatment_path_XX[idx[1L]] else NA
    },
    by = group_path_int_XX]

  # Per-group first-switch period F_g (smallest int-time with
  # treatment != d_sq). NA if no switch.
  d[, F_g_path_XX := {
      ok <- !is.na(treatment_path_XX) & !is.na(d_sq_path_XX) &
            treatment_path_XX != d_sq_path_XX
      if (any(ok)) as.integer(min(time_path_int_XX[ok])) else NA_integer_
    },
    by = group_path_int_XX]

  T_max_int <- max(d$time_path_int_XX)

  # Per-group path string: D at F_g - 1 (baseline), F_g, F_g + 1, ...,
  # F_g + effects - 1 (all in int-time space). Never-switchers get
  # "no_switch".
  paths_per_g <- d[, {
    fg <- F_g_path_XX[1L]
    g_orig <- group_path_XX[1L]
    if (is.na(fg)) {
      list(group_orig = g_orig, path = "no_switch")
    } else {
      idx <- (fg - 1L):(fg + effects - 1L)
      idx <- idx[idx >= 1L & idx <= T_max_int]
      tvec <- treatment_path_XX[match(idx, time_path_int_XX)]
      if (length(tvec) == 0L) {
        list(group_orig = g_orig, path = "no_obs")
      } else {
        list(group_orig = g_orig,
             path = paste(.fmt_path_value(tvec), collapse = ","))
      }
    }
  }, by = group_path_int_XX]

  # Optional: keep only the top_n most common paths. top_n = 0 means
  # "drop all" — all paths become NA, leaving zero usable rows.
  if (!is.null(top_n)) {
    if (top_n == 0L) {
      paths_per_g[, path := NA_character_]
    } else {
      counts <- sort(table(paths_per_g$path), decreasing = TRUE)
      keep <- names(counts)[seq_len(min(top_n, length(counts)))]
      paths_per_g[, path := ifelse(path %in% keep, path, NA_character_)]
    }
  }

  # Clean up internal scratch on d.
  d[, c("group_path_XX", "time_path_XX", "treatment_path_XX",
        "group_path_int_XX", "time_path_int_XX",
        "d_sq_path_XX", "F_g_path_XX") := NULL]

  # Build a small lookup keyed by the ORIGINAL group column for a
  # type-safe merge.
  lookup <- data.table::data.table(
    grp = paths_per_g$group_orig,
    path = paths_per_g$path
  )
  data.table::setnames(lookup, "grp", group)
  out <- merge(d, lookup, by = group, all.x = TRUE)
  out <- out[order(out[[group]], out[[time]])]

  if (original_class == "data.frame") {
    as.data.frame(out)
  } else {
    out
  }
}


#' Run didgpu separately per treatment trajectory
#'
#' Wrapper: calls `didgpu_compute_paths()` to add a `path` column to the
#' panel, then runs `didgpu_by(df, by_var = "path", ...)` on the result.
#' Subgroups (rows with `path = NA` because they're not in `top_n`) are
#' dropped before estimation.
#'
#' @inheritParams didgpu_compute_paths
#' @param ... Additional arguments forwarded to `didgpu()` (e.g.
#'   `placebo`, `bootstrap_reps`, `normalized`).
#' @param checkpoint_dir If non-NULL, each subgroup writes into a
#'   subdirectory of this path named after the path string.
#' @param verbose Logical. Print one line per subgroup as it starts.
#' @return A `didgpu_by_result`: named list of `didgpu_result` objects,
#'   one per path.
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 100L, n_periods = 12L,
#'                             tau_profile = c(0.5, 1.0),
#'                             seed = 17L)
#' fit_by_path <- didgpu_by_path(
#'   p, outcome = "Y", group = "unit", time = "period", treatment = "D",
#'   effects = 2L, top_n = 3L,
#'   bootstrap_reps = 0L, backend = "r", verbose = FALSE
#' )
#' print(fit_by_path)
#' }
#' @export
didgpu_by_path <- function(df, outcome, group, time, treatment,
                             effects = 1L, top_n = NULL,
                             ...,
                             checkpoint_dir = NULL, verbose = TRUE) {
  augmented <- didgpu_compute_paths(df, outcome = outcome, group = group,
                                      time = time, treatment = treatment,
                                      effects = effects, top_n = top_n)
  augmented <- augmented[!is.na(augmented$path), , drop = FALSE]
  if (nrow(augmented) == 0L) {
    stop("After path computation, no rows remain (try a larger `top_n` ",
         "or NULL).")
  }
  didgpu_by(augmented, by_var = "path",
             outcome = outcome, group = group,
             time = time, treatment = treatment,
             effects = effects,
             ...,
             checkpoint_dir = checkpoint_dir, verbose = verbose)
}


# Format a treatment value for the path string. Integers stay as
# integers ("0", "1"); floats get 3 significant digits to avoid noise.
#' @keywords internal
#' @noRd
.fmt_path_value <- function(x) {
  if (is.numeric(x)) {
    is_int <- !is.na(x) & abs(x - round(x)) < 1e-9
    out <- character(length(x))
    out[is.na(x)] <- "NA"
    out[is_int & !is.na(x)] <- as.character(as.integer(round(x[is_int & !is.na(x)])))
    other <- !is_int & !is.na(x)
    out[other] <- formatC(x[other], digits = 3, format = "g")
    out
  } else {
    as.character(x)
  }
}
