# ============================================================================
# didgpu_summarize_panel: pre-flight check of a panel's shape and feasibility
# before fitting. Lets users see at a glance:
#   - Is the panel balanced?
#   - What are the baseline treatment cohorts and their sizes?
#   - How many units actually switch, and in which direction?
#   - When do they switch (F_g distribution)?
#   - What are the feasible effect / placebo horizons?
#
# The output is also useful for catching data issues: e.g., a panel where
# everyone switches at the same period (estimator fails), or a panel with
# almost no never-treated units (low control mass).
# ============================================================================


#' Estimate the runtime of a planned didgpu() call
#'
#' Times a couple of single-iter fits and extrapolates to the
#' bootstrap_reps total, accounting for parallel workers and any
#' already-completed cells in the checkpoint dir. Useful for deciding
#' whether to grab coffee or to scale down `bootstrap_reps`.
#'
#' @inheritParams didgpu
#' @param probes Integer. Number of timing probes to average (default 2).
#'
#' @return Invisibly a list with `wall_per_iter` (median seconds per
#'   point-estimate fit), `n_remaining` (cells still to do), `n_workers`,
#'   `total_seconds`, `total_human` (formatted string).
#' @examples
#' p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L,
#'                             tau_profile = c(0.5, 1.0),
#'                             seed = 17L)
#' info <- didgpu_estimate_runtime(p, "Y", "unit", "period", "D",
#'                                  effects = 2L, bootstrap_reps = 100L,
#'                                  backend = "r", probes = 1L)
#' info$total_human
#' @export
didgpu_estimate_runtime <- function(
    df, outcome, group, time, treatment,
    effects = 1L, placebo = 0L,
    cluster = NULL, controls = NULL, weight = NULL,
    trends_nonparam = NULL,
    only_never_switchers = FALSE, same_switchers = FALSE,
    dont_drop_larger_lower = FALSE, switchers = "",
    bootstrap_reps = 100L,
    checkpoint_dir = NULL,
    backend = "auto",
    n_workers = 1L, probes = 2L) {

  # How many cells remain to do? Same accounting as didgpu() itself.
  n_total <- bootstrap_reps + 1L
  n_done <- 0L
  if (!is.null(checkpoint_dir) &&
      file.exists(file.path(checkpoint_dir, "manifest.csv"))) {
    chk <- didgpu_load_checkpoint(checkpoint_dir)
    n_done <- nrow(chk$manifest)
  }
  n_remaining <- n_total - n_done

  if (n_remaining == 0L) {
    if (n_workers != 0L)
      message("All cells already done -- nothing to estimate.")
    return(invisible(list(wall_per_iter = 0, n_remaining = 0L,
                           n_workers = n_workers,
                           total_seconds = 0, total_human = "0s")))
  }

  # Time `probes` point-estimate fits (no bootstrap, no checkpoint).
  message(sprintf("[didgpu] timing %d probe fit(s)...", probes))
  times <- numeric(probes)
  for (i in seq_len(probes)) {
    t0 <- Sys.time()
    invisible(didgpu(df = df, outcome = outcome, group = group, time = time,
                      treatment = treatment, effects = effects, placebo = placebo,
                      cluster = cluster, controls = controls, weight = weight,
                      trends_nonparam = trends_nonparam,
                      only_never_switchers = only_never_switchers,
                      same_switchers = same_switchers,
                      dont_drop_larger_lower = dont_drop_larger_lower,
                      switchers = switchers, bootstrap_reps = 0L,
                      backend = backend, verbose = FALSE))
    times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  wall_per_iter <- stats::median(times)

  # Total = remaining iters / n_workers * wall_per_iter (+ ~1s cluster overhead).
  effective_workers <- max(1L, n_workers)
  total_seconds <- (n_remaining / effective_workers) * wall_per_iter
  if (effective_workers > 1L) total_seconds <- total_seconds + 1.5  # cluster setup

  total_human <- if (total_seconds < 60) {
    sprintf("%.1fs", total_seconds)
  } else if (total_seconds < 3600) {
    sprintf("%.1f min", total_seconds / 60)
  } else {
    sprintf("%.1f hours", total_seconds / 3600)
  }

  message(sprintf("[didgpu] %d cells remaining; %.2fs/iter; %d worker%s -> ~%s total",
                  n_remaining, wall_per_iter, effective_workers,
                  if (effective_workers == 1L) "" else "s", total_human))

  invisible(list(
    wall_per_iter = wall_per_iter,
    n_remaining = as.integer(n_remaining),
    n_workers = as.integer(effective_workers),
    total_seconds = total_seconds,
    total_human = total_human
  ))
}


#' Summarize a panel before fitting
#'
#' Computes summary statistics that let you sanity-check a panel and
#' know what the estimator can recover. Prints a human-readable
#' summary and returns the underlying numbers invisibly.
#'
#' @param df A panel.
#' @param outcome,group,time,treatment Column names.
#' @param verbose Print the summary. Default TRUE.
#'
#' @return Invisibly, a list with `n_units`, `n_periods`, `n_rows`,
#'   `is_balanced`, `d_sq_dist` (table), `n_never_change`,
#'   `n_switchers_in`, `n_switchers_out`, `f_g_dist` (table),
#'   `max_effects`, `max_placebo`, `na_outcome`, `na_treatment`.
#' @examples
#' p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L,
#'                             tau_profile = c(0.5, 1.0),
#'                             seed = 17L)
#' s <- didgpu_summarize_panel(p, "Y", "unit", "period", "D")
#' s$n_switchers_in
#' s$max_effects
#' @export
didgpu_summarize_panel <- function(df, outcome, group, time, treatment,
                                     verbose = TRUE) {
  for (nm in c("outcome", "group", "time", "treatment")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v)) {
      stop("`", nm, "` must be a single non-empty character column name.")
    }
    if (!v %in% names(df)) stop("column not in df: ", v)
  }

  d <- data.table::as.data.table(df)
  data.table::setnames(d, c(outcome, group, time, treatment),
                       c("Y_", "G_", "T_", "D_"), skip_absent = FALSE)

  n_units   <- length(unique(d$G_))
  n_periods <- length(unique(d$T_))
  n_rows    <- nrow(d)
  is_balanced <- n_rows == n_units * n_periods

  t_min <- min(d$T_, na.rm = TRUE)
  d_sq  <- d[T_ == t_min, list(d_sq = D_[1L]), by = G_]
  d_sq_dist <- table(d_sq$d_sq, useNA = "ifany")

  # F_g: first period at which D != d_sq for each unit.
  d <- merge(d, d_sq, by = "G_", all.x = TRUE)
  switches <- d[!is.na(D_) & !is.na(d_sq) & D_ != d_sq,
                list(F_g = min(T_)), by = G_]
  T_max <- max(d$T_, na.rm = TRUE)
  fg_all <- merge(d_sq, switches, by = "G_", all.x = TRUE)
  fg_all[is.na(F_g), F_g := as.integer(T_max + 1L)]

  # Switcher direction: avg post-switch D minus d_sq.
  d2 <- merge(d, fg_all[, list(G_, F_g, d_sq2 = d_sq)], by = "G_", all.x = TRUE)
  avg_post <- d2[T_ >= F_g & !is.na(D_),
                  list(avg_post = mean(D_)), by = G_]
  fg_all <- merge(fg_all, avg_post, by = "G_", all.x = TRUE)
  # Never-switchers have F_g = T_max + 1 (the sentinel we wrote earlier).
  is_switcher <- fg_all$F_g <= T_max
  fg_all[, direction := ifelse(is.na(avg_post) | !is_switcher, NA_integer_,
                       ifelse(avg_post > d_sq, 1L,
                       ifelse(avg_post < d_sq, 0L, NA_integer_)))]

  n_never_change   <- sum(!is_switcher)
  n_switchers_in   <- sum(fg_all$direction == 1L, na.rm = TRUE)
  n_switchers_out  <- sum(fg_all$direction == 0L, na.rm = TRUE)

  f_g_dist <- table(fg_all$F_g[is_switcher])

  # Feasible horizons (only switcher units contribute).
  L_g <- ifelse(is_switcher,
                pmax(0L, as.integer(T_max - fg_all$F_g + 1L)), 0L)
  max_effects <- max(L_g, na.rm = TRUE)

  L_g_placebo <- ifelse(is_switcher & fg_all$F_g >= 3L,
                         pmin(L_g, as.integer(fg_all$F_g - 2L)), NA_integer_)
  max_placebo <- suppressWarnings(max(L_g_placebo, na.rm = TRUE))
  if (!is.finite(max_placebo)) max_placebo <- 0L

  na_outcome   <- sum(is.na(d$Y_))
  na_treatment <- sum(is.na(d$D_))

  if (verbose) {
    cat(sprintf("Panel: %d units x %d periods = %d cells (observed: %d)\n",
                n_units, n_periods, n_units * n_periods, n_rows))
    if (!is_balanced) {
      cat(sprintf("  WARNING: panel is UNBALANCED (missing %d cells).\n",
                  n_units * n_periods - n_rows))
    } else {
      cat("  Panel is balanced.\n")
    }
    cat(sprintf("\nBaseline (t = %s) treatment d_sq distribution:\n", t_min))
    print(d_sq_dist)
    cat(sprintf("\nSwitcher composition (%d units total):\n", n_units))
    cat(sprintf("  Never-switchers (controls): %d\n", n_never_change))
    cat(sprintf("  Switchers in   (S_g = 1):   %d\n", n_switchers_in))
    cat(sprintf("  Switchers out  (S_g = 0):   %d\n", n_switchers_out))
    if (length(f_g_dist) > 0L) {
      cat("\nF_g (first-switch period) among switchers:\n")
      print(f_g_dist)
    }
    cat(sprintf("\nFeasible horizons:\n"))
    cat(sprintf("  max effects = %d  (longest post-switch window)\n", max_effects))
    cat(sprintf("  max placebo = %d  (longest pre-switch window for placebos)\n",
                max_placebo))
    if (na_outcome > 0 || na_treatment > 0) {
      cat(sprintf("\nMissing values: %d outcome NAs, %d treatment NAs\n",
                  na_outcome, na_treatment))
    }
    if (n_never_change == 0L && (n_switchers_in == 0L || n_switchers_out == 0L)) {
      cat("\n  WARNING: no never-switchers AND only one switcher direction --\n",
          "  estimation will rely on cross-cohort comparisons within d_sq.\n",
          sep = "")
    }
    if (length(f_g_dist) == 1L) {
      cat("\n  WARNING: all switchers switch at the same period --\n",
          "  the DIDmultiplegtDYN estimator cannot be used in this case.\n",
          sep = "")
    }
  }

  invisible(list(
    n_units = n_units, n_periods = n_periods, n_rows = n_rows,
    is_balanced = is_balanced,
    d_sq_dist = d_sq_dist,
    n_never_change = n_never_change,
    n_switchers_in = n_switchers_in,
    n_switchers_out = n_switchers_out,
    f_g_dist = f_g_dist,
    max_effects = as.integer(max_effects),
    max_placebo = as.integer(max_placebo),
    na_outcome = na_outcome,
    na_treatment = na_treatment
  ))
}
