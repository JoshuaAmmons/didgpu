# ============================================================================
# Public entry point: didgpu()
#
# The user-facing function. Validates args, initialises (or loads) the
# checkpoint directory, runs the point estimate + bootstrap loop via the
# chosen backend, saving each cell as it completes, and aggregates the
# committed cells into a result object compatible in structure with
# DIDmultiplegtDYN.
#
# Re-invoking with the same checkpoint_dir and resume = TRUE skips any
# bootstrap iter already in the manifest. If the manifest's recorded
# config differs from the current call (panel hash, effects, placebo,
# bootstrap_reps, seed), the call errors rather than silently producing
# mixed-config output.
# ============================================================================


#' Run a checkpointed, GPU-capable dynamic DiD estimation
#'
#' Reimplementation entry point. See `vignette("reference_internals")`
#' for the design spec and `didgpu_backend_info()` for what's available
#' on this machine.
#'
#' @param df A data.frame (or data.table) panel. Must contain the
#'   columns named by `outcome`, `group`, `time`, `treatment`.
#' @param outcome,group,time,treatment Character. Column names.
#' @param effects Integer. Number of post-treatment event-times to
#'   estimate (e = 1..effects). Must be >= 1.
#' @param placebo Integer. Number of pre-treatment placebos. 0 to skip.
#' @param cluster Character or NULL. Column name for cluster bootstrap;
#'   default NULL clusters by `group`.
#' @param controls Character vector or NULL. Names of covariate columns
#'   to control for. The point estimate adjusts diff_y by the FWL
#'   projection on these covariates' first differences, restricted to
#'   never-switcher rows within each baseline-treatment cohort.
#' @param weight Character or NULL. Name of a per-row weight column.
#'   If NULL, every observation gets weight 1 (the binary case).
#' @param continuous Integer or NULL. If set (typically `1` for linear
#'   continuous treatment), collapses cohorts to a single one (`d_sq = 0`
#'   for all units) and adds polynomial features
#'   `(baseline_D)^1, ..., (baseline_D)^continuous` as auto-controls in
#'   the FWL projection. Use when the treatment column is genuinely
#'   continuous (e.g., dosage, tax rate) so that no two units share an
#'   identical baseline value.
#' @param trends_nonparam Character or NULL. Name of a categorical
#'   column to extend the cohort grouping by. With this set, cohort
#'   averages are computed per `(time, d_sq, trends_nonparam)`
#'   instead of per `(time, d_sq)`, allowing time trends to vary by
#'   the extra grouping (e.g., industry).
#' @param trends_lin Logical. If TRUE, allow group-specific linear
#'   trends: the estimator runs on the first-difference of the outcome
#'   (and any user controls), then for each k = 1..effects returns the
#'   cumulative sum of per-event-time first-difference DIDs (mapping
#'   back to a level effect). Forces `same_switchers = TRUE` and
#'   suppresses ATE. Reference: Section 1.3 of the Web Appendix of
#'   de Chaisemartin and D'Haultfoeuille (2024).
#' @param predict_het Optional list of length 2:
#'   `list(covariates, event_times)`, where `covariates` is a character
#'   vector of time-invariant column names and `event_times` is an integer
#'   vector of which event-times to do the heterogeneity regression for
#'   (use `-1` for all). Regresses the per-group ATE contribution at each
#'   event-time on the covariates (with cohort interaction dummies and
#'   HC1 robust SEs), to characterise effect heterogeneity. Returns the
#'   regression in `result$results$predict_het`. Not compatible with
#'   `normalized = TRUE` (the regression target is the unnormalised DID).
#' @param only_never_switchers Logical. If TRUE, restrict controls to
#'   strictly never-switched units (drop pre-switch rows of units that
#'   eventually switch).
#' @param same_switchers Logical. If TRUE, restrict switcher rows to
#'   units that have valid controls at every event-time `q` in
#'   `1..effects`. Mirrors the reference's same_switchers option.
#' @param same_switchers_pl Logical. If TRUE, additionally restrict the
#'   placebo computations to units that have valid pre-period diff_y at
#'   every placebo horizon `q` in `1..placebo`. Mirrors the reference's
#'   same_switchers_pl option.
#' @param dont_drop_larger_lower Logical. By default (FALSE), drop the
#'   post-non-monotone rows of any unit whose treatment both increases
#'   strictly above and decreases strictly below the baseline. Set TRUE
#'   to keep those rows.
#' @param switchers One of `""` (default — both directions),
#'   `"in"` (only switcher-in units), or `"out"` (only switcher-out
#'   units). Mirrors the `switchers` arg in DIDmultiplegtDYN.
#' @param normalized Logical. If TRUE, divide each per-event-time DID
#'   estimate by its pooled cumulative treatment-change magnitude
#'   `delta_D_k`, so the reported number is a per-unit-of-treatment
#'   effect. For binary on/off treatment this divides the k-th effect
#'   by `k` (so cumulative ATTs become per-period). For continuous or
#'   multivalued treatment, divides by the average per-switcher
#'   cumulative change in actual treatment magnitude over event-times
#'   `1..k`. Mirrors `normalized` in DIDmultiplegtDYN.
#' @param bootstrap_reps Integer. Number of bootstrap iterations.
#' @param ci_level Numeric in (0, 100). Confidence level for CIs.
#' @param seed Integer. RNG seed for bootstrap iter 1 onward (iter 0 is
#'   the deterministic point estimate). Per-iter seed is `seed + iter`.
#' @param checkpoint_dir Character or NULL. If non-NULL, every cell is
#'   saved here and the run is resumable. If NULL, all work is in-memory
#'   and is lost on crash.
#' @param resume Logical. If TRUE and the manifest exists, skip cells
#'   already in it. If FALSE, error rather than overwrite.
#' @param backend One of `"auto"`, `"reference"`, `"r"`, `"cpu"`, `"cuda"`.
#'   See `didgpu_backend_info()`.
#' @param n_workers Integer >= 1. Number of parallel worker processes
#'   for the bootstrap loop. 1 = sequential (default). With n_workers > 1,
#'   uses `parallel::makeCluster()` to distribute bootstrap iters across
#'   cores. Each cell is saved to disk atomically so resume still works.
#'   The point estimate (cell 0) always runs sequentially first.
#' @param verbose Logical. Print one line per completed cell.
#' @param on_iter Optional function called as `on_iter(iter, value)`
#'   after each cell save. Useful for emailing/logging.
#'
#' @return An object of class `didgpu_result` (see [print.didgpu_result()]).
#' @examples
#' # Simulate a small panel with a known event-time profile.
#' p <- didgpu_simulate_panel(
#'   n_units = 40L, n_periods = 10L,
#'   frac_treated = 0.6,
#'   tau_profile = c(0.5, 1.0, 1.2),
#'   sigma = 0.4, seed = 17L
#' )
#'
#' # Point estimate only (no bootstrap), pure-R backend.
#' fit <- didgpu(p, "Y", "unit", "period", "D",
#'                effects = 3L, placebo = 1L,
#'                bootstrap_reps = 0L, backend = "r",
#'                verbose = FALSE)
#' coef(fit)
#'
#' # With a small bootstrap for SEs / CIs.
#' \donttest{
#' fit_boot <- didgpu(p, "Y", "unit", "period", "D",
#'                     effects = 3L, placebo = 1L,
#'                     bootstrap_reps = 20L, seed = 1L,
#'                     backend = "r", verbose = FALSE)
#' print(fit_boot)
#' confint(fit_boot)
#' }
#' @export
didgpu <- function(
    df,
    outcome,
    group,
    time,
    treatment,
    effects        = 1L,
    placebo        = 0L,
    cluster        = NULL,
    controls       = NULL,
    weight         = NULL,
    continuous     = NULL,
    trends_nonparam = NULL,
    trends_lin     = FALSE,
    predict_het    = NULL,
    only_never_switchers = FALSE,
    same_switchers = FALSE,
    same_switchers_pl = FALSE,
    dont_drop_larger_lower = FALSE,
    switchers      = "",
    normalized     = FALSE,
    bootstrap_reps = 100L,
    ci_level       = 95,
    seed           = 1L,
    checkpoint_dir = NULL,
    resume         = TRUE,
    backend        = "auto",
    n_workers      = 1L,
    verbose        = TRUE,
    on_iter        = NULL) {

  # ---- validate ----
  stopifnot(is.data.frame(df) || data.table::is.data.table(df))
  for (nm in c("outcome", "group", "time", "treatment")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v)) {
      stop("`", nm, "` must be a single non-empty character column name.")
    }
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  if (!is.null(cluster) && !cluster %in% names(df)) {
    stop("cluster column not in df: ", cluster)
  }
  if (!is.null(controls)) {
    missing_ctrl <- setdiff(controls, names(df))
    if (length(missing_ctrl)) {
      stop("controls not in df: ", paste(missing_ctrl, collapse = ", "))
    }
  }
  if (!is.null(weight) && !weight %in% names(df)) {
    stop("weight column not in df: ", weight)
  }
  if (!is.null(trends_nonparam) && !trends_nonparam %in% names(df)) {
    stop("trends_nonparam column not in df: ", trends_nonparam)
  }
  effects        <- .as_int1(effects,        "effects",        min = 1L)
  placebo        <- .as_int1(placebo,        "placebo",        min = 0L)
  bootstrap_reps <- .as_int1(bootstrap_reps, "bootstrap_reps", min = 0L)
  seed           <- .as_int1(seed,           "seed",           min = 0L)
  n_workers      <- .as_int1(n_workers,      "n_workers",      min = 1L)
  stopifnot(is.numeric(ci_level), length(ci_level) == 1L,
            ci_level > 0, ci_level < 100)
  if (!is.character(switchers) || length(switchers) != 1L ||
      !switchers %in% c("", "in", "out")) {
    stop('`switchers` must be one of "", "in", "out".')
  }
  if (isTRUE(same_switchers_pl)) {
    if (!isTRUE(same_switchers)) {
      stop("`same_switchers_pl = TRUE` requires `same_switchers = TRUE` ",
           "as well (the placebo-side gate is meaningful only when the ",
           "effects-side gate is also active). Mirrors DIDmultiplegtDYN.")
    }
    if (placebo == 0L) {
      stop("`same_switchers_pl = TRUE` requires `placebo > 0`.")
    }
  }
  if (!is.null(predict_het)) {
    if (!is.list(predict_het) || length(predict_het) != 2L) {
      stop("`predict_het` must be a list of length 2: ",
           "list(covariates, event_times).")
    }
    het_vars <- unlist(predict_het[[1L]])
    if (!is.character(het_vars) || length(het_vars) == 0L) {
      stop("`predict_het[[1]]` must be a non-empty character vector of ",
           "covariate names.")
    }
    miss_het <- setdiff(het_vars, names(df))
    if (length(miss_het)) {
      stop("predict_het covariates not in df: ",
           paste(miss_het, collapse = ", "))
    }
    if (isTRUE(normalized)) {
      message("`normalized = TRUE` together with `predict_het` is not ",
              "supported (the heterogeneity regression is on the ",
              "unnormalised DID). `predict_het` will be ignored.")
      predict_het <- NULL
    }
  }
  if (!is.null(on_iter)) stopifnot(is.function(on_iter))

  fit_one <- .resolve_backend(backend)

  # ---- canonical args bundle ----
  args <- list(
    outcome = outcome, group = group, time = time, treatment = treatment,
    effects = effects, placebo = placebo,
    cluster = cluster, controls = controls, weight = weight,
    continuous = continuous,
    trends_nonparam = trends_nonparam,
    trends_lin = trends_lin,
    predict_het = predict_het,
    only_never_switchers = only_never_switchers,
    same_switchers = same_switchers,
    same_switchers_pl = same_switchers_pl,
    dont_drop_larger_lower = dont_drop_larger_lower,
    switchers = switchers, normalized = normalized,
    ci_level = ci_level,
    bootstrap_reps = bootstrap_reps, seed = seed,
    backend = backend
  )

  ph <- .panel_hash(df, outcome, group, time, treatment)

  # ---- checkpoint init / load ----
  if (!is.null(checkpoint_dir)) {
    manifest_path <- file.path(checkpoint_dir, "manifest.csv")
    if (file.exists(manifest_path) && resume) {
      chk <- didgpu_load_checkpoint(checkpoint_dir)
      .check_meta_compatibility(chk$meta, ph, args)
      manifest <- chk$manifest
      if (verbose) {
        message(sprintf("[didgpu] resuming %s: %d/%d cells already done",
                        checkpoint_dir, nrow(manifest), bootstrap_reps + 1L))
      }
    } else {
      pkgv <- tryCatch(as.character(utils::packageVersion("didgpu")),
                       error = function(e) "0.0.0.dev")
      meta <- c(args, list(panel_hash = ph, package_version = pkgv))
      didgpu_init_checkpoint(checkpoint_dir, meta, force = !resume)
      manifest <- .empty_manifest()
    }
  } else {
    manifest <- .empty_manifest()
  }

  # ---- iter plan ----
  todo <- .cells_todo(manifest, bootstrap_reps)
  n_total <- bootstrap_reps + 1L
  n_done0 <- n_total - length(todo)

  # ---- main loop ----
  in_memory <- list()  # used only when checkpoint_dir is NULL
  # Sequential vs parallel dispatch. Always run cell 0 (the point
  # estimate) sequentially first; the bootstrap cells go to workers.
  use_parallel <- n_workers > 1L && length(setdiff(todo, 0L)) > 1L

  # Helper: run one iter and save it.
  .run_iter <- function(iter) {
    iter_seed <- if (iter == 0L) 0L else (seed + iter)
    t0 <- Sys.time()
    value <- fit_one(df, args, iter_seed)
    value$wall_seconds_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    value
  }

  if (use_parallel) {
    # Sequential pass for cell 0 (if not already done).
    seq_iters <- intersect(todo, 0L)
    par_iters <- setdiff(todo, 0L)

    for (iter in seq_iters) {
      value <- .run_iter(iter)
      if (!is.null(checkpoint_dir)) {
        .save_cell(checkpoint_dir, b = iter, value = value,
                   wall_seconds = value$wall_seconds_total)
      } else {
        in_memory[[as.character(iter)]] <- value
      }
      if (!is.null(on_iter)) on_iter(iter, value)
      if (verbose) message(sprintf("[didgpu] cell b=%-5d  %.2fs  (seq)",
                                    iter, value$wall_seconds_total))
    }

    if (length(par_iters) > 0L) {
      if (verbose) message(sprintf("[didgpu] dispatching %d bootstrap iters across %d workers",
                                    length(par_iters), n_workers))
      cl <- parallel::makeCluster(n_workers)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterEvalQ(cl, suppressMessages(library(didgpu)))
      parallel::clusterExport(cl, c("df", "args", "seed"),
                               envir = environment())
      # We don't need clusterSetRNGStream — each iter sets its own seed
      # explicitly via .cluster_resample(set.seed(iter_seed, "Mersenne-Twister")).
      worker_run <- function(iter) {
        iter_seed <- if (iter == 0L) 0L else (seed + iter)
        # `.resolve_backend` is internal; on a worker we need
        # utils::getFromNamespace because the worker only sees exported names.
        resolve <- utils::getFromNamespace(".resolve_backend", "didgpu")
        fit_one_local <- resolve(args$backend %||% "auto")
        t0 <- Sys.time()
        v <- fit_one_local(df, args, iter_seed)
        v$wall_seconds_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        list(iter = iter, value = v)
      }
      results <- parallel::parLapplyLB(cl, par_iters, worker_run)

      # Save cells sequentially after workers return.
      for (r in results) {
        if (!is.null(checkpoint_dir)) {
          .save_cell(checkpoint_dir, b = r$iter, value = r$value,
                     wall_seconds = r$value$wall_seconds_total)
        } else {
          in_memory[[as.character(r$iter)]] <- r$value
        }
        if (!is.null(on_iter)) on_iter(r$iter, r$value)
      }
      if (verbose) {
        total <- sum(sapply(results, function(r) r$value$wall_seconds_total))
        message(sprintf("[didgpu] %d parallel iters done (worker-sum %.2fs)",
                        length(results), total))
      }
    }
    return(.finalize_result(checkpoint_dir, in_memory, args, ph))
  }

  # ---- sequential path ----
  for (idx in seq_along(todo)) {
    iter <- todo[idx]
    value <- .run_iter(iter)

    if (!is.null(checkpoint_dir)) {
      .save_cell(checkpoint_dir, b = iter, value = value,
                 wall_seconds = value$wall_seconds_total)
    } else {
      in_memory[[as.character(iter)]] <- value
    }

    if (!is.null(on_iter)) on_iter(iter, value)

    if (verbose) {
      message(sprintf("[didgpu] cell b=%-5d  %.2fs   (%d/%d total)",
                      iter, value$wall_seconds_total,
                      n_done0 + idx, n_total))
    }
  }

  # ---- aggregate ----
  cells <- if (!is.null(checkpoint_dir)) {
    didgpu_aggregate_cells(checkpoint_dir)
  } else {
    in_memory
  }

  result <- .aggregate_to_result(cells, args, ph)
  result$checkpoint_dir <- if (!is.null(checkpoint_dir)) {
    normalizePath(checkpoint_dir, winslash = "/", mustWork = FALSE)
  } else NA_character_
  class(result) <- c("didgpu_result", "list")
  result
}


# -------- helpers --------

# Shared finalization for sequential and parallel paths: collect cells,
# aggregate into a result object, attach checkpoint_dir, set class.
.finalize_result <- function(checkpoint_dir, in_memory, args, panel_hash) {
  cells <- if (!is.null(checkpoint_dir)) {
    didgpu_aggregate_cells(checkpoint_dir)
  } else {
    in_memory
  }
  result <- .aggregate_to_result(cells, args, panel_hash)
  result$checkpoint_dir <- if (!is.null(checkpoint_dir)) {
    normalizePath(checkpoint_dir, winslash = "/", mustWork = FALSE)
  } else NA_character_
  class(result) <- c("didgpu_result", "list")
  result
}

.as_int1 <- function(x, name, min = NULL) {
  if (length(x) != 1L || is.na(x)) stop("`", name, "` must be a length-1 integer.")
  v <- suppressWarnings(as.integer(x))
  if (is.na(v) || v != x) stop("`", name, "` must be coercible to integer; got ", x)
  if (!is.null(min) && v < min) stop("`", name, "` must be >= ", min)
  v
}

.check_meta_compatibility <- function(meta, panel_hash, args) {
  must_match <- list(
    panel_hash     = panel_hash,
    outcome        = args$outcome,
    group          = args$group,
    time           = args$time,
    treatment      = args$treatment,
    effects        = args$effects,
    placebo        = args$placebo,
    switchers      = args$switchers %||% "",
    bootstrap_reps = args$bootstrap_reps,
    seed           = args$seed
  )
  for (k in names(must_match)) {
    have <- meta[[k]]
    want <- must_match[[k]]
    # jsonlite restores ints as numeric — compare loosely on numbers
    if (is.numeric(have) || is.numeric(want)) {
      if (!isTRUE(all.equal(unname(have), unname(want)))) {
        stop("Checkpoint config mismatch on `", k, "`: manifest has ",
             format(have), ", call has ", format(want),
             ". Use force = TRUE to discard the old checkpoint, or fix the call.")
      }
    } else if (!identical(as.character(have), as.character(want))) {
      stop("Checkpoint config mismatch on `", k, "`: manifest has '",
           have, "', call has '", want, "'.")
    }
  }
  invisible(TRUE)
}
