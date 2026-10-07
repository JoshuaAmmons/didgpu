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
#' @param bootstrap_reps Integer. Number of bootstrap iterations, default
#'   `0`. Standard errors, confidence intervals and the joint nullity
#'   tests are computed analytically from the estimator's asymptotic
#'   linear representation, exactly as `DIDmultiplegtDYN` computes them,
#'   so no resampling is needed for them. Set this above zero only if you
#'   specifically want bootstrap quantities; it does not change the
#'   reported `SE` column.
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
#' @param verbose Logical. Print one line per completed cell. Default
#'   `FALSE`, so that a fit prints what `DIDmultiplegtDYN` prints and
#'   nothing else.
#' @param on_iter Optional function called as `on_iter(iter, value)`
#'   after each cell save. Useful for emailing/logging.
#' @param bootstrap `DIDmultiplegtDYN`'s way of asking for a bootstrap:
#'   the number of replications, `c(reps, seed)`, `c(reps = , seed = )`
#'   or `list(reps, seed)`. When given it sets `bootstrap_reps` (and
#'   `seed`).
#' @param graph_off Logical. As in `DIDmultiplegtDYN`, the event-study
#'   graph is drawn after the fit unless this is `TRUE`; either way it is
#'   returned in `$plot` (a ggplot, when ggplot2 and cowplot are
#'   installed).
#' @param ggplot_args A list of ggplot2 layers added to that graph.
#' @param less_conservative_se Logical. As in `DIDmultiplegtDYN`: demean each switcher's outcome change within switchers sharing its full treatment path, rather than its (baseline, switch date, dose) cell, which gives smaller standard errors.
#' @param more_granular_demeaning Logical. Turns `less_conservative_se` on, as in `DIDmultiplegtDYN`.
#' @param drop_if_d_miss_before_first_switch Logical. Drop a group's outcomes from the first period, before its first switch, at which its treatment is missing.
#' @param effects_equal `TRUE`, `"all"` or `"lb, ub"`: test that all effects, or effects `lb` to `ub`, are equal. Reported in `$results$p_equality_effects`.
#' @param predict_het_hc2bm Logical. With `predict_het`, cluster its HC2 standard errors (by `cluster`, or by group) with the Bell-McCaffrey small-sample adjustment.
#' @param normalized_weights Logical, requires `normalized = TRUE`. Report the weight each normalized effect puts on the current treatment and on each of its lags, in `$normalized_weights`.
#' @param save_results Path of a CSV file to which the estimates, standard errors, confidence intervals and counts are written.
#' @param save_sample Logical. Return, in `$save_sample`, the data with `did_sample` (control, switcher-in or switcher-out) and `did_effect` (the effect a switcher's cell is used for).
#' @param design `c(p, path)`: tabulate the switchers' treatment paths over the effects' window, the most common paths that together cover share `p` of switchers, in `$design`; `path` is `"console"` or an Excel file.
#' @param date_first_switch `c(option, path)`: tabulate the switchers' first switch dates, overall (`option = ""`) or by baseline treatment (`"by_baseline_treat"`), in `$date_first_switch`; `path` is `"console"` or an Excel file.
#' @param by Name of a group-level, time-invariant variable: the estimation is run separately for each of its levels, in `$by_level_1`, `$by_level_2`, ..., with a combined plot.
#' @param by_path A positive integer, `-1` or `"all"`: the estimation is run separately for the switchers following each of the most common `by_path` treatment paths (all of them with `-1` or `"all"`), with the not-yet-switched groups of the same baseline as controls.
#' @param reset Non-negative integer, as in `DIDmultiplegtDYN` 2.4.0. With `reset = k > 0`, the panel is restricted to the groups with the most non-missing treatment observations, and a group is split into a new sub-group whenever its treatment has stayed unchanged for `k` periods after a change; standard errors are then clustered by the original groups unless `cluster` is given.
#' @param avg_time_periods Logical, as in `DIDmultiplegtDYN` 2.4.0: report the average number of periods over which the effect of a dose is accumulated, in `$avg_time_periods` and `$avg_cumul`.
#'
#' @details
#' `didgpu()` behaves as `DIDmultiplegtDYN::did_multiplegt_dyn()` does
#' around the estimates too: it stops with the reference's error when no
#' effect can be estimated, prints the reference's messages (horizons cut
#' to what the data allow, effects or placebos that cannot be estimated,
#' fixest's note on rows dropped from the control regressions), records
#' the joint tests' caveats in `$results$vcov_warnings`, and prints the
#' same tables.
#'
#' @return An object of class `didgpu_result` (see [print.didgpu_result()]).
#'   Its `$results` holds three matrices, each laid out exactly as
#'   `DIDmultiplegtDYN::did_multiplegt_dyn()` lays out its own:
#'   \describe{
#'     \item{`Effects`}{one row per event-time `k`: the cumulative DID
#'       estimate, its SE, CI bounds and switcher / observation counts.}
#'     \item{`ATE`}{a single row, the reference's `Av_tot_eff`: the average
#'       total effect **per unit of treatment**, namely
#'       `sum_k (N_k * DID_k) / sum_k (N_k * delta_k)`, where `N_k` is the
#'       switcher mass at event-time `k` and `delta_k` the average
#'       treatment change among those switchers. That denominator is what
#'       makes it a per-unit quantity, so the ATE is **not** in general the
#'       switcher-weighted mean of the `Effects` column: the two coincide
#'       only when treatment is binary and absorbing, where every `delta_k`
#'       is 1. `normalized` does not affect it, and it is suppressed (`NA`)
#'       under `trends_lin`.}
#'     \item{`Placebos`}{one row per placebo, same columns as `Effects`.}
#'   }
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
    bootstrap_reps = 0L,
    ci_level       = 95,
    seed           = 1L,
    checkpoint_dir = NULL,
    resume         = TRUE,
    backend        = "auto",
    n_workers      = 1L,
    verbose        = FALSE,
    on_iter        = NULL,
    bootstrap      = NULL,
    graph_off      = FALSE,
    ggplot_args    = NULL,
    less_conservative_se = FALSE,
    more_granular_demeaning = FALSE,
    drop_if_d_miss_before_first_switch = FALSE,
    effects_equal  = FALSE,
    predict_het_hc2bm = FALSE,
    normalized_weights = FALSE,
    save_results   = NULL,
    save_sample    = FALSE,
    design         = NULL,
    date_first_switch = NULL,
    by             = NULL,
    by_path        = NULL,
    reset          = 0,
    avg_time_periods = FALSE) {

  # ---- validate ----
  # DIDmultiplegtDYN's own checks and wording first
  # (did_multiplegt_dyn.R, "General syntax checks").
  if (!is.null(cluster) && identical(cluster, group)) cluster <- NULL
  .dcdh_check_args(environment())
  if (!(is.numeric(reset) && length(reset) == 1L && !is.na(reset) &&
        reset %% 1 == 0 && reset >= 0)) {
    stop("Syntax error in reset option. Non-negative integer required.")
  }
  # reset (did_multiplegt_dyn.R:393-399 in 2.4.0): split the groups, and
  # cluster by the original ones unless a cluster was given. Done once,
  # before any by / by_path split.
  if (reset > 0) {
    rs <- .dcdh_reset(df, group, time, treatment, cluster, reset)
    df <- rs$df
    cluster <- rs$cluster
    reset <- 0
  }
  # by / by_path: one run per level (R/by_option.R).
  if (!is.null(by) || !is.null(by_path)) {
    if (!is.null(by_path)) {
      ok <- (is.character(by_path) && length(by_path) == 1L && tolower(by_path) == "all") ||
            (is.numeric(by_path) && length(by_path) == 1L && !is.na(by_path) &&
               by_path %% 1 == 0 && (by_path > 0 || by_path == -1))
      if (!ok) stop("Syntax error in by_path option. Positive integer, -1, or \"all\" required.")
    }
    a <- as.list(environment())
    a <- a[intersect(names(a), names(formals(didgpu)))]
    if (!is.null(bootstrap)) {
      bs <- .dcdh_parse_bootstrap(bootstrap)
      a$bootstrap_reps <- bs$reps
    }
    a$bootstrap_reps <- as.integer(a$bootstrap_reps)
    return(.didgpu_by_levels(a))
  }
  seed_given <- !missing(seed)
  if (!is.null(bootstrap)) {
    bs <- .dcdh_parse_bootstrap(bootstrap)
    bootstrap_reps <- bs$reps
    if (!is.null(bs$seed)) { seed <- bs$seed; seed_given <- TRUE }
  }
  if (!is.null(predict_het)) {
    if (isTRUE(normalized)) {
      stop("The options predict_het and normalized cannot be specified together!")
    }
    if (!is.null(controls)) {
      stop("The options predict_het and controls cannot be specified together!")
    }
  }
  # effects_equal: TRUE, "all" or "lb, ub" (did_multiplegt_dyn.R, the
  # effects_equal block), with the reference's errors.
  if (!is.null(effects_equal) &&
      !(inherits(effects_equal, "logical") || inherits(effects_equal, "character"))) {
    stop("Syntax error in effects_equal option. Logical or string (e.g., 'all' or 'lb, ub') required.")
  }
  effects_equal_lb <- NULL; effects_equal_ub <- NULL
  if (inherits(effects_equal, "character")) {
    if (effects_equal == "all") {
      effects_equal <- TRUE
    } else {
      parts <- strsplit(effects_equal, ",")[[1]]
      if (length(parts) != 2) {
        stop("Syntax error in effects_equal option. Use TRUE, 'all', or 'lb, ub' format (e.g., '2, 5').")
      }
      effects_equal_lb <- suppressWarnings(as.integer(trimws(parts[1])))
      effects_equal_ub <- suppressWarnings(as.integer(trimws(parts[2])))
      if (is.na(effects_equal_lb) | is.na(effects_equal_ub)) {
        stop("Syntax error in effects_equal option. Bounds must be integers.")
      }
      if (effects_equal_ub <= effects_equal_lb | effects_equal_lb < 1) {
        stop("Syntax error in effects_equal option: The bounds specified are out of range.")
      }
      effects_equal <- TRUE
    }
  }
  for (v in c("design", "date_first_switch")) {
    x <- get(v)
    if (!is.null(x) && !(length(x) == 2)) {
      stop(sprintf("Syntax error in %s option. Array with two arguments required.", v))
    }
  }
  if (!is.null(design) && !is.null(continuous)) {
    stop("The design option can not be specified together with the continuous option!")
  }
  if (isTRUE(normalized_weights) && isFALSE(normalized)) {
    stop("normalized option required to compute normalized_weights")
  }
  if (isTRUE(predict_het_hc2bm) && is.null(predict_het)) {
    stop("Option predict_het_hc2bm only available when predict_het is specified.")
  }
  if (isTRUE(same_switchers_pl) && !isTRUE(same_switchers)) {
    stop("The same_switchers_pl option only works if same_switchers is specified as well!")
  }
  for (nm in c("outcome", "group", "time", "treatment")) {
    v <- get(nm)
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  if (!is.null(cluster) && !cluster %in% names(df)) {
    stop("cluster column not in df: ", cluster)
  }
  # Rows without a cluster are dropped, as in the reference
  # (did_multiplegt_main.R:96-98).
  if (!is.null(cluster) && anyNA(df[[cluster]])) {
    df <- df[!is.na(df[[cluster]]), , drop = FALSE]
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
  if (!switchers %in% c("", "in", "out")) {
    stop('`switchers` must be one of "", "in", "out".')
  }
  if (!is.null(predict_het)) {
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
    # A covariate that varies within a group is dropped, with the
    # reference's message (did_multiplegt_main.R:113-126).
    good <- .het_time_invariant(df, het_vars, group, outcome, treatment)
    predict_het <- if (length(good)) list(good, predict_het[[2L]]) else NULL
  }
  if (!is.null(on_iter)) stopifnot(is.function(on_iter))

  # What DIDmultiplegtDYN says about the standard errors before it starts
  # (did_multiplegt_dyn.R:186-191); said once, before the levels, under by.
  in_by <- isTRUE(getOption("didgpu.in_by"))
  if (!in_by && bootstrap_reps > 0L && is.null(continuous)) {
    message("did_multiplegt_dyn computes by default analytical standard ",
            "errors - in most cases, there is no need to use the bootstrap ",
            "option.\nBootstrapping is a much slower alternative and we ",
            "recommend it only in combination with the continuous option.")
  }
  if (!in_by && bootstrap_reps == 0L && !is.null(continuous)) {
    message("You specified the continuous option without the bootstrap ",
            "option. \nPlease be aware that we recommend to compute ",
            "bootstraped standard errors when you are using the continuous ",
            "option as the analytical standard errors can be liberal in that ",
            "case.")
  }

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
    less_conservative_se = less_conservative_se,
    more_granular_demeaning = more_granular_demeaning,
    drop_if_d_miss_before_first_switch = drop_if_d_miss_before_first_switch,
    effects_equal = isTRUE(effects_equal),
    predict_het_hc2bm = isTRUE(predict_het_hc2bm),
    normalized_weights = isTRUE(normalized_weights),
    save_sample = isTRUE(save_sample),
    design = design, date_first_switch = date_first_switch,
    avg_time_periods = isTRUE(avg_time_periods),
    effects_equal_lb = effects_equal_lb, effects_equal_ub = effects_equal_ub,
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
      meta <- c(args, list(panel_hash = ph, package_version = pkgv,
                           cell_rev = .didgpu_cell_rev))
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

  # Helper: run one iter.
  .run_iter <- function(iter) {
    iter_seed <- if (iter == 0L) 0L else (seed + iter)
    t0 <- Sys.time()
    value <- .fit_or_empty(fit_one, df, args, iter_seed)
    value$wall_seconds_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    value
  }
  # Helper: keep one finished cell.
  .store <- function(iter, value) {
    if (!is.null(checkpoint_dir)) {
      .save_cell(checkpoint_dir, b = iter, value = value,
                 wall_seconds = value$wall_seconds_total)
    } else {
      in_memory[[as.character(iter)]] <<- value
    }
    if (!is.null(on_iter)) on_iter(iter, value)
  }

  # The point estimate always comes first. Everything DIDmultiplegtDYN
  # says while estimating is said here, from cell 0 -- also when that
  # cell is read back from a checkpoint -- and before any bootstrap work,
  # which is where the reference says it.
  if (0L %in% todo) {
    value <- .run_iter(0L)
    .store(0L, value)
    if (verbose) message(sprintf("[didgpu] cell b=%-5d  %.2fs   (%d/%d total)",
                                  0L, value$wall_seconds_total,
                                  n_done0 + 1L, n_total))
    n_done0 <- n_done0 + 1L
  }
  cell0 <- if (!is.null(checkpoint_dir)) .load_cell0(checkpoint_dir) else in_memory[["0"]]
  .emit_notes(attr(.aggregate_to_result(list("0" = cell0), args, ph),
                   "didgpu_notes"))

  boot_todo <- setdiff(todo, 0L)
  if (bootstrap_reps > 0L) {
    if (seed_given) {
      message(sprintf("\nBootstrap, %.0f reps (seed = %.0f):", bootstrap_reps, seed))
    } else {
      message(sprintf("\nBootstrap, %.0f reps:", bootstrap_reps))
    }
  }
  # The reference's progress line: a dot per rep, the count every 5th,
  # a line break every 70th and at the end.
  .progress <- function(j) {
    cat(".")
    if (j %% 5 == 0) cat(sprintf("%.0f", j))
    if (j %% 70 == 0) cat("\n")
    if (j == bootstrap_reps) cat("\n")
  }

  if (n_workers > 1L && length(boot_todo) > 1L) {
    if (verbose) message(sprintf("[didgpu] dispatching %d bootstrap iters across %d workers",
                                  length(boot_todo), n_workers))
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Workers load didgpu from the library this session loaded it from,
    # not whatever copy their default library path finds first. The
    # loader must not be a closure over didgpu's namespace: unserializing
    # one on a worker would load the namespace from the default path
    # before the loader runs.
    load_here <- function(lib) {
      suppressMessages(library("didgpu", lib.loc = lib, character.only = TRUE))
      NULL
    }
    environment(load_here) <- globalenv()
    parallel::clusterCall(cl, load_here, dirname(find.package("didgpu")))
    parallel::clusterExport(cl, c("df", "args", "seed"),
                             envir = environment())
    # We don't need clusterSetRNGStream — each iter sets its own seed
    # explicitly via .cluster_resample(set.seed(iter_seed, "Mersenne-Twister")).
    worker_run <- function(iter) {
      iter_seed <- seed + iter
      # `.resolve_backend` is internal; on a worker we need
      # utils::getFromNamespace because the worker only sees exported names.
      resolve <- utils::getFromNamespace(".resolve_backend", "didgpu")
      fit_one_local <- resolve(args$backend %||% "auto")
      t0 <- Sys.time()
      fit_or_empty <- utils::getFromNamespace(".fit_or_empty", "didgpu")
      v <- fit_or_empty(fit_one_local, df, args, iter_seed)
      v$wall_seconds_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      list(iter = iter, value = v)
    }
    results <- parallel::parLapplyLB(cl, boot_todo, worker_run)
    # Save cells sequentially after workers return.
    for (r in results) {
      .store(r$iter, r$value)
      .progress(r$iter)
    }
    if (verbose) {
      total <- sum(sapply(results, function(r) r$value$wall_seconds_total))
      message(sprintf("[didgpu] %d parallel iters done (worker-sum %.2fs)",
                      length(results), total))
    }
  } else {
    for (idx in seq_along(boot_todo)) {
      iter <- boot_todo[idx]
      value <- .run_iter(iter)
      .store(iter, value)
      if (verbose) {
        message(sprintf("[didgpu] cell b=%-5d  %.2fs   (%d/%d total)",
                        iter, value$wall_seconds_total,
                        n_done0 + idx, n_total))
      }
      .progress(iter)
    }
  }

  result <- .finalize_result(checkpoint_dir, in_memory, args, ph)
  # The event-study graph, kept in $plot and drawn unless graph_off, as
  # DIDmultiplegtDYN does (did_multiplegt_dyn.R:540-583).
  result$plot <- .dcdh_graph(result$results, ggplot_args)
  # save_sample: the user's rows with the reference's two tags
  # (did_multiplegt_dyn.R:543-550), and save_results: its results matrix
  # as a CSV (did_multiplegt_main.R:2301-2304).
  if (isTRUE(save_sample) && !is.null(cell0$save_sample_df)) {
    result$save_sample <- cell0$save_sample_df
  } else if (isTRUE(save_sample) && !is.null(cell0$save_sample)) {
    ss <- cell0$save_sample
    names(ss) <- c(group, time, "did_sample", "did_effect")
    df_m <- merge(as.data.frame(df), ss, by = c(group, time))
    result$save_sample <- df_m[order(df_m[[group]], df_m[[time]]), ]
  }
  if (!is.null(save_results)) .dcdh_save_results(result$results, save_results)
  # design / date_first_switch tables, and their Excel file if one was named.
  if (!is.null(cell0$desc_tables$design)) {
    result$design <- cell0$desc_tables$design
    result$design$design_path <- design[2]
  }
  if (!is.null(cell0$desc_tables$dfs)) {
    result$date_first_switch <- cell0$desc_tables$dfs
    result$date_first_switch$dfs_path <- date_first_switch[2]
  }
  .dcdh_write_xlsx(result$design, result$date_first_switch, design, date_first_switch)
  if (isFALSE(graph_off) && !is.null(result$plot)) print(result$plot)
  result
}


# -------- helpers --------

# Shared finalization: collect cells, aggregate into a result object,
# attach checkpoint_dir, set class. The notes were said after cell 0.
.finalize_result <- function(checkpoint_dir, in_memory, args, panel_hash) {
  cells <- if (!is.null(checkpoint_dir)) {
    didgpu_aggregate_cells(checkpoint_dir)
  } else {
    in_memory
  }
  result <- .aggregate_to_result(cells, args, panel_hash)
  attr(result, "didgpu_notes") <- NULL
  result$checkpoint_dir <- if (!is.null(checkpoint_dir)) {
    normalizePath(checkpoint_dir, winslash = "/", mustWork = FALSE)
  } else NA_character_
  class(result) <- c("didgpu_result", "list")
  result
}

# One cell. A resample in which nothing is estimable is an empty cell,
# which the aggregator drops; on the point estimate the reference's
# error stands.
.fit_or_empty <- function(fit_one, df, args, iter_seed) {
  if (iter_seed == 0L) return(fit_one(df, args, iter_seed))
  tryCatch(fit_one(df, args, iter_seed),
           didgpu_no_effect = function(e) {
             list(effects = numeric(0), placebos = numeric(0), ate = NA_real_,
                  iter_seed = as.integer(iter_seed))
           })
}

# predict_het covariates that are constant within every group. The
# reference averages each row's within-group SD over the rows of groups
# with some outcome and treatment, and keeps a variable only if that
# average is 0 (or NA).
.het_time_invariant <- function(df, het_vars, group, outcome, treatment) {
  d <- data.table::as.data.table(df)[, c(group, outcome, treatment, het_vars),
                                     with = FALSE]
  data.table::setnames(d, c(group, outcome, treatment), c("g_", "y_", "d_"))
  d <- d[, if (any(!is.na(y_)) && any(!is.na(d_))) .SD, by = g_]
  good <- character(0)
  for (v in het_vars) {
    s <- d[, list(s = stats::sd(get(v), na.rm = TRUE), n = .N), by = g_]
    s$s[is.na(s$s)] <- 0
    m <- sum(s$s * s$n) / sum(s$n)
    if (is.na(m) || m == 0) {
      good <- c(good, v)
    } else {
      message(sprintf(paste0("The variable %s specified in the option ",
        "predict_het is time-varying, the command will therefore ignore it."),
        v))
    }
  }
  good
}

# The reference's results matrix -- effects, the ATE row, placebos, with
# a Time column -- written where save_results says.
.dcdh_save_results <- function(res, path) {
  m <- rbind(cbind(res$Effects, Time = seq_len(nrow(res$Effects))),
             cbind(res$ATE, Time = 0))
  if (!is.null(res$Placebos)) {
    m <- rbind(m, cbind(res$Placebos, Time = -seq_len(nrow(res$Placebos))))
  }
  suppressWarnings(utils::write.csv(m, path, row.names = TRUE, col.names = TRUE))
  invisible(path)
}

# The point-estimate cell from a checkpoint directory.
.load_cell0 <- function(checkpoint_dir) {
  chk <- didgpu_load_checkpoint(checkpoint_dir)
  i <- which(as.integer(chk$manifest$b) == 0L)[1L]
  if (is.na(i)) stop("Cannot aggregate: cell b=0 (point estimate) is missing ",
                     "from the checkpoint manifest.")
  readRDS(file.path(checkpoint_dir, "cells", chk$manifest$cell_file[i]))
}

# DIDmultiplegtDYN's argument checks for the options didgpu shares with
# it, with its messages (did_multiplegt_dyn.R, "General syntax checks").
.dcdh_check_args <- function(env) {
  get1 <- function(v) get(v, envir = env)
  if (!inherits(get1("df"), "data.frame")) {
    stop("Syntax error in df option. Dataframe object required.")
  }
  for (v in c("outcome", "group", "time", "treatment", "by", "cluster", "weight",
              "switchers", "save_results")) {
    x <- get1(v)
    if (!is.null(x) && !(length(x) == 1 && inherits(x, "character"))) {
      stop(sprintf("Syntax error in %s option. Only one string allowed.", v))
    }
  }
  for (v in c("effects", "ci_level", "continuous")) {
    x <- get1(v)
    if (!is.null(x) && !(is.numeric(x) && length(x) == 1L && !is.na(x) &&
                           x %% 1 == 0 && x > 0)) {
      stop(sprintf("Syntax error in %s option. Positive integer required.", v))
    }
  }
  x <- get1("placebo")
  if (!(is.numeric(x) && length(x) == 1L && !is.na(x) && x %% 1 == 0 && x >= 0)) {
    stop("Syntax error in placebo option. Non-negative integer required.")
  }
  x <- get1("predict_het")
  if (!is.null(x) && !(inherits(x, "list") && length(x) == 2)) {
    stop("Syntax error in predict_het option. List with two arguments required.")
  }
  for (v in c("controls", "trends_nonparam")) {
    x <- get1(v)
    if (!is.null(x) && !inherits(x, "character")) {
      stop(sprintf("Syntax error in %s option. String or string array required.", v))
    }
  }
  for (v in c("normalized", "trends_lin", "same_switchers", "same_switchers_pl",
              "dont_drop_larger_lower", "only_never_switchers", "graph_off",
              "less_conservative_se", "more_granular_demeaning",
              "drop_if_d_miss_before_first_switch", "predict_het_hc2bm",
              "normalized_weights", "save_sample", "avg_time_periods")) {
    x <- get1(v)
    if (!is.null(x) && !inherits(x, "logical")) {
      stop(sprintf("Syntax error in %s option. Logical required.", v))
    }
  }
  invisible(NULL)
}

# DIDmultiplegtDYN's `bootstrap` argument: reps, c(reps, seed),
# c(reps = , seed = ) or list(reps, seed).
.dcdh_parse_bootstrap <- function(x) {
  if (is.list(x)) {
    nm <- names(x)
    if (length(x) < 1L || length(x) > 2L) {
      stop("Syntax error in bootstrap option. Provide reps (and optional seed).")
    }
    reps <- if (!is.null(nm) && "reps" %in% nm) x[["reps"]] else x[[1]]
    seed <- if (!is.null(nm) && "seed" %in% nm) x[["seed"]] else if (length(x) == 2L) x[[2]] else NULL
  } else if (is.numeric(x)) {
    nm <- names(x)
    if (length(x) < 1L || length(x) > 2L) {
      stop("Syntax error in bootstrap option. Provide reps as a scalar or c(reps, seed) of length 2.")
    }
    reps <- if (!is.null(nm) && "reps" %in% nm) unname(x[["reps"]]) else unname(x[1])
    seed <- if (!is.null(nm) && "seed" %in% nm) unname(x[["seed"]]) else if (length(x) == 2L) unname(x[2]) else NULL
  } else {
    stop("Syntax error in bootstrap option. Provide a numeric scalar, c(reps, seed), or named c(reps = ..., seed = ...).")
  }
  if (!(is.numeric(reps) && length(reps) == 1L && !is.na(reps) && reps %% 1 == 0 && reps > 1)) {
    stop("Syntax error in bootstrap option: reps must be an integer greater than 1.")
  }
  if (!is.null(seed) && !(is.numeric(seed) && length(seed) == 1L && !is.na(seed) && seed %% 1 == 0)) {
    stop("Syntax error in bootstrap option: seed must be an integer.")
  }
  list(reps = as.integer(reps), seed = if (is.null(seed)) NULL else as.integer(seed))
}

.as_int1 <- function(x, name, min = NULL) {
  if (length(x) != 1L || is.na(x)) stop("`", name, "` must be a length-1 integer.")
  v <- suppressWarnings(as.integer(x))
  if (is.na(v) || v != x) stop("`", name, "` must be coercible to integer; got ", x)
  if (!is.null(min) && v < min) stop("`", name, "` must be >= ", min)
  v
}

# Revision of what a saved cell MEANS. Bump it whenever an estimator
# change alters the numbers a cell holds, so checkpoints written before
# the change are refused on resume rather than silently mixed in.
#   1  original
#   2  ATE is Av_tot_eff, per unit of treatment
#   3  analytic SEs and influence vectors carried in the point-estimate cell
#   4  horizon counts exactly as the reference clamps them; the point-
#      estimate cell carries the ATE-row counts, delta_D_avg_total,
#      max_pl / max_pl_gap and the reference's horizon messages
#' @keywords internal
#' @noRd
.didgpu_cell_rev <- 4L

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
  # A checkpoint written by a different version of didgpu may hold cells
  # computed by a different estimator. That has happened twice: the ATE
  # changed when it was corrected to Av_tot_eff, and the SEs changed when
  # they moved from bootstrap to analytic. Resuming into such a directory
  # silently mixes old and new numbers in one result, which is worse than
  # either, so refuse rather than warn.
  # The package version alone cannot do this job: several estimator
  # changes shipped under the same version string (0.1.2), so a
  # version-only check would wave stale cells through. Every checkpoint
  # therefore also records .didgpu_cell_rev, bumped whenever a cell's
  # contents change meaning, and a checkpoint without one predates the
  # stamp and is treated as stale.
  have_r <- meta[["cell_rev"]]
  have_v <- meta[["package_version"]]
  if (is.null(have_r) || !isTRUE(as.integer(have_r) == .didgpu_cell_rev)) {
    stop("Checkpoint was written by a different didgpu build (didgpu ",
         have_v %||% "unknown", ", estimator revision ",
         have_r %||% "none", "; this build is revision ", .didgpu_cell_rev,
         "). Its cells may come from a different estimator, so resuming ",
         "would mix old and new numbers. Re-run into a fresh ",
         "checkpoint_dir, or pass resume = FALSE to discard and recompute.",
         call. = FALSE)
  }
  # On/off options that change what a cell holds. A checkpoint that does
  # not record one was written with it off.
  for (k in c("normalized", "trends_lin", "same_switchers", "same_switchers_pl",
              "only_never_switchers", "dont_drop_larger_lower",
              "less_conservative_se", "more_granular_demeaning",
              "drop_if_d_miss_before_first_switch")) {
    if (!identical(isTRUE(as.logical(meta[[k]])), isTRUE(args[[k]]))) {
      stop("Checkpoint config mismatch on `", k, "`: manifest has ",
           isTRUE(as.logical(meta[[k]])), ", call has ", isTRUE(args[[k]]), ".")
    }
  }
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
