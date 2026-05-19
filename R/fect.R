# ============================================================================
# didgpu_fect: counterfactual-prediction DiD estimators
#
# Companion family of estimators to the de Chaisemartin / D'Haultfoeuille
# core in `didgpu()`. fect (Liu, Wang, Xu 2024) estimates Y(0) for
# treated units using a model fit on controls only, then reports
# ATT = mean over (group, time) treated cells of (Y_observed − Y_predicted).
#
# Three estimation methods:
#   "fe"  — two-way fixed effects (iterative demeaning until convergence)
#   "ife" — interactive fixed effects (Bai 2009): low-rank factor model
#           Y = u_i + v_t + lambda_i' F_t + epsilon  with r factors,
#           fit by alternating between demeaning and a rank-r SVD
#   "mc"  — matrix completion (Athey et al. 2021): low-rank recovery
#           via nuclear-norm soft-thresholding on the controls-only matrix
#
# All three benefit from GPU acceleration:
#   - "fe":  iterated weighted sums (group/time fixed-effect updates)
#   - "ife": dense SVD inside each iteration → cuSOLVER cusolverDnDgesvdj
#   - "mc":  soft-thresholded SVD per iteration → same kernel as ife
#
# Status: scaffolded. The R-side public API is stable; the actual
# estimation logic is not yet implemented — each method raises a
# NotImplementedError with a roadmap message. Once implemented, the
# CUDA / Rcpp+Eigen backends compose with the same checkpoint /
# resume / parallel-bootstrap infrastructure that didgpu() already uses.
#
# Design notes for the eventual implementation:
#   - Panel preparation reuses .prep_panel() with adjustments (need the
#     control-only matrix Y_c, the treatment mask M, and a balanced
#     time index — same primitives).
#   - The bootstrap loop reuses .cluster_resample() and the per-cell
#     checkpoint format (one .rds per bootstrap rep, manifest.csv).
#   - The orchestrator reuses didgpu()'s outer loop almost verbatim;
#     only the per-cell fit function differs.
#   - The S3 result class is `didgpu_fect_result`, parallel to
#     `didgpu_result`. The print / coef / confint / vcov / plot methods
#     inherit from a shared base.
# ============================================================================


#' Counterfactual-prediction DiD estimators (fect family)
#'
#' Fits one of three counterfactual-prediction estimators (Liu, Wang, Xu
#' 2024) on a panel and returns an estimate of the average treatment
#' effect on the treated (ATT) along with placebo-style robustness
#' diagnostics. The three methods differ in how they model the
#' counterfactual outcome Y(0):
#'
#' \itemize{
#'   \item `"fe"`  — two-way fixed effects (unit + time FE), no factor
#'         loadings. Fast; requires the strict parallel-trends
#'         assumption.
#'   \item `"ife"` — interactive fixed effects (Bai 2009): unit FE +
#'         time FE + `r` latent factors `lambda_i' F_t`. Relaxes
#'         parallel trends to allow unit-specific time-varying
#'         confounders.
#'   \item `"mc"`  — matrix completion (Athey et al. 2021): low-rank
#'         recovery of the control-only outcome matrix via
#'         nuclear-norm soft-thresholding. No need to choose `r`
#'         in advance.
#' }
#'
#' All three share the per-cell checkpoint + resume + parallel
#' bootstrap infrastructure of [didgpu()].
#'
#' **Status: SCAFFOLDED.** The R-side public API is stable; the actual
#' estimation logic is not yet implemented. Each method currently
#' raises a `NotImplementedError`. The CUDA + Rcpp+Eigen backends for
#' these will land in subsequent releases.
#'
#' @inheritParams didgpu
#' @param method One of `"fe"`, `"ife"`, `"mc"`. Default `"fe"`.
#' @param r Integer. For method `"ife"`, the number of latent factors.
#'   Default `2L`. Ignored for `"fe"` and `"mc"`.
#' @param lambda Numeric or NULL. For method `"mc"`, the nuclear-norm
#'   penalty parameter. If NULL (default), chosen by cross-validation
#'   over a default grid. Ignored for `"fe"` and `"ife"`.
#' @param tol Numeric. Convergence tolerance for the iterative methods
#'   (`"fe"`, `"ife"`, `"mc"`). Default `1e-5`.
#' @param max_iter Integer. Maximum iterations for the iterative
#'   methods. Default `500L`.
#'
#' @return An object of class `didgpu_fect_result` — same general shape
#'   as `didgpu_result` (effects table, placebo table, ATE, S3 methods
#'   work the same way).
#'
#' @references
#' - Liu, L., Wang, Y., and Xu, Y. (2024). "A practical guide to
#'   counterfactual estimators for causal inference with time-series
#'   cross-sectional data." *American Journal of Political Science*.
#' - Bai, J. (2009). "Panel data models with interactive fixed
#'   effects." *Econometrica* 77 (4): 1229-1279.
#' - Athey, S., Bayati, M., Doudchenko, N., Imbens, G., and Khosravi, K.
#'   (2021). "Matrix completion methods for causal panel data models."
#'   *JASA* 116 (536): 1716-1730.
#'
#' @examples
#' \dontrun{
#' # NOT YET IMPLEMENTED -- the example shows the planned interface only.
#' p <- didgpu_simulate_panel(n_units = 100L, n_periods = 20L, seed = 17L)
#' fit <- didgpu_fect(p, outcome = "Y", group = "unit",
#'                     time = "period", treatment = "D",
#'                     method = "ife", r = 2L,
#'                     bootstrap_reps = 100L,
#'                     checkpoint_dir = "checkpoints/fect_run1")
#' print(fit)
#' plot(fit)
#' }
#'
#' @export
didgpu_fect <- function(
    df,
    outcome,
    group,
    time,
    treatment,
    method         = c("fe", "ife", "mc"),
    effects        = 1L,
    r              = 2L,
    lambda         = NULL,
    tol            = 1e-5,
    max_iter       = 500L,
    bootstrap_reps = 100L,
    seed           = 1L,
    checkpoint_dir = NULL,
    backend        = "auto",
    n_workers      = 1L,
    verbose        = TRUE) {

  method <- match.arg(method)

  # Common arg validation (mirrors didgpu()).
  stopifnot(is.data.frame(df) || data.table::is.data.table(df))
  for (nm in c("outcome", "group", "time", "treatment")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v)) {
      stop("`", nm, "` must be a single non-empty character column name.")
    }
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  effects <- as.integer(effects %||% 1L)
  bootstrap_reps <- as.integer(bootstrap_reps)
  seed <- as.integer(seed)
  n_workers <- as.integer(n_workers)

  # All three methods (fe, ife, mc) are now implemented.

  # Resolve backend: cuda > r. CUDA fect kernel is in src/cuda_fect_fe.cu;
  # if not built, fall back to the r-side reference implementation.
  cuda_ok <- isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE))
  resolved_backend <- if (backend %in% c("auto", "cuda") && cuda_ok) "cuda"
                       else "r"

  # Canonical args bundle. `placebo = 0L` is fixed for fect: placebos
  # aren't part of the fect_fe estimator (no pre-period DID concept).
  # Kept in the bundle so didgpu_init_checkpoint can use the same meta
  # schema as didgpu().
  args <- list(
    outcome = outcome, group = group, time = time, treatment = treatment,
    method = method, effects = effects, placebo = 0L,
    r = r, lambda = lambda,
    tol = tol, max_iter = max_iter,
    bootstrap_reps = bootstrap_reps, seed = seed,
    backend = resolved_backend
  )

  # Per-iter fit closure: pick the right implementation.
  fit_one <- switch(method,
    fe  = .fect_fe_one_iter,
    ife = .fect_ife_one_iter,
    mc  = .fect_mc_one_iter,
    stop("unknown method")
  )

  ph <- .panel_hash(df, outcome, group, time, treatment)

  # ---- checkpoint init / load ----
  if (!is.null(checkpoint_dir)) {
    manifest_path <- file.path(checkpoint_dir, "manifest.csv")
    if (file.exists(manifest_path)) {
      chk <- didgpu_load_checkpoint(checkpoint_dir)
      manifest <- chk$manifest
      if (verbose) {
        message(sprintf("[didgpu_fect] resuming %s: %d/%d cells already done",
                        checkpoint_dir, nrow(manifest), bootstrap_reps + 1L))
      }
    } else {
      pkgv <- tryCatch(as.character(utils::packageVersion("didgpu")),
                       error = function(e) "0.0.0.dev")
      meta <- c(args, list(panel_hash = ph, package_version = pkgv))
      didgpu_init_checkpoint(checkpoint_dir, meta)
      manifest <- .empty_manifest()
    }
  } else {
    manifest <- .empty_manifest()
  }

  # ---- iter plan ----
  todo <- .cells_todo(manifest, bootstrap_reps)
  n_total <- bootstrap_reps + 1L
  n_done0 <- n_total - length(todo)

  in_memory <- list()
  for (idx in seq_along(todo)) {
    iter <- todo[idx]
    iter_seed <- if (iter == 0L) 0L else (seed + iter)
    t0 <- Sys.time()
    value <- fit_one(df, args, iter_seed)
    value$wall_seconds_total <- as.numeric(difftime(Sys.time(), t0,
                                                     units = "secs"))
    if (!is.null(checkpoint_dir)) {
      .save_cell(checkpoint_dir, b = iter, value = value,
                 wall_seconds = value$wall_seconds_total)
    } else {
      in_memory[[as.character(iter)]] <- value
    }
    if (verbose) {
      message(sprintf("[didgpu_fect] cell b=%-5d  %.2fs   (%d/%d total)",
                      iter, value$wall_seconds_total,
                      n_done0 + idx, n_total))
    }
  }

  # ---- aggregate (re-use didgpu's aggregator) ----
  cells <- if (!is.null(checkpoint_dir)) {
    didgpu_aggregate_cells(checkpoint_dir)
  } else {
    in_memory
  }
  result <- .aggregate_to_result(cells, args, ph)
  result$checkpoint_dir <- if (!is.null(checkpoint_dir)) {
    normalizePath(checkpoint_dir, winslash = "/", mustWork = FALSE)
  } else NA_character_
  result$method <- method
  class(result) <- c("didgpu_fect_result", "didgpu_result", "list")
  result
}


# .fect_ife_one_iter and .fect_mc_one_iter live in R/fect_ife.R and R/fect_mc.R.


# ----------------------------------------------------------------------------
# Per-method stubs. Each raises a clear NotImplemented with a roadmap
# pointer. Once implementation lands, these become the per-cell fit
# functions (one bootstrap rep) and the orchestrator (didgpu_fect)
# wires them through the existing checkpoint / parallel infrastructure.
# ----------------------------------------------------------------------------

# Retained as a sentinel for any future method names that get added to
# the public dispatcher before their implementation lands. All three
# current methods (fe, ife, mc) ARE implemented.
.not_implemented_fect <- function(method) {
  stop(sprintf(
    "didgpu_fect(method = '%s') is not yet implemented.", method),
    call. = FALSE)
}
