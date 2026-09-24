# ============================================================================
# Backend dispatch
#
# didgpu supports multiple compute backends with the same numerical contract:
#
#   "reference" : delegate to DIDmultiplegtDYN::did_multiplegt_main. Used as
#                 the parity oracle and as the working backend while our
#                 native implementations are under development. Always
#                 available if DIDmultiplegtDYN is installed.
#
#   "r"         : pure-R reimplementation of the binary non-absorbing core
#                 (see vignettes/reference_internals.md, section 3). Has no
#                 external dependency on the reference. Slower than the
#                 reference where the reference is itself in pure R, faster
#                 where the reference does redundant work in the bootstrap.
#                 STATUS: skeleton only; raises NotImplemented.
#
#   "cpu"       : C++ port via Rcpp + Eigen. Drop-in for "r" with the same
#                 numerical contract but ~10-30x faster on real panels.
#                 STATUS: stub.
#
#   "cuda"      : CUDA port via nvcc + cuBLAS + cuSOLVER for the bootstrap
#                 loop. Same numerical contract within tolerance documented
#                 in tests/testthat/test-backends-agree.R.
#                 STATUS: stub.
#
#   "auto"      : pick the best available at runtime, preferring cuda > cpu
#                 > r > reference.
#
# Every backend exposes the same per-iter fit function:
#
#   fit_one(df, args, iter_seed) -> list(
#     effects = numeric(n_effects),   # ATT_e for e = 1..n_effects
#     ate     = numeric(1),
#     placebos = numeric(n_placebos), # NULL if no placebos
#     N_inc   = numeric(n_effects),   # incidence counts per event time
#     N_inc_pl = numeric(n_placebos), # placebo incidence counts
#     direction = character(),        # "in", "out", or "both"
#     wall_seconds = numeric(1)
#   )
#
# `iter_seed` controls bootstrap resampling: iter_seed == 0 means "compute
# the point estimate on the original panel (no resampling)". iter_seed > 0
# means "set the seed, then cluster-resample with replacement".
# ============================================================================


#' Report installed backends and their availability
#'
#' @return A data.frame with columns `backend`, `available`, `notes`.
#' @examples
#' didgpu_backend_info()
#' @export
didgpu_backend_info <- function() {
  ref_ok <- requireNamespace("DIDmultiplegtDYN", quietly = TRUE)
  cuda_ok <- isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE))
  data.frame(
    backend   = c("reference", "r", "cpu", "cuda"),
    available = c(ref_ok, TRUE, TRUE, cuda_ok),
    notes     = c(
      if (ref_ok) paste0("DIDmultiplegtDYN ",
                         as.character(utils::packageVersion("DIDmultiplegtDYN")))
      else        "install.packages('DIDmultiplegtDYN')",
      "binary; full reference parity (controls, weight, trends_*, normalized, predict_het, ...)",
      "binary-no-controls effects in C++; placebos + complex options fall back to r",
      if (cuda_ok) "effects only (binary, no controls); placebos use r-backend"
      else        "stub; install CUDA Toolkit with full headers + nvcc"
    ),
    stringsAsFactors = FALSE
  )
}


#' Resolve a backend name to an actual fit_one closure
#'
#' @param backend Character: "auto", "reference", "r", "cpu", "cuda".
#' @return A function with signature `function(df, args, iter_seed)`.
#' @keywords internal
#' @noRd
.resolve_backend <- function(backend) {
  backend <- match.arg(backend, c("auto", "reference", "r", "cpu", "cuda"))

  if (backend == "auto") {
    info <- didgpu_backend_info()
    available <- info$backend[info$available]
    pref <- c("cuda", "cpu", "r", "reference")
    chosen <- pref[pref %in% available]
    if (length(chosen) == 0L) {
      stop("No backend available. Install DIDmultiplegtDYN for the reference ",
           "backend, or build didgpu against Rcpp/CUDA for native backends.")
    }
    backend <- chosen[1L]
  }

  switch(
    backend,
    reference = .backend_reference(),
    r         = .backend_r(),
    cpu       = .backend_cpu(),
    cuda      = .backend_cuda(),
    stop("unknown backend: ", backend)
  )
}


# -------- backend: reference (delegates to DIDmultiplegtDYN) --------

.backend_reference <- function() {
  if (!requireNamespace("DIDmultiplegtDYN", quietly = TRUE)) {
    stop("'reference' backend requires the DIDmultiplegtDYN package. ",
         "Install with install.packages('DIDmultiplegtDYN').")
  }

  function(df, args, iter_seed) {
    # iter_seed 0 = point estimate on original panel.
    # iter_seed > 0 = cluster-resample with that seed.
    df_use <- if (iter_seed == 0L) df else .cluster_resample(df, args, iter_seed)

    t0 <- Sys.time()
    suppressMessages(suppressWarnings({
      # Call the orchestrator. We force graph_off=TRUE and bootstrap=NULL
      # because we are providing our own bootstrap loop on the outside.
      # NOTE: DIDmultiplegtDYN's arg validator uses inherits(x, "numeric")
      # which is FALSE for integer storage. We coerce numerics to double
      # to keep the reference's brittle check happy.
      res <- DIDmultiplegtDYN::did_multiplegt_dyn(
        df         = as.data.frame(df_use),
        outcome    = args$outcome,
        group      = args$group,
        time       = args$time,
        treatment  = args$treatment,
        effects    = as.double(args$effects),
        placebo    = as.double(args$placebo),
        cluster    = args$cluster,
        controls   = args$controls,
        weight     = args$weight,
        continuous = args$continuous,
        trends_nonparam = args$trends_nonparam,
        trends_lin = isTRUE(args$trends_lin),
        only_never_switchers = isTRUE(args$only_never_switchers),
        same_switchers = isTRUE(args$same_switchers),
        same_switchers_pl = isTRUE(args$same_switchers_pl),
        dont_drop_larger_lower = isTRUE(args$dont_drop_larger_lower),
        switchers  = args$switchers %||% "",
        normalized = isTRUE(args$normalized),
        predict_het = args$predict_het,
        ci_level   = as.double(args$ci_level %||% 95),
        graph_off  = TRUE
      )
    }))
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    # Pull the per-event-time point estimates. With bootstrap=NULL the
    # reference still computes SEs analytically; we only consume the
    # point estimates here (col 1 of Effects/Placebos) because the
    # outer bootstrap loop is providing iteration replicates from which
    # we will compute SEs ourselves at aggregate time.
    effects <- as.numeric(res$results$Effects[, 1])
    ate     <- if (!is.null(res$results$ATE)) as.numeric(res$results$ATE[1, 1])
               else NA_real_
    placebos <- if (!is.null(res$results$Placebos) &&
                    nrow(res$results$Placebos) > 0L) {
      as.numeric(res$results$Placebos[, 1])
    } else numeric(0)

    # Sample sizes from reference's results matrix. The reference reports
    # FOUR count columns: 5 = N (unweighted obs), 6 = Switchers (unweighted
    # switchers), 7 = N.w (weighted obs), 8 = Switchers.w (weighted
    # switchers). Forward all four so weighted reference runs stay faithful;
    # on unweighted panels cols 7/8 coincide with 5/6.
    col_or_na <- function(m, j, n) {
      if (!is.null(m) && ncol(m) >= j) as.numeric(m[, j]) else rep(NA_real_, n)
    }
    n_eff_effects     <- col_or_na(res$results$Effects, 5L, length(effects))
    n_sw_unw_effects  <- col_or_na(res$results$Effects, 6L, length(effects))
    n_eff_w_effects   <- col_or_na(res$results$Effects, 7L, length(effects))
    n_sw_w_effects    <- col_or_na(res$results$Effects, 8L, length(effects))
    # n_inc_effects keeps the reference's Switchers column (col 6) -- it is
    # only used downstream for ATE weighting / fallback, not the estimate.
    n_inc_effects  <- n_sw_unw_effects
    pl_mat <- if (length(placebos) > 0L) res$results$Placebos else NULL
    n_eff_placebos    <- col_or_na(pl_mat, 5L, length(placebos))
    n_sw_unw_placebos <- col_or_na(pl_mat, 6L, length(placebos))
    n_eff_w_placebos  <- col_or_na(pl_mat, 7L, length(placebos))
    n_sw_w_placebos   <- col_or_na(pl_mat, 8L, length(placebos))
    n_inc_placebos <- n_sw_unw_placebos

    list(
      effects        = effects,
      ate            = ate,
      placebos       = placebos,
      n_effects      = length(effects),
      n_placebos     = length(placebos),
      n_inc_effects  = n_inc_effects,
      n_inc_placebos = n_inc_placebos,
      n_eff_effects  = n_eff_effects,
      n_eff_placebos = n_eff_placebos,
      n_eff_w_effects   = n_eff_w_effects,
      n_eff_w_placebos  = n_eff_w_placebos,
      n_sw_unw_effects  = n_sw_unw_effects,
      n_sw_unw_placebos = n_sw_unw_placebos,
      n_sw_w_effects    = n_sw_w_effects,
      n_sw_w_placebos   = n_sw_w_placebos,
      iter_seed      = as.integer(iter_seed),
      wall_seconds   = wall,
      backend        = "reference"
    )
  }
}


# -------- backend: pure R (partial) --------

.backend_r <- function() {
  .backend_r_impl()
}


# -------- backend: CPU C++ (Rcpp) --------

# The CPU backend uses `didgpu_cpp_core_one_event_time` (in
# src/cpu_core.cpp) for the per-(event-time, direction) inner kernel.
# It currently supports the binary-no-controls case. Any unsupported
# feature combination (controls, normalized, trends_lin, continuous,
# same_switchers, ...) falls back to the r-backend transparently.
.backend_cpu <- function() {
  function(df, args, iter_seed) {
    # Feature compatibility check. Fall back to r-backend if any
    # unsupported feature is requested.
    unsupported <- list(
      controls               = !is.null(args$controls),
      weight                 = !is.null(args$weight),
      continuous             = !is.null(args$continuous),
      trends_nonparam        = !is.null(args$trends_nonparam),
      trends_lin             = isTRUE(args$trends_lin),
      normalized             = isTRUE(args$normalized),
      same_switchers         = isTRUE(args$same_switchers),
      same_switchers_pl      = isTRUE(args$same_switchers_pl),
      only_never_switchers   = isTRUE(args$only_never_switchers),
      predict_het            = !is.null(args$predict_het)
    )
    if (any(vapply(unsupported, isTRUE, logical(1)))) {
      return(.backend_r_impl()(df, args, iter_seed))
    }
    df_use <- if (iter_seed == 0L) df else .cluster_resample(df, args, iter_seed)

    t0 <- Sys.time()
    prepped <- .prep_panel(df_use, args$outcome, args$group, args$time,
                            args$treatment,
                            dont_drop_larger_lower = isTRUE(args$dont_drop_larger_lower))
    sw <- args$switchers %||% ""
    h <- .clamp_horizons(prepped, args$effects, args$placebo, switchers = sw)
    ce <- .compute_effects_cpp(prepped, h$l_eff, switchers = sw,
                                want_se = (iter_seed == 0L),
                                cluster_col = args$cluster)
    cp <- .compute_placebos_cpp(prepped, h$l_pl, switchers = sw,
                            want_se = (iter_seed == 0L),
                            cluster_col = args$cluster)
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    # Av_tot_eff -- see .ate_weighted in core_r.R. Neither kernel backend
    # accepts `normalized`, so ce$effects is already the raw DID here.
    ate <- .ate_weighted(ce$effects, ce$n_inc, ce$delta_ate)
    se_ate <- if (is.null(ce$u_mat)) NA_real_ else {
      ok <- !is.na(ce$effects) & !is.na(ce$delta_ate) &
            is.finite(ce$n_inc) & ce$n_inc > 0
      den <- if (any(ok)) sum(ce$n_inc[ok] * ce$delta_ate[ok]) else 0
      if (!any(ok) || !is.finite(den) || den == 0) NA_real_ else {
        .se_from_u(as.numeric(ce$u_mat[, ok, drop = FALSE] %*% ce$n_inc[ok]) / den,
                   ce$G, ce$cluster_of_group)
      }
    }

    list(
      effects        = ce$effects,
      ate            = ate,
      placebos       = cp$placebos,
      n_effects      = h$l_eff,
      n_placebos     = h$l_pl,
      se_effects     = ce$se,
      se_placebos    = cp$se,
      se_ate         = se_ate,
      u_mat_effects  = ce$u_mat,
      u_mat_placebos = cp$u_mat,
      se_G           = ce$G,
      se_cluster_of_group = ce$cluster_of_group,
      n_inc_effects  = ce$n_inc,
      n_inc_placebos = cp$n_inc,
      n_eff_effects  = ce$n_eff,
      n_eff_placebos = cp$n_eff,
      iter_seed      = as.integer(iter_seed),
      wall_seconds   = wall,
      backend        = "cpu"
    )
  }
}


# Wire the C++ kernel through .compute_effects / .compute_placebos shape.
# We borrow .core_one_event_time's per-direction Neyman-pooling logic
# from the R code; only the per-(k, direction) compute is in C++.
#' @keywords internal
#' @noRd
.compute_effects_cpp <- function(prepped, effects, switchers = "",
                                  want_se = FALSE, cluster_col = NULL) {
  out <- numeric(effects); n_inc <- integer(effects); n_eff <- integer(effects)
  delta_ate <- numeric(effects)
  G_all <- length(unique(prepped$group_XX))
  cog <- if (!is.null(cluster_col) && nzchar(cluster_col) &&
               cluster_col %in% names(prepped)) {
    prepped[, list(cl = .SD[[1L]][1L]), by = group_XX,
            .SDcols = cluster_col]$cl
  } else NULL
  u_mat <- if (isTRUE(want_se)) matrix(0, nrow = G_all, ncol = effects) else NULL
  se_vec <- rep(NA_real_, effects)
  for (k in seq_len(effects)) {
    res_in  <- if (switchers != "out") .cpu_one_event_time(prepped, k = k, direction = 1L)
               else list(att = NA_real_, N_inc = 0L, N_eff = 0L)
    res_out <- if (switchers != "in")  .cpu_one_event_time(prepped, k = k, direction = 0L)
               else list(att = NA_real_, N_inc = 0L, N_eff = 0L)
    n_in <- res_in$N_inc; n_out <- res_out$N_inc
    if (n_in + n_out == 0L) {
      out[k] <- NA_real_; n_inc[k] <- 0L; n_eff[k] <- 0L
      delta_ate[k] <- NA_real_; next
    }
    att_in       <- if (n_in  > 0L)  res_in$att      else 0
    att_out_pool <- if (n_out > 0L) -res_out$att     else 0
    w_in <- n_in / (n_in + n_out)
    out[k]   <- w_in * att_in + (1 - w_in) * att_out_pool
    n_inc[k] <- n_in + n_out
    n_eff[k] <- (res_in$N_eff %||% 0L) + (res_out$N_eff %||% 0L)
    # The C++ kernel returns only att and N_inc, so the switcher mask is
    # rebuilt here to get Av_tot_eff's denominator and, when asked, the
    # analytic-SE influence contribution -- one mask build for both.
    ex_in  <- if (n_in  > 0L) .delta_ate_kernel_path(prepped, k, 1L, want_se,
                                                      G_all, cluster_col) else NULL
    ex_out <- if (n_out > 0L) .delta_ate_kernel_path(prepped, k, 0L, want_se,
                                                      G_all, cluster_col) else NULL
    da_in  <- if (!is.null(ex_in))  ex_in$delta_ate  else NA_real_
    da_out <- if (!is.null(ex_out)) ex_out$delta_ate else NA_real_
    delta_ate[k] <- if (n_in == 0L) da_out
                    else if (n_out == 0L) da_in
                    else w_in * da_in + (1 - w_in) * da_out
    if (isTRUE(want_se)) {
      ucol <- numeric(G_all)
      if (!is.null(ex_in)  && !is.null(ex_in$u_var))  ucol <- ucol + w_in * ex_in$u_var
      if (!is.null(ex_out) && !is.null(ex_out$u_var)) ucol <- ucol - (1 - w_in) * ex_out$u_var
      u_mat[, k] <- ucol
      se_vec[k] <- .se_from_u(ucol, G_all, cog)
    }
  }
  list(effects = out, n_inc = n_inc, n_eff = n_eff, delta_ate = delta_ate,
       se = se_vec, u_mat = u_mat, G = G_all, cluster_of_group = cog)
}

# Placebos still use the R backend (the C++ port is only the effects
# kernel for now; placebo kernel is identical structurally but uses a
# different diff_y formula and we haven't ported it).
#' @keywords internal
#' @noRd
.compute_placebos_cpp <- function(prepped, placebo, switchers = "",
                                   want_se = FALSE, cluster_col = NULL) {
  .compute_placebos(prepped, placebo, switchers = switchers,
                     want_se = want_se, cluster_col = cluster_col)
}

# Build the prepped panel into the contiguous columnar shape the C++
# kernel needs, then call it. Cache the layout on `prepped` so the
# per-k call only sees the kernel cost, not the layout cost.
#' @keywords internal
#' @noRd
.cpu_one_event_time <- function(prepped, k, direction) {
  layout <- attr(prepped, "didgpu_cpu_layout")
  if (is.null(layout)) {
    layout <- .cpu_build_layout(prepped)
    data.table::setattr(prepped, "didgpu_cpu_layout", layout)
  }
  didgpu_cpp_core_one_event_time(
    outcome      = layout$outcome,
    N_gt         = layout$N_gt,
    group_id     = layout$group_id,
    time_id      = layout$time_id,
    cohort_id    = layout$cohort_id,
    F_g          = layout$F_g,
    S_g          = layout$S_g,
    T_g          = layout$T_g,
    L_g          = layout$L_g,
    group_offset = layout$group_offset,
    n_cohorts    = layout$n_cohorts,
    k            = as.integer(k),
    direction    = as.integer(direction))
}

# Build the columnar layout once per prep. The prepped panel is
# guaranteed sorted by (group_XX, time_XX) already, so no setorder.
#' @keywords internal
#' @noRd
.cpu_build_layout <- function(prepped) {
  d <- prepped
  # 0-based int keys (already consecutive 1..n in prep).
  g0 <- as.integer(d$group_XX) - 1L
  t0 <- as.integer(d$time_XX) - 1L
  # Cohort = (time, d_sq) pair. Encode as time * (max_dsq+1) + dsq for
  # a quick integer key — faster than factor(paste(...)).
  d_sq_vals <- d$d_sq_XX
  d_sq_int <- as.integer(d_sq_vals)
  # If d_sq has non-integer values (e.g. NaN from no-baseline groups),
  # fall back to factor-based encoding.
  if (anyNA(d_sq_int) || any(d_sq_int != d_sq_vals, na.rm = TRUE)) {
    c0 <- as.integer(factor(paste(t0, d_sq_vals, sep = "_"))) - 1L
  } else {
    span <- max(d_sq_int, na.rm = TRUE) - min(d_sq_int, na.rm = TRUE) + 1L
    if (!is.finite(span) || span <= 0L) span <- 1L
    shift_int <- if (any(d_sq_int < 0L)) -min(d_sq_int) else 0L
    c0 <- t0 * span + (d_sq_int + shift_int)
    # Re-pack to consecutive 0..n_cohorts-1.
    c0 <- as.integer(factor(c0)) - 1L
  }
  n_cohorts <- max(c0) + 1L

  # Per-group summaries (one row per group, in 0..n_groups-1 order).
  g_first <- !duplicated(g0)
  F_g <- as.integer(d$F_g_XX[g_first]) - 1L
  Sg_raw <- d$S_g_XX[g_first]
  S_g <- ifelse(is.na(Sg_raw), -1L, as.integer(Sg_raw))
  T_g <- as.integer(d$T_g_XX[g_first]) - 1L
  L_g <- as.integer(d$L_g_XX[g_first])
  T_g[is.na(T_g)] <- -1L
  L_g[is.na(L_g)] <- 0L
  n_groups <- length(F_g)

  # group_offset = prefix-sum of per-group row counts. Use tabulate
  # (faster than data.table grouped count for this).
  counts <- tabulate(g0 + 1L, nbins = n_groups)
  group_offset <- c(0L, cumsum(counts))

  list(
    outcome      = as.numeric(d$outcome_XX),
    N_gt         = as.numeric(d$N_gt_XX),
    group_id     = g0,
    time_id      = t0,
    cohort_id    = c0,
    F_g          = F_g,
    S_g          = S_g,
    T_g          = T_g,
    L_g          = L_g,
    group_offset = group_offset,
    n_cohorts    = as.integer(n_cohorts)
  )
}


# -------- backend: CUDA --------

.backend_cuda <- function() {
  if (!isTRUE(didgpu_has_cuda_support())) {
    return(function(df, args, iter_seed) {
      stop("Backend 'cuda' is not available. Install the NVIDIA CUDA Toolkit ",
           "with the full include headers (`include/crt/` must contain ",
           "host_config.h etc.), set CUDA_HOME or CUDA_PATH, and reinstall ",
           "didgpu so nvcc compiles src/cuda_*.cu. See inst/doc/cuda_setup.md.")
    })
  }
  # CUDA is built. Run the binary-no-controls path through .cuda_one_event_time
  # (R/cuda_glue.R), which has the same return shape as .core_one_event_time.
  function(df, args, iter_seed) {
    # Feature compatibility check. The CUDA kernel implements the plain
    # effects path ONLY; every other option must fall back to r-backend.
    #
    # This list used to cover just controls/weight/trends_nonparam, so
    # anything else passed straight through to a kernel that does not
    # implement it and a WRONG NUMBER came back silently. `normalized`
    # was the worst case: with a multivalued treatment, backend "cuda"
    # returned the UNnormalised effects, which do not vary with dose --
    # so a dose-response analysis looked as though the treatment had
    # been binarised. Measured against DIDmultiplegtDYN with
    # normalized = TRUE on a 3-level dose:
    #     r     |diff| 5.551e-17
    #     cpu   |diff| 5.551e-17
    #     cuda  |diff| 8.465e-01   <- silently unnormalised
    # backend = "auto" resolves to cuda whenever a GPU is present, so
    # this was the default path on CUDA machines.
    #
    # Kept deliberately identical to .backend_cpu()'s list: a backend
    # must never answer a question it cannot compute.
    unsupported <- list(
      controls               = !is.null(args$controls),
      weight                 = !is.null(args$weight),
      continuous             = !is.null(args$continuous),
      trends_nonparam        = !is.null(args$trends_nonparam),
      trends_lin             = isTRUE(args$trends_lin),
      normalized             = isTRUE(args$normalized),
      same_switchers         = isTRUE(args$same_switchers),
      same_switchers_pl      = isTRUE(args$same_switchers_pl),
      only_never_switchers   = isTRUE(args$only_never_switchers),
      predict_het            = !is.null(args$predict_het)
    )
    if (any(vapply(unsupported, isTRUE, logical(1)))) {
      if (iter_seed == 0L) {
        hit <- names(unsupported)[vapply(unsupported, isTRUE, logical(1))]
        message("[didgpu] CUDA backend does not implement: ",
                paste(hit, collapse = ", "),
                ". Falling back to r-backend for this call.")
      }
      return(.backend_r_impl()(df, args, iter_seed))
    }
    df_use <- if (iter_seed == 0L) df else .cluster_resample(df, args, iter_seed)
    t0 <- Sys.time()
    # dont_drop_larger_lower must reach .prep_panel or it is silently
    # ignored on this backend (it is honoured on r and cpu).
    prepped <- .prep_panel(df_use, args$outcome, args$group, args$time,
                            args$treatment,
                            dont_drop_larger_lower = isTRUE(args$dont_drop_larger_lower))
    sw <- args$switchers %||% ""
    h <- .clamp_horizons(prepped, args$effects, args$placebo, switchers = sw)
    # Run the per-event-time CUDA kernel. We mimic .compute_effects's
    # direction loop here so the result shape stays consistent.
    ce <- .compute_effects_cuda(prepped, h$l_eff, switchers = sw,
                                 want_se = (iter_seed == 0L),
                                 cluster_col = args$cluster)
    cp <- .compute_placebos(prepped, h$l_pl, switchers = sw,
                            want_se = (iter_seed == 0L),
                            cluster_col = args$cluster)  # placebos still r-side
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    # This branch has had the ATE wrong twice. It first read
    #     ate <- if (h$l_eff == 1L) ce$effects[1] else NA_real_
    # so backend "cuda" returned NA whenever effects > 1 while the
    # per-horizon effects themselves were fine; then it returned the
    # switcher-weighted mean, which is Av_tot_eff only for binary
    # absorbing treatment. Both were default-path bugs, because
    # backend = "auto" resolves to "cuda" on any machine with a GPU.
    # There is now one definition, in .ate_weighted() (core_r.R). This
    # backend rejects `normalized`, so ce$effects is already the raw DID.
    ate <- .ate_weighted(ce$effects, ce$n_inc, ce$delta_ate)
    se_ate <- if (is.null(ce$u_mat)) NA_real_ else {
      ok <- !is.na(ce$effects) & !is.na(ce$delta_ate) &
            is.finite(ce$n_inc) & ce$n_inc > 0
      den <- if (any(ok)) sum(ce$n_inc[ok] * ce$delta_ate[ok]) else 0
      if (!any(ok) || !is.finite(den) || den == 0) NA_real_ else {
        .se_from_u(as.numeric(ce$u_mat[, ok, drop = FALSE] %*% ce$n_inc[ok]) / den,
                   ce$G, ce$cluster_of_group)
      }
    }
    list(
      effects        = ce$effects,
      ate            = ate,
      placebos       = cp$placebos,
      n_effects      = h$l_eff,
      n_placebos     = h$l_pl,
      se_effects     = ce$se,
      se_placebos    = cp$se,
      se_ate         = se_ate,
      u_mat_effects  = ce$u_mat,
      u_mat_placebos = cp$u_mat,
      se_G           = ce$G,
      se_cluster_of_group = ce$cluster_of_group,
      n_inc_effects  = ce$n_inc,
      n_inc_placebos = cp$n_inc,
      n_eff_effects  = ce$n_inc,
      n_eff_placebos = cp$n_eff,
      iter_seed      = as.integer(iter_seed),
      wall_seconds   = wall,
      backend        = "cuda"
    )
  }
}

# .compute_effects but using .cuda_one_event_time (R/cuda_glue.R) instead
# of .core_one_event_time. Same Neyman pooling logic.
.compute_effects_cuda <- function(prepped, effects, switchers = "",
                                   want_se = FALSE, cluster_col = NULL) {
  out <- numeric(effects); n_inc <- integer(effects)
  delta_ate <- numeric(effects)
  G_all <- length(unique(prepped$group_XX))
  cog <- if (!is.null(cluster_col) && nzchar(cluster_col) &&
               cluster_col %in% names(prepped)) {
    prepped[, list(cl = .SD[[1L]][1L]), by = group_XX,
            .SDcols = cluster_col]$cl
  } else NULL
  u_mat <- if (isTRUE(want_se)) matrix(0, nrow = G_all, ncol = effects) else NULL
  se_vec <- rep(NA_real_, effects)
  for (k in seq_len(effects)) {
    res_in <- if (switchers != "out") .cuda_one_event_time(prepped, k = k, direction = 1L, want_se = want_se, cluster_col = cluster_col)
              else list(att = NA_real_, N_inc = 0L, delta_ate = NA_real_)
    res_out <- if (switchers != "in") .cuda_one_event_time(prepped, k = k, direction = 0L, want_se = want_se, cluster_col = cluster_col)
               else list(att = NA_real_, N_inc = 0L, delta_ate = NA_real_)
    n_in <- res_in$N_inc; n_out <- res_out$N_inc
    if (n_in + n_out == 0L) {
      out[k] <- NA_real_; n_inc[k] <- 0L; delta_ate[k] <- NA_real_; next
    }
    att_in       <- if (n_in  > 0L)  res_in$att  else 0
    att_out_pool <- if (n_out > 0L) -res_out$att else 0
    w_in <- n_in / (n_in + n_out)
    out[k] <- w_in * att_in + (1 - w_in) * att_out_pool
    n_inc[k] <- n_in + n_out
    # Av_tot_eff's denominator, pooled with the same Neyman weights.
    da_in  <- if (n_in  > 0L) res_in$delta_ate  else NA_real_
    da_out <- if (n_out > 0L) res_out$delta_ate else NA_real_
    delta_ate[k] <- if (n_in == 0L) da_out
                    else if (n_out == 0L) da_in
                    else w_in * da_in + (1 - w_in) * da_out
    if (isTRUE(want_se)) {
      ucol <- numeric(G_all)
      if (n_in  > 0L && !is.null(res_in$u_var))  ucol <- ucol + w_in * res_in$u_var
      if (n_out > 0L && !is.null(res_out$u_var)) ucol <- ucol - (1 - w_in) * res_out$u_var
      u_mat[, k] <- ucol
      se_vec[k] <- .se_from_u(ucol, G_all, cog)
    }
  }
  list(effects = out, n_inc = n_inc, delta_ate = delta_ate,
       se = se_vec, u_mat = u_mat, G = G_all, cluster_of_group = cog)
}


# -------- shared: cluster-resample one bootstrap iter --------

#' @keywords internal
#' @noRd
.cluster_resample <- function(df, args, iter_seed) {
  # Force Mersenne-Twister so the sample is identical regardless of
  # whether we're in a parallel worker (which defaults to L'Ecuyer-CMRG)
  # or the main R session.
  set.seed(as.integer(iter_seed), kind = "Mersenne-Twister")
  cluster_col <- args$cluster %||% args$group
  ids <- df[[cluster_col]]
  # Build index list per cluster level, then sample clusters with replacement.
  by_clust <- split(seq_len(nrow(df)), ids)
  picks <- sample(length(by_clust), length(by_clust), replace = TRUE)

  # When the same cluster is picked more than once, each pick must be
  # treated as a distinct unit for the bootstrap. Otherwise the
  # downstream code sees duplicate (group, time) keys and either
  # cartesian-explodes (didgpu's r-backend) or silently collapses them
  # (changing semantics).
  #
  # Relabel strategy: assign each pick a *fresh integer ID* in a space
  # that doesn't collide with the originals. This keeps both group and
  # cluster columns numeric so the reference backend (which coerces
  # cluster to numeric internally) keeps working. We use:
  #   new_id = (pick_count - 1) * OFFSET + original_id
  # where OFFSET is large enough to never collide with originals.
  # OFFSET must be larger than the max ID in BOTH the cluster column
  # and the group column, so relabeled IDs don't collide with existing
  # original IDs in either column.
  orig_ids <- if (is.numeric(ids)) ids else suppressWarnings(as.numeric(ids))
  if (any(is.na(orig_ids))) orig_ids <- as.integer(factor(ids))
  orig_groups <- if (is.numeric(df[[args$group]])) df[[args$group]]
                 else suppressWarnings(as.numeric(df[[args$group]]))
  if (any(is.na(orig_groups))) orig_groups <- as.integer(factor(df[[args$group]]))
  OFFSET <- max(orig_ids, orig_groups, na.rm = TRUE) + 1L

  pick_counts <- integer(length(by_clust))
  out_blocks <- vector("list", length(picks))
  for (i in seq_along(picks)) {
    j <- picks[i]
    pick_counts[j] <- pick_counts[j] + 1L
    block <- df[by_clust[[j]], , drop = FALSE]
    if (pick_counts[j] > 1L) {
      shift <- (pick_counts[j] - 1L) * OFFSET
      block[[args$group]] <- as.numeric(block[[args$group]]) + shift
      if (!is.null(args$cluster) && args$cluster != args$group) {
        block[[args$cluster]] <- as.numeric(block[[args$cluster]]) + shift
      }
    }
    out_blocks[[i]] <- block
  }
  out <- do.call(rbind, out_blocks)
  # Stable order for downstream determinism.
  out <- out[order(out[[args$group]], out[[args$time]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a
