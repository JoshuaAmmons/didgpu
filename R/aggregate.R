# ============================================================================
# Aggregate committed cells into a result object.
#
# Cell 0 is the point estimate (no resampling). Cells 1..bootstrap_reps are
# the bootstrap iters. The aggregator:
#
#   - Pulls the per-iter effects/placebos vectors into matrices.
#   - SE per coefficient = sd across bootstrap iters.
#   - CI per coefficient = point +/- z * SE  (normal approx).
#   - Joint chi-square p-values for the effects vector and the placebos
#     vector (using the bootstrap empirical covariance).
#
# Output structure mirrors DIDmultiplegtDYN where it can:
#   $coef     : list(b = named numeric, vcov = matrix)
#   $results  : list(N_Effects, N_Placebos, Effects, ATE, Placebos,
#                    p_jointeffects, p_jointplacebo)
#   $args     : the canonical args bundle
#   $manifest : the cell manifest data.frame
#
# Differences from the reference's output:
#   - Standard errors are bootstrap-derived (not the reference's
#     analytical SE). With enough bootstrap reps this converges, but
#     for small bootstrap_reps the SEs WILL differ from the reference.
#     That is intentional: we are checkpointing the bootstrap, which
#     means SEs come out of that bootstrap rather than out of analytic
#     formulas. The Effects[, 1] column (the point estimate) IS expected
#     to match the reference to floating-point precision under the
#     'reference' backend.
#   - Coefficient names are NOT padded with trailing spaces (the
#     reference's "Effect_1    " naming is preserved as a tolerance in
#     compat.R for downstream code that depends on it; the canonical
#     name in didgpu is "Effect_1" with no padding).
# ============================================================================


#' @keywords internal
#' @noRd
.aggregate_to_result <- function(cells, args, panel_hash) {
  if (length(cells) == 0L) {
    stop("No cells available to aggregate. ",
         "Either bootstrap_reps == 0 with no point estimate yet, or the ",
         "checkpoint manifest is empty.")
  }

  # Pull cells in iter order. Convert names from "0","1",... to integers
  # to be sure the point estimate ends up at the right position.
  iters <- sort(as.integer(names(cells)))
  cells <- cells[as.character(iters)]
  has_point <- 0L %in% iters
  if (!has_point) {
    stop("Cannot aggregate: cell b=0 (point estimate) is missing from the ",
         "checkpoint manifest.")
  }
  boot_iters <- setdiff(iters, 0L)

  # Effects matrix : boot_iters x n_effects (point goes in $coef$b).
  e0 <- cells[["0"]]$effects
  n_e <- length(e0)
  p0 <- cells[["0"]]$placebos
  n_p <- length(p0)

  # Drop degenerate bootstrap iterations instead of crashing. Under
  # sparse-switching treatments a resample can contain no valid switcher
  # cell at some horizon, yielding an effects/placebos vector of the wrong
  # length (typically length 0). The vapply() calls below hard-require
  # exact lengths, so a single such iteration used to abort the whole
  # aggregation ("values must be length K ... result is length 0") -- and
  # the failure probability grows with bootstrap_reps, so exactly the
  # large-rep runs users want for final inference were the ones crashing.
  # Dropping failed resamples is standard bootstrap practice; SEs and the
  # bootstrap covariance are computed from the surviving iterations, and a
  # warning reports how many were dropped.
  n_boot_dropped <- 0L
  if (length(boot_iters) > 0L) {
    ok <- vapply(as.character(boot_iters), function(i) {
      ce <- cells[[i]]$effects
      cp <- cells[[i]]$placebos
      length(ce) == n_e && (n_p == 0L || length(cp) == n_p) && !all(is.na(ce))
    }, logical(1))
    if (any(!ok)) {
      n_boot_dropped <- sum(!ok)
      warning(sprintf(
        paste0("didgpu: dropped %d of %d bootstrap iteration(s) whose resample ",
               "produced degenerate cells (no valid switchers at some horizon); ",
               "SEs/covariance use the remaining %d iterations."),
        n_boot_dropped, length(boot_iters), sum(ok)), call. = FALSE)
      boot_iters <- boot_iters[ok]
    }
  }

  e_mat <- if (length(boot_iters) > 0L) {
    m <- vapply(as.character(boot_iters),
                function(i) cells[[i]]$effects,
                numeric(n_e))
    if (is.null(dim(m))) matrix(m, nrow = n_e, ncol = length(boot_iters))
    else m
  } else matrix(numeric(0), nrow = n_e, ncol = 0L)
  e_mat <- t(e_mat)   # iter x effect

  p_mat <- if (length(boot_iters) > 0L && n_p > 0L) {
    m <- vapply(as.character(boot_iters),
                function(i) cells[[i]]$placebos,
                numeric(n_p))
    if (is.null(dim(m))) matrix(m, nrow = n_p, ncol = length(boot_iters))
    else m
  } else matrix(numeric(0), nrow = n_p, ncol = 0L)
  p_mat <- t(p_mat)

  # The horizon counts are the reference's own l_XX / l_placebo_XX (see
  # .clamp_horizons), so every row is reported, as the reference reports
  # it: a horizon nobody reaches is an NA row with a message, not a
  # missing row.

  ate_vec <- if (length(boot_iters) > 0L) {
    vapply(as.character(boot_iters),
           function(i) cells[[i]]$ate %||% NA_real_,
           numeric(1))
  } else numeric(0)

  # SEs and CIs. A kept bootstrap iteration can still carry NA at SOME
  # horizons (its resample has switchers overall but none reaching horizon
  # j). Plain sd()/cov() would then return NA for every affected column --
  # with few switchers this wiped out ALL SEs. Compute per-horizon SEs
  # from the finite draws, requiring a minimum bootstrap support of
  # MIN_BOOT_SUPPORT finite draws per horizon (below that the SE is
  # genuinely not estimable and stays NA).
  # The support floor scales with the requested rep count: users running
  # quick small-rep fits (e.g. 16 reps, all finite) still get SEs; large
  # runs require at least 30 finite draws per horizon.
  .col_sd <- function(m) {
    if (nrow(m) < 2L) return(rep(NA_real_, ncol(m)))
    floor_n <- max(2L, min(30L, as.integer(nrow(m) %/% 2L)))
    apply(m, 2L, function(col) {
      v <- col[is.finite(col)]
      if (length(v) >= floor_n) stats::sd(v) else NA_real_
    })
  }
  z <- stats::qnorm(0.5 + (args$ci_level %||% 95) / 200)
  # Analytic SEs, from the estimator's asymptotic linear representation,
  # are what DIDmultiplegtDYN reports and are preferred whenever the
  # point-estimate cell carries them. The bootstrap SD is the fallback
  # for backends and option combinations that do not produce them.
  # See R/analytic_se.R.
  cell0 <- cells[["0"]]
  se_an_e <- cell0$se_effects  %||% NULL
  se_an_p <- cell0$se_placebos %||% NULL
  e_se <- if (nrow(e_mat) >= 2L) .col_sd(e_mat) else rep(NA_real_, n_e)
  p_se <- if (nrow(p_mat) >= 2L) .col_sd(p_mat) else rep(NA_real_, n_p)
  if (!is.null(se_an_e) && length(se_an_e) == n_e) {
    e_se <- ifelse(is.finite(se_an_e), se_an_e, e_se)
  }
  if (!is.null(se_an_p) && length(se_an_p) == n_p) {
    p_se <- ifelse(is.finite(se_an_p), se_an_p, p_se)
  }
  ate_ok <- ate_vec[is.finite(ate_vec)]
  ate_floor <- max(2L, min(30L, as.integer(length(ate_vec) %/% 2L)))
  ate_se <- if (length(ate_ok) >= ate_floor) stats::sd(ate_ok) else NA_real_
  se_an_ate <- cell0$se_ate %||% NA_real_
  if (is.finite(se_an_ate)) ate_se <- se_an_ate

  e_ci_lo <- e0 - z * e_se;  e_ci_hi <- e0 + z * e_se
  p_ci_lo <- p0 - z * p_se;  p_ci_hi <- p0 + z * p_se
  ate0 <- cells[["0"]]$ate %||% NA_real_
  ate_ci_lo <- ate0 - z * ate_se
  ate_ci_hi <- ate0 + z * ate_se

  # Row names as the reference writes them: padded to 12 characters
  # ("Effect_1    ", "Av_tot_eff  "), so code written against
  # DIDmultiplegtDYN's matrices indexes didgpu's the same way.
  # NB: paste0("Effect_", seq_len(0)) is the length-1 "Effect_", hence
  # the guards.
  .pad12 <- function(x) sprintf("%-12s", x)
  effect_names  <- if (n_e > 0L) paste0("Effect_",  seq_len(n_e)) else character(0)
  placebo_names <- if (n_p > 0L) paste0("Placebo_", seq_len(n_p)) else character(0)

  # Count columns come from cell b=0 (the point estimate). The reference
  # reports FOUR distinct count columns:
  #   N           = unweighted contributing observations
  #   Switchers   = unweighted switcher cells
  #   N.w         = weighted   contributing observations (sum of weights)
  #   Switchers.w = weighted   switcher cells (sum of switcher weights)
  # The r/reference backends populate dedicated fields for each; the
  # kernel backends only run on UNWEIGHTED panels where N_gt is 0/1, so
  # the weighted and unweighted columns coincide there.
  c0 <- cells[["0"]]
  n_eff_e    <- (c0$n_eff_effects %||% rep(NA_integer_, n_e))[seq_len(n_e)]              # N
  n_eff_p    <- (c0$n_eff_placebos %||% rep(NA_integer_, n_p))[seq_len(n_p)]
  n_sw_unw_e <- (c0$n_sw_unw_effects %||% c0$n_inc_effects %||% rep(NA_integer_, n_e))[seq_len(n_e)]   # Switchers
  n_sw_unw_p <- (c0$n_sw_unw_placebos %||% c0$n_inc_placebos %||% rep(NA_integer_, n_p))[seq_len(n_p)]
  n_eff_w_e  <- (c0$n_eff_w_effects %||% n_eff_e)[seq_len(n_e)]                          # N.w
  n_eff_w_p  <- (c0$n_eff_w_placebos %||% n_eff_p)[seq_len(n_p)]
  n_sw_w_e   <- (c0$n_sw_w_effects %||% n_sw_unw_e)[seq_len(n_e)]                        # Switchers.w
  n_sw_w_p   <- (c0$n_sw_w_placebos %||% n_sw_unw_p)[seq_len(n_p)]
  # The point estimate's own per-event-time observation counts, which
  # count a row once across switching directions (see .dcdh_ate_extras).
  if (length(c0$n_eff_union) == n_e && n_e > 0L) {
    n_eff_e   <- c0$n_eff_union
    n_eff_w_e <- c0$n_eff_w_union
  }

  cols8 <- c("Estimate", "SE", "LB CI", "UB CI",
             "N", "Switchers", "N.w", "Switchers.w")
  .mat8 <- function(est, se, lo, hi, n, sw, nw, sww, rn) {
    m <- matrix(c(est, se, lo, hi, n, sw, nw, sww), ncol = 8L,
                nrow = length(est), dimnames = list(.pad12(rn), cols8))
    if (!length(est)) rownames(m) <- NULL
    m
  }
  Effects  <- .mat8(e0, e_se, e_ci_lo, e_ci_hi,
                    n_eff_e, n_sw_unw_e, n_eff_w_e, n_sw_w_e, effect_names)
  Placebos <- .mat8(p0, p_se, p_ci_lo, p_ci_hi,
                    n_eff_p, n_sw_unw_p, n_eff_w_p, n_sw_w_p, placebo_names)

  # The ATE row (did_multiplegt_main.R:1494-1531): Switchers columns sum
  # the per-event-time switchers, N columns count the observations used
  # by any event-time. Under trends_lin the reference leaves it all NA.
  tl <- isTRUE(args$trends_lin)
  ate_sw   <- if (tl || !n_e) NA_real_ else sum(n_sw_unw_e, na.rm = TRUE)
  ate_sw_w <- if (tl || !n_e) NA_real_ else sum(n_sw_w_e,   na.rm = TRUE)
  ATE <- .mat8(ate0, ate_se, ate_ci_lo, ate_ci_hi,
               c0$ate_N %||% NA_real_, ate_sw,
               c0$ate_N_w %||% NA_real_, ate_sw_w, "Av_tot_eff")

  # Coefficient vector and its covariance, effects then placebos.
  #
  # When the point-estimate cell carries the influence vectors, the
  # covariance is ANALYTIC: se^2 on the diagonal, polarisation off it,
  # the way DIDmultiplegtDYN builds its joint tests. That one matrix then
  # feeds vcov(), fit$results$p_jointeffects / p_jointplacebo, and
  # didgpu_joint_placebo(), so the three can never disagree with each
  # other or with the reported SE column. The bootstrap covariance is the
  # fallback when no influence vectors exist -- controls, continuous,
  # trends_lin.
  #
  # This is deliberately NOT the reference's coef$vcov. That matrix
  # (did_multiplegt_main.R:2259-2296) wraps each influence column in
  # ifelse(is.null(col), NA, col), which returns only the column's FIRST
  # element, so its off-diagonal terms are built from one group's value
  # recycled over all groups: they change when the groups are relabelled
  # and the matrix is not positive semi-definite. Its joint tests use the
  # full columns and agree with this one.
  b <- c(e0, p0)
  names(b) <- .pad12(c(effect_names, placebo_names))
  vnames <- c(effect_names, placebo_names)
  G_se <- cell0$se_G %||% NA_real_
  cog  <- cell0$se_cluster_of_group
  U_e <- .se_scale_u(cell0$u_mat_effects,  cell0$u_scale_effects)
  U_p <- .se_scale_u(cell0$u_mat_placebos, cell0$u_scale_placebos)
  have_an <- !is.null(U_e) && ncol(U_e) == n_e &&
             (n_p == 0L || (!is.null(U_p) && ncol(U_p) == n_p)) &&
             all(is.finite(c(e_se, p_se)))
  V <- NULL
  if (have_an) {
    U_all <- if (n_p > 0L) cbind(U_e, U_p) else U_e
    V <- .se_vcov_from_u(U_all, c(e_se, p_se), G_se, cog)
    if (!all(is.finite(V))) V <- NULL
  }
  if (is.null(V)) {
    V <- if (nrow(e_mat) >= 2L) {
      full <- if (n_p > 0L) cbind(e_mat, p_mat) else e_mat
      stats::cov(full, use = "pairwise.complete.obs")
    } else matrix(NA_real_, nrow = length(b), ncol = length(b))
  }
  dimnames(V) <- list(vnames, vnames)

  # ---- what the reference says while estimating, in its order ----
  # (horizon counts, unestimable rows, then the two joint tests)
  ref_notes <- c0$ref_notes
  notes <- lapply(c0$horizon_notes %||% character(0),
                  function(s) list(type = "message", text = s))
  .say <- function(type, text) {
    notes[[length(notes) + 1L]] <<- list(type = type, text = text)
  }
  .zero <- function(x) is.na(x) | x == 0
  for (i in seq_len(n_e)) {
    if (.zero(n_sw_w_e[i]) || .zero(n_eff_w_e[i])) {
      .say("message", paste0("Effect_", i, " cannot be estimated. There is ",
                             "no switcher or no control for this effect."))
    }
  }
  for (i in seq_len(n_p)) {
    if (.zero(n_sw_w_p[i]) || .zero(n_eff_w_p[i])) {
      .say("message", paste0("Placebo_", i, " cannot be estimated. There is ",
                             "no switcher or no control for this placebo."))
    }
  }

  # Joint nullity tests (did_multiplegt_main.R:1714-1806, 1810-1904): run
  # only with two or more horizons, all of them estimated; NA when the
  # covariance is not invertible, and a caveat when it is close to it.
  # The reference raises both as warnings inside suppressWarnings(), so
  # they reach the user only through results$vcov_warnings and the
  # "Warnings" block of print(); they are recorded here the same way.
  vcov_warnings <- character(0)
  nrm <- isTRUE(args$normalized)
  .joint <- function(est, Vb, sw, scale, what, boot) {
    l <- length(est)
    ok_n <- sum(!.zero(sw)) == l
    ok_d <- !nrm || (!is.null(scale) && length(scale) == l &&
                       all(is.finite(scale) & scale != 0))
    if (!(ok_n && ok_d)) {
      .say("message", sprintf(paste0("Some %s could not be estimated. ",
        "Therefore, the test of joint nullity of the %s could not be ",
        "computed."), what, what))
      return(NA_real_)
    }
    if (!have_an) return(.joint_pvalue(est, boot))
    ev <- eigen(Vb, only.values = TRUE)$values
    ev <- Re(ev[abs(Im(ev)) < 1e-10])
    ev <- ev[ev > 1e-10]
    one <- if (what == "effects") "effect" else "placebo"
    if (length(ev) < l) {
      w <- sprintf(paste0("The F-test that all %s are equal to zero is not ",
        "computed because the variance of %s is not invertible. This can for ",
        "instance happen if you cluster standard errors and you have more %s ",
        "estimators than clusters."), what, what, one)
      vcov_warnings <<- c(vcov_warnings, w)
      return(NA_real_)
    }
    if (max(ev) / min(ev) >= 1000) {
      w <- sprintf(paste0("The F-test that all %s are equal to zero may not ",
        "be reliable, because the variance of the %s is close to not being ",
        "invertible (the ratio of its largest and smallest eigenvalues is ",
        "larger than 1000). This can for instance happen when you compute ",
        "many %s estimators, or when your %s are very strongly correlated."),
        # The reference's wording: "many effects estimators", but "many
        # placebo estimators".
        what, what, if (what == "effects") "effects" else "placebo", what)
      vcov_warnings <<- c(vcov_warnings, w)
    }
    chi2 <- as.numeric(t(est) %*% MASS::ginv(Vb) %*% est)
    1 - stats::pchisq(chi2, df = l)
  }
  ie <- seq_len(n_e); ip <- n_e + seq_len(n_p)
  p_joint_e <- if (n_e > 1L) {
    .joint(e0, V[ie, ie, drop = FALSE], n_sw_w_e, cell0$u_scale_effects,
           "effects", e_mat)
  } else NULL
  p_joint_p <- if (n_p > 1L) {
    .joint(p0, V[ip, ip, drop = FALSE], n_sw_w_p, cell0$u_scale_placebos,
           "placebos", p_mat)
  } else NULL

  # Test that the effects in [lb, ub] are equal (did_multiplegt_main.R:
  # 2125-2253): a chi-square on their deviations from their mean.
  p_equal <- NULL
  if (isTRUE(args$effects_equal) && n_e > 1L) {
    lb <- args$effects_equal_lb %||% 1L
    ub <- args$effects_equal_ub %||% n_e
    if (ub > n_e) {
      .say("message", sprintf("Upper bound %d exceeds number of effects %d. Using %d as upper bound.",
                              as.integer(ub), as.integer(n_e), as.integer(n_e)))
      ub <- n_e
    }
    whole <- lb == 1L && ub == n_e
    rng <- lb:ub; L <- length(rng)
    p_equal <- NA_real_
    if (sum(!.zero(n_sw_w_e[rng])) == L) {
      Dm <- cbind(diag(L - 1L), 0) - matrix(1 / L, L - 1L, L)
      te <- Dm %*% e0[rng]
      tv <- Dm %*% V[rng, rng, drop = FALSE] %*% t(Dm)
      tv <- (tv + t(tv)) / 2
      ev <- eigen(tv, only.values = TRUE)$values
      ev <- Re(ev[abs(Im(ev)) < 1e-10]); ev <- ev[ev > 1e-10]
      if (!all(is.finite(tv))) {
        p_equal <- NA_real_
      } else if (length(ev) < L - 1L) {
        vcov_warnings <- c(vcov_warnings, if (whole)
          "The F-test that all effects are equal is not computed because the variance of effects is not invertible. This may be due to perfect multicollinearity among the effects. Consider reducing the number of effects estimated."
          else sprintf("The F-test that effects %d to %d are equal is not computed because the variance of effects is not invertible. This may be due to perfect multicollinearity among the effects.", as.integer(lb), as.integer(ub)))
      } else {
        if (max(ev) / min(ev) >= 1000) {
          vcov_warnings <- c(vcov_warnings, if (whole)
            "The F-test that all effects are equal may not be reliable, because the variance of the effects is close to not being invertible (the ratio of its largest and smallest eigenvalues is larger than 1000). This may be due to strong multicollinearity among the effects. Consider reducing the number of effects estimated."
            else sprintf("The F-test that effects %d to %d are equal may not be reliable, because the variance of the effects is close to not being invertible (the ratio of its largest and smallest eigenvalues is larger than 1000).", as.integer(lb), as.integer(ub)))
        }
        chi2 <- as.numeric(t(te) %*% MASS::ginv(tv) %*% te)
        p_equal <- 1 - stats::pchisq(chi2, df = L - 1L)
      }
    } else {
      .say("message", if (whole)
        "Some effects could not be estimated. Therefore, the test of equality of effects could not be computed."
        else sprintf("Some effects in range %d to %d could not be estimated. Therefore, the test of equality of effects could not be computed.", as.integer(lb), as.integer(ub)))
    }
  }

  # backend = "reference" hands back DIDmultiplegtDYN's own tests and its
  # own words; nothing is recomputed or re-said for it.
  if (!is.null(ref_notes)) {
    notes <- ref_notes
    vcov_warnings <- c0$vcov_warnings_ref %||% character(0)
    if (n_e > 1L) p_joint_e <- c0$p_joint_effects_ref %||% NA_real_
    if (n_p > 1L) p_joint_p <- c0$p_joint_placebo_ref %||% NA_real_
    if (!is.null(p_equal)) p_equal <- c0$p_equality_ref %||% NA_real_
  }

  # predict_het: carry the iter-0 cell's block (a data.frame) through
  # into results$predict_het. If absent, omit the field.
  het_block <- cells[["0"]]$predict_het

  # $results in DIDmultiplegtDYN's order and with its presence rules
  # (did_multiplegt_main.R:2316-2365): no Placebos block without
  # placebos, a joint test only where the reference reports one.
  results_list <- list(
    N_Effects         = as.numeric(n_e),
    N_Placebos        = as.numeric(n_p),
    Effects           = Effects,
    ATE               = ATE,
    delta_D_avg_total = c0$delta_D_avg_total %||% NA_real_,
    max_pl            = c0$max_pl %||% NA_real_,
    max_pl_gap        = c0$max_pl_gap %||% NA_real_
  )
  if (!is.null(p_joint_e)) results_list$p_jointeffects <- p_joint_e
  if (!is.null(p_equal)) results_list$p_equality_effects <- p_equal
  if (n_p > 0L) {
    results_list$Placebos <- Placebos
    if ((args$placebo %||% 0L) > 1L && n_p > 1L) {
      results_list$p_jointplacebo <- p_joint_p
    }
  }
  if (!is.null(het_block)) results_list$predict_het <- het_block
  if (length(vcov_warnings)) results_list$vcov_warnings <- vcov_warnings
  results_list$n_boot         <- length(boot_iters)
  results_list$n_boot_dropped <- n_boot_dropped

  out <- list(
    coef = list(b = b, vcov = V),
    results = results_list,
    args = c(args, list(panel_hash = panel_hash)),
    cells_used = length(cells)
  )
  # normalized_weights, formatted as the reference returns it.
  if (!is.null(c0$norm_weights_ref)) out$normalized_weights <- c0$norm_weights_ref
  if (!is.null(c0$norm_weights)) {
    W <- c0$norm_weights
    tot <- matrix(1, 1, ncol(W)) %*% ifelse(is.na(W), 0, W)
    W <- rbind(W, tot)
    dimnames(W) <- list(c(paste0("k=", seq_len(ncol(W)) - 1L), "Total"),
                        paste0("\u2113", "=", seq_len(ncol(W))))
    W[, ] <- sprintf("%s", format(round(W[, ], 3), big.mark = ",",
                                  scientific = FALSE, trim = TRUE))
    out$normalized_weights <- list(norm_weight_mat = noquote(W))
  }
  # DIDmultiplegtDYN 2.4.0 also exposes the per-effect switcher counts
  # at the top level, as Stata's e(N_switchers_effect_k).
  for (k in seq_len(n_e)) {
    out[[paste0("N_switchers_effect_", k)]] <- as.numeric(Effects[k, "Switchers"])
  }
  # avg_time_periods (did_multiplegt_dyn.R:492-508 in 2.4.0).
  if (!is.null(c0$avg_cumul)) {
    av <- c0$avg_cumul
    out$avg_time_periods <- av
    out$avg_cumul <- av$avg_cumul
    for (k in seq_along(av$nswitch)) out[[paste0("N_switch_avg_", k)]] <- av$nswitch[k]
    if (is.null(ref_notes)) {
      notes[[length(notes) + 1L]] <- list(type = "message", text = sprintf(
        "Average number of time periods over which a treatment's effect is accumulated = %s",
        format(av$avg_cumul, nsmall = 4)))
    }
  }
  attr(out, "didgpu_notes") <- notes
  out
}


# Say the notes .aggregate_to_result collected, as DIDmultiplegtDYN says
# them: messages as messages, warnings as warnings.
.emit_notes <- function(notes) {
  for (n in notes %||% list()) {
    if (identical(n$type, "warning")) warning(n$text, call. = FALSE)
    else message(n$text)
  }
  invisible(NULL)
}


# Joint chi-square p-value: theta0' V^-1 theta0 ~ chi2(k) under H0:
# theta = 0. V is the bootstrap covariance.
#
# Two robustness rules (0.1.2):
#  - Horizons with fewer than 30 finite bootstrap draws are excluded from
#    the joint test (their variance is not estimable), and NA draws at
#    kept horizons are handled with a pairwise-complete covariance.
#  - When the covariance of a high-dimensional coefficient block is
#    near-singular (rcond < 1e-10, common with many horizons and few
#    switchers), the chi-square statistic is numerically unstable: a
#    warning tells the user to prefer a low-dimensional prespecified
#    test (e.g. the leads nearest treatment) over this omnibus p.
.joint_pvalue <- function(theta0, boot_mat) {
  if (nrow(boot_mat) < 2L) return(NA_real_)
  support <- colSums(is.finite(boot_mat))
  floor_n <- max(2L, min(30L, as.integer(nrow(boot_mat) %/% 2L)))
  keep <- is.finite(theta0) & support >= floor_n
  if (!any(keep)) return(NA_real_)
  th <- theta0[keep]
  bm <- boot_mat[, keep, drop = FALSE]
  if (nrow(bm) < length(th) + 1L) return(NA_real_)
  V <- stats::cov(bm, use = "pairwise.complete.obs")
  if (any(!is.finite(V))) return(NA_real_)
  rc <- suppressWarnings(tryCatch(rcond(V), error = function(e) NA_real_))
  if (is.finite(rc) && rc < 1e-10) {
    warning(sprintf(
      paste0("didgpu: the %d-dimensional bootstrap covariance behind a joint ",
             "test is near-singular (rcond = %.1e); the omnibus chi-square ",
             "p-value is numerically unreliable. Prefer a low-dimensional ",
             "prespecified test (e.g. the leads nearest treatment)."),
      length(th), rc), call. = FALSE)
  }
  inv <- try(solve(V), silent = TRUE)
  if (inherits(inv, "try-error")) {
    inv <- MASS::ginv(V)
  }
  q <- as.numeric(t(th) %*% inv %*% th)
  stats::pchisq(q, df = length(th), lower.tail = FALSE)
}


# -------- printing, as DIDmultiplegtDYN prints --------
#
# print() and summary() reproduce print.did_multiplegt_dyn and mat_print
# from DIDmultiplegtDYN 2.4.0 (R/print.R), line for line, so a didgpu
# result reads exactly as the same model does under the reference.
# DIDmultiplegtDYN is MIT-licensed, Copyright (c) 2024 Diego Ciccia,
# Felix Knau, Melitine Malezieux, Doulo Sow, Clement de Chaisemartin.
# Two lines are left out: the closing acknowledgement of the European
# Union grant that funded DIDmultiplegtDYN, which did not fund didgpu.

.dcdh_mat_print <- function(mat) {
  if (inherits(mat, "matrix")) {
    dis <- matrix(data = 0, nrow = nrow(mat), ncol = ncol(mat))
    dis[, 1:4] <- sprintf("%s", format(round(mat[, 1:4], 5), big.mark = ",",
                                       scientific = FALSE, trim = TRUE))
    dis[, 5:ncol(dis)] <- sprintf("%s", format(round(mat[, 5:ncol(dis)], 0),
                                               big.mark = ",",
                                               scientific = FALSE, trim = TRUE))
    rownames(dis) <- rownames(mat)
    colnames(dis) <- colnames(mat)
    print(noquote(dis[, , drop = FALSE]))
  } else {
    dis <- vector(length = length(mat))
    dis[1:4] <- sprintf("%s", format(round(mat[1:4], 5), big.mark = ",",
                                     scientific = FALSE, trim = TRUE))
    dis[5:length(mat)] <- sprintf("%s", format(round(mat[5:length(mat)], 0),
                                               big.mark = ",",
                                               scientific = FALSE, trim = TRUE))
    names(dis) <- names(mat)
    print(noquote(dis[, drop = FALSE]))
  }
}

#' Print method for didgpu_result
#'
#' Prints the estimation tables exactly as
#' `DIDmultiplegtDYN::did_multiplegt_dyn()` prints its own: the event-study
#' effects and their joint test, the average total effect, and the
#' placebos and their joint test.
#'
#' @param x A `didgpu_result` object.
#' @param ... Unused (for S3 method compatibility).
#' @return The input invisibly.
#' @export
print.didgpu_result <- function(x, ...) {
  cat("\n")
  by_levels <- x$by_levels %||% "_no_by"
  for (b in seq_along(by_levels)) {
    if (by_levels[b] == "_no_by") {
      ref <- x
    } else {
      ref <- x[[paste0("by_level_", b)]]
      section <- if (!is.null(x$args[["by"]])) {
        paste(" By", x$args$by, "=", by_levels[b], "###")
      } else {
        paste0(" By treatment path: (", by_levels[b], ") ", "###")
      }
      cat(noquote(strrep("#", 70 - nchar(section) - 1))); cat(section)
      cat("\n"); cat("\n")
    }
    .print_level(ref, x)
  }
  cat("\n")
  invisible(x)
}

# One level's tables (the body of the reference's per-level loop).
.print_level <- function(ref, x) {
  ncol_show <- 6 + ((!is.null(x$args$weight)) * 2)
  rule <- function(n = 70) { cat(noquote(strrep("-", n))); cat("\n") }
  boot <- (x$args$bootstrap_reps %||% 0L) > 0L

  rule()
  cat(strrep(" ", 7)); cat("Estimation of treatment effects: Event-study effects"); cat("\n")
  rule()
  .dcdh_mat_print(ref$results$Effects[, 1:ncol_show])
  cat("\n")
  if (!is.null(ref$results$p_jointeffects)) {
    if (is.na(ref$results$p_jointeffects)) {
      cat("Test of joint nullity of the effects : p-value = not computed (see warnings)")
    } else {
      cat(sprintf("Test of joint nullity of the effects : p-value = %.4f",
                  ref$results$p_jointeffects))
    }
    cat("\n")
  }
  if (!is.null(ref$results$p_equality_effects)) {
    if (is.na(ref$results$p_equality_effects)) {
      cat("Test of equality of the effects : p-value = not computed (see warnings)")
    } else {
      cat(sprintf("Test of equality of the effects : p-value = %.4f",
                  ref$results$p_equality_effects))
    }
    cat("\n"); cat("\n")
  }

  if (isTRUE(x$args$trends_lin)) {
    rule()
    cat(strrep(" ", 4)); cat("When the trends_lin is specified no average effects are reported"); cat("\n")
    rule()
  } else {
    rule()
    cat(strrep(" ", 4)); cat("Average cumulative (total) effect per treatment unit"); cat("\n")
    rule()
    .dcdh_mat_print(ref$results$ATE[, 1:ncol_show])
    cat(sprintf("Average number of time periods over which a treatment effect is accumulated: %.4f",
                ref$results$delta_D_avg_total))
    cat("\n")
  }
  cat("\n")

  if (ref$results$N_Placebos != 0) {
    rule()
    cat(strrep(" ", 4)); cat(" Testing the parallel trends and no anticipation assumptions"); cat("\n")
    rule()
    .dcdh_mat_print(ref$results$Placebos[, 1:ncol_show])
    if (!boot) {
      cat("\n")
      if (!is.null(ref$results$p_jointplacebo) && is.na(ref$results$p_jointplacebo)) {
        cat("Test of joint nullity of the placebos : p-value = not computed (see warnings)")
      } else if (!is.null(ref$results$p_jointplacebo)) {
        cat(sprintf("Test of joint nullity of the placebos : p-value = %.4f",
                    ref$results$p_jointplacebo))
      }
      cat("\n")
    }
    cat("\n")
  }

  if (!is.null(ref$design)) {
    if (ref$design$design_path == "console") {
      cat("\n")
      rule()
      cat(strrep(" ", 4)); cat(sprintf("Detection of treatment paths - %.0f periods after first switch", ref$design$design_const[1])); cat("\n")
      rule()
      print(ref$design$design_mat); cat("\n")
      cat(sprintf("Treatment paths detected in at least %.2f%% of the %.0f switching groups for which %.0f effects could be estimated",
                  ref$design$design_const[2], ref$design$design_const[3], ref$design$design_const[1]))
      cat(sprintf(" (Total %% = %.2f%%)", ref$design$design_const[4])); cat("\n"); cat("\n")
      cat("Design interpretation (first row):"); cat("\n")
      n_groups <- ref$design$design_mat[1, 1]
      d_start <- ref$design$design_mat[1, 3]
      d_vec <- "("
      for (i in 1:ref$design$design_const[1]) {
        d_vec <- paste0(d_vec, ref$design$design_mat[1, 3 + i], ",")
      }
      d_vec <- paste0(substr(d_vec, 1, nchar(d_vec) - 1), ")")
      cat(sprintf("%s groups started with treatment %s and experienced treatment path %s", n_groups, d_start, d_vec))
      cat("\n")
    } else {
      cat(sprintf("Design exported to %s", ref$design$design_path)); cat("\n")
    }
  }

  if (!is.null(ref$date_first_switch)) {
    dfs <- ref$date_first_switch
    if (dfs$dfs_opt != "by_baseline_treat") {
      if (dfs$dfs_path == "console") {
        cat("\n")
        rule(40)
        cat(strrep(" ", 7)); cat("Switching dates"); cat("\n")
        rule(40)
        cat("By any status quo treatment"); cat("\n")
        print(dfs$dfs_mat)
        cat("\n")
      } else {
        cat(sprintf("Switching dates exported to %s", dfs$dfs_path)); cat("\n")
      }
    } else {
      if (dfs$dfs_path == "console") {
        cat("\n")
        rule(40)
        cat(strrep(" ", 7)); cat("Switching dates"); cat("\n")
        rule(40)
        for (l in 1:dfs$levels_baseline_treat) {
          cat(sprintf("Status quo treatment = %s", dfs[[paste0("level", l)]])); cat("\n")
          print(dfs[[paste0("dfs_mat", l)]])
          cat("\n")
        }
      }
      if (dfs$dfs_path != "console") {
        cat(sprintf("Switching dates exported to %s", dfs$dfs_path)); cat("\n")
      }
    }
  }

  if (!is.null(ref$normalized_weights)) {
    cat("\n")
    rule(60)
    cat(strrep(" ", 13)); cat("Weights on treatment lags"); cat("\n")
    rule(60)
    print(ref$normalized_weights$norm_weight_mat)
    cat("\n")
  }

  if (!is.null(ref$results$predict_het)) {
    cat("\n")
    rule(60)
    cat(strrep(" ", 13)); cat("Predicting effect heterogeneity"); cat("\n")
    rule(60)
    ph <- ref$results$predict_het
    .het_block <- function(tab, label) {
      for (l in levels(factor(tab$effect))) {
        het_tab <- subset(tab, tab$effect == l)
        het_mat <- as.matrix(het_tab[, c(3, 4, 6, 7, 8)])
        rownames(het_mat) <- het_tab$covariate
        colnames(het_mat) <- c("Estimate", "SE", "LB CI", "UB CI", "N")
        cat(sprintf("%s %s:\n", label, l))
        .dcdh_mat_print(het_mat)
        cat(sprintf("Test of joint nullity of the estimates : p-value = %.4f\n",
                    mean(het_tab$pF)))
        cat("\n")
      }
    }
    .het_block(subset(ph, ph$effect > 0), "Effect")
    pl <- subset(ph, ph$effect < 0)
    if (nrow(pl) > 0L) {
      pl$effect <- -pl$effect
      .het_block(pl, "Placebo")
    }
  }

  if (!is.null(ref$results$vcov_warnings)) {
    cat("\n")
    rule()
    cat(strrep(" ", 4)); cat("Warnings"); cat("\n")
    rule()
    for (w in ref$results$vcov_warnings) {
      cat(paste0("- ", w)); cat("\n")
    }
  }
  invisible(NULL)
}

#' Summary method for didgpu_result
#'
#' The same display as [print.didgpu_result()], as in DIDmultiplegtDYN.
#'
#' @param object A `didgpu_result` object.
#' @param ... Unused (for S3 method compatibility).
#' @return The input invisibly.
#' @export
summary.didgpu_result <- function(object, ...) {
  print.didgpu_result(object, ...)
}


#' Coefficient extractor for didgpu_result
#'
#' Returns a named numeric vector of point estimates. Names are
#' "Effect_1", ..., "Effect_n_effects" optionally followed by
#' "Placebo_1", ..., "Placebo_n_placebos" and (if available) "ATE".
#'
#' @param object A `didgpu_result` object.
#' @param which One of `"effects"`, `"placebos"`, `"ate"`, or `"all"`
#'   (default). Controls which subset is returned.
#' @param ... Unused.
#' @return Named numeric vector.
#' @export
coef.didgpu_result <- function(object, which = "all", ...) {
  which <- match.arg(which, c("all", "effects", "placebos", "ate"))
  pieces <- list()
  if (which %in% c("all", "effects")) {
    e <- object$results$Effects
    if (!is.null(e) && nrow(e) > 0L) {
      v <- as.numeric(e[, "Estimate"])
      names(v) <- trimws(rownames(e))
      pieces$effects <- v
    }
  }
  if (which %in% c("all", "placebos") && object$results$N_Placebos > 0L) {
    p <- object$results$Placebos
    if (!is.null(p) && nrow(p) > 0L) {
      v <- as.numeric(p[, "Estimate"])
      names(v) <- trimws(rownames(p))
      pieces$placebos <- v
    }
  }
  if (which %in% c("all", "ate") && !is.null(object$results$ATE)) {
    a <- object$results$ATE
    if (!is.na(a[1, "Estimate"])) {
      v <- as.numeric(a[, "Estimate"])
      names(v) <- "ATE"
      pieces$ate <- v
    }
  }
  v <- unlist(pieces, use.names = TRUE)
  # Strip the leading list-name prefix ("effects.", "placebos.", "ate.")
  # that unlist() injects when the parent list is named.
  names(v) <- sub("^(effects|placebos|ate)\\.", "", names(v))
  v
}


#' Confidence intervals for didgpu_result
#'
#' Returns the percentile-based CI matrix already computed during
#' aggregation (from the bootstrap distribution). Re-computing at a
#' different level requires a fresh run because the cell-level
#' percentiles are not stored on disk.
#'
#' @param object A `didgpu_result` object.
#' @param parm Optional character vector of coefficient names to subset.
#' @param level Confidence level. Must match the level used at fit time;
#'   otherwise a warning is issued and the stored CI is returned anyway.
#' @param ... Unused.
#' @return A 2-column matrix with the lower and upper bounds.
#' @export
confint.didgpu_result <- function(object, parm = NULL, level = NULL, ...) {
  fit_level <- object$args$ci_level %||% 95
  if (!is.null(level) && abs(level * 100 - fit_level) > 1e-9 &&
      abs(level     - fit_level) > 1e-9) {
    warning(sprintf("Stored CIs are at level %.1f%% (from the fit); ",
                    fit_level),
            "requested level=", level, " is ignored. Re-run didgpu() ",
            "with ci_level = ", level, " to change.")
  }
  # Pick lower/upper-bound columns from the printed tables.
  lo_hi_from <- function(m) {
    if (is.null(m) || nrow(m) == 0L) return(NULL)
    cn <- colnames(m)
    lo <- grep("^LB", cn)[1L]; hi <- grep("^UB", cn)[1L]
    out <- cbind(as.numeric(m[, lo]), as.numeric(m[, hi]))
    rownames(out) <- trimws(rownames(m))
    out
  }
  pieces <- list()
  pieces$effects <- lo_hi_from(object$results$Effects)
  if (object$results$N_Placebos > 0L) pieces$placebos <- lo_hi_from(object$results$Placebos)
  if (!is.null(object$results$ATE) &&
      !is.na(object$results$ATE[1, "Estimate"])) {
    pieces$ate <- lo_hi_from(object$results$ATE)
    if (!is.null(pieces$ate)) rownames(pieces$ate) <- "ATE"
  }
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) == 0L) return(matrix(numeric(0), nrow = 0L, ncol = 2L,
                                           dimnames = list(NULL,
                                                            c("LB", "UB"))))
  out <- do.call(rbind, pieces)
  colnames(out) <- c(sprintf("%g%% LB", fit_level),
                     sprintf("%g%% UB", fit_level))
  if (!is.null(parm)) {
    miss <- setdiff(parm, rownames(out))
    if (length(miss)) {
      warning("parm not in result: ", paste(miss, collapse = ", "))
    }
    out <- out[intersect(parm, rownames(out)), , drop = FALSE]
  }
  out
}


#' Event-study plot of a didgpu result
#'
#' Plots estimates against event-time horizon: pre-treatment placebos
#' at negative horizons, post-treatment effects at positive horizons,
#' with the stored CIs as vertical error bars and a horizontal dashed
#' line at 0 for reference. Uses base R graphics — no ggplot2
#' dependency. Returns the input invisibly so calls can be chained.
#'
#' @param x A `didgpu_result` object.
#' @param ... Extra graphical parameters passed to the underlying
#'   `plot()` call (e.g. `main`, `xlab`, `ylab`, `col`, `pch`, `lwd`,
#'   `xlim`, `ylim`).
#' @param show_zero_line Logical. Draw a dashed line at y = 0
#'   (default TRUE).
#' @param show_zero_horizon Logical. Mark the boundary between placebo
#'   and effect horizons with a vertical dashed line (default TRUE).
#' @param ci Logical. Draw error bars at the stored CI level
#'   (default TRUE; suppressed when bootstrap_reps = 0 because the
#'   CIs are NA).
#' @return The input invisibly.
#' @export
plot.didgpu_result <- function(x, ...,
                                show_zero_line = TRUE,
                                show_zero_horizon = TRUE,
                                ci = TRUE) {
  e  <- x$results$Effects
  pl <- if (x$results$N_Placebos > 0L) x$results$Placebos else NULL
  if ((is.null(e) || nrow(e) == 0L) && is.null(pl)) {
    stop("Nothing to plot: result has no effects or placebos.")
  }

  # Build the (horizon, estimate, LB, UB) frame.
  rows_eff <- if (!is.null(e) && nrow(e) > 0L) {
    data.frame(h  = seq_len(nrow(e)),
               y  = as.numeric(e[, "Estimate"]),
               lb = as.numeric(e[, grep("^LB", colnames(e))[1L]]),
               ub = as.numeric(e[, grep("^UB", colnames(e))[1L]]))
  } else NULL
  rows_pl <- if (!is.null(pl) && nrow(pl) > 0L) {
    data.frame(h  = -seq_len(nrow(pl)),
               y  = as.numeric(pl[, "Estimate"]),
               lb = as.numeric(pl[, grep("^LB", colnames(pl))[1L]]),
               ub = as.numeric(pl[, grep("^UB", colnames(pl))[1L]]))
  } else NULL
  df <- rbind(rows_pl, rows_eff)
  df <- df[order(df$h), , drop = FALSE]

  # Suppress CIs if no bootstrap was run.
  if (isTRUE(ci) && (is.null(x$args$bootstrap_reps) ||
                      x$args$bootstrap_reps == 0L)) ci <- FALSE

  # Defaults for graphical params, overridable via ...
  dots <- list(...)
  if (is.null(dots$xlab)) dots$xlab <- "Event-time horizon"
  if (is.null(dots$ylab)) dots$ylab <- "Estimate"
  if (is.null(dots$main)) dots$main <- "didgpu event-study"
  if (is.null(dots$pch))  dots$pch  <- 19L
  if (is.null(dots$col))  dots$col  <- "black"
  if (is.null(dots$xlim)) dots$xlim <- range(df$h) + c(-0.5, 0.5)
  if (is.null(dots$ylim)) {
    yvals <- if (isTRUE(ci)) c(df$y, df$lb, df$ub) else df$y
    yvals <- yvals[is.finite(yvals)]
    if (length(yvals) == 0L) yvals <- c(-1, 1)
    pad <- 0.05 * diff(range(yvals))
    if (pad == 0) pad <- 0.1 * max(1, abs(yvals[1L]))
    dots$ylim <- range(yvals) + c(-pad, pad)
  }

  do.call(plot,
          c(list(df$h, df$y, type = "p"), dots))

  if (isTRUE(show_zero_line)) {
    graphics::abline(h = 0, lty = 2L, col = "grey50")
  }
  if (isTRUE(show_zero_horizon)) {
    graphics::abline(v = 0.5, lty = 2L, col = "grey50")
  }
  if (isTRUE(ci)) {
    graphics::segments(df$h, df$lb, df$h, df$ub,
                       col = dots$col, lwd = max(1L, dots$lwd %||% 1L))
    # Endcaps.
    w <- 0.1
    graphics::segments(df$h - w, df$lb, df$h + w, df$lb,
                       col = dots$col, lwd = max(1L, dots$lwd %||% 1L))
    graphics::segments(df$h - w, df$ub, df$h + w, df$ub,
                       col = dots$col, lwd = max(1L, dots$lwd %||% 1L))
  }

  invisible(x)
}


#' Variance-covariance matrix of estimates
#'
#' Returns the covariance matrix of `(Effects, Placebos)`, computed at fit
#' time and stored on the result. It is analytic -- built from the same
#' influence functions as the reported SEs, so its diagonal equals
#' `SE^2` exactly and it needs no bootstrap. Where analytic SEs are not
#' available (`controls`, `continuous`, `trends_lin`) it is the empirical
#' covariance of the bootstrap replicates instead, and a square NA matrix
#' if `bootstrap_reps` is 0 or 1.
#'
#' The ordering matches `coef(object)`'s default (effects first, then
#' placebos). The ATE row/column is NOT included — it is a linear
#' combination of the per-event-time effects, so its variance can be
#' recovered as `t(w) %*% vcov(object) %*% w` where `w` is the
#' incidence-weighted vector.
#'
#' @param object A `didgpu_result` object.
#' @param ... Unused.
#' @return A square matrix.
#' @export
vcov.didgpu_result <- function(object, ...) {
  V <- object$coef$vcov
  if (is.null(V)) {
    nm <- names(object$coef$b)
    return(matrix(NA_real_, length(nm), length(nm),
                  dimnames = list(nm, nm)))
  }
  V
}


# The event-study graph DIDmultiplegtDYN draws (did_multiplegt_dyn_graph,
# R/did_multiplegt_dyn_graph.R in 2.4.0): estimates against time relative
# to the last period before the switch, the ATE row pinned at (0, 0),
# red CI bars except at t = 0. Built only when ggplot2 and cowplot are
# installed, as they are wherever DIDmultiplegtDYN is.
.dcdh_graph <- function(results, ggplot_args = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("cowplot", quietly = TRUE)) return(NULL)
  grmat <- rbind(cbind(results$Effects, seq_len(nrow(results$Effects))),
                 cbind(results$ATE, 0))
  if (!is.null(results$Placebos)) {
    grmat <- rbind(grmat, cbind(results$Placebos, -seq_len(nrow(results$Placebos))))
  }
  colnames(grmat)[ncol(grmat)] <- "Time"
  grmat[nrow(results$Effects) + 1, c(1, 3, 4)] <- 0
  grmat <- data.frame(grmat[, c(1, 3, 4, 9)])
  keep <- grmat$Estimate != 0
  p <- ggplot2::ggplot(grmat, ggplot2::aes(x = Time, y = Estimate, group = 1)) +
    ggplot2::geom_line(colour = "blue") +
    ggplot2::geom_errorbar(data = function(x) x[keep, , drop = FALSE],
                           ggplot2::aes(ymin = LB.CI, ymax = UB.CI),
                           position = ggplot2::position_dodge(0.05),
                           width = 0.2, colour = "red") +
    ggplot2::geom_point(colour = "blue") +
    ggplot2::ggtitle("DID, from last period before treatment changes (t=0) to t") +
    ggplot2::xlab("Relative time to last period before treatment changes (t=0)") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    cowplot::theme_minimal_grid()
  for (layer in ggplot_args) p <- p + layer
  p
}
