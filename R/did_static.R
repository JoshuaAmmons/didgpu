# ============================================================================
# didgpu_did_static(): de Chaisemartin & D'Haultfoeuille (2020) DID_M.
#
# The instantaneous DiD estimator that ALLOWS treatment to turn on AND off
# (non-absorbing), unlike the staggered-adoption methods (didgpu_cs(),
# didgpu_bacon()). For each pair of consecutive periods it compares the
# t-1 -> t outcome change of units that SWITCH treatment to the change of
# "stayers" with the same period-(t-1) treatment, then averages over all
# switch events:
#
#   DID_M = (1 / N_S) * sum_t [ N+_t * DID+_t + N-_t * DID-_t ]
#   DID+_t = dYbar(switch 0->1 at t) - dYbar(stay 0->0 at t)     (gaining D)
#   DID-_t = dYbar(stay 1->1 at t)   - dYbar(switch 1->0 at t)   (holding D)
#
# where dYbar() is the (weighted) mean of Y_t - Y_{t-1}. DID+_t is the
# effect of gaining treatment; DID-_t recovers the effect of treatment from
# the switch-OUT direction (stayers-treated minus switchers-out). N+_t / N-_t
# are the (weighted) switcher counts and N_S the total switchers that have a
# valid same-baseline counterfactual. Under common trends + no anticipation
# DID_M consistently estimates the ATT among switchers at the switch period,
# even with heterogeneous effects and dynamic/reversing treatment.
#
# Reference: de Chaisemartin & D'Haultfoeuille (2020), "Two-Way Fixed Effects
# Estimators with Heterogeneous Treatment Effects", AER 110(9): 2964-2996.
# Cross-validated against DIDmultiplegt::did_multiplegt (the DID_M estimand).
# ============================================================================


# Point estimator on a prepped data.table (cols G, T, D in {0,1}, Y, W).
# Returns list(did, n_switchers, per_period). Uses the previous OBSERVED
# period within each unit as t-1 (exact for panels with consecutive-integer
# periods; the standard DID_M setting).
#' @keywords internal
#' @noRd
.didm_point <- function(d) {
  data.table::setkey(d, G, T)
  d[, `:=`(D_lag = data.table::shift(D, 1L),
           Y_lag = data.table::shift(Y, 1L)), by = G]
  dd <- d[!is.na(D_lag)]
  if (nrow(dd) == 0L) return(list(did = NA_real_, n_switchers = 0, per_period = NULL))
  dd[, dY := Y - Y_lag]
  dd[, trans := data.table::fcase(
        D_lag == 0 & D == 1, "in",
        D_lag == 1 & D == 0, "out",
        D_lag == 0 & D == 0, "stay0",
        default = "stay1")]

  wmean <- function(x, w) if (length(x) == 0L || sum(w) == 0) NA_real_ else sum(x * w) / sum(w)
  num <- 0; den <- 0; per <- list()
  for (tt in sort(unique(dd$T))) {
    sub <- dd[T == tt]
    s_in  <- sub[trans == "in"];    s_out <- sub[trans == "out"]
    st0   <- sub[trans == "stay0"]; st1   <- sub[trans == "stay1"]
    n_in  <- sum(s_in$W);  n_out <- sum(s_out$W)
    didp <- if (n_in  > 0 && nrow(st0) > 0)
              wmean(s_in$dY,  s_in$W)  - wmean(st0$dY, st0$W) else NA_real_
    didm <- if (n_out > 0 && nrow(st1) > 0)
              wmean(st1$dY, st1$W) - wmean(s_out$dY, s_out$W) else NA_real_
    if (!is.na(didp)) { num <- num + n_in  * didp; den <- den + n_in }
    if (!is.na(didm)) { num <- num + n_out * didm; den <- den + n_out }
    per[[length(per) + 1L]] <- data.frame(
      time = tt, n_switch_in = n_in, n_switch_out = n_out,
      did_in = didp, did_out = didm, stringsAsFactors = FALSE)
  }
  list(did = if (den > 0) num / den else NA_real_,
       n_switchers = den,
       per_period = do.call(rbind, per))
}


#' de Chaisemartin-D'Haultfoeuille (2020) DID_M instantaneous estimator
#'
#' The DID_M estimator for a binary treatment that may switch on AND off
#' (non-absorbing). It compares each switching unit's period-over-period
#' outcome change to that of stayers with the same prior-period treatment,
#' and averages over all switch events. Consistent for the ATT among
#' switchers under common trends + no anticipation, even with heterogeneous
#' and dynamic treatment effects (where TWFE is biased). Standard errors are
#' from a cluster (by default unit) bootstrap, matching the package's other
#' estimators.
#'
#' @param df A data.frame / data.table panel.
#' @param outcome,group,time,treatment Character column names: outcome, unit
#'   id, time id, and a binary 0/1 (possibly non-absorbing) treatment.
#' @param weight Optional character column name of observation weights.
#' @param bootstrap_reps Integer; cluster-bootstrap replicates for the SE
#'   (0 to skip and return only the point estimate). Default 100.
#' @param cluster Optional character column name to resample in the
#'   bootstrap. Defaults to `group`.
#' @param ci_level Confidence level in percent (default 95).
#' @param seed RNG seed for the bootstrap.
#' @param verbose Logical.
#' @return A `didgpu_did_static_result`: a list with `did` (the DID_M point
#'   estimate), `se`, `ci`, `n_switchers`, `per_period` (per-period switcher
#'   counts and directional DiDs), and `args`.
#' @references de Chaisemartin, C. & D'Haultfoeuille, X. (2020). Two-Way
#'   Fixed Effects Estimators with Heterogeneous Treatment Effects.
#'   \emph{American Economic Review} 110(9): 2964-2996.
#' @seealso [didgpu()] for the dynamic generalization, [didgpu_twfe()] /
#'   [didgpu_bacon()] for the biased TWFE baseline and its decomposition.
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 80L, n_periods = 8L,
#'                            tau_profile = c(0.5, 1.0), seed = 5L)
#' p$D <- as.integer(p$D >= 0.5)
#' didgpu_did_static(p, "Y", "unit", "period", "D", bootstrap_reps = 0L)
#' }
#' @export
didgpu_did_static <- function(df, outcome, group, time, treatment,
                              weight = NULL, bootstrap_reps = 100L,
                              cluster = NULL, ci_level = 95, seed = 1L,
                              verbose = TRUE) {
  stopifnot(is.data.frame(df) || data.table::is.data.table(df))
  for (nm in c("outcome", "group", "time", "treatment")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v))
      stop("`", nm, "` must be a single non-empty column name.")
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  cluster_col <- cluster %||% group
  if (!cluster_col %in% names(df)) stop("cluster column not in df: ", cluster_col)
  bootstrap_reps <- as.integer(bootstrap_reps)

  src <- data.table::as.data.table(df)
  d <- data.table::data.table(
    Y  = as.numeric(src[[outcome]]),
    G  = src[[group]],
    T  = src[[time]],
    D  = as.numeric(src[[treatment]]),
    W  = if (is.null(weight)) 1.0 else as.numeric(src[[weight]]),
    CL = src[[cluster_col]])
  d <- d[!is.na(Y) & !is.na(D) & !is.na(G) & !is.na(T)]
  if (!all(d$D %in% c(0, 1)))
    stop("didgpu_did_static requires a binary 0/1 treatment.")

  pt <- .didm_point(data.table::copy(d))
  if (is.na(pt$did))
    warning("DID_M is NA: no switchers with a valid same-baseline ",
            "counterfactual were found.")

  se <- NA_real_; ci <- c(NA_real_, NA_real_)
  if (bootstrap_reps > 0L && !is.na(pt$did)) {
    clusters <- unique(d$CL)
    nc <- length(clusters)
    # Pre-split row indices by cluster ONCE. The naive approach filters
    # `d[CL == draw[i]]` for every drawn cluster on every replicate, which is
    # O(n_units^2 * T) per rep and dominates the runtime (a B=1000 bootstrap at
    # 1000 units took ~21 min). Splitting once turns each replicate into an O(n)
    # gather + vectorized relabel. Results are bit-identical: the RNG draw
    # sequence is unchanged and each drawn copy still gets a distinct group id.
    by_cl <- split(seq_len(nrow(d)), factor(d$CL, levels = clusters))
    set.seed(seed)
    boot <- numeric(bootstrap_reps)
    for (b in seq_len(bootstrap_reps)) {
      draw <- sample(clusters, nc, replace = TRUE)
      idx_list <- by_cl[as.character(draw)]
      all_idx  <- unlist(idx_list, use.names = FALSE)
      block    <- rep.int(seq_len(nc), lengths(idx_list))  # distinct id per draw copy
      db <- d[all_idx]
      db[, G := paste0(G, "__b", block)]                   # relabel drawn copies uniquely
      boot[b] <- .didm_point(db)$did
    }
    boot <- boot[is.finite(boot)]
    if (length(boot) >= 2L) {
      se <- stats::sd(boot)
      a <- (1 - ci_level / 100) / 2
      ci <- pt$did + c(stats::qnorm(a), stats::qnorm(1 - a)) * se
    }
  }

  out <- structure(list(
    did = pt$did, se = se, ci = ci,
    n_switchers = pt$n_switchers,
    per_period = pt$per_period,
    args = list(outcome = outcome, group = group, time = time,
                treatment = treatment, weight = weight,
                cluster = cluster_col, bootstrap_reps = bootstrap_reps,
                ci_level = ci_level, n_obs = nrow(d))
  ), class = c("didgpu_did_static_result", "list"))
  if (verbose) print(out)
  invisible(out)
}


#' Print method for didgpu_did_static_result
#' @param x A `didgpu_did_static_result`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.didgpu_did_static_result <- function(x, ...) {
  cat("de Chaisemartin-D'Haultfoeuille (2020) DID_M (instantaneous, reversible)\n")
  cat(sprintf("  DID_M = %.6f", x$did))
  if (!is.na(x$se))
    cat(sprintf("   SE = %.6f   %d%% CI [%.4f, %.4f]",
                x$se, x$args$ci_level, x$ci[1], x$ci[2]))
  cat(sprintf("\n  switchers (weighted): %.4g   over %d obs\n",
              x$n_switchers, x$args$n_obs))
  cat("  Allows treatment to turn on AND off; robust to heterogeneous,\n")
  cat("  dynamic effects (unlike TWFE). See didgpu() for dynamic effects.\n")
  invisible(x)
}
