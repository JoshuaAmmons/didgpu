# ============================================================================
# Benchmark the 8 newer estimators (#98-#105) to verify the CPU-only decision
# empirically (Caveat 2 of the GPU dispatch audit). For each, measure CPU cost
# and scaling at realistic + large sizes, decompose point-estimate vs
# bootstrap-loop time, and report per-replicate cost. The verdict on
# GPU-worthiness is rendered against the package's documented GPU economics:
#   - cluster bootstrap that refits an EXPENSIVE per-rep estimate -> GPU ~179x
#     (CS), IF-shortcut + batched matmul.
#   - IF-based multiplier bootstrap                                -> GPU ~1.5x.
#   - per-rep dominated by data.table reshaping / tiny ops         -> GPU LOSES
#     (LOO lesson); algorithmic/parallel R wins instead.
#   - sub-second absolute time                                     -> nothing to
#     accelerate (GPU launch floor ~0.1s).
# ============================================================================
suppressWarnings(suppressMessages(library(didgpu)))
suppressWarnings(suppressMessages(library(data.table)))

timeit <- function(fn, reps = 3L) {
  ts <- numeric(reps)
  for (i in seq_len(reps)) { gc(FALSE); t0 <- Sys.time(); fn(); ts[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs")) }
  stats::median(ts)
}

# ---- DGPs -----------------------------------------------------------------
make_binary_nonabsorbing <- function(nU, Tn, seed = 1L) {  # did_static (on/off)
  set.seed(seed)
  id <- rep(seq_len(nU), each = Tn); t <- rep(seq_len(Tn), nU)
  a <- rnorm(nU)
  D <- rbinom(nU * Tn, 1L, 0.4)                            # switches on AND off
  Y <- a[id] + 0.1 * t + 0.5 * D + rnorm(nU * Tn, 0, 0.5)
  data.frame(id = id, t = t, D = D, Y = Y)
}
make_staggered_absorbing <- function(nU, Tn, seed = 1L) {  # twfe / bacon
  set.seed(seed)
  id <- rep(seq_len(nU), each = Tn); t <- rep(seq_len(Tn), nU)
  a <- rnorm(nU)
  cohort <- sample(c(Inf, seq(2L, Tn - 1L)), nU, replace = TRUE)
  D <- as.integer(t >= cohort[id])
  Y <- a[id] + 0.1 * t + 0.5 * D + rnorm(nU * Tn, 0, 0.5)
  data.frame(id = id, t = t, D = D, Y = Y)
}
make_fhs <- function(nU, Tn = 24L, seed = 1L) {            # freyaldenhoven
  set.seed(seed)
  d <- data.frame(id = rep(1:nU, each = Tn), t = rep(1:Tn, times = nU))
  adopt <- sample(c(rep(Inf, nU %/% 4L), sample(6:20, nU - nU %/% 4L, replace = TRUE)), nU)
  d$z <- as.integer(d$t >= adopt[d$id])
  eta <- rnorm(nU)[d$id] + 0.2 * d$t
  d$x <- eta + rnorm(nrow(d), 0, 0.3)
  d$y <- rnorm(nU)[d$id] + 0.1 * d$t + 0.5 * d$z + 0.4 * eta + rnorm(nrow(d), 0, 0.3)
  d
}
make_cont_dose <- function(nU, seed = 1L) {                # cs_continuous
  set.seed(seed); half <- nU %/% 2L
  G <- rep(c(0, 2), times = c(nU - half, half))
  D <- ifelse(G == 2, runif(nU, 0.1, 1), 0)
  a <- rnorm(nU)
  d <- data.frame(id = rep(seq_len(nU), each = 2L), time_period = rep(1:2, times = nU),
                  G = rep(G, each = 2L), D = rep(D, each = 2L))
  d$Y <- a[d$id] + 0.3 * d$time_period +
         ifelse(d$time_period == 2L, 2 * d$D - 1.2 * d$D^2, 0) + rnorm(nrow(d), 0, 0.3)
  d
}
make_nostayer <- function(nU, seed = 1L) {                 # did_continuous
  set.seed(seed)
  dD <- rnorm(nU, 0, 1); dY <- 0.3 + 2 * dD - 0.5 * dD^2 + rnorm(nU, 0, 0.5)
  data.frame(id = rep(seq_len(nU), each = 2L), t = rep(1:2, nU),
             D = as.numeric(rbind(0, dD)), Y = as.numeric(rbind(0, dY)))
}

rows <- list()
add <- function(est, scn, nU, Tn, B, point_s, total_s, note = "") {
  boot_s <- if (is.na(point_s) || is.na(total_s)) NA_real_ else max(total_s - point_s, 0)
  per_rep_ms <- if (is.na(boot_s) || is.na(B) || B == 0) NA_real_ else boot_s / B * 1000
  row <- data.frame(
    estimator = est, scenario = scn, n_units = nU, n_periods = Tn, B = B,
    point_s = round(point_s, 4), total_s = round(total_s, 4),
    boot_s = round(boot_s, 4), per_rep_ms = round(per_rep_ms, 3),
    note = note, stringsAsFactors = FALSE)
  rows[[length(rows) + 1L]] <<- row
  print(row, row.names = FALSE)                            # progress
  try(saveRDS(do.call(rbind, rows), "C:/Users/ammonsj/bench_new_estimators.rds"), silent = TRUE)  # incremental
}

cat("Running benchmarks (median of 3, fewer for the largest configs)...\n\n")

## ---- did_static (cluster bootstrap, refit .didm_point per rep) ----
# n=5000 uses B=100 (per-rep cost is B-independent; extrapolate total to B=1000).
for (cfg in list(c(1000, 12, 200), c(1000, 12, 1000), c(5000, 20, 100))) {
  nU <- cfg[1]; Tn <- cfg[2]; B <- cfg[3]
  d <- make_binary_nonabsorbing(nU, Tn)
  reps <- if (nU >= 5000) 1L else 3L
  pt <- timeit(function() didgpu_did_static(d, "Y", "id", "t", "D", bootstrap_reps = 0L, verbose = FALSE), reps)
  tt <- timeit(function() didgpu_did_static(d, "Y", "id", "t", "D", bootstrap_reps = B, seed = 1L, verbose = FALSE), reps)
  add("did_static", sprintf("n%d T%d B%d", nU, Tn, B), nU, Tn, B, pt, tt, "refit per rep")
}

## ---- cs_continuous (multiplier/Rademacher bootstrap, IF-based) ----
for (cfg in list(c(1000, 200), c(1000, 1000), c(5000, 1000))) {
  nU <- cfg[1]; B <- cfg[2]
  d <- make_cont_dose(nU)
  reps <- if (nU >= 5000) 2L else 3L
  pt <- timeit(function() didgpu_cs_continuous(d, "Y", "D", "G", "time_period", "id", degree = 3, num_knots = 2, bootstrap_reps = 0L, verbose = FALSE), reps)
  tt <- timeit(function() didgpu_cs_continuous(d, "Y", "D", "G", "time_period", "id", degree = 3, num_knots = 2, bootstrap_reps = B, seed = 1L, verbose = FALSE), reps)
  add("cs_continuous", sprintf("n%d B%d", nU, B), nU, 2L, B, pt, tt, "multiplier IF")
}

## ---- did_continuous (resample bootstrap, refit per rep) ----
for (cfg in list(c(2000, 200), c(2000, 1000), c(10000, 1000))) {
  nU <- cfg[1]; B <- cfg[2]
  d <- make_nostayer(nU)
  reps <- if (nU >= 10000) 1L else 3L
  pt <- timeit(function() didgpu_did_continuous(d, "Y", "D", "id", "t", estimator = "parametric", degree = 2, bootstrap_reps = 0L, verbose = FALSE), reps)
  tt <- timeit(function() didgpu_did_continuous(d, "Y", "D", "id", "t", estimator = "parametric", degree = 2, bootstrap_reps = B, seed = 1L, verbose = FALSE), reps)
  add("did_continuous(par)", sprintf("n%d B%d", nU, B), nU, 2L, B, pt, tt, "refit per rep (poly OLS)")
}
# nonparametric local-linear (heavier per rep)
{
  d <- make_nostayer(4000)
  pt <- timeit(function() didgpu_did_continuous(d, "Y", "D", "id", "t", estimator = "nonparametric", bootstrap_reps = 0L, verbose = FALSE), 3L)
  tt <- timeit(function() didgpu_did_continuous(d, "Y", "D", "id", "t", estimator = "nonparametric", bootstrap_reps = 200L, seed = 1L, verbose = FALSE), 3L)
  add("did_continuous(np)", "n4000 B200", 4000L, 2L, 200L, pt, tt, "refit per rep (local-linear)")
}

## ---- twfe (analytic CR1, no bootstrap) ----
for (cfg in list(c(1000, 12), c(5000, 20))) {
  nU <- cfg[1]; Tn <- cfg[2]
  d <- make_staggered_absorbing(nU, Tn)
  tt <- timeit(function() didgpu_twfe(d, "Y", "id", "t", "D", effects = 4L, placebo = 2L, verbose = FALSE), 3L)
  add("twfe", sprintf("n%d T%d", nU, Tn), nU, Tn, 0L, NA_real_, tt, "analytic CR1, no boot")
}

## ---- freyaldenhoven (analytic CR1/CR2, no bootstrap) ----
for (cfg in list(c(200, 24), c(1000, 24))) {
  nU <- cfg[1]; Tn <- cfg[2]
  d <- make_fhs(nU, Tn)
  ttO <- timeit(function() didgpu_freyaldenhoven(d, "y", "z", "id", "t", estimator = "OLS", pre = 0, post = 3, verbose = FALSE), 3L)
  add("freyaldenhoven(OLS)", sprintf("n%d T%d", nU, Tn), nU, Tn, 0L, NA_real_, ttO, "analytic, no boot")
  ttF <- timeit(function() didgpu_freyaldenhoven(d, "y", "z", "id", "t", estimator = "FHS", proxy = "x", pre = 0, post = 3, verbose = FALSE), 3L)
  add("freyaldenhoven(FHS)", sprintf("n%d T%d", nU, Tn), nU, Tn, 0L, NA_real_, ttF, "analytic 2SLS, no boot")
}

## ---- bacon (closed-form 2x2, no SE) ----
for (cfg in list(c(1000, 12), c(5000, 20))) {
  nU <- cfg[1]; Tn <- cfg[2]
  d <- make_staggered_absorbing(nU, Tn)
  tt <- timeit(function() didgpu_bacon(d, "Y", "id", "t", "D"), 3L)
  add("bacon", sprintf("n%d T%d", nU, Tn), nU, Tn, 0L, NA_real_, tt, "closed-form")
}

## ---- equivalence + joint_placebo (post-hoc on stored vcov) ----
{
  p <- didgpu_simulate_panel(n_units = 200L, n_periods = 10L, tau_profile = c(0.5, 1), seed = 1L)
  fit <- didgpu(p, "Y", "unit", "period", "D", effects = 3L, placebo = 3L,
                bootstrap_reps = 100L, seed = 1L, backend = "r", verbose = FALSE)
  teq <- timeit(function() didgpu_equivalence(fit, delta = 0.5), 5L)
  add("equivalence", "posthoc(fit B100)", 200L, 10L, 0L, NA_real_, teq, "post-hoc on stored vcov")
  tjp <- timeit(function() didgpu_joint_placebo(fit, horizons = 1:3), 5L)
  add("joint_placebo", "posthoc(fit B100)", 200L, 10L, 0L, NA_real_, tjp, "post-hoc on stored vcov")
}

res <- do.call(rbind, rows)
cat("\n================= RESULTS =================\n")
print(res, row.names = FALSE)
saveRDS(res, "C:/Users/ammonsj/bench_new_estimators.rds")
cat("\nSaved to C:/Users/ammonsj/bench_new_estimators.rds\n")
