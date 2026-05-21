# Is didgpu_twfe() fast? Compare against the obvious base-R alternative
# for a two-way FE event study: lm(Y ~ leads/lags + factor(unit) +
# factor(time)). (fixest would be the specialized C++ benchmark, but
# it's not installed in this env; see notes in the printed output.)
suppressPackageStartupMessages({ library(didgpu); library(data.table) })

make_panel <- function(n_units, n_periods, seed = 5L) {
  set.seed(seed)
  ufe <- rnorm(n_units); tfe <- rnorm(n_periods, sd = 0.5)
  g <- rep(seq_len(n_units), each = n_periods)
  t <- rep(seq_len(n_periods), times = n_units)
  D <- as.integer(runif(n_units * n_periods) < 0.25)
  Y <- ufe[g] + tfe[t] + 1.0 * D + rnorm(n_units * n_periods, sd = 0.3)
  data.frame(unit = g, period = t, D = D, Y = Y)
}

# Base-R lm with explicit lead/lag + FE dummies (point estimates only).
lm_twfe <- function(p, effects, placebo) {
  d <- as.data.table(p); setkey(d, unit, period)
  regs <- character(0)
  for (k in seq_len(effects)) {
    nm <- paste0("E", k); src <- d[, list(unit, period = period + (k-1L), vv = D)]
    d[src, (nm) := i.vv, on = c("unit","period")]; d[is.na(get(nm)), (nm) := 0]
    regs <- c(regs, nm)
  }
  for (j in seq_len(placebo)) {
    nm <- paste0("P", j); src <- d[, list(unit, period = period - j, vv = D)]
    d[src, (nm) := i.vv, on = c("unit","period")]; d[is.na(get(nm)), (nm) := 0]
    regs <- c(regs, nm)
  }
  f <- as.formula(paste0("Y ~ ", paste(regs, collapse = "+"),
                         " + factor(unit) + factor(period)"))
  stats::lm(f, data = as.data.frame(d))
}

tt <- function(expr) { t0 <- Sys.time(); force(expr)
  as.numeric(difftime(Sys.time(), t0, units = "secs")) }

cat(sprintf("%-6s %-4s  %10s  %10s  %8s\n",
            "units","per","didgpu_twfe","lm+dummies","speedup"))
for (nu in c(100L, 500L, 2000L, 5000L)) {
  np <- 12L
  p <- make_panel(nu, np)
  t_dg <- tt(didgpu_twfe(p, "Y","unit","period","D",
                         effects = 4L, placebo = 3L, verbose = FALSE))
  # lm with thousands of unit dummies gets very slow / memory-heavy;
  # cap it so the benchmark itself doesn't take forever.
  t_lm <- if (nu <= 2000L) tt(lm_twfe(p, 4L, 3L)) else NA_real_
  sp <- if (is.na(t_lm)) NA else t_lm / t_dg
  cat(sprintf("%-6d %-4d  %9.3fs  %9s  %7s\n", nu, np, t_dg,
              if (is.na(t_lm)) "  (skip)" else sprintf("%8.3fs", t_lm),
              if (is.na(sp)) "  -" else sprintf("%6.1fx", sp)))
}
cat("\nfixest::feols would be the specialized-tool benchmark (optimized\n")
cat("C++ demeaning, multithreaded); not installed here. didgpu_twfe's\n")
cat("value is integration with the robust estimators + same output\n")
cat("shape, not beating fixest on raw TWFE speed.\n")
