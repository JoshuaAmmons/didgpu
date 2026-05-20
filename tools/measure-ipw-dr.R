# Measure CUDA-vs-R agreement for the new IPW/DR GPU kernel, across
# several panel sizes WITH covariates (the path that exercises the
# logistic regression). Prints max abs ATT diff per method/size so we
# can choose a defensible equivalence tolerance.
suppressPackageStartupMessages({ library(didgpu) })
if (!didgpu_has_cuda_support()) { cat("no CUDA\n"); quit(status = 0) }

make_panel <- function(n_units, n_periods, seed = 17L) {
  p <- didgpu_simulate_panel(n_units = n_units, n_periods = n_periods,
                              tau_profile = c(0.5, 1.0), seed = seed)
  p$D <- as.integer(p$D >= 0.5)
  set.seed(seed + 1L)
  uvals <- stats::rnorm(length(unique(p$unit)))
  names(uvals) <- as.character(sort(unique(p$unit)))
  p$x1 <- uvals[as.character(p$unit)]
  set.seed(seed + 2L)
  uvals2 <- stats::rnorm(length(unique(p$unit)))
  names(uvals2) <- as.character(sort(unique(p$unit)))
  p$x2 <- uvals2[as.character(p$unit)]
  p
}

cat(sprintf("%-5s %-4s %-4s  %14s  %14s\n", "meth", "nu", "np",
            "max|att diff|", "max|IF diff|"))
for (method in c("IPW", "DR")) {
  for (nu in c(60L, 120L, 240L)) {
    for (np in c(8L, 12L)) {
      p <- make_panel(nu, np)
      fr <- didgpu_cs(p, "Y", "unit", "period", "D", covariates = c("x1","x2"),
                       est_method = method, bootstrap_reps = 0L,
                       backend = "r", verbose = FALSE)
      fc <- didgpu_cs(p, "Y", "unit", "period", "D", covariates = c("x1","x2"),
                       est_method = method, bootstrap_reps = 0L,
                       backend = "cuda", verbose = FALSE)
      att_diff <- max(abs(fr$att_gt$att - fc$att_gt$att), na.rm = TRUE)
      # IF diff: compare per-cell IF vectors.
      IFr <- attr(fr$att_gt, "IF_per_cell"); IFc <- attr(fc$att_gt, "IF_per_cell")
      if_diff <- 0
      for (i in seq_along(IFr)) {
        d <- suppressWarnings(max(abs(IFr[[i]]$IF - IFc[[i]]$IF)))
        if (is.finite(d) && d > if_diff) if_diff <- d
      }
      cat(sprintf("%-5s %-4d %-4d  %14.3e  %14.3e\n",
                  method, nu, np, att_diff, if_diff))
    }
  }
}
