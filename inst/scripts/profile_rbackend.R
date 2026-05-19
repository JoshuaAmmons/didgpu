# Profile the r-backend to find hotspots.
library(didgpu)
library(profvis)

p <- didgpu_simulate_panel(
  n_units = 500L, n_periods = 100L, frac_treated = 0.6,
  tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2),
  sigma = 0.4, seed = 17L
)

# Warmup.
invisible(didgpu(p, "Y", "unit", "period", "D",
                  effects = 5L, placebo = 2L,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE))

# Profile a single point-estimate fit.
pf <- profvis::profvis({
  for (rep in 1:5) {
    didgpu(p, "Y", "unit", "period", "D",
            effects = 5L, placebo = 2L,
            bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  }
}, interval = 0.005)

# Save the html so we can inspect it.
htmlwidgets::saveWidget(pf, "inst/scripts/profile.html", selfcontained = TRUE)
cat("Profile saved to inst/scripts/profile.html\n")
