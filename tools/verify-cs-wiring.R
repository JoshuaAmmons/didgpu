library(didgpu)
cat("=== CUDA / R fallback wiring test for didgpu_cs ===\n")
cat("CUDA compiled in: ", didgpu_has_cuda_support(), "\n", sep = "")

# Simulate a tiny CS-friendly panel: 20 units, 6 periods, staggered adoption.
p <- didgpu_simulate_panel(n_units = 20L, n_periods = 6L,
                           tau_profile = c(0.5, 1.0), seed = 17L)

# Drop continuous treatment to integer for CS.
p$D <- as.integer(p$D >= 0.5)

cat("\n[1] backend = 'r' (baseline)...\n")
fit_r <- didgpu_cs(p, "Y", "unit", "period", "D",
                   est_method = "OR", aggregation = "event",
                   bootstrap_reps = 0L, backend = "r", verbose = FALSE)
cat("    ATT(g,t) rows: ", nrow(fit_r$att_gt), "\n", sep = "")
cat("    First ATT:    ", round(fit_r$att_gt$att[1], 6), "\n", sep = "")

cat("\n[2] backend = 'cuda' (should fall back to R since kernel returns -1)...\n")
fit_c <- didgpu_cs(p, "Y", "unit", "period", "D",
                   est_method = "OR", aggregation = "event",
                   bootstrap_reps = 0L, backend = "cuda", verbose = FALSE)
cat("    ATT(g,t) rows: ", nrow(fit_c$att_gt), "\n", sep = "")
cat("    First ATT:    ", round(fit_c$att_gt$att[1], 6), "\n", sep = "")

cat("\n[3] Equivalence check (CUDA-fallback path == R baseline)...\n")
stopifnot(nrow(fit_r$att_gt) == nrow(fit_c$att_gt))
stopifnot(all.equal(fit_r$att_gt$att, fit_c$att_gt$att,
                    tolerance = 1e-10))
cat("    PASS\n")

cat("\n[4] backend = 'auto' (auto-resolves to cuda if available)...\n")
fit_a <- didgpu_cs(p, "Y", "unit", "period", "D",
                   est_method = "OR", aggregation = "event",
                   bootstrap_reps = 0L, backend = "auto", verbose = FALSE)
cat("    Resolved backend: ", fit_a$args$backend, "\n", sep = "")
stopifnot(all.equal(fit_r$att_gt$att, fit_a$att_gt$att,
                    tolerance = 1e-10))
cat("    PASS\n")

cat("\n=== All Phase-1 #79 wiring checks PASS ===\n")
