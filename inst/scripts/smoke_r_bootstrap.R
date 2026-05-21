# Run a full bootstrap through the r-backend orchestrator. Verifies:
#   - Per-cell checkpointing works with the r-backend
#   - Aggregator handles r-backend cell shape
#   - Resume works
#   - Point estimates match reference-backend exactly (same panel)
#   - Bootstrap SEs are close to reference-backend SEs (Monte Carlo noise)
library(didgpu)
library(DIDmultiplegtDYN)

p <- didgpu_simulate_panel(
  n_units = 80L, n_periods = 18L, frac_treated = 0.6,
  min_treat_period = 4L, max_treat_period = 9L,
  tau_profile = c(0.5, 1.0, 1.2, 1.0),
  sigma = 0.4, seed = 17L
)

cat("=== r-backend, full bootstrap (effects=3, placebos=1, 20 reps) ===\n")
cdir_r <- tempfile("didgpu_r_boot_")
t0 <- Sys.time()
fit_r <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 3L, placebo = 1L,
                 bootstrap_reps = 20L, seed = 1L,
                 checkpoint_dir = cdir_r,
                 backend = "r", verbose = FALSE)
cat(sprintf("  wall: %.2fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("  cells written: %d\n", nrow(didgpu_load_checkpoint(cdir_r)$manifest)))
print(fit_r)

cat("\n=== reference-backend on same panel + seed for comparison ===\n")
cdir_ref <- tempfile("didgpu_ref_boot_")
t0 <- Sys.time()
fit_ref <- didgpu(p, "Y", "unit", "period", "D",
                   effects = 3L, placebo = 1L,
                   bootstrap_reps = 20L, seed = 1L,
                   checkpoint_dir = cdir_ref,
                   backend = "reference", verbose = FALSE)
cat(sprintf("  wall: %.2fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("\n=== compare point estimates ===\n")
diff_eff <- max(abs(fit_r$results$Effects[, "Estimate"] -
                    fit_ref$results$Effects[, "Estimate"]))
diff_pl  <- max(abs(fit_r$results$Placebos[, "Estimate"] -
                    fit_ref$results$Placebos[, "Estimate"]))
diff_ate <- abs(fit_r$results$ATE[1, "Estimate"] -
                fit_ref$results$ATE[1, "Estimate"])
cat(sprintf("  max diff Effects: %.2e\n", diff_eff))
cat(sprintf("  max diff Placebos: %.2e\n", diff_pl))
cat(sprintf("  diff ATE: %.2e\n", diff_ate))

cat("\n=== compare SEs (should be Monte Carlo close) ===\n")
r_se <- fit_r$results$Effects[, "SE"]
ref_se <- fit_ref$results$Effects[, "SE"]
cat("  r-backend SEs:    ", paste(sprintf("%.4f", r_se), collapse = "  "), "\n")
cat("  reference SEs:    ", paste(sprintf("%.4f", ref_se), collapse = "  "), "\n")
cat("  relative diffs:   ", paste(sprintf("%.1f%%", 100 * abs(r_se - ref_se) / ref_se),
                                  collapse = "  "), "\n")

cat("\n=== resume idempotency ===\n")
fit_r2 <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 1L,
                  bootstrap_reps = 20L, seed = 1L,
                  checkpoint_dir = cdir_r,
                  backend = "r", verbose = FALSE)
cat(sprintf("  effects match? %s\n",
            all.equal(fit_r$results$Effects[, "Estimate"],
                      fit_r2$results$Effects[, "Estimate"])))
