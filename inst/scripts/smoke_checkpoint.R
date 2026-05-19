# Smoke-test the checkpoint layer in isolation.
suppressPackageStartupMessages({ library(data.table) })
source("R/checkpoint.R")
source("R/simulate.R")

cdir <- tempfile("didgpu_chk_")

# Build a fake panel and hash it. Both functions are in global env after source().
p <- didgpu_simulate_panel(n_units = 30L, n_periods = 8L, seed = 11L)
ph <- .panel_hash(p, "Y", "unit", "period", "D")
cat("panel hash:", ph, "\n")

meta <- list(
  panel_hash = ph,
  seed = 7L,
  bootstrap_reps = 10L,
  effects = 5L,
  placebo = 3L,
  outcome = "Y", group = "unit", time = "period", treatment = "D",
  package_version = "0.0.0.9000"
)

didgpu_init_checkpoint(cdir, meta)
cat("init OK ->", cdir, "\n")
cat("contents:\n"); print(list.files(cdir, recursive = TRUE))
cat("meta.json:\n"); cat(readLines(file.path(cdir, "meta.json")), sep = "\n")

# Simulate writing 3 cells.
.save_cell(cdir, b = 0L, value = list(coef = c(Effect_1 = 0.5, Effect_2 = 1.0, Placebo_1 = 0.1)),
           wall_seconds = 1.2)
.save_cell(cdir, b = 1L, value = list(coef = c(Effect_1 = 0.48, Effect_2 = 1.05, Placebo_1 = 0.12)),
           wall_seconds = 1.3)
.save_cell(cdir, b = 3L, value = list(coef = c(Effect_1 = 0.51, Effect_2 = 0.98, Placebo_1 = 0.08)),
           wall_seconds = 1.1)
cat("\nsaved 3 cells (b = 0, 1, 3 -- note 2 is skipped)\n")

chk <- didgpu_load_checkpoint(cdir)
cat("\n--- manifest after 3 saves ---\n")
print(chk$manifest)
cat("\n--- todo (bootstrap_reps = 10) ---\n")
cells_todo <- get(".cells_todo")
print(cells_todo(chk$manifest, 10L))

cat("\n--- aggregate cells ---\n")
agg <- didgpu_aggregate_cells(cdir)
print(agg)
cat("\nall good\n")
